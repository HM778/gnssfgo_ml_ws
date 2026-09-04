#include "ProcessingNodeExtra.hpp"

#include "../include/models/SelfAdjustedFactorGraph.hpp"
#include "../include/tools/tic_toc.h"
#include "../include/tools/transformer_bridge.hpp"

#include <atomic>
#include <algorithm>
#include <chrono>
#include <mutex>

class SelfAdjustedTR : public ProcessingNodeExtra
{
    TicToc OptTime;

    std::vector<gnss_comm::ObsPtr> fgo_input_raw;
    Eigen::Vector3d cov_enu;
    
    std::thread optimizationThread;

    double last_ingested_time_frame = -1.0;
    std::map<int, sv_info> last_sv_info_map;    //更新前的卫星信息
    double max_running_time_ms;                 //单次优化最大运行时间
    int max_psr_factors_per_epoch = 24;         //每轮优化最大伪距因子数
    int min_tr_factor_num = 20;                 //每轮优化最少TR因子数，避免约束过弱
    bool SAME_RELIABLE,SAME_TIME_WEIGHT;        //对照组开关参数
    
    // ===== OSQA Transformer 质量评估桥接 相关参数初始化=====
#ifdef ENABLE_TRANSFORMER_BRIDGE
    bool osqa_enabled = false;
    std::string osqa_input_path;
    std::string osqa_output_path;
    double osqa_quality_threshold = 0.3;
    int osqa_epoch_counter = 0;
    double last_osqa_export_time_frame = -1.0;
#endif

    std::mutex m_factor_graph_mux;
    
private:
    SelfAdjustedFactorGraph factor_graph;


public:
    SelfAdjustedTR() 
    {   
        // Initialization
        current_sys_time.sec = -1;

        // params setting
        loadParams();
        nh.param<double>("max_running_time_ms",max_running_time_ms,50.0);
        nh.param<int>("max_psr_factors_per_epoch", max_psr_factors_per_epoch, 12);
        nh.param<int>("min_tr_factor_num", min_tr_factor_num, 20);

        nh.param<bool>("same_time_weight", SAME_TIME_WEIGHT, false);
        nh.param<bool>("same_reliable", SAME_RELIABLE, false);

        // 均匀频率输出: launch 参数 output_rate_hz (Hz), 上限 20
        double output_rate_hz = 1.0;
        nh.param<double>("output_rate_hz", output_rate_hz, 1.0);
        StartOutputThread(output_rate_hz);
        
        ROS_ERROR("Set time_dicount/reliable method: %d / %d", SAME_TIME_WEIGHT,SAME_RELIABLE);

        // ===== OSQA Transformer 质量评估桥接 相关参数加载=====
#ifdef ENABLE_TRANSFORMER_BRIDGE
        nh.param<bool>("osqa_enabled", osqa_enabled, false);
        if (osqa_enabled)
        {
            nh.param<std::string>("osqa_input_path", osqa_input_path,
                ros::package::getPath("gnssfgo") + "/../osqa_input.jsonl");
            nh.param<std::string>("osqa_output_path", osqa_output_path,
                ros::package::getPath("gnssfgo") + "/../osqa_output.jsonl");
            nh.param<double>("osqa_quality_threshold", osqa_quality_threshold, 0.3);
            ROS_INFO("[OSQA Bridge] Enabled. Input: %s, Output: %s, Threshold: %.2f",
                     osqa_input_path.c_str(), osqa_output_path.c_str(), osqa_quality_threshold);
        }
        else
        {
            ROS_INFO("[OSQA Bridge] Disabled. Running without transformer quality evaluation.");
        }
#endif
        factor_graph.windowSize = windowSize;
        factor_graph.MARGINAL_ENABLE = marginal_enable;
        factor_graph.SAME_RELIABLE = SAME_RELIABLE;
        factor_graph.SAME_TIME_WEIGHT = SAME_TIME_WEIGHT;
        factor_graph.max_psr_factors_per_epoch_ = std::max(4, max_psr_factors_per_epoch);
        factor_graph.MinTRFactorNum = std::max(10, min_tr_factor_num);
        
        InitialSubTopics();
        InitialPubTopics();
        // callback队列
        StartSpinners(); 
        optimizationThread = std::thread(&SelfAdjustedTR::Optimization, this);
    }

    void Optimization()
    {
        {
            std::lock_guard<std::mutex> lk(m_factor_graph_mux);
            factor_graph.getOptimizationSettingInfo();
            factor_graph.setupSolverOptions();
            factor_graph.setupLossFunction("Huber");
        }
        
        ROS_INFO("\033[1;32mFactor graph running....\033[0m");
        
        bool first_opt = true;
        Eigen::Matrix<double, 3, 1> last_pos_ecef = Eigen::Vector3d::Zero();

        while(ros::ok())
        {
            if(waitForNewData(50))  // blocking wait instead of busy-spin
            {
                synchronizeAndInput();
            }

            int data_stream_size = 0;
            bool has_new_data = false;
            {
                std::lock_guard<std::mutex> lk(m_factor_graph_mux);
                data_stream_size = factor_graph.getDataStreamSize();
                has_new_data = factor_graph.has_new_data;
            }

            if(data_stream_size>0) //滑动窗口内有数据
            {
                ROS_DEBUG_THROTTLE(1.0, "Current valid data group size in slip window: %d", data_stream_size);
                if(has_new_data)
                {
                    OptTime.tic();
                    double solve_move_m = 0.0;  // 本轮优化对最新状态的实际改动量
                    double spp_gap_m = 0.0;     // 发布解与 SPP 初值的偏差

                    Eigen::Matrix<double, 3, 1> init_pos_ecef;
                    Eigen::Vector3d ref_llh;
                    {
                        std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
                        init_pos_ecef = latest_pos_ecef;
                        ref_llh = enu_ref_llh;
                    }

                    {
                        std::lock_guard<std::mutex> lk(m_factor_graph_mux);

                        if(factor_graph.time_frame_now < 0)
                        {
                            continue;
                        }

                        // Rebuild the factor graph every iteration to avoid stale residual blocks
                        // constraining reused state memory after the sliding window moves.
                        factor_graph.resetProblem();

                        // ===== OSQA Transformer 质量评分导入 =====
                        // 在构建因子图之前，读取 OSQA 对之前 epoch 的卫星信号质量评估结果
                        // 质量评分将用于调整伪距因子和 TR 双差因子的权重
#ifdef ENABLE_TRANSFORMER_BRIDGE
                        if (osqa_enabled)
                        {
                            // 使用带 time_frame 匹配的读取重载:
                            // OSQA 处理存在滞后, 盲取最后一行会把评分套到错误的历元上;
                            // 精确匹配失败时回退到 60s 内最近的历元。
                            std::map<int, transformer_bridge::SatQualityInfo> quality_scores =
                                transformer_bridge::readLatestQualityScores(
                                    osqa_output_path, factor_graph.time_frame_now);
                            if (!quality_scores.empty())
                            {
                                std::map<int, double> simple_scores;
                                for (const auto &kv : quality_scores)
                                {
                                    simple_scores[kv.first] = kv.second.quality;  //只使用了融合后的质量评分
                                }
                                factor_graph.setQualityScores(simple_scores);
                                factor_graph.osqa_quality_threshold_ = osqa_quality_threshold;
                            }
                        }
#endif

                        if(factor_graph.first_run)
                        {
                            //为待估计的位置分配内存，保存窗口内估计状态
                            factor_graph.setupStateMemory();
                            factor_graph.first_run = false;
                        }

                        // Populate initial guesses for states already optimized before.
                        factor_graph.initializeOldGraph();

                        // 为窗口内新增观测组待估计变量设置初值为伪距定位结果
                        // 在外部进行伪距定位，避免传入ephem_array
                        factor_graph.initializeGraph(init_pos_ecef);
                        
                        // Add parameter blocks for the current sliding window
                        factor_graph.addAllParameterBlocksToGraph();
                        factor_graph.configureStateParameterizationForWindow();
                        factor_graph.anchorFirstStateInWindow(factor_graph.measSize > 1);

                        // 纯伪距 + 多普勒 因子图
                        factor_graph.addDopplerFactors();
                        factor_graph.addPsrFactors();

                        // factor_graph.addTRDDCPFactors();

                        // 求解前最新状态初值, 用于衡量本轮优化对解的实际改动量
                        const Eigen::Vector3d pre_solve_pos(factor_graph.state_array[factor_graph.measSize - 1][0],
                                                            factor_graph.state_array[factor_graph.measSize - 1][1],
                                                            factor_graph.state_array[factor_graph.measSize - 1][2]);

                        ceres::Solver::Options local_options = factor_graph.options;
                        if (factor_graph.measSize <= 2)
                        {
                            local_options.linear_solver_type = ceres::DENSE_QR;
                            local_options.num_threads = 1;
                            local_options.num_linear_solver_threads = 1;
                        }
                        ceres::Solve(local_options, &factor_graph.problem, &factor_graph.summary);
                        const Eigen::Vector3d solved_pos(factor_graph.state_array[factor_graph.measSize - 1][0],
                                                         factor_graph.state_array[factor_graph.measSize - 1][1],
                                                         factor_graph.state_array[factor_graph.measSize - 1][2]);
                        solve_move_m = (solved_pos - pre_solve_pos).norm();
                        // 发散防护: 个别历元求解器会收敛到公里级错误盆地
                        // (观测到与接收机钟跳步相关的约束重构), 回退为 SPP 初值
                        factor_graph.guardLatestStateAgainstDivergence(init_pos_ecef);
                        spp_gap_m = (factor_graph.getLatestPosECEF() - init_pos_ecef).norm();
                        // factor_graph.solveFloatAmbiguity();
                        factor_graph.saveGraphStateToVector(true);
                        factor_graph.CPresidualsUpdate();
                        // 真实后验残差汇总 (double 评估, 区别于因子内 Jet 误打印)
                        factor_graph.logPostFitResidualSummary();

                        // ===== OSQA Transformer 数据导出 =====
                        // 在优化完成后，将所有预处理数据和因子残差导出到 JSONL 文件
                        // 供外部 OSQA Transformer 程序读取和评估卫星信号质量
#ifdef ENABLE_TRANSFORMER_BRIDGE
                        if (osqa_enabled)
                        {
                            exportEpochToOSQA(factor_graph);
                        }
#endif

                        if (factor_graph.getCovarianceMatrixOfLatestEpoch())
                        {
                            factor_graph.printLatestPosCovarianceENU(cov_enu);
                        }

                        if(first_opt)
                        {
                            //设置enu参考点，检查位置有效性
                            if(posValid(ref_llh))
                            {
                                factor_graph.setupReferencePoint(ref_llh(0), ref_llh(1), ref_llh(2));
                                first_opt = false;
                            }
                        }
                        

                        const double curr_time_sec = factor_graph.time_frame_now / 10.0;
                        // fgo_enu坐标输出
                        update_fgo_enu(factor_graph.getLatestPosENU()(0), factor_graph.getLatestPosENU()(1), factor_graph.getLatestPosENU()(2));
                        update_fgo_llh(factor_graph.getLatestPosLLH()(0), factor_graph.getLatestPosLLH()(1), factor_graph.getLatestPosLLH()(2));
                        
                        last_pos_ecef(0)=factor_graph.getLatestPosECEF()(0);
                        last_pos_ecef(1)=factor_graph.getLatestPosECEF()(1);
                        last_pos_ecef(2)=factor_graph.getLatestPosECEF()(2);
                        // 移除滑动窗口外的状态和观测数据（输入时已维护，这里保留为兜底）
                        factor_graph.removeStatesOutsideSlidingWindow();

                        factor_graph.time_frame_last = factor_graph.time_frame_now;
                        factor_graph.has_new_data = false;
                    }
                    double run_time = OptTime.toc();
                    
                    factor_graph.resizeMaxTRFactorNum(run_time, max_running_time_ms);

                    savetofile(factor_graph.getLatestPosLLH(), factor_graph.getLatestPosENU(), cov_enu, factor_graph.time_frame_now, factor_graph.TRFactorCount, run_time);
                    const int num_vars = factor_graph.measSize * state_size;
                    printf("[ OPT SUMMARY] Window=%d epochs | Variables=%d | TRDDCP=%d PSR=%d DOPPLER=%d | Time=%.1f ms\n",
                           factor_graph.measSize, num_vars,
                           factor_graph.TRFactorCount,
                           factor_graph.last_added_psr_factor_count,
                           factor_graph.last_added_doppler_factor_count,
                           run_time);
                    printf("[ OPT EFFECT] cost %.4g -> %.4g | |Δpos_solve|=%.4f m | FGO-SPP gap=%.4f m\n",
                           factor_graph.summary.initial_cost, factor_graph.summary.final_cost,
                           solve_move_m, spp_gap_m);
                }
            }
            // ros::Duration(0.005).sleep();
        }
        {
            std::lock_guard<std::mutex> lk(m_factor_graph_mux);
            factor_graph.gnss_raw_map.clear();
            factor_graph.doppler_map.clear();
        }
        ROS_INFO("\033[1;32mFactor graph optimization thread finished.\033[0m");
    }

    void synchronizeAndInput()
    {
        // 伪距、速度均可用
        if (!can_solve())
        {
            hasNewData.store(false, std::memory_order_release);
            return;
        }

        // 观测数据的局部保存变量
        std::vector<gnss_comm::ObsPtr> local_fgo_input_raw;
        std::vector<gnss_comm::EphemBasePtr> local_ephem_array;
        gnss_comm::gtime_t local_sys_time;
        double local_gpst_sec = -1.0;
        Eigen::Vector3d local_latest_pos_ecef;
        Eigen::Vector3d local_latest_pos_llh;
        Eigen::Vector3d local_enu_ref_llh;
        bool local_enu_ref_set = false;
        nav_msgs::Odometry local_dop_meas;

        // 转存数据到局部变量，减少锁的持有时间
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            local_fgo_input_raw = meas;  // all systems, all frequencies
            local_ephem_array = ephem_array;

            local_sys_time = current_sys_time;
            local_gpst_sec = current_gpst_sec;
            local_latest_pos_ecef = latest_pos_ecef;  //最小二乘结果
            local_latest_pos_llh = latest_pos_llh;
            local_enu_ref_llh = enu_ref_llh;
            local_enu_ref_set = enu_ref_set;
            local_dop_meas = dop_vel_meas;
           
        }

        if (!gpsTimeValid(local_gpst_sec) && !sysTimeValid(local_sys_time))
        {
            hasNewData.store(false, std::memory_order_release);
            return;
        }

        // Incrementally update satellite states for this epoch:
        // start from the last available result, then refresh/add current satellites.
        // 增量更新卫星信息
        std::map<int, sv_info> local_sv_info_map = last_sv_info_map;
        
        if (gpsTimeValid(local_gpst_sec) || sysTimeValid(local_sys_time))
        {
            for (const auto &ephem_base : local_ephem_array)
            {
                // Try Keplerian ephemeris first (GPS, BDS, Galileo)
                gnss_comm::EphemPtr ephem = std::dynamic_pointer_cast<gnss_comm::Ephem>(ephem_base);
                if (ephem)
                {
                    const int sat_i = ephem->sat;
                    sv_info sv;
                    sv.sat_prn = sat_i;
                    sv.sys = satsys(sat_i, NULL);
                    sv.ddt = eph2svdt(local_sys_time, ephem);
                    sv.pos = eph2pos(local_sys_time, ephem, &sv.dt);
                    sv.vel = eph2vel(local_sys_time, ephem, &sv.ddt);

                    // Set L1/L2 frequencies and wavelengths from observation data
                    for (const auto& obs : local_fgo_input_raw)
                    {
                        if (!obs || obs->sat != sat_i) continue;
                        int l1_idx = -1, l2_idx = -1;
                        sv.freq_l1 = L1_freq(obs, &l1_idx);
                        sv.lamda_l1 = (sv.freq_l1 > 0.0) ? LIGHT_SPEED / sv.freq_l1 : 0.0;
                        sv.freq_l2 = L2_freq(obs, &l2_idx);
                        sv.lamda_l2 = (sv.freq_l2 > 0.0) ? LIGHT_SPEED / sv.freq_l2 : 0.0;
                        break;
                    }

                    if(!satPosValid(sv.pos[0], sv.pos[1], sv.pos[2]))
                    {
                        sv.avaliable = false;
                        // Mark the satellite as unavailable,discard updating
                        continue;
                    }
                                        
                    if (local_enu_ref_set && posValid(local_latest_pos_llh))
                    {
                        const double d_sr = (local_latest_pos_ecef - sv.pos).norm();
                        const double tau = d_sr / LIGHT_SPEED;
                        const gnss_comm::gtime_t transmit_time = gnss_comm::time_add(local_sys_time, -tau);
                        sv.pos = eph2pos(transmit_time, ephem, &sv.dt);

                        const Eigen::Vector3d rev2sat_ecef = (sv.pos - local_latest_pos_ecef).normalized();
                        const Eigen::Vector3d rev2sat_enu = gnss_comm::ecef2enu(local_enu_ref_llh, rev2sat_ecef);
                        sv.azimuth = (rev2sat_enu.head<2>().norm() < 1e-12) ? 0.0 : atan2(rev2sat_enu.x(), rev2sat_enu.y());
                        sv.elevation = asin(std::max(-1.0, std::min(1.0, rev2sat_enu.z())));

                        sv.avaliable = true;
                        local_sv_info_map[sat_i] = sv;
                    }

                    continue;
                }

                // Try GLONASS ephemeris
                gnss_comm::GloEphemPtr gephem = std::dynamic_pointer_cast<gnss_comm::GloEphem>(ephem_base);
                if (gephem)
                {
                    const int sat_i = gephem->sat;
                    sv_info sv;
                    sv.sat_prn = sat_i;
                    sv.sys = satsys(sat_i, NULL);
                    sv.ddt = geph2svdt(local_sys_time, gephem);
                    sv.pos = geph2pos(local_sys_time, gephem, &sv.dt);
                    sv.vel = geph2vel(local_sys_time, gephem, &sv.ddt);

                    // GLONASS FDMA: frequencies depend on freqo channel number
                    sv.freq_l1 = FREQ1_GLO + gephem->freqo * DFRQ1_GLO;
                    sv.lamda_l1 = LIGHT_SPEED / sv.freq_l1;
                    sv.freq_l2 = FREQ2_GLO + gephem->freqo * DFRQ2_GLO;
                    sv.lamda_l2 = LIGHT_SPEED / sv.freq_l2;

                    if (!std::isfinite(sv.pos[0]) || sv.pos[0] == 0.0)
                    {
                        sv.avaliable = false;
                        // Mark the satellite as unavailable, discard updating
                        continue;
                    }

                    if (local_enu_ref_set && posValid(local_latest_pos_llh))
                    {
                        const double d_sr = (local_latest_pos_ecef - sv.pos).norm();
                        const double tau = d_sr / LIGHT_SPEED;
                        const gnss_comm::gtime_t transmit_time = gnss_comm::time_add(local_sys_time, -tau);
                        sv.pos = geph2pos(transmit_time, gephem, &sv.dt);

                        const Eigen::Vector3d rev2sat_ecef = (sv.pos - local_latest_pos_ecef).normalized();
                        const Eigen::Vector3d rev2sat_enu = gnss_comm::ecef2enu(local_enu_ref_llh, rev2sat_ecef);
                        sv.azimuth = (rev2sat_enu.head<2>().norm() < 1e-12) ? 0.0 : atan2(rev2sat_enu.x(), rev2sat_enu.y());
                        sv.elevation = asin(std::max(-1.0, std::min(1.0, rev2sat_enu.z())));

                        sv.avaliable = true;
                        local_sv_info_map[sat_i] = sv;
                    }                    
                }
            }
        }

        last_sv_info_map = local_sv_info_map;

        const double base_gpst_sec = (gpsTimeValid(local_gpst_sec))? local_gpst_sec: (local_sys_time.time + local_sys_time.sec);
        const double time_frame = base_gpst_sec * 10.0;  //0.1s为单位的时间帧
        
        // 防止重复输入同一时刻的数据（如多次收到相同的观测数据或星历数据更新时）
        if (std::isfinite(last_ingested_time_frame) && std::abs(time_frame - last_ingested_time_frame) < 1e-6)
        {
            hasNewData.store(false, std::memory_order_release);
            return;
        }

        {
            std::lock_guard<std::mutex> lk(m_factor_graph_mux);
            // Drop stale/out-of-order epochs to keep the optimizer on the newest stream.
            //  丢弃过时/乱序的观测数据，保持优化器处理最新的数据流
            if (factor_graph.time_frame_now > 0.0 && time_frame <= factor_graph.time_frame_now)
            {
                hasNewData.store(false, std::memory_order_release);
                return;
            }
            
            // Different input standards:
            // - Doppler is used for continuous propagation (looser gate).
            // - DD pseudorange/TDCP needs at least 4 common satellites.

            const bool gnss_obs_ok = ((int)local_fgo_input_raw.size() >= 4) && ((int)local_ephem_array.size() >= 4);

            fgo_input_raw = gnss_obs_ok ? local_fgo_input_raw : std::vector<gnss_comm::ObsPtr>();

            factor_graph.input_sv_info_window_map(local_sv_info_map, time_frame);
            // factor_graph.sv_info_window_map[time_frame] = local_sv_info_map;

            // Always create an epoch in gnss_raw_map (even if empty) so the state stream is continuous.
            factor_graph.input_gnss_raw_data(fgo_input_raw, time_frame);
            printf("Ingested new data for time frame %.2f with %lu observations.\n", time_frame, fgo_input_raw.size());
            {
                // Count per-system and per-frequency statistics
                int cnt_gps=0, cnt_glo=0, cnt_gal=0, cnt_bds=0;
                int cnt_l1=0, cnt_l2=0, cnt_dual=0;
                for (const auto& o : fgo_input_raw) {
                    if (!o) continue;
                    const int sys = satsys(o->sat, nullptr);
                    if (sys==SYS_GPS) cnt_gps++; else if (sys==SYS_GLO) cnt_glo++;
                    else if (sys==SYS_GAL) cnt_gal++; else if (sys==SYS_BDS) cnt_bds++;
                    int l1=-1,l2=-1; L1_freq(o,&l1); L2_freq(o,&l2);
                    if (l1>=0) cnt_l1++; if (l2>=0) cnt_l2++;
                    if (l1>=0 && l2>=0) cnt_dual++;
                }
                printf("  Systems: G=%d R=%d E=%d C=%d | L1=%d L2=%d Dual=%d\n",
                       cnt_gps, cnt_glo, cnt_gal, cnt_bds, cnt_l1, cnt_l2, cnt_dual);
            }
            factor_graph.input_ephem_data(local_ephem_array, time_frame);
            factor_graph.input_doppler_data(local_dop_meas, time_frame);

            factor_graph.has_new_data = true;

            last_ingested_time_frame = time_frame;

            // Keep sliding window size bounded before Optimization() reads measSize/state_array.
            factor_graph.removeStatesOutsideSlidingWindow();
            
        }

        hasNewData.store(false, std::memory_order_release);
    }

#ifdef ENABLE_TRANSFORMER_BRIDGE
    // ===== OSQA Transformer 数据导出 =====
    // 导出优化后的后验残差, 而非原始伪距减几何距离。
    // 使用与 addPsrFactors() 完全相同的校正量:
    //   est = geo_range + sagnac − sat_clock + TGD + iono + tropo + rx_clock_opt
    //   post_fit_residual = est − raw_psr   (单位: 米)
    void exportEpochToOSQA(SelfAdjustedFactorGraph &fg)
    {
        osqa_epoch_counter++;

        double time_frame = fg.time_frame_now;
        double gpst_sec = time_frame / 10.0;

        if (std::abs(time_frame - last_osqa_export_time_frame) < 5e-2)
        {
            return;
        }
        last_osqa_export_time_frame = time_frame;

        int latest_idx = fg.measSize - 1;
        if (latest_idx < 0 || latest_idx >= static_cast<int>(fg.state_array.size()))
        {
            return;
        }

        // ---- 优化后的接收机状态 ----
        const double rx_x = fg.state_array[latest_idx][0];
        const double rx_y = fg.state_array[latest_idx][1];
        const double rx_z = fg.state_array[latest_idx][2];
        const double clk_gps = fg.state_array[latest_idx][3];
        const double clk_glo = fg.state_array[latest_idx][4];
        const double clk_gal = fg.state_array[latest_idx][5];
        const double clk_bds = fg.state_array[latest_idx][6];
        Eigen::Vector3d receiver_ecef(rx_x, rx_y, rx_z);
        Eigen::Vector3d receiver_enu(0, 0, 0);
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            receiver_enu = latest_pos_enu;
        }

        // ---- 最新 epoch 数据 ----
        auto raw_iter = fg.gnss_raw_map.rbegin();
        if (raw_iter == fg.gnss_raw_map.rend()) return;
        const auto &observations = raw_iter->second;
        double epoch_time = raw_iter->first;

        // ---- 获取星历 ----
        std::vector<gnss_comm::EphemBasePtr> epoch_ephems;
        {
            auto eph_it = fg.ephems_map.find(epoch_time);
            if (eph_it != fg.ephems_map.end()) epoch_ephems = eph_it->second;
        }

        // ---- 卫星信息 ----
        std::map<int, sv_info> epoch_sv_info;
        auto sv_iter = fg.sv_info_window_map.find(epoch_time);
        if (sv_iter != fg.sv_info_window_map.end()) epoch_sv_info = sv_iter->second;

        // ---- 用优化后状态调用 buildCorrectedPseudorangeMeasurements ----
        // 获取与 addPsrFactors 完全相同的校正量 (sat_clock, TGD, iono, tropo)
        std::vector<double> iono_params(8, 0.0);
        std::vector<gnss_comm_extra::CorrectedPseudorangeMeasurement> corrected =
            gnss_comm_extra::buildCorrectedPseudorangeMeasurements(
                observations, epoch_ephems, receiver_ecef, true, iono_params);

        // ---- 计算后验残差 ----
        // 只保留 L1: buildCorrectedPseudorangeMeasurements 对每颗卫星依次 push L1/L2
        // 两条测量, 若不做过滤 map 会被 L2 残差覆盖, 与 "psr_residual_l1" 标签不符。
        // 所有残差统一量纲为米。
        std::map<int, double> psr_residual_map;
        for (const auto &m : corrected)
        {
            if (!m.valid || m.sat_sys == "Unknown") continue;
            if (m.freq != 1) continue;

            int sat = m.sat;
            double geo_range = (receiver_ecef - m.sat_pos).norm();
            double sagnac = OMGE_ / CLIGHT_ * (m.sat_pos.x() * rx_y - m.sat_pos.y() * rx_x);
            double sv_clk_m = m.sv_dt_sec * CLIGHT_;
            double tgd_m = m.tgd_sec * CLIGHT_;
            double iono_m = m.ion_delay_m;
            double tropo_m = m.tro_delay_m;

            // 接收机钟差 (米) — 与 pseudorangeFactor 中的 state[3..6] 对应
            double rx_clk = 0.0;
            if (m.sat_sys == "GPS")      rx_clk = clk_gps;
            else if (m.sat_sys == "GLONASS") rx_clk = clk_glo;
            else if (m.sat_sys == "Galileo")  rx_clk = clk_gal;
            else if (m.sat_sys == "BeiDou")   rx_clk = clk_bds;

            // est = geo + sagnac - sv_clk + tgd + iono + tropo + rx_clk
            double est = geo_range + sagnac - sv_clk_m + tgd_m + iono_m + tropo_m + rx_clk;
            double post_fit_residual = est - m.pseudorange;  // 米, 未加权

            // Guard: 把钟差和大气延迟主导的粗差截断, 保留真实测量噪声尺度
            if (std::isfinite(post_fit_residual) && std::abs(post_fit_residual) < 1e4)
            {
                psr_residual_map[sat] = post_fit_residual;
            }
        }

        // ---- 因子残差 (量纲统一为米: psr 伪距后验残差 / dop_cp 载波相位变化残差) ----
        // TR DD PR/CP 双差残差无法计算, 已在导出与分析侧停用。
        std::map<int, std::map<std::string, double>> factor_residuals;
        for (const auto &obs : observations)
        {
            if (!obs) continue;
            const int sat = static_cast<int>(obs->sat);
            std::map<std::string, double> res;
            auto fg_res_it = factor_graph.residuals.find(sat);
            if (fg_res_it != factor_graph.residuals.end())
            {
                res = fg_res_it->second;
            }
            auto psr_it = psr_residual_map.find(sat);
            if (psr_it != psr_residual_map.end())
            {
                res["psr"] = psr_it->second;
            }
            factor_residuals[sat] = res;
        }

        bool ok = transformer_bridge::exportEpochData(
            osqa_input_path, gpst_sec, time_frame,
            receiver_ecef, receiver_enu, observations,
            epoch_sv_info,
            fg.sat_lock_count_l1, fg.sat_lock_count_l2,
            psr_residual_map, factor_residuals);

        if (ok)
        {
            ROS_DEBUG("[OSQA Bridge] Exported epoch %d (time_frame=%.1f, %zu satellites) to %s",
                      osqa_epoch_counter, time_frame, observations.size(), osqa_input_path.c_str());
        }
        else
        {
            ROS_WARN("[OSQA Bridge] Failed to export epoch %d to %s",
                     osqa_epoch_counter, osqa_input_path.c_str());
        }
    }
#endif


    ~SelfAdjustedTR()
    {
        StopOutputThread();
        freeSpinner();
        if (optimizationThread.joinable()) optimizationThread.join();
        saveEphems();
    }

    bool savetofile(Eigen::Vector3d lla, Eigen::Vector3d enu, Eigen::Vector3d cov_enu, double gpst_sec, int factor_count,double opt_time)
    {
        std::ofstream outfile(output_filename, std::ios::out | std::ios::app);
        if (!outfile.is_open())
        {
            ROS_ERROR("Failed to open output file: %s", output_filename.c_str());
            return false;
        }

        // Write header if file is new
        if (outfile.tellp() == 0)
        {
            outfile << "timestamp,latitude,longitude,altitude,enu_x,enu_y,enu_z,cov_xx,cov_yy,cov_zz,trddcp_factor_count,opt_time\n";
        }

        // Write data
        outfile << std::fixed << std::setprecision(6) << gpst_sec << ","
                << lla[0] << ","
                << lla[1] << ","
                << lla[2] << ","
                << enu[0] << ","
                << enu[1] << ","
                << enu[2] << ","
                << cov_enu[0] << "," 
                << cov_enu[1] << "," 
                << cov_enu[2] << ","
                << factor_count << ","
                << opt_time << "\n";

        outfile.close();
        return true;
    }
   
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "self_adjusted_node"); 
    ROS_INFO("\033[1;32m----> trbinfgo_node Started (4-system L1+L2).\033[0m"); 
    // ...
    SelfAdjustedTR sa_node;
    ros::waitForShutdown();
    return 0;
}
