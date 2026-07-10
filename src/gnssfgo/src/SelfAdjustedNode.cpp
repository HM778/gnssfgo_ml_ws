#include "ProcessingNodeExtra.hpp"

#include "../include/models/SelfAdjustedFactorGraph.hpp"
#include "../include/tools/tic_toc.h"

#include <atomic>
#include <chrono>
#include <mutex>

class SelfAdjustedTR : public ProcessingNodeExtra
{
    TicToc OptTime;

    std::vector<gnss_comm::ObsPtr> fgo_input_raw;
    Eigen::Vector3d cov_enu;
    
    std::thread optimizationThread;

    double last_ingested_time_frame = -1.0;
    std::map<int, sv_info> last_sv_info_map;
    double max_running_time_ms;
    double min_output_dt = 0.0;
    bool SAME_RELIABLE,SAME_TIME_WEIGHT;
    double last_output_time_sec = -1.0;

    std::mutex m_factor_graph_mux;
    std::atomic<bool> history_worker_running{false};
    std::atomic<bool> history_update_pending{false};

    
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

        nh.param<bool>("same_time_weight", SAME_TIME_WEIGHT, false);
        nh.param<bool>("same_reliable", SAME_RELIABLE, false);
        
        ROS_ERROR("Set time_dicount/reliable method: %d / %d", SAME_TIME_WEIGHT,SAME_RELIABLE);
        factor_graph.windowSize = windowSize;
        factor_graph.MARGINAL_ENABLE = marginal_enable;
        factor_graph.SAME_RELIABLE = SAME_RELIABLE;
        factor_graph.SAME_TIME_WEIGHT = SAME_TIME_WEIGHT;
        
        InitialSubTopics();
        InitialPubTopics();
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

                        factor_graph.addTRDDCPFactors();

                        ceres::Solver::Options local_options = factor_graph.options;
                        if (factor_graph.measSize <= 2)
                        {
                            local_options.linear_solver_type = ceres::DENSE_QR;
                            local_options.num_threads = 1;
                            local_options.num_linear_solver_threads = 1;
                        }
                        ceres::Solve(local_options, &factor_graph.problem, &factor_graph.summary);
                        // factor_graph.solveFloatAmbiguity();
                        factor_graph.saveGraphStateToVector(true);

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
                        last_output_time_sec = curr_time_sec;
                        
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
        // 伪距、速度至少一者可用
        if (!can_solve())
        {
            hasNewData.store(false, std::memory_order_release);
            return;
        }

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
                        sv.lamda = (sv.freq_l1 > 0.0) ? LIGHT_SPEED / sv.freq_l1 : 0.0;
                        sv.freq_l2 = L2_freq(obs, &l2_idx);
                        sv.lamda_l2 = (sv.freq_l2 > 0.0) ? LIGHT_SPEED / sv.freq_l2 : 0.0;
                        break;
                    }

                    if (!std::isfinite(sv.pos[0]) || sv.pos[0] == 0.0)
                    {
                        sv.avaliable = false;
                        local_sv_info_map[sat_i] = sv;
                        continue;
                    }
                    sv.avaliable = true;

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
                    }

                    local_sv_info_map[sat_i] = sv;
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
                    sv.lamda = LIGHT_SPEED / sv.freq_l1;
                    sv.freq_l2 = FREQ2_GLO + gephem->freqo * DFRQ2_GLO;
                    sv.lamda_l2 = LIGHT_SPEED / sv.freq_l2;

                    if (!std::isfinite(sv.pos[0]) || sv.pos[0] == 0.0)
                    {
                        sv.avaliable = false;
                        local_sv_info_map[sat_i] = sv;
                        continue;
                    }
                    sv.avaliable = true;

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
                    }

                    local_sv_info_map[sat_i] = sv;
                }
            }
        }

        last_sv_info_map = local_sv_info_map;

        const double base_gpst_sec = (gpsTimeValid(local_gpst_sec))? local_gpst_sec: (local_sys_time.time + local_sys_time.sec);
        const double time_frame = base_gpst_sec * 10.0;
        // ROS_WARN("------------------end-------------------------Processing time %f", time_frame);

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
            static constexpr int kMinGnssSatsForTdcp = 4;
            const bool gnss_obs_ok = ((int)local_fgo_input_raw.size() >= kMinGnssSatsForTdcp) && ((int)local_ephem_array.size() >= kMinGnssSatsForTdcp);

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


    ~SelfAdjustedTR()
    {
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
