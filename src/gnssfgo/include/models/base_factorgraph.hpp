/*******************************************************
 *
 *******************************************************/
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/GnssEphemMsg.h>
#include <gnss_comm/GnssGloEphemMsg.h>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_constant.hpp>
#include <ceres/ceres.h>

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <ros/ros.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <new>

#include "../factors/pseudorange_factor.hpp"
#include "../factors/doppler_factor.hpp"
#include "../factors/carrier_phase_factor.hpp"
#include "../factors/prior_factor.hpp"
#include "../tools/gnss_comm_extra.h"
#include "../tools/checking.hpp"
#include "../tools/differential_functions.hpp"
#include "../datatype.h"

#define state_size 7 // x,y,z, clk_gps, clk_glo, clk_gal, clk_bds
#define ar_size 100

#define enable_static_MM 0 // enable zero velocity motion model?
#define SHOW_COVMAX 0 

class FactorGraph{
public:
    /* continuous data stream */
    std::map<double, std::vector<gnss_comm::ObsPtr>> gnss_raw_map;
    std::map<double,std::vector<gnss_comm::EphemBasePtr>> ephems_map;
    std::map<double, nav_msgs::Odometry> doppler_map;
    std::map<double, std::vector<gnss_comm::ObsPtr>> station_gnss_raw_map;
    std::map<int,sv_info> sv_info_map;
    std::map<double, std::map<int, sv_info>> sv_info_window_map;
    double station_x, station_y, station_z;
    /* Ceres solver object */
    ceres::Problem problem;
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    ceres::LossFunction *loss_function;
    // timestamp, sat_no, available or not
    std::map<double,std::vector<int>> available_satlist; // available satellites for each epoch

    /* size of factor graph */
    int windowSize = 10;
    
    /* position state array of factor graph */
    std::vector<double*> state_array;

    /* ambiguity state array of factor graph */
    std::vector<double*> ar_state_array;

    /* array save the num of ambiguity unknowns for each epoch */
    std::vector<int*> ar_state_num;

    /* gps second array */
    // 这个数组的作用是记录每个时刻对应的GPS秒数，方便在添加观测因子时进行时间同步和匹配
    std::vector<int> gps_sec_array;
    
    /* position state saved in vector pair */
    std::vector<std::pair<double, Eigen::Vector3d>> Ps;
    std::vector<std::pair<double, Eigen::Vector4d>> Clocks;

    /* ambiguity state saved in vector pair <time, <PRN, ambiguity>> */
    std::vector<std::pair<double, std::vector<std::pair<int, double>>>> AR;

    /* ambiguity index -> satellite id mapping for each epoch */
    std::vector<std::vector<int>> epoch_ar_sat_ids;

    /* size of the last factor graph optimization */
    // 在初始窗口未填满阶段，用于初始化窗口内新增观测的状态
    int lastFactorGraphSize = 0;

    /* fixed variance of doppler measurements */
    double var = 0.6;

    /* reference point for ENU calculation */
    Eigen::MatrixXd enu_ref_llh, enu_ref_ecef;

    /* measurements size */
    int measSize = 0;

    /* factor counts for per-iteration summary */
    int last_added_psr_factor_count = 0;
    int last_added_doppler_factor_count = 0;

    /* latest GNSS-RTK solution with LAMBDA */
    Eigen::Matrix<double, 3,1> fixedStateGNSSRTK;
    Eigen::MatrixXd covMatrix;

    int fixed_cnt = 0;
    
    bool first_run = true;
    bool has_new_data = false;
    double time_frame_last = -1.0;
    double time_frame_now = -1.0;

    // Sliding-window marginalization prior for the first state in the next window.
    bool MARGINAL_ENABLE = true;
    bool has_marginalization_prior_ = false;
    double marginalization_prior_time_ = -1.0;
    std::array<double, state_size> marginalization_prior_mean_{};
    Eigen::Matrix<double, state_size, state_size> marginalization_prior_sqrt_info_ = Eigen::Matrix<double, state_size, state_size>::Identity();

protected:
    ceres::Problem::Options problem_options;

    void freeStateMemory()
    {
        for (double *&ptr : state_array)
        {
            delete[] ptr;
            ptr = nullptr;
        }
        state_array.clear();
    }

    void freeARStateMemory()
    {
        for (double *&ptr : ar_state_array)
        {
            delete[] ptr;
            ptr = nullptr;
        }
        ar_state_array.clear();
    }

    void freeARStateNumMemory()
    {
        for (int *&ptr : ar_state_num)
        {
            delete[] ptr;
            ptr = nullptr;
        }
        ar_state_num.clear();
    }


public:
    FactorGraph()
    {
        loss_function = nullptr;
        problem_options.cost_function_ownership = ceres::TAKE_OWNERSHIP;
        problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
        problem_options.local_parameterization_ownership = ceres::TAKE_OWNERSHIP;
        problem.~Problem();
        new (&problem) ceres::Problem(problem_options);
    }

    ~FactorGraph()
    {
        freeStateMemory();
        freeARStateMemory();
        freeARStateNumMemory();
        delete loss_function;
        loss_function = nullptr;
    }

     /* input gnss raw (pseudorange/carrier-phase) data  */
    bool input_gnss_raw_data(std::vector<gnss_comm::ObsPtr> GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            gnss_raw_map[timestamp] = GNSS_data;
            time_frame_now = timestamp;
            return true;
        }
    }

    /* input Doppler data  */
    bool input_doppler_data(nav_msgs::Odometry dopplerData, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            doppler_map[timestamp] = dopplerData;
            return true;
        }
    }

    /* input GNSS data from station  */
    bool input_station_data(std::vector<gnss_comm::ObsPtr> GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            station_gnss_raw_map[timestamp] = GNSS_data;
            return true;
        }
    }

    bool input_ephem_data(std::vector<gnss_comm::EphemBasePtr> ephem_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            ephems_map[timestamp] = ephem_data;
            return true;
        }
    }


    bool input_sv_info_window_map(std::map<int, sv_info> sv_info_map, double timestamp)
    {
        this->sv_info_window_map[timestamp] = sv_info_map;
        auto current_time = timestamp;
        this->sv_info_map = sv_info_window_map[current_time];
        return true;
    }

    /* set up ceres-solver options */
    bool setupSolverOptions()
    {
        options.use_nonmonotonic_steps = true;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.trust_region_strategy_type = ceres::TrustRegionStrategyType::DOGLEG;
        options.dogleg_type = ceres::DoglegType::SUBSPACE_DOGLEG;

        // Keep solver deterministic and avoid any thread-safety surprises.
        options.num_threads = 1;
        options.num_linear_solver_threads = 1;
        options.max_num_iterations = 50;
        // options.max_solver_time_in_seconds = 0.04;      // 绝对禁止超过40ms
        return true;
    }

    /* set up Loss functions options */
    bool setupLossFunction(std::string loss)
    {
        if (loss_function)
        {
            delete loss_function;
            loss_function = nullptr;
        }
        if(loss=="Huber")
            loss_function = new ceres::HuberLoss(1.0);
        else 
        {
            loss_function = new ceres::CauchyLoss(1.0);
        }
        return true;
    }

    /* get data stream size */  
    int getDataStreamSize()
    {
        measSize = gnss_raw_map.size();
        return measSize;
    }

    bool resetProblem() 
    {
        problem.~Problem();
        new (&problem) ceres::Problem(problem_options);
        summary = ceres::Solver::Summary();
        covMatrix.resize(0, 0);
        return true;
    }

    // 优化窗口内每个时刻对应一个状态向量，每个状态向量包含位置和钟差
    bool addAllParameterBlocksToGraph()
    {
        int length = measSize;
        for (int i = 0; i < length; ++i)
        {
            problem.AddParameterBlock(state_array[i], state_size);
        }
        return true;
    }

    // 根据当前窗口内的观测数据类型，配置状态参数化（如果某个系统没有观测，则对应钟差维度固定）
    bool configureStateParameterizationForWindow()
    {
        int length = measSize;
        auto iter_pr = gnss_raw_map.begin();
        for (int i = 0; i < length && iter_pr != gnss_raw_map.end(); ++i, ++iter_pr)
        {
            std::vector<int> fixed_dims;
            bool has_gps = false;
            bool has_glo = false;
            bool has_gal = false;
            bool has_bds = false;
            for (const auto &obs : iter_pr->second)
            {
                if (!obs)
                {
                    continue;
                }
                const int sys = satsys(obs->sat, NULL);
                if (sys == SYS_GPS)       has_gps = true;
                else if (sys == SYS_GLO)  has_glo = true;
                else if (sys == SYS_GAL)  has_gal = true;
                else if (sys == SYS_BDS)  has_bds = true;
            }

            // 如果当前时刻没有某系统观测，则固定对应钟差维度
            // state: [x, y, z, clk_gps, clk_glo, clk_gal, clk_bds]
            if (!has_gps) fixed_dims.push_back(3);
            if (!has_glo) fixed_dims.push_back(4);
            if (!has_gal) fixed_dims.push_back(5);
            if (!has_bds) fixed_dims.push_back(6);

            if (!fixed_dims.empty())
            {
                problem.SetParameterization(state_array[i], new ceres::SubsetParameterization(state_size, fixed_dims));
            }
        }
        return true;
    }

    // 固定窗口内第一个状态（绝对位置）以消除全局平移模糊度
    bool anchorFirstStateInWindow(bool flag)
    {
        if (!flag || measSize <= 0)
        {
            return true;
        }

        // 边缘化历史状态为先验约束，不固定第一个状态
        if (MARGINAL_ENABLE && addMarginalizationPriorFactorForFirstState())
        {
            return true;
        }
        // 未开启边缘化，直接固定第一个状态
        problem.SetParameterBlockConstant(state_array[0]);
        return true;
    }

    // 将将被滑窗移除的历史历元影响，近似压缩为“下一窗口首历元”的高斯先验因子。 todo：可以考虑更精细的边缘化策略（例如保留更多状态或使用更复杂的先验分布）以减少信息损失。
    // 把“下一窗口首历元”状态保存为先验均值，并设置弱信息矩阵（位置 3 m 钟差 30 m 标准差）
    bool marginalizeHistoryEpochToPriorFactor()
    {
        if (gnss_raw_map.empty())
        {
            return false;
        }

        const int numElementsToRemove = static_cast<int>(gnss_raw_map.size()) - windowSize;
        if (numElementsToRemove <= 0)
        {
            return true;
        }

        auto next_first_it = gnss_raw_map.begin();
        std::advance(next_first_it, numElementsToRemove);
        if (next_first_it == gnss_raw_map.end())
        {
            return false;
        }

        const int next_first_index = numElementsToRemove;
        std::array<double, state_size> mean{};
        bool got_mean = false;

        // Preferred: use current optimized state buffer.
        if (next_first_index >= 0 && next_first_index < static_cast<int>(state_array.size()) && state_array[next_first_index])
        {
            for (int i = 0; i < state_size; ++i)
            {
                mean[i] = state_array[next_first_index][i];
            }
            got_mean = true;
        }

        // Fallback: read from the persisted state map.
        if (!got_mean)
        {
            const double t_keep = next_first_it->first;
            for (const auto &kv : Ps)
            {
                if (std::abs(kv.first - t_keep) < 1e-6)
                {
                    mean[0] = kv.second[0];
                    mean[1] = kv.second[1];
                    mean[2] = kv.second[2];
                    mean[3] = 0.0;
                    mean[4] = 0.0;
                    got_mean = true;
                    break;
                }
            }
        }

        if (!got_mean)
        {
            return false;
        }

        Eigen::Matrix<double, state_size, state_size> sqrt_info = Eigen::Matrix<double, state_size, state_size>::Zero();

        // A practical weak prior:
        // - position: 3.0 m std
        // - clocks:   30.0 m std (kept loose to avoid over-constraining)
        const std::array<double, state_size> sigma{{3.0, 3.0, 3.0, 30.0, 30.0, 30.0, 30.0}};
        for (int i = 0; i < state_size; ++i)
        {
            const double s = std::max(1e-6, sigma[i]);
            sqrt_info(i, i) = 1.0 / s;
        }

        has_marginalization_prior_ = true;
        marginalization_prior_time_ = next_first_it->first;
        marginalization_prior_mean_ = mean;
        marginalization_prior_sqrt_info_ = sqrt_info;
        return true;
    }

    bool fixHistoryStateInWindow(bool flag)
    {
        if (!flag)
        {
            return true;
        }

        // Optimize all states except the first one in the window.
        return anchorFirstStateInWindow(measSize > 1);
    }

    /* setup state size */
    // 设置状态参数内存
    bool setupStateMemory()
    {
        state_array.resize(windowSize);
        ar_state_num.resize(windowSize);

        for(int i = 0; i < windowSize;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            if (state_array[i])
            {
                delete[] state_array[i];
                state_array[i] = nullptr;
            }
            if (ar_state_num[i])
            {
                delete[] ar_state_num[i];
                ar_state_num[i] = nullptr;
            }

            state_array[i]  = new double[state_size];
            ar_state_num[i] = new int[1];
            ar_state_num[i][0] = 0;
        }

        return true;
    }

    /* initialize the previous optimzied states */
    // 加载上一次优化的状态结果作为当前窗口内状态的初值
    bool initializeOldGraph()
    {
        std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator iter_pr;
        iter_pr = gnss_raw_map.begin();
        int length = measSize;

        /* 加载历史数据，meassize应当与窗口大小保持一致 */
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            double time = iter_pr->first;

            for(int i = 0; i < Ps.size(); i++)
            {
                if(time == Ps[i].first)
                {
                    state_array[m][0] = Ps[i].second[0];
                    state_array[m][1] = Ps[i].second[1];
                    state_array[m][2] = Ps[i].second[2];
                    state_array[m][3] = Clocks[i].second[0];
                    state_array[m][4] = Clocks[i].second[1];
                    state_array[m][5] = Clocks[i].second[2];
                    state_array[m][6] = Clocks[i].second[3];
                }
            } 
            
        }
        return true;
    }

    // 对于当前窗口内新增的状态，设置其初值为伪距定位结果（初始窗口未填满阶段）或上一个时刻的优化结果（窗口填满后）
    bool initializeGraph(Eigen::MatrixXd psr_ECEF)
    {
        int length = measSize;
        std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator iter;
        iter = gnss_raw_map.begin();
        bool add_tr_factor = true;
        for(int i = 0; i <length; i++,iter++)
        {
            //初始窗口未填满阶段
            if(i >= (lastFactorGraphSize-1))  
            {
                std::vector<gnss_comm::ObsPtr> gnss_data = (iter->second);
                double index = iter->first;
                // Eigen::MatrixXd psr_ECEF = gnss_comm::psr_pos(gnss_data,ephems_map[index],{8,0.0});
                if( i <= 0)
                {
                    state_array[i][0] = psr_ECEF(0);
                    state_array[i][1] = psr_ECEF(1);
                    state_array[i][2] = psr_ECEF(2);
                    state_array[i][3] = 0.0;
                    state_array[i][4] = 0.0;
                    state_array[i][5] = 0.0;
                    state_array[i][6] = 0.0;
                }
                else
                {
                    // 先用上一历元状态进行预测（多普勒积分），避免每个历元都被伪距高度噪声“拉回”。
                    Eigen::Vector3d pred_state(state_array[i-1][0], state_array[i-1][1], state_array[i-1][2]);
                    double dt_sec = 0.0;
                    if (iter != gnss_raw_map.begin())
                    {
                        auto iter_prev = std::prev(iter);
                        dt_sec = (iter->first - iter_prev->first) / 10.0;
                    }
                    if (!std::isfinite(dt_sec) || dt_sec <= 0.0)
                    {
                        dt_sec = 0.1; // fallback to 10Hz
                    }

                    const auto doppler_it = doppler_map.find(iter->first);
                    if (doppler_it != doppler_map.end())
                    {
                        const auto &v = doppler_it->second.twist.twist.linear;
                        pred_state.x() += v.x * dt_sec;
                        pred_state.y() += v.y * dt_sec;
                        pred_state.z() += v.z * dt_sec;
                    }

                    state_array[i][0] = pred_state.x();
                    state_array[i][1] = pred_state.y();
                    state_array[i][2] = pred_state.z();
                    state_array[i][3] = state_array[i-1][3];
                    state_array[i][4] = state_array[i-1][4];
                    state_array[i][5] = state_array[i-1][5];
                    state_array[i][6] = state_array[i-1][6];

                    // 伪距结果可用时，仅在偏差过大且伪距质量良好时重置，避免高度周期性抖动。
                    int valid_psr_cnt = 0;
                    double psr_std_sum = 0.0;
                    for (const auto& obs : gnss_data)
                    {
                        if (!obs)
                        {
                            continue;
                        }
                        int l1_idx = -1;
                        L1_freq(obs, &l1_idx);
                        if (l1_idx < 0 || l1_idx >= static_cast<int>(obs->psr_std.size()))
                        {
                            continue;
                        }
                        const double std_m = obs->psr_std[l1_idx];
                        if (std::isfinite(std_m) && std_m > 0.0)
                        {
                            psr_std_sum += std_m;
                            ++valid_psr_cnt;
                        }
                    }
                    const double avg_psr_std = (valid_psr_cnt > 0) ? (psr_std_sum / valid_psr_cnt) : 1e9;
                    const bool psr_quality_good = (valid_psr_cnt >= 6) && (avg_psr_std <= 8.0);

                    if (std::isfinite(psr_ECEF.norm()) && psr_ECEF.norm() > 1.0)
                    {
                        const double delta_norm = (pred_state - psr_ECEF).norm();
                        if (delta_norm > 100.0 && psr_quality_good)
                        {
                            ROS_INFO("RE_INITIAL BY PSR_ECEF");
                            state_array[i][0] = psr_ECEF(0);
                            state_array[i][1] = psr_ECEF(1);
                            state_array[i][2] = psr_ECEF(2);
                            state_array[i][3] = 0.0;
                            state_array[i][4] = 0.0;
                            state_array[i][5] = 0.0;
                            state_array[i][6] = 0.0;
                        }
                    }
                }
                
            }
        }
        return 1;
    }

    bool initializeGraphByPSR(Eigen::MatrixXd psr_ECEF)
    {
        int length = measSize;
        std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator iter;
        iter = gnss_raw_map.begin();
        bool add_tr_factor = true;
        for(int i = 0; i <length; i++,iter++)
        {
            //初始窗口未填满阶段
            if(i >= (lastFactorGraphSize-1))  
            {
                std::vector<gnss_comm::ObsPtr> gnss_data = (iter->second);
                double index = iter->first;
                // Eigen::MatrixXd psr_ECEF = gnss_comm::psr_pos(gnss_data,ephems_map[index],{8,0.0});
                if( i <= 0)
                {
                    state_array[i][0] = psr_ECEF(0);
                    state_array[i][1] = psr_ECEF(1);
                    state_array[i][2] = psr_ECEF(2);
                    state_array[i][3] = 0.0;
                    state_array[i][4] = 0.0;
                    state_array[i][5] = 0.0;
                    state_array[i][6] = 0.0;
                }
                else
                {
                    state_array[i][0] = psr_ECEF(0);
                    state_array[i][1] = psr_ECEF(1);
                    state_array[i][2] = psr_ECEF(2);
                    state_array[i][3] = 0.0;
                    state_array[i][4] = 0.0;
                    state_array[i][5] = 0.0;
                    state_array[i][6] = 0.0;
                }

            }
        }
        return add_tr_factor;
    }

    bool addMarginalizationPriorFactorForFirstState()
    {
        if (!has_marginalization_prior_ || measSize <= 0 || gnss_raw_map.empty())
        {
            return false;
        }

        const double first_time = gnss_raw_map.begin()->first;
        if (std::abs(first_time - marginalization_prior_time_) > 1e-6)
        {
            return false;
        }

        ceres::CostFunction *prior_function =
            new ceres::AutoDiffCostFunction<FirstStatePriorFactor, state_size, state_size>(
                new FirstStatePriorFactor(marginalization_prior_mean_, marginalization_prior_sqrt_info_));
        problem.AddResidualBlock(prior_function, nullptr, state_array[0]);
        return true;
    }

    bool clearMarginalizationPriorFactor()
    {
        has_marginalization_prior_ = false;
        marginalization_prior_time_ = -1.0;
        marginalization_prior_mean_.fill(0.0);
        marginalization_prior_sqrt_info_.setIdentity();
        return true;
    }

    /* 固定窗口内除最新状态以外的历史状态，以增强优化收敛性（类似于滑动窗口中的Marginalization） */
    bool fixHistoryState(bool flag)
    {
        /* fixed the first state only after first optimization */
        if(flag)
        {
            for(int i = 0; i < (lastFactorGraphSize-1); i++)
            {
                // printf("Fixing state at window %d, x: %.2f, y: %.2f, z: %.2f\n", i, state_array[i][0], state_array[i][1], state_array[i][2]);
                problem.SetParameterBlockConstant(state_array[i]);
                
            }   
        }    
        return true;
    }

    /* add Doppler FACTORS */
    bool addDopplerFactors()
    {
        /* process doppler measurements */
        std::map<double, nav_msgs::Odometry>::iterator iterdopp, iterdoppNext;
        
        const double kMaxCov = 300.0;
        const double kMinConfidence = 0.5;

        // 将协方差映射为 (0,1] 的置信度；数值越大代表测量越可靠、约束越强。
        auto confidence_from_cov = [&](double cov_diag) -> double {
            if (!std::isfinite(cov_diag) || cov_diag <= 0.0)
                return kMinConfidence;
            const double normalized = cov_diag / kMaxCov;
            return std::max(kMinConfidence, std::min(1.0, normalized));
        };
        int added_doppler_factor_count = 0;
        int i = 0;
        for(iterdopp = doppler_map.begin(); iterdopp != doppler_map.end()&& (i + 1) < measSize;iterdopp++, i++)
        {
            /* add doppler measurements */
            iterdoppNext = iterdopp;
            iterdoppNext ++;
            if(iterdoppNext != doppler_map.end())
            {
                double delta_t = (iterdoppNext->first - iterdopp->first) / 10.0;
                if (!std::isfinite(delta_t) || std::abs(delta_t) < 1e-6)
                {
                    continue;
                }
                double v_x_i = iterdopp->second.twist.twist.linear.x;
                double v_y_i = iterdopp->second.twist.twist.linear.y;
                double v_z_i = iterdopp->second.twist.twist.linear.z;
                #if enable_static_MM
                v_x_i = 0.0;
                v_y_i = 0.0;
                v_z_i = 0.0;
                #endif
                double var_x = confidence_from_cov(iterdopp->second.twist.covariance[0]);
                double var_y = confidence_from_cov(iterdopp->second.twist.covariance[1]);
                double var_z = confidence_from_cov(iterdopp->second.twist.covariance[2]);
                Eigen::Vector3d var_vec(var_x,var_y,var_z);
                ceres::CostFunction* doppler_function = new ceres::AutoDiffCostFunction<dopplerFactor, 3 
                                                        , state_size,state_size>(new 
                                                        dopplerFactor(v_x_i, v_y_i, v_z_i, delta_t, var_vec));
                problem.AddResidualBlock(doppler_function, loss_function, state_array[i],state_array[i+1]);
                ++added_doppler_factor_count;
            }
        }
        printf("[DOPPLER FACTOR] Added %d .\n", added_doppler_factor_count);
        last_added_doppler_factor_count = added_doppler_factor_count;
        return true;
    }

    bool addDopplerFactorsForWindow()
    {
        return addDopplerFactors();
    }

    bool addPsrFactors()
    {
        int length = measSize;
        auto iter_pr = gnss_raw_map.begin();
        auto iter_ephem = ephems_map.begin();
        int added_pseudorange_factor_count = 0;
        int if_psr_count = 0;
        for (int epoch_idx = 0; epoch_idx < length && iter_pr != gnss_raw_map.end(); ++epoch_idx, ++iter_pr)
        {
            const double epoch_time = iter_pr->first;
            const std::vector<gnss_comm::ObsPtr> &epoch_gnss_data = iter_pr->second;
            static constexpr int kMinPsrObs = 4;
            if ((int)epoch_gnss_data.size() < kMinPsrObs)
            {
                continue;
            }

            while (iter_ephem != ephems_map.end() && iter_ephem->first < epoch_time)
            {
                ++iter_ephem;
            }
            if (iter_ephem == ephems_map.end() || iter_ephem->first != epoch_time)
            {
                continue;
            }
            const std::vector<gnss_comm::EphemBasePtr> &epoch_ephems = iter_ephem->second;
            if ((int)epoch_ephems.size() < kMinPsrObs)
            {
                continue;
            }

            // 使用初始窗口内状态的当前时刻位置作为卫星位置校正的参考，如果没有有效初值则不进行状态校正（即不使用状态先验）
            Eigen::Vector3d state_guess(state_array[epoch_idx][0], state_array[epoch_idx][1], state_array[epoch_idx][2]);
            const bool has_state_guess = std::isfinite(state_guess.norm()) && state_guess.norm() > 1e3;
            std::vector<double> iono_params(8, 0.0);

            // 获取当前历元的校正伪距测量，包括卫星位置校正、钟差校正、离子层和对流层延迟校正等，作为后续添加伪距因子的输入
            const std::vector<gnss_comm_extra::CorrectedPseudorangeMeasurement> corrected_measurements =
                gnss_comm_extra::buildCorrectedPseudorangeMeasurements(epoch_gnss_data, epoch_ephems, state_guess, has_state_guess, iono_params);

            for (const auto &measurement : corrected_measurements)
            {
                if (!measurement.valid || measurement.sat_sys == "Unknown")
                {
                    continue;
                }

                ceres::CostFunction* ps_function = new ceres::AutoDiffCostFunction<pseudorangeFactor, 1
                                                                , state_size>(new pseudorangeFactor(
                                                                    measurement.sat_sys,
                                                                    measurement.sat_pos.x(),
                                                                    measurement.sat_pos.y(),
                                                                    measurement.sat_pos.z(),
                                                                    measurement.pseudorange,
                                                                    measurement.sigma,
                                                                    measurement.sv_dt_sec,
                                                                    measurement.tgd_sec,
                                                                    measurement.ion_delay_m,
                                                                    measurement.tro_delay_m
                                                                ));
                problem.AddResidualBlock(ps_function, loss_function, state_array[epoch_idx]);
                ++added_pseudorange_factor_count;
            }
        }
        printf("[  PSR   FACTOR] Added %d (L1=%d IF=%d).\n", added_pseudorange_factor_count,
                 added_pseudorange_factor_count - if_psr_count, if_psr_count);
        last_added_psr_factor_count = added_pseudorange_factor_count;
        return true;
    }

    bool getCovarianceMatrixOfLatestEpoch()
    {
        /* setup the covariance matrix */
        int length = measSize;
        if (length <= 0)
        {
            ROS_WARN("[FactorGraph::getCovarianceMatrixOfLatestEpoch] empty window");
            return false;
        }
        const int AR_cnt = (length > 0 && (length - 1) < (int)ar_state_num.size() && ar_state_num[length-1])
            ? ar_state_num[length-1][0]
            : 0;

        // Current implementation only requests the covariance of the latest position-state block.
        // When AR_cnt == 0 (e.g. pseudorange-only), avoid resizing to a huge/garbage size.
        const int covSize = state_size + std::max(0, AR_cnt);
        covMatrix.resize(covSize, covSize);

        /* initialize covariance computation options */
        ceres::Covariance::Options cov_options;
        cov_options.algorithm_type = ceres::DENSE_SVD;
        cov_options.min_reciprocal_condition_number = 1.95146e-40;
        // NOTE:
        // `null_space_rank` should not be set larger than the effective rank of this
        // local block. A fixed value 10 with `state_size == 5` can over-truncate SVD
        // and collapse the returned covariance (often to near/all-zero values).
        // Keep it at 0 here and let Ceres estimate covariance from the observable rank.
        cov_options.null_space_rank = 0;
        ceres::Covariance covariance(cov_options);

        /* declare covariance blocks */
        std::vector<std::pair<const double*, const double*> > cov_block;
        cov_block.push_back(std::make_pair(state_array[length-1], state_array[length-1]));

        /* compute covariance matrix*/
        if (!covariance.Compute(cov_block, &problem))
        {
            ROS_WARN("[FactorGraph::getCovarianceMatrixOfLatestEpoch] Covariance::Compute failed (length=%d)", length);
            return false;
        }

        if (!covariance.GetCovarianceBlock(state_array[length-1], state_array[length-1], covMatrix.data()))
        {
            ROS_WARN("[FactorGraph::getCovarianceMatrixOfLatestEpoch] GetCovarianceBlock failed (length=%d)", length);
            return false;
        }

        // std::cout << "covMatrix-> \n" <<std::setprecision(4)<< covMatrix<<"\n";

        return true;
    }

    bool solveAmbiguityResolutionViaCovCeresSolver()
    {
        /* construct the ambiguity resolution problem */
        int length = measSize;
        int AR_cnt = ar_state_num[length-1][0];

        if (AR_cnt <= 0)
        {
            // No ambiguity states to fix (e.g. pseudorange-only mode)
            return true;
        }

        /* assign covariance */
        Eigen::MatrixXd Qa, Qb; // Qa: covariance of phase bias
        Eigen::MatrixXd jacobian_mat = covMatrix;
        Qa.resize(AR_cnt, AR_cnt);
        Qb.resize(3, 3);
        for(int i = 3; i < AR_cnt+3; i++)
            for(int j = 3; j < AR_cnt+3; j++)
            {
                Qa(i-3,j-3) = jacobian_mat(i,j);
            }
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            {
                Qb(i,j) = jacobian_mat(i,j);
            }
        // std::cout<<"Qa-> \n"<<std::setprecision(4)<<Qa<<std::endl;

        /* assign covariance */
        // int covSize = pose_state_size + ar_state_num[length-1][0];
        // Qb = covMatrix.block<3,3>(0,0); // cov of posiiton state
        // Qa = covMatrix.block<ar_state_num[length-1][0],ar_state_num[length-1][0]>(pose_state_size-1,pose_state_size-1); // cov of AR state

        double *a;
        double *Q;
        double *F;
        double S[2];
        double n  = AR_cnt;
        double m = 2;
        F = new double[static_cast<size_t>(n) * 2];
        // Q = new double[static_cast<size_t>(n) * n];
        size_t n_size = static_cast<size_t>(n);
        size_t total = n_size * n_size;  // 两个 size_t 相乘
        Q = new double[total];
        a= new double[static_cast<size_t>(n) * 1];

        for(int i = 0; i < n; i++)
        {
            // a[i] = ar_state_array[length-1][i];
            a[i] = ar_state_array[length-1][i];
        }
        std::vector<double> cov_vector;
        // fi_cov_ar = fi_cov_ar.inverse();
        for(int i = 0; i < n; i++) // cols
            for(int j = 0; j < n; j++) // rows
            {
                cov_vector.push_back(Qa(j,i));
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
            
        for(int i=0; i<cov_vector.size();i++)
        {
            Q[i] = cov_vector[i];
        }

        lambda(n, m, a, Q, F, S);
        // F中存储了整数解，S中存储了残差平方和
        Eigen::MatrixXd Q_ba, Q_ab;
        Q_ba.resize(3,AR_cnt);
        Q_ab.resize(AR_cnt,3);
        for(int i = 0; i < 3; i++) // rows
            for(int j = 3; j < (3+ AR_cnt); j++) // cols
            {
                Q_ba(i,j-3) = jacobian_mat(i,j);
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
        for(int i = 3; i < (3 + AR_cnt); i++) // rows
            for(int j = 0; j < 3; j++) // cols
            {
                Q_ab(i-3,j) = jacobian_mat(i,j);
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
        Eigen::MatrixXd integer_b, float_b;
        integer_b.resize(AR_cnt, 1);
        float_b.resize(AR_cnt, 1);
        for(int i=0; i<AR_cnt; i++)
        {
            // float_b(i) = ar_state_array[length-1][i];
            float_b(i) = ar_state_array[length-1][i];
            integer_b(i) = F[i];
        }
        // std::cout<<"Q_ba-> \n"<<std::setprecision(4)<<Q_ba<<std::endl;
        Eigen::Matrix<double, 3,1> a_state;
        a_state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]; 
        a_state = a_state + Q_ba * Qa.inverse() * (float_b - integer_b);
        // a_state = a_state - Q_ab * Qb.inverse() * (float_b - integer_b);

        Eigen::Matrix<double ,3,1> ENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2];
        double fix_flag = 0;
        //s[0]：固定解的残差平方和。固定解是指通过整数最小二乘估计（LAMBDA算法）得到的整数模糊度解。
        //s[1]：浮点解的残差平方和。浮点解是指通过卡尔曼滤波得到的双差模糊度的浮点解
        //s[0]/s[1] > 1.5时，认为固定解的可靠性较高，可以接受为固定状态
        if((S[1]/S[0])> 1.5)
        {
            state = a_state; // 
            fix_flag = 1;
            ROS_INFO("<<<<<Fixed!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>");
            fixed_cnt++;
        }
        fixedStateGNSSRTK = state;

        // Save latest epoch ambiguities with satellite IDs when mapping is available.
        if (length > 0 && (length - 1) < static_cast<int>(epoch_ar_sat_ids.size()))
        {
            const std::vector<int> &sat_ids = epoch_ar_sat_ids[length - 1];
            std::vector<std::pair<int, double>> ar_snapshot;
            const int save_cnt = std::min(AR_cnt, static_cast<int>(sat_ids.size()));
            ar_snapshot.reserve(save_cnt);
            for (int i = 0; i < save_cnt; ++i)
            {
                ar_snapshot.emplace_back(sat_ids[i], integer_b(i));
                printf("sat_id: %d, integer_ambiguity: %.2f\n",sat_ids[i], integer_b(i));
            }

            if (!ar_snapshot.empty() && !gnss_raw_map.empty())
            {
                const double epoch_time = std::prev(gnss_raw_map.end())->first;
                bool updated = false;
                for (auto &entry : AR)
                {
                    if (std::abs(entry.first - epoch_time) < 1e-6)
                    {
                        entry.second = ar_snapshot;
                        updated = true;
                        break;
                    }
                }
                if (!updated)
                {
                    AR.emplace_back(epoch_time, ar_snapshot);
                }
            }
        }

        return true;
    }

    /*
     * @brief LAMBDA integer least squares with decorrelation and search
     * @param n number of ambiguities
     * @param m number of candidates (expects m>=2 in current usage)
     * @param a float ambiguities (n x 1)
     * @param Q covariance (n x n), column-major
     * @param F integer candidates (n x m), column-major
     * @param S residuals (m x 1), S[0]=best, S[1]=second-best
     */
    bool lambda(int n, int m, const double* a, const double* Q, double* F, double* S)
    {
        if (n <= 0 || m <= 0 || !a || !Q || !F || !S)
        {
            return false;
        }

        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> Qm(Q, n, n);
        Eigen::Map<const Eigen::VectorXd> a_vec(a, n);

        Eigen::LLT<Eigen::MatrixXd> llt(Qm);
        if (llt.info() != Eigen::Success)
        {
            return false;
        }

        Eigen::MatrixXd Lc = llt.matrixL();
        Eigen::VectorXd d = Lc.diagonal();
        if ((d.array() <= 0.0).any())
        {
            return false;
        }

        Eigen::VectorXd D = d.array().square();
        Eigen::MatrixXd L = Lc * d.cwiseInverse().asDiagonal();
        Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(n, n);

        int k = n - 2;
        while (k >= 0)
        {
            for (int j = k + 1; j < n; ++j)
            {
                double mu = std::round(L(j, k));
                if (mu != 0.0)
                {
                    L.row(j).head(k + 1) -= mu * L.row(k).head(k + 1);
                    Z.col(j) -= mu * Z.col(k);
                }
            }

            double delta = D(k) + L(k + 1, k) * L(k + 1, k) * D(k + 1);
            if (delta < D(k + 1))
            {
                double lam = D(k) / delta;
                double eta = D(k + 1) * L(k + 1, k) / delta;
                D(k) = lam * D(k + 1);
                D(k + 1) = delta;

                for (int i = 0; i < k; ++i)
                {
                    std::swap(L(k, i), L(k + 1, i));
                }
                L(k + 1, k) = eta;
                for (int i = k + 2; i < n; ++i)
                {
                    std::swap(L(i, k), L(i, k + 1));
                }
                Z.col(k).swap(Z.col(k + 1));
                k = std::min(k + 1, n - 2);
            }
            else
            {
                --k;
            }
        }

        Eigen::VectorXd z = Z.transpose() * a_vec;

        Eigen::MatrixXd Fz = Eigen::MatrixXd::Zero(n, m);
        Eigen::VectorXd Svec = Eigen::VectorXd::Constant(m, std::numeric_limits<double>::infinity());

        Eigen::VectorXd zcond = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd zhat = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd dist = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd dist_sum = Eigen::VectorXd::Zero(n);

        auto insert_candidate = [&](const Eigen::VectorXd& cand, double s) {
            if (s >= Svec(m - 1))
            {
                return;
            }
            int idx = m - 1;
            while (idx > 0 && s < Svec(idx - 1))
            {
                Svec(idx) = Svec(idx - 1);
                Fz.col(idx) = Fz.col(idx - 1);
                --idx;
            }
            Svec(idx) = s;
            Fz.col(idx) = cand;
        };

        int kk = n - 1;
        zcond(kk) = z(kk);
        zhat(kk) = std::round(zcond(kk));
        step(kk) = (zcond(kk) - zhat(kk) >= 0.0) ? 1.0 : -1.0;
        dist(kk) = (zcond(kk) - zhat(kk)) * (zcond(kk) - zhat(kk)) / D(kk);
        dist_sum(kk) = dist(kk);

        double maxdist = std::numeric_limits<double>::infinity();

        while (true)
        {
            double newdist = (kk == n - 1) ? dist(kk) : dist(kk) + dist_sum(kk + 1);
            if (newdist < maxdist)
            {
                if (kk == 0)
                {
                    insert_candidate(zhat, newdist);
                    maxdist = Svec(m - 1);

                    zhat(0) += step(0);
                    double y = zcond(0) - zhat(0);
                    step(0) = (step(0) > 0.0) ? -step(0) - 1.0 : -step(0) + 1.0;
                    dist(0) = y * y / D(0);
                }
                else
                {
                    --kk;
                    zcond(kk) = z(kk);
                    for (int i = kk + 1; i < n; ++i)
                    {
                        zcond(kk) += L(i, kk) * (zcond(i) - zhat(i));
                    }
                    zhat(kk) = std::round(zcond(kk));
                    double y = zcond(kk) - zhat(kk);
                    step(kk) = (y >= 0.0) ? 1.0 : -1.0;
                    dist(kk) = y * y / D(kk);
                }
            }
            else
            {
                if (kk == n - 1)
                {
                    break;
                }
                ++kk;
                zhat(kk) += step(kk);
                double y = zcond(kk) - zhat(kk);
                step(kk) = (step(kk) > 0.0) ? -step(kk) - 1.0 : -step(kk) + 1.0;
                dist(kk) = y * y / D(kk);
            }

            if (kk < n - 1)
            {
                dist_sum(kk) = dist(kk) + dist_sum(kk + 1);
            }
        }

        Eigen::MatrixXd ZinvT = Z.transpose().inverse();
        Eigen::MatrixXd F_a = ZinvT * Fz;

        for (int j = 0; j < m; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                F[i + n * j] = F_a(i, j);
            }
            S[j] = Svec(j);
        }

        return true;
    }

    /* solve the ambiguity resolution */
    // bool solveAmbiguityResolutionFixedSolution()
    // {
    //     /* construct the ambiguity resolution problem*/
    //     int length = measSize;
    //     /* get the Jacobian matrix for covariance propogation */
    //     Eigen::MatrixXd jacobian_matrix; 
    //     Eigen::MatrixXd weighting_matrix;

    //     /* number of ambiguity variables */
    //     int AR_cnt = 0;
    //     for(int i = 0; i < ar_size; i++)
    //     {
    //         if(fabs(ar_state_array[length-1][i])>0)
    //         AR_cnt++;
    //     }
    //     ROS_INFO("AR_cnt-> %d", AR_cnt);
    //     ROS_INFO("ar_state_num[k][0]-> %d", ar_state_num[length-1][0]);
    //     AR_cnt = ar_state_num[length-1][0];
    //     jacobian_matrix.resize(40, 3 + AR_cnt);
    //     jacobian_matrix.setIdentity();
    //     weighting_matrix.resize(40, 40);
    //     weighting_matrix.setIdentity();

    //     std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator iter_cov; // station gnss measurements map iterator
    //     iter_cov = station_gnss_raw_map.end();
    //     iter_cov--;

    //     std::vector<gnss_comm::ObsPtr> st_gnss_data = (iter_cov->second);
    //     int sv_cnt = st_gnss_data.size();
    //     double t = iter_cov->first;

    //     /* find user end gnss data with closest time */
    //     std::vector<gnss_comm::ObsPtr> closest_gnss_data;
    //     findClosestEpoch(t, gnss_raw_map, closest_gnss_data);

    //     /* get the dd measurements between
    //     * 1. st_gnss_data from station
    //     * 2. closest_gnss_data from user end
    //     */

    //     /* tranverse station gnss map */
    //     int carrier_phase_index = 0;
    //     int jac_row = 0;
    //     Eigen::Vector3d u_pose(state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]);
    //     for(int q = 0; q < sv_cnt; q++)
    //     {
    //         double sat_id = st_gnss_data[q]->sat;
    //         gnss_comm::ObsPtr u_master_sv, u_iSV, r_master_sv;
            
    //         /* find the master satellite from the user end */
    //         if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV,sv_info_map))
    //         {
    //             /* find the satellite from station gnss iwth the same id with master satellite */
    //             findSatellitewithSameId(u_master_sv->sat, st_gnss_data, r_master_sv);

    //             DDMeasurement DD_measurement;
    //             DD_measurement.u_master_SV = u_master_sv;
    //             DD_measurement.u_iSV = u_iSV;

    //             DD_measurement.r_master_SV = r_master_sv;
    //             DD_measurement.r_iSV = st_gnss_data[q];

    //             Eigen::Vector3d base_pose(station_x, station_y, station_z);

    //             /* get the row for jocabian of pr DD*/
                

    //             gnss_comm_extra::getPrDDJacobian(sv_info_map, u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, weighting_matrix);
                
    //             jac_row++;
    //             if(checkCarrierPhaseConsistency(DD_measurement) && (carrier_phase_index<AR_cnt)) 
    //             // (carrier_phase_index<10))
    //             {
    //                     /* get the row for jocabian of cp DD*/
    //                     gnss_comm_extra::getCpDDJacobian(sv_info_map, u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, carrier_phase_index, weighting_matrix);

    //                 carrier_phase_index++;
    //                 // LOG(INFO)<<"add double-difference carrier phase factor";
    //                 jac_row++;
    //             }
    //             else
    //             {
    //                 // LOG(INFO)<<"no carrier-phase measurement";
    //             }
    //         }
    //     }

    //     Eigen::MatrixXd jacobian_matrix_;
    //     Eigen::MatrixXd weighting_matrix_;
    //     weighting_matrix_.resize(jac_row,jac_row);
    //     jacobian_matrix_.resize(jac_row, jacobian_matrix.cols());
    //     for(int i = 0; i < jac_row; i++)
    //         for(int j = 0; j < jacobian_matrix.cols(); j++)
    //         {
    //             jacobian_matrix_(i,j) = jacobian_matrix(i,j);
    //         }

    //     for(int i = 0; i <jac_row; i++)
    //         for(int j = 0; j <jac_row; j++)
    //         {
    //             weighting_matrix_(i,j) = weighting_matrix(i,j);
    //         }

    //     Eigen::MatrixXd jacobian_mat;
    //     jacobian_mat = (jacobian_matrix_.transpose() * weighting_matrix_ *jacobian_matrix_).inverse();
    //     // std::cout<<"jacobian_mat-> \n"<<std::setprecision(4)<<jacobian_mat<<std::endl;
    //     // std::cout<<"weighting_matrix_-> \n"<<std::setprecision(4)<<weighting_matrix_<<std::endl;
    //     Eigen::MatrixXd Qa, Qb; // Qa: covariance of phase bias
    //     Qa.resize(AR_cnt, AR_cnt);
    //     Qb.resize(3, 3);
    //     for(int i = 3; i < AR_cnt+3; i++)
    //         for(int j = 3; j < AR_cnt+3; j++)
    //         {
    //             Qa(i-3,j-3) = jacobian_mat(i,j);
    //         }
    //     for(int i = 0; i < 3; i++)
    //         for(int j = 0; j < 3; j++)
    //         {
    //             Qb(i,j) = jacobian_mat(i,j);
    //         }
    //     // std::cout<<"Qa-> \n"<<std::setprecision(4)<<Qa<<std::endl;

    //     double *a;
    //     double *Q;
    //     double *F;
    //     double S[2];
    //     double n  = AR_cnt;
    //     double m = 2;
    //     F = new double[static_cast<size_t>(n) * 2];
    //     size_t n_size = static_cast<size_t>(n);
    //     size_t total = n_size * n_size;  // 两个 size_t 相乘
    //     Q = new double[total];
    //     a= new double[static_cast<size_t>(n) * 1];

    //     for(int i = 0; i < n; i++)
    //     {
    //         a[i] = ar_state_array[length-1][i];
    //     }
    //     std::vector<double> cov_vector;
    //     // fi_cov_ar = fi_cov_ar.inverse();
    //     for(int i = 0; i < n; i++) // cols
    //         for(int j = 0; j < n; j++) // rows
    //         {
    //             cov_vector.push_back(Qa(j,i));
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
            
    //     for(int i=0; i<cov_vector.size();i++)
    //     {
    //         Q[i] = cov_vector[i];
    //     }

    //     lambda(n, m, a, Q, F, S);
    //     Eigen::MatrixXd Q_ba, Q_ab;
    //     Q_ba.resize(3,AR_cnt);
    //     Q_ab.resize(AR_cnt,3);
    //     for(int i = 0; i < 3; i++) // rows
    //         for(int j = 3; j < (3+ AR_cnt); j++) // cols
    //         {
    //             Q_ba(i,j-3) = jacobian_mat(i,j);
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
    //     for(int i = 3; i < (3 + AR_cnt); i++) // rows
    //         for(int j = 0; j < 3; j++) // cols
    //         {
    //             Q_ab(i-3,j) = jacobian_mat(i,j);
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
    //     Eigen::MatrixXd integer_b, float_b;
    //     integer_b.resize(AR_cnt, 1);
    //     float_b.resize(AR_cnt, 1);
    //     for(int i=0; i<AR_cnt; i++)
    //     {
    //         float_b(i) = ar_state_array[length-1][i];
    //         integer_b(i) = F[i];
    //     }
    //     // std::cout<<"Q_ba-> \n"<<std::setprecision(4)<<Q_ba<<std::endl;
    //     Eigen::Matrix<double, 3,1> a_state;
    //     a_state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]; 
    //     a_state = a_state + Q_ba * Qa.inverse() * (float_b - integer_b);
    //     // a_state = a_state - Q_ab * Qb.inverse() * (float_b - integer_b);

    //     Eigen::Matrix<double ,3,1> ENU;
    //     Eigen::Matrix<double, 3,1> state;
    //     state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2];
    //     double fix_flag = 0;
    //     if((S[1]/S[0])> 3)
    //     {
    //         state = a_state; // 
    //         fix_flag = 1;
    //         ROS_INFO("<<<<<Fixed!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>");
    //         fixed_cnt++;
    //     }
    //     fixedStateGNSSRTK = state;

    //     return true;
    // }

    /* solve the ambiguity resolution */
    /* 求解整周模糊度*/
    /* solve the ambiguity resolution of current epoch  使用协方差传播（Covariance Propagation）的方式进行模糊度确定
    * 这种方法通常基于卡尔曼滤波或其他递推滤波方法，通过动态模型和观测模型逐步更新状态估计。
    * 该方法没有动态模型，而是通过数学推导将动态信息转化为先验因子 */
    // bool solveAmbiguityResolutionViaCovPropogation()
    // {
    //     /* construct the ambiguity resolution problem*/
    //     int length = measSize;

    //     /* get the Jacobian matrix for covariance propogation */
    //     /* 给定随机变量 x 及其协方差矩阵 Pₓₓ，以及函数 y = f(x)，协方差传播计算输出 y 的协方差矩阵 Pᵧᵧ。*/
    //     Eigen::MatrixXd jacobian_matrix; 
    //     Eigen::MatrixXd weighting_matrix;

    //     /* number of ambiguity variables */
    //     int AR_cnt = 0;
    //     AR_cnt = ar_state_num[length-1][0];
    //     // LOG(INFO) << "AR_cnt-> " << AR_cnt;
    //     // LOG(INFO) << "ar_state_num[k][0]-> "<< ar_state_num[length-1][0];
    //     AR_cnt = ar_state_num[length-1][0];
    //     jacobian_matrix.resize(40, 3 + AR_cnt);
    //     jacobian_matrix.setIdentity();
    //     weighting_matrix.resize(40, 40);
    //     weighting_matrix.setIdentity();

    //     std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator iter_cov; // station gnss measurements map iterator
    //     iter_cov = station_gnss_raw_map.end();
    //     iter_cov--;

    //     std::vector<gnss_comm::ObsPtr> st_gnss_data = (iter_cov->second);
    //     int sv_cnt = st_gnss_data.size();
    //     double t = iter_cov->first; //最近一个历元的时间

    //     /* find user end gnss data with closest time */
    //     /* 使用移动端最近历元的观测数据*/
    //     std::vector<gnss_comm::ObsPtr> closest_gnss_data;
    //     findClosestEpoch(t, gnss_raw_map, closest_gnss_data);

    //     /* get the dd measurements between
    //     * 1. st_gnss_data from station
    //     * 2. closest_gnss_data from user end
    //     */

    //     /* tranverse station gnss map */
    //     int carrier_phase_index = 0;
    //     int jac_row = 0;
    //     //最新的位置状态信息
    //     Eigen::Vector3d u_pose(state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]);
    //     for(int q = 0; q < sv_cnt; q++)
    //     {
    //         double sat_id = st_gnss_data[q]->sat;
    //         gnss_comm::ObsPtr u_master_sv, u_iSV, r_master_sv;
            
    //         /* find the master satellite from the user end */
    //         //构建双差观测方程（又一次重复操作???）
    //         if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV,sv_info_map))
    //         {
    //             /* find the satellite from station gnss iwth the same id with master satellite */
    //             findSatellitewithSameId(u_master_sv->sat, st_gnss_data, r_master_sv);

    //             DDMeasurement DD_measurement;
    //             DD_measurement.u_master_SV = u_master_sv;
    //             DD_measurement.u_iSV = u_iSV;

    //             DD_measurement.r_master_SV = r_master_sv;
    //             DD_measurement.r_iSV = st_gnss_data[q];

    //             Eigen::Vector3d base_pose(station_x, station_y, station_z);

    //             /* get the row for jocabian of pr DD*/
    //             gnss_comm_extra::getPrDDJacobian(sv_info_map, u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, weighting_matrix);
                
    //             jac_row++;
    //             if(checkCarrierPhaseConsistency(DD_measurement) && (carrier_phase_index<AR_cnt)) 
    //             // (carrier_phase_index<10))
    //             {
    //                     /* get the row for jocabian of cp DD*/
    //                     gnss_comm_extra::getCpDDJacobian(sv_info_map, u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, carrier_phase_index, weighting_matrix);

    //                 carrier_phase_index++;
    //                 // LOG(INFO)<<"add double-difference carrier phase factor";
    //                 jac_row++;
    //             }
    //             else
    //             {
    //                 // LOG(INFO)<<"no carrier-phase measurement";
    //             }
    //         }
    //     }

    //     Eigen::MatrixXd jacobian_matrix_;
    //     Eigen::MatrixXd weighting_matrix_;
    //     weighting_matrix_.resize(jac_row,jac_row);
    //     jacobian_matrix_.resize(jac_row, jacobian_matrix.cols());
    //     for(int i = 0; i < jac_row; i++)
    //         for(int j = 0; j < jacobian_matrix.cols(); j++)
    //         {
    //             jacobian_matrix_(i,j) = jacobian_matrix(i,j);
    //         }

    //     for(int i = 0; i <jac_row; i++)
    //         for(int j = 0; j <jac_row; j++)
    //         {
    //             weighting_matrix_(i,j) = weighting_matrix(i,j);
    //         }

    //     Eigen::MatrixXd jacobian_mat;
    //     jacobian_mat = (jacobian_matrix_.transpose() * weighting_matrix_ *jacobian_matrix_).inverse();

    //     // std::cout<<"jacobian_mat-> \n"<<std::setprecision(4)<<jacobian_mat<<std::endl;
    //     // std::cout<<"weighting_matrix_-> \n"<<std::setprecision(4)<<weighting_matrix_<<std::endl;
        
    //     Eigen::MatrixXd Qa, Qb; // Qa: covariance of phase bias
    //     Qa.resize(AR_cnt, AR_cnt);
    //     Qb.resize(3, 3);
    //     for(int i = 3; i < AR_cnt+3; i++)
    //         for(int j = 3; j < AR_cnt+3; j++)
    //         {
    //             Qa(i-3,j-3) = jacobian_mat(i,j);
    //         }
    //     for(int i = 0; i < 3; i++)
    //         for(int j = 0; j < 3; j++)
    //         {
    //             Qb(i,j) = jacobian_mat(i,j);
    //         }
    //     // std::cout<<"Qa-> \n"<<std::setprecision(4)<<Qa<<std::endl;

    //     double *a;
    //     double *Q;
    //     double *F;
    //     double S[2];
    //     double n  = AR_cnt;
    //     double m = 2;
    //     F = new double[static_cast<size_t>(n) * 2];
    //     size_t n_size = static_cast<size_t>(n);
    //     size_t total = n_size * n_size;  // 两个 size_t 相乘
    //     Q = new double[total];
    //     a= new double[static_cast<size_t>(n) * 1];

    //     for(int i = 0; i < n; i++)
    //     {
    //         // a[i] = ar_state_array[length-1][i];
    //         a[i] = ar_state_array[length-1][i];
    //     }
    //     std::vector<double> cov_vector;
    //     // fi_cov_ar = fi_cov_ar.inverse();
    //     for(int i = 0; i < n; i++) // cols
    //         for(int j = 0; j < n; j++) // rows
    //         {
    //             cov_vector.push_back(Qa(j,i));
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
            
    //     for(int i=0; i<cov_vector.size();i++)
    //     {
    //         Q[i] = cov_vector[i];
    //     }

    //     lambda(n, m, a, Q, F, S);
    //     Eigen::MatrixXd Q_ba, Q_ab;
    //     Q_ba.resize(3,AR_cnt);
    //     Q_ab.resize(AR_cnt,3);
    //     for(int i = 0; i < 3; i++) // rows
    //         for(int j = 3; j < (3+ AR_cnt); j++) // cols
    //         {
    //             Q_ba(i,j-3) = jacobian_mat(i,j);
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
    //     for(int i = 3; i < (3 + AR_cnt); i++) // rows
    //         for(int j = 0; j < 3; j++) // cols
    //         {
    //             Q_ab(i-3,j) = jacobian_mat(i,j);
    //             // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
    //             // if(i=i) 
    //         }
    //     Eigen::MatrixXd integer_b, float_b;
    //     integer_b.resize(AR_cnt, 1);
    //     float_b.resize(AR_cnt, 1);
    //     for(int i=0; i<AR_cnt; i++)
    //     {
    //         // float_b(i) = ar_state_array[length-1][i];
    //         float_b(i) = ar_state_array[length-1][i];
    //         integer_b(i) = F[i];
    //     }
    //     // std::cout<<"Q_ba-> \n"<<std::setprecision(4)<<Q_ba<<std::endl;
    //     Eigen::Matrix<double, 3,1> a_state;
    //     a_state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]; 
    //     a_state = a_state + Q_ba * Qa.inverse() * (float_b - integer_b);
    //     // a_state = a_state - Q_ab * Qb.inverse() * (float_b - integer_b);

    //     Eigen::Matrix<double ,3,1> ENU;
    //     Eigen::Matrix<double, 3,1> state;
    //     state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2];
    //     double fix_flag = 0;
    //     if((S[1]/S[0])> 1.5)
    //     {
    //         state = a_state; // 
    //         fix_flag = 1;
    //         ROS_INFO("<<<<<Fixed!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>");
    //         fixed_cnt++;
    //     }
    //     fixedStateGNSSRTK = state;

    //     return true;
    // }

    /* save graph state to vector for next solving */
    bool saveGraphStateToVector(bool allsave = false)
    {
        lastFactorGraphSize = measSize;

        int length = measSize;
        int m;

        if (!allsave) { m = length -1; } // 只保存最新的状态估计
        else { m = 0; } // 保存全部状态估计
       
        while(m < length)
        {
            auto find_time_iter = gnss_raw_map.end();
            std::advance(find_time_iter, windowSize-m);
            double time = find_time_iter->first; 
            
            /* if the state vector is empty, override */
            // printf("time-> %.2f, state-> [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n", time, state_array[m][0], state_array[m][1], state_array[m][2], state_array[m][3], state_array[m][4], state_array[m][5], state_array[m][6]);
            if(Ps.size()==0)
            {
                Ps.push_back(std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2])));
                Clocks.push_back(std::make_pair(time, Eigen::Vector4d(state_array[m][3],state_array[m][4],state_array[m][5],state_array[m][6])));
            }
            /* if the state vector is NOT empty, update */
            else
            {
                bool findTimeKey = false;
                for(int i = 0; i < Ps.size(); i++)
                {
                    if(time == Ps[i].first)
                    {
                        Ps[i] = std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2]));
                        Clocks[i] = std::make_pair(time, Eigen::Vector4d(state_array[m][3],state_array[m][4],state_array[m][5],state_array[m][6]));
                        findTimeKey = true;
                    }
                }
                /* new time frame, add to state vector*/
                if(findTimeKey==false)
                {
                    Ps.push_back(std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2])));
                    Clocks.push_back(std::make_pair(time, Eigen::Vector4d(state_array[m][3],state_array[m][4],state_array[m][5],state_array[m][6])));
                } 
            }
            m++;
        }
        return true;
    }

    /**
   * @brief maintain sliding window slidingWindowSize
   * @param gnss raw msg and doppler msg
   * @return void
   @ 
   */

    void removeStatesOutsideSlidingWindow()
    {
        // Before erasing old epochs, compress their effect into a prior on
        // the next-window first state.
        // 边缘化历史数据到先验因子
        if(MARGINAL_ENABLE)
        {
            marginalizeHistoryEpochToPriorFactor();
        }


        int numElementsToRemove = 0;

        /* sliding window gnss raw pseudorange */
        numElementsToRemove = static_cast<int>(gnss_raw_map.size()) - windowSize;
        if (numElementsToRemove > 0)
        {
            auto it = gnss_raw_map.begin();
            while (it != gnss_raw_map.end() && numElementsToRemove-- > 0)
            {
                sv_info_window_map.erase(it->first);
                it = gnss_raw_map.erase(it);
            }
        }

        /* sliding window ephemeris */
        numElementsToRemove = static_cast<int>(ephems_map.size()) - windowSize;
        if (numElementsToRemove > 0)
        {
            auto it = ephems_map.begin();
            while (it != ephems_map.end() && numElementsToRemove-- > 0)
            {
                it = ephems_map.erase(it);
            }
        }

        /* sliding window station gnss raw pseudorange */
        numElementsToRemove = static_cast<int>(station_gnss_raw_map.size()) - windowSize;
        if (numElementsToRemove > 0)
        {
            auto it = station_gnss_raw_map.begin();
            while (it != station_gnss_raw_map.end() && numElementsToRemove-- > 0)
            {
                it = station_gnss_raw_map.erase(it);
            }
        }

        /* sliding window gnss raw doppler */
        numElementsToRemove = static_cast<int>(doppler_map.size()) - windowSize;
        if (numElementsToRemove > 0)
        {
            auto it = doppler_map.begin();
            while (it != doppler_map.end() && numElementsToRemove-- > 0)
            {
                it = doppler_map.erase(it);
            }
        }

        // Keep this quiet in real-time; use ROS_DEBUG if needed.
    }

    void getOptimizationSettingInfo()
    {
        #if TIME_WEIGHT
            ROS_INFO("Time weight enabled.[only valid for TRRTK method]");
        #else
            ROS_ERROR("Time weight disabled.[only valid for TRRTK method]");
        #endif

        #if SHOW_COVMAX
            ROS_INFO("Covariance matrix showing enabled");
        #else
            ROS_ERROR("Covariance matrix showing disabled.");
        #endif

        if(MARGINAL_ENABLE)
            ROS_INFO("Marginalization enabled.");
        else 
            ROS_ERROR("Marginalization disabled.");

        ROS_INFO("Current window size: %d", windowSize);

        return;
    }

    /* get the latest state in ENU */
    Eigen::Matrix<double ,3,1> getLatestPosENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> fgo_enu;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        Eigen::Matrix<double, 3,1> v_ecef = state - enu_ref_ecef;
        fgo_enu =  gnss_comm::ecef2enu(enu_ref_llh, v_ecef);    
        return fgo_enu;
    }

        /* set up the reference point for ENU calculation */
    bool setupReferencePoint(double ref_lon, double ref_lat, double ref_alt)
    {
        /* reference point for ENU calculation */
        enu_ref_llh.resize(3,1);
        enu_ref_ecef.resize(3,1);
        enu_ref_llh<< ref_lon, ref_lat, ref_alt;
        enu_ref_ecef = gnss_comm::geo2ecef(enu_ref_llh);
        return true;
    }


    Eigen::Matrix<double ,3,1> getLatestPosECEF()
    {
        int length = measSize;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        return state;
    }

    Eigen::Matrix<double ,3,1> getLatestPosLLH()
    {
        int length = measSize;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        Eigen::Matrix<double ,3,1> llh_state = gnss_comm::ecef2geo(state);
        return llh_state;
    }

     /* print the latest state in ENU */
    bool printLatestFloatStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> ENUstate;
        ENUstate<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        Eigen::Matrix<double, 3,1> v_ecef = ENUstate - enu_ref_ecef;
        FGOENU =  gnss_comm::ecef2enu(enu_ref_llh, v_ecef);    
        std::cout << "FGOENU-> "<< FGOENU<< std::endl;  
        return true;
    }

    /* print the latest state in ENU */
    bool printLatestFixedStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> fixedFGOENU;
        Eigen::Matrix<double, 3,1> v_ecef = fixedStateGNSSRTK - enu_ref_ecef;
        fixedFGOENU = gnss_comm::ecef2enu(enu_ref_llh, v_ecef);    
        std::cout << "fixedFGOENU-> "<< fixedFGOENU<< std::endl;  
        return true;
    }

    /* get the path of FGO in ENU */
    nav_msgs::Path getPathENU(nav_msgs::Path& fgo_path)
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        fgo_path.poses.clear();
        fgo_path.header.frame_id = "map";
        for(int i = 0; i < length;i++)
        {
            state<< state_array[i][0], 
                    state_array[i][1], 
                    state_array[i][2];
            Eigen::Matrix<double, 3,1> v_ecef = state - enu_ref_ecef;
            FGOENU = gnss_comm::ecef2enu(enu_ref_llh, v_ecef);  
            geometry_msgs::PoseStamped pose_stamped;
            pose_stamped.header.stamp = ros::Time::now();
            pose_stamped.header.frame_id = "map";
            pose_stamped.pose.position.x = FGOENU(0);
            pose_stamped.pose.position.y = FGOENU(1);
            pose_stamped.pose.position.z = 10;
            fgo_path.poses.push_back(pose_stamped);
            // std::cout << "pose_stamped- FGO-> "<< std::endl<< pose_stamped;
        }
              
        return fgo_path;
    }

};


