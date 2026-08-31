/*******************************************************
 * Copyright (C) 2026
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

#include  <ros/ros.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <map>
#include <set>

#include "base_factorgraph.hpp"
#include "../tools/checking.hpp"

#define state_size 7 // x,y,z, clk_gps, clk_glo, clk_gal, clk_bds

class SelfAdjustedFactorGraph: public FactorGraph{
    
public:
    std::vector<TRDDMeasurement> tr_measurments; // 当前窗口内的TR测量列表

    std::map<double, Eigen::Matrix<double, 7, 1>> Ps_full_state; // 保存完整状态（位置+钟差）

    std::map<int,int> sat_lock_count_l1,sat_lock_count_l2;  // 卫星连续观测计数器L1,L2频段分开统计
    std::map<int,bool> sat_lock_check,sat_lock_check_l1,sat_lock_check_l2; // 卫星连续观测检查标志

    int MaxTRFactorNum = 100;   // TR因子数量上限，动态调整以控制优化时间
    int MinTRFactorNum = 20;    // TR因子数量下限，避免约束过弱
    int TRFactorCount;          // 统计优化规模

    bool SAME_RELIABLE = false;
    bool SAME_TIME_WEIGHT = false;

    std::map<int, std::map<std::string, double>> residuals;  // prns -> factor_name -> residual_value
    // TDEC 周跳检测的共模项缓存 (按历元+频段): 接收机钟跳步会使所有卫星的
    // Δcp 同步偏移, 不扣除会触发大规模"周跳"误报、锁定集体清零、TDCP 约束瞬间真空
    double tdec_common_l1_ = 0.0, tdec_common_l2_ = 0.0;
    double tdec_common_epoch_ = -1.0;
    bool tdec_common_valid_l1_ = false, tdec_common_valid_l2_ = false;

public:
    SelfAdjustedFactorGraph()
    {
        loss_function = nullptr;
        problem_options.cost_function_ownership = ceres::TAKE_OWNERSHIP;
        problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
        problem_options.local_parameterization_ownership = ceres::TAKE_OWNERSHIP;
        problem.~Problem();
        new (&problem) ceres::Problem(problem_options);
    }


    ~SelfAdjustedFactorGraph()
    {
        freeStateMemory();
        delete loss_function;
        loss_function = nullptr;
    }

    bool initializeOldGraph()
    {
        auto iter_pr = gnss_raw_map.begin();
        for (int m = 0; m < measSize && iter_pr != gnss_raw_map.end(); ++m, ++iter_pr)
        {
            const double time = iter_pr->first;
            const auto full_state_it = Ps_full_state.find(time);
            // 加载已优化状态
            if (full_state_it != Ps_full_state.end())
            {
                state_array[m][0] = full_state_it->second(0);
                state_array[m][1] = full_state_it->second(1);
                state_array[m][2] = full_state_it->second(2);
                state_array[m][3] = full_state_it->second(3);
                state_array[m][4] = full_state_it->second(4);
                state_array[m][5] = full_state_it->second(5);
                state_array[m][6] = full_state_it->second(6);
                continue;
            }

        }
        return true;
    }

    // 保存本轮优化结果
    bool saveGraphStateToVector(bool allsave = false)
    {
        // FactorGraph::saveGraphStateToVector(allsave);

        int length = measSize;
        int m = allsave ? 0 : (length - 1);
        Ps_full_state.clear();

        auto iter_pr = gnss_raw_map.begin();
        int idx = 0;
        while (iter_pr != gnss_raw_map.end())
        {
            if (idx >= m && idx < length)
            {
                const double time = iter_pr->first;
                Eigen::Matrix<double, 7, 1> full_state;
                full_state << state_array[idx][0], state_array[idx][1], state_array[idx][2],
                              state_array[idx][3], state_array[idx][4], state_array[idx][5], state_array[idx][6];
                Ps_full_state[time] = full_state;
            }
            ++iter_pr;
            ++idx;
        }
        return true;
    }

    

    bool addTRDDCPFactors()
    {
        tr_measurments.clear();
        TRFactorCount = 0;

        if(measSize < 2)
        {
            return false; // 需要至少两个历元的数据才能构建TR因子
        }

        auto iter = gnss_raw_map.end();
        --iter;
        // ROS_INFO("CURRENT TIME: %f  ---- LATEST ONS TIME: %f", iter->first, time_frame_now);
        
        ObservationCheck();

        // for(const auto& p : sat_lock_count_l1) { 
        //     printf("[TR DEBUG] sat=%d L1_lock=%d L2_lock=%d\n", 
        //            p.first, p.second, sat_lock_count_l2[p.first]); 
        // }

        //遍历窗口，挑选候选TR因子
        auto curr_iter = iter;
        auto gnss_data = curr_iter->second;
        for(const auto& obs : gnss_data)
        {
            if(!obs)
            {
                continue;
            }

            // 提前确定可用频段 (L1=1, L2=2), 避免后续循环覆盖
            int l1_idx=-1, l2_idx=-1;
            L1_freq(obs,&l1_idx);
            L2_freq(obs,&l2_idx);
            std::vector<int> freq_list;
            if(l1_idx >= 0) freq_list.push_back(1);
            if(l2_idx >= 0) freq_list.push_back(2);

            auto ref_iter = curr_iter;
            int i = 1;
            while (i <= windowSize && ref_iter != gnss_raw_map.begin())
            {
                --ref_iter; // move to previous epoch; i means epoch-gap (>=1)

                for(int freq_idx : freq_list)
                {
                    // 不满足观测连续性 → 跳过当前频段, 继续尝试下一个频段或gap
                    if( ((freq_idx == 1) && (sat_lock_count_l1[obs->sat] <= i)) ||
                        ((freq_idx == 2) && (sat_lock_count_l2[obs->sat] <= i)) )
                    {
                        continue;
                    }

                    TRDDMeasurement tr_meas;
                    tr_meas.freq_idx = freq_idx;
                    tr_meas.u_master_SV = obs;

                    // 在连续锁定的卫星中寻找副卫星
                    std::map<int,int> sat_lock_count;
                    if(freq_idx == 1) sat_lock_count = sat_lock_count_l1;
                    else if (freq_idx == 2) sat_lock_count = sat_lock_count_l2;
                    for(auto pair: sat_lock_count)
                    {
                        if(pair.first == obs->sat)
                        {
                            continue; // 跳过主卫星
                        }

                        if (satsys(pair.first, nullptr) != satsys(obs->sat, nullptr))
                        {
                            continue; // 仅在同一星座/系统内构建TR双差
                        }

                        if(pair.second > (i + 1)) // 连续观测至少i+1个历元（含当前）
                        {
                            findSatellitewithSameId(pair.first, gnss_data, tr_meas.u_iSV, tr_meas.freq_idx);
                            findSatellitewithSameId(pair.first, ref_iter->second, tr_meas.r_iSV,tr_meas.freq_idx);

                            // 挑出参考历元的主卫星观测
                            findSatellitewithSameId(obs->sat, ref_iter->second, tr_meas.r_master_SV,tr_meas.freq_idx);

                            tr_meas.prev_epoch_index = measSize - 1 - i;
                            tr_meas.curr_epoch_index = measSize - 1;
                            tr_meas.prev_time = ref_iter->first;
                            tr_meas.curr_time = curr_iter->first;

                            if (tr_meas.prev_epoch_index < 0 || tr_meas.prev_epoch_index >= measSize ||
                                tr_meas.curr_epoch_index < 0 || tr_meas.curr_epoch_index >= measSize ||
                                tr_meas.prev_epoch_index == tr_meas.curr_epoch_index)
                            {
                                continue;
                            }

                            // 检查星历信息
                            auto reference_sv_info = sv_info_window_map[tr_meas.prev_time];
                            auto current_sv_info = sv_info_window_map[tr_meas.curr_time];
                            if (!hasSatelliteInfo(current_sv_info, tr_meas.u_master_SV) ||
                                !hasSatelliteInfo(current_sv_info, tr_meas.u_iSV) ||
                                !hasSatelliteInfo(reference_sv_info, tr_meas.r_master_SV) ||
                                !hasSatelliteInfo(reference_sv_info, tr_meas.r_iSV))
                            {
                                continue;
                            }

                            // 将预测的移动距离关联到TR测量中
                            Eigen::Vector3d prev_pos(state_array[tr_meas.prev_epoch_index][0], state_array[tr_meas.prev_epoch_index][1], state_array[tr_meas.prev_epoch_index][2]);
                            Eigen::Vector3d pred_pos(state_array[tr_meas.curr_epoch_index-1][0], state_array[tr_meas.curr_epoch_index-1][1], state_array[tr_meas.curr_epoch_index-1][2]);
                            Eigen::Vector3d dop_vel ;
                            dop_vel.x() = doppler_map[tr_meas.curr_time].twist.twist.linear.x;
                            dop_vel.y() = doppler_map[tr_meas.curr_time].twist.twist.linear.y;
                            dop_vel.z() = doppler_map[tr_meas.curr_time].twist.twist.linear.z;
                            pred_pos = pred_pos + dop_vel*((time_frame_now - time_frame_last)/10.0);
                            double pred_move = (pred_pos - prev_pos).norm();
                            tr_meas.pred_move = pred_move;

                            if(checkRedundantPair(tr_meas))
                            {
                                continue; // 已存在相同卫星组合的TR因子
                            }

                            // 因子权重自适应调整
                            double time_discount = 1.0;
                            double p_rel = 1.0;
                            if(!SAME_RELIABLE)
                            {
                                p_rel = 1.0;
                            }

                            if(!SAME_TIME_WEIGHT)
                            {
                                time_discount = 1.0;
                            }

                            // OSQA Transformer 质量评分集成
                            // 可信卫星 (q≈1) 的 TDCP 因子获得全额权重, 低可信卫星按比例降权;
                            // 保留 0.1 的下限, 避免极端情况下优化退化为纯伪距解。
                            // TDCP 精度 (σ≈3-8cm) 远高于伪距 (σ≈1.4m), 配合 σ 模型统一后
                            // 可信卫星的载波相位约束在解中占主导。
                            double quality_score = 1.0;
#ifdef ENABLE_TRANSFORMER_BRIDGE
                            quality_score = getQualityScore(obs->sat);
#endif

                            tr_meas.tr_score = std::max(0.1, pow(p_rel,1.0/6.0) * time_discount * quality_score);
                            tr_measurments.push_back(tr_meas);

                            TRFactorCount++;
                            if(TRFactorCount >= MaxTRFactorNum)
                            {
                                break;
                            }
                        }
                    }
                    if(TRFactorCount >= MaxTRFactorNum) break;
                }
                ++i;
                if(TRFactorCount >= MaxTRFactorNum) break;
            }
            if(TRFactorCount >= MaxTRFactorNum) break;
        }

        for(int i =0; i<tr_measurments.size(); i++)
        {
            // 双差伪距因子
            auto reference_sv_info = sv_info_window_map[tr_measurments[i].prev_time];
            auto current_sv_info = sv_info_window_map[tr_measurments[i].curr_time];
            int prev_epoch_index = tr_measurments[i].prev_epoch_index;
            int epoch_index = tr_measurments[i].curr_epoch_index;

            // OSQA Transformer 质量评分集成: DD 因子质量取两颗卫星的几何平均
            double dd_pr_conf = 1.0;
#ifdef ENABLE_TRANSFORMER_BRIDGE
            if (tr_measurments[i].u_master_SV && tr_measurments[i].u_iSV)
            {
                double q_master = getQualityScore(static_cast<int>(tr_measurments[i].u_master_SV->sat));
                double q_iSV = getQualityScore(static_cast<int>(tr_measurments[i].u_iSV->sat));
                dd_pr_conf = std::sqrt(std::max(q_master, 0.01) * std::max(q_iSV, 0.01));
            }
#endif

            ceres::CostFunction* dd_pr_function =
            new ceres::AutoDiffCostFunction<DDPseudorangeFactor, 1, state_size, state_size>(
                new DDPseudorangeFactor(tr_measurments[i], current_sv_info, reference_sv_info, dd_pr_conf, tr_measurments[i].freq_idx));

            problem.AddResidualBlock(dd_pr_function, loss_function, state_array[prev_epoch_index], state_array[epoch_index]);
            
            
            ceres::CostFunction* dd_cp_function =
                new ceres::AutoDiffCostFunction<TRDDCPFactor, 1, state_size, state_size>(
                    new TRDDCPFactor(tr_measurments[i], current_sv_info, reference_sv_info, tr_measurments[i].tr_score, tr_measurments[i].freq_idx));
            problem.AddResidualBlock(dd_cp_function, loss_function, state_array[prev_epoch_index], state_array[epoch_index]);
        }

        // 因子信息统计
        int cnt_gps=0, cnt_glo=0, cnt_gal=0, cnt_bds=0;
        int cnt_l1=0, cnt_l2=0;
        for (const auto& tr_m : tr_measurments) {
            if (!tr_m.u_master_SV) continue;
            const int sys = satsys(tr_m.u_master_SV->sat, nullptr);
            if (sys==SYS_GPS)       cnt_gps++; 
            else if (sys==SYS_GLO)  cnt_glo++;
            else if (sys==SYS_GAL)  cnt_gal++; 
            else if (sys==SYS_BDS)  cnt_bds++;
           
            int l1=-1,l2=-1; 
            L1_freq(tr_m.u_master_SV,&l1); 
            L2_freq(tr_m.u_master_SV,&l2);
            if (l1>=0) cnt_l1++; 
            if (l2>=0) cnt_l2++;
        }
        printf("[ TRDDCP FACTOR] Added %d (G=%d R=%d E=%d C=%d | L1=%d L2=%d)\n",
               (int)tr_measurments.size(), cnt_gps, cnt_glo, cnt_gal, cnt_bds, cnt_l1, cnt_l2);
        return true;
    }

    void resizeMaxTRFactorNum(double run_time,double max_running_time_ms)
    {
        const int min_tr = std::max(10, MinTRFactorNum);
        if(run_time < 10.0)
        {
            MaxTRFactorNum = 100;
            return; // 优化时间在可接受范围内，无需调整
        }
        MaxTRFactorNum = (MaxTRFactorNum / (run_time / max_running_time_ms)); // 如果优化时间过长，减少TR因子数量上限
        if (MaxTRFactorNum < min_tr)
        {
            MaxTRFactorNum = min_tr;
        }
        printf("Adjusted MaxTRNum : %d. \n", MaxTRFactorNum);
    }

    bool checkRedundantPair(TRDDMeasurement tr_meas)
    {
        for(int i=0;i<tr_measurments.size();i++)
        {
            // 主主 == 副副
            if(tr_meas.r_master_SV->sat == tr_measurments[i].r_master_SV->sat 
                && tr_meas.r_iSV->sat == tr_measurments[i].r_iSV->sat 
                && tr_meas.curr_epoch_index == tr_measurments[i].curr_epoch_index 
                && tr_meas.prev_epoch_index == tr_measurments[i].prev_epoch_index
                && tr_meas.freq_idx == tr_measurments[i].freq_idx)
            {
                // ROS_INFO("Redundant pair found: Master SV %d, iSV %d, Epoch gap %d", tr_meas.r_master_SV->sat, tr_meas.r_iSV->sat, tr_meas.curr_epoch_index - tr_meas.prev_epoch_index);
                return true;
            }
            // 主副 == 副主
            else if(tr_meas.r_master_SV->sat == tr_measurments[i].r_iSV->sat 
                && tr_meas.r_iSV->sat == tr_measurments[i].r_master_SV->sat
                && tr_meas.curr_epoch_index == tr_measurments[i].curr_epoch_index
                && tr_meas.prev_epoch_index == tr_measurments[i].prev_epoch_index
                && tr_meas.freq_idx == tr_measurments[i].freq_idx)
            {
                // ROS_INFO("Redundant pair found (reversed): Master SV %d, iSV %d, Epoch gap %d", tr_meas.r_master_SV->sat, tr_meas.r_iSV->sat, tr_meas.curr_epoch_index - tr_meas.prev_epoch_index);
                return true;
            }
        }
        return false;
    }

    void ObservationCheck()
    {
        // 重置卫星连续观测检查标志
        for (auto& pair : sat_lock_check) 
        {
            pair.second = false;
        }
        
        // 由现在向过去检查
        auto iter = gnss_raw_map.end();
        --iter;
        if(iter!=gnss_raw_map.begin())
        {
            for(int i=0;i<iter->second.size();i++)
            {
                int sat_id = iter->second[i]->sat;
                if(sat_lock_check.find(sat_id) == sat_lock_check.end() )
                {
                    // 没有该卫星的检查记录--首次发现
                    sat_lock_check[sat_id] = true;
                    sat_lock_count_l1[sat_id] = 1;
                    sat_lock_count_l2[sat_id] = 1;
                    continue;
                }
                
                int l1_idx=-1,l2_idx=-1;
                L1_freq(iter->second[i],&l1_idx);
                L2_freq(iter->second[i],&l2_idx);
                if(l1_idx >= 0)
                {
                    if(!cycleSlipDetect(sat_id, 1))
                        sat_lock_count_l1[sat_id]++;
                }
                if(l2_idx >= 0)
                {
                    if(!cycleSlipDetect(sat_id, 2))
                        sat_lock_count_l2[sat_id]++;
                }
                if(l1_idx < 0 && l2_idx < 0)
                {
                    continue;
                }
                // 无论有没有周跳都标记为检查过了
                sat_lock_check[sat_id] = true;

            }
            // 未检查过的失锁判定
            for(auto& pair : sat_lock_check)
            {
                if(pair.second == false)
                {
                    sat_lock_count_l1[pair.first] = 0;
                    sat_lock_count_l2[pair.first] = 0;
                    // ROS_INFO("Satellite %d lock LOSE: %d", pair.first, sat_lock_count_l1[pair.first]);
                }
            }
        }
    }

    // 计算当前历元对 TDEC 残差 (Δcp + doppler_mean·Δt) 的全体卫星中位数, 按历元+频段缓存。
    // 样本 <4 时返回 0 (早期历元统计不足, 不做共模修正)。
    double tdecCommonTerm(int freq)
    {
        if (tdec_common_epoch_ != time_frame_now)
        {
            tdec_common_l1_ = 0.0;
            tdec_common_l2_ = 0.0;
            tdec_common_epoch_ = time_frame_now;
        }
        double &cached = (freq == 1) ? tdec_common_l1_ : tdec_common_l2_;
        bool &cached_valid = (freq == 1) ? tdec_common_valid_l1_ : tdec_common_valid_l2_;
        if (cached_valid) return cached;                        // 本历元已计算

        if (gnss_raw_map.size() < 2) return 0.0;
        const auto allobs_now = gnss_raw_map.rbegin();
        const auto allobs_prev = std::next(allobs_now);
        const auto prev_sv_it = sv_info_window_map.find(allobs_prev->first);
        const auto curr_sv_it = sv_info_window_map.find(allobs_now->first);
        if (prev_sv_it == sv_info_window_map.end() || curr_sv_it == sv_info_window_map.end())
        {
            return 0.0;
        }

        std::vector<double> samples;
        for (const auto &obs : allobs_now->second)
        {
            if (!obs) continue;
            const int sat = static_cast<int>(obs->sat);
            int ci = -1, pi = -1;
            if (freq == 1) { L1_freq(obs, &ci); } else { L2_freq(obs, &ci); }
            if (ci < 0 || ci >= static_cast<int>(obs->cp.size()) ||
                ci >= static_cast<int>(obs->dopp.size())) continue;
            const double curr_cp = obs->cp[ci];
            const double dopp_curr = obs->dopp[ci];
            double prev_cp = 0.0, dopp_prev = 0.0, prev_sec = 0.0, curr_sec = 0.0;
            bool found = false;
            for (const auto &prev_obs : allobs_prev->second)
            {
                if (prev_obs && static_cast<int>(prev_obs->sat) == sat)
                {
                    if (freq == 1) { L1_freq(prev_obs, &pi); } else { L2_freq(prev_obs, &pi); }
                    if (pi >= 0 && pi < static_cast<int>(prev_obs->cp.size()) &&
                        pi < static_cast<int>(prev_obs->dopp.size()))
                    {
                        prev_cp = prev_obs->cp[pi];
                        dopp_prev = prev_obs->dopp[pi];
                        prev_sec = prev_obs->time.time + prev_obs->time.sec;
                        curr_sec = obs->time.time + obs->time.sec;
                        found = true;
                    }
                    break;
                }
            }
            if (!found || prev_cp == 0.0 || curr_cp == 0.0) continue;
            const double dt = curr_sec - prev_sec;
            if (!std::isfinite(dt) || dt <= 0.0) continue;
            samples.push_back((curr_cp - prev_cp) + 0.5 * (dopp_curr + dopp_prev) * dt);
        }
        if (samples.size() < 4) return 0.0;   // 样本不足: 本历元保持共模 0
        cached_valid = true;
        std::sort(samples.begin(), samples.end());
        const size_t mid = samples.size() / 2;
        cached = (samples.size() % 2 == 1) ? samples[mid]
                                           : 0.5 * (samples[mid - 1] + samples[mid]);
        return cached;
    }

    // for single freq check
    bool cycleSlipDetect(int sat_id, int freq)
    {
        auto mark_slip = [this, sat_id, freq](int reason) {
            if(freq == 1)
            {
                sat_lock_count_l1[sat_id] = 0;
                sat_lock_check_l1[sat_id] = true;
            }
            else
            {
                sat_lock_count_l2[sat_id] = 0;
                sat_lock_check_l2[sat_id] = true;
            }
            ROS_INFO("SAT: %d [%d] --- SLIP REASON: %d",sat_id, freq, reason);
        };

        if (gnss_raw_map.size() < 2 || measSize < 2)
        {
            mark_slip(1);
            return true;
        }

        // 只与上一个历元的观测对比
        const auto allobs_now = std::prev(gnss_raw_map.end());
        const auto allobs_prev = std::prev(allobs_now);

        gnss_comm::ObsPtr prev_obs;
        gnss_comm::ObsPtr curr_obs;
        findSatellitewithSameId(sat_id, allobs_prev->second, prev_obs,freq);
        findSatellitewithSameId(sat_id, allobs_now->second, curr_obs,freq);
        if (!prev_obs || !curr_obs)
        {
            // 没有该卫星的观测
            mark_slip(2);
            return true;
        }

        double prev_cp_cycle = 0.0;
        double curr_cp_cycle = 0.0;
        if (!getValidCarrierPhase(prev_obs, prev_cp_cycle,freq) ||
            !getValidCarrierPhase(curr_obs, curr_cp_cycle,freq))
        {
            // 载波相位、观测值异常
            mark_slip(3);
            return true;
        }

        const double prev_time_sec = prev_obs->time.time + prev_obs->time.sec;
        const double curr_time_sec = curr_obs->time.time + curr_obs->time.sec;
        const double dt = curr_time_sec - prev_time_sec;
        if (!std::isfinite(dt) || dt <= 0.0)
        {  
            // 时间顺序错误
            mark_slip(4);
            return true;
        }
       
        const auto prev_sv_it = sv_info_window_map.find(allobs_prev->first);
        const auto curr_sv_it = sv_info_window_map.find(allobs_now->first);
        if (prev_sv_it == sv_info_window_map.end() || curr_sv_it == sv_info_window_map.end())
        {
            // 没有最新的卫星星历
            mark_slip(5);
            return true;
        }

        const auto prev_sat_it = prev_sv_it->second.find(sat_id);
        const auto curr_sat_it = curr_sv_it->second.find(sat_id);
        if (prev_sat_it == prev_sv_it->second.end() || curr_sat_it == curr_sv_it->second.end())
        {
            // 没有该卫星星历
            mark_slip(6);
            return true;
        }


        // TDEC 周跳检测: 载波相位历元间差分 vs 该星自身多普勒积分
        // 预测 Δcp_pred = -doppler·Δt (相邻历元多普勒取均值, 梯形积分)。
        // 该检验只依赖观测流内部的 cp↔doppler 相干性, 不依赖星历几何:
        // - 卫星运动与接收机钟漂均已体现在 doppler 中, 无需额外建模;
        // - 周跳使 cp 跳变 k 个整周而 doppler 连续, 残差出现 λ·k 量级跳变;
        // - 旧实现用锁定时刻快照 sat_cp_const 抵消卫星运动项, 快照随时间失稳,
        //   且预测式漏掉卫星运动项, 导致每 3~5 历元必然误报周跳、锁定计数反复清零。
        double dopp_curr = 0.0, dopp_prev = 0.0;
        int dopp_idx = -1;
        if (freq == 1) { L1_freq(curr_obs, &dopp_idx); }
        else           { L2_freq(curr_obs, &dopp_idx); }
        if (dopp_idx >= 0 && dopp_idx < static_cast<int>(curr_obs->dopp.size()))
        {
            dopp_curr = curr_obs->dopp[dopp_idx];
        }
        dopp_idx = -1;
        if (freq == 1) { L1_freq(prev_obs, &dopp_idx); }
        else           { L2_freq(prev_obs, &dopp_idx); }
        if (dopp_idx >= 0 && dopp_idx < static_cast<int>(prev_obs->dopp.size()))
        {
            dopp_prev = prev_obs->dopp[dopp_idx];
        }

        const double delta_cp_pred = -0.5 * (dopp_curr + dopp_prev) * dt;
        const double delta_cp_obs = curr_cp_cycle - prev_cp_cycle;

        // 扣除本历元对所有卫星共同的 TDEC 偏移 (中位数):
        // 接收机钟跳步/钟瞬态会同步影响全部卫星, 属接收机行为而非周跳;
        // 单颗卫星的真实周跳仍会相对中位数偏离而被捕获
        const double tdec_common = tdecCommonTerm(freq);

        // 阈值 2 周: 实测数据 cp↔doppler 相干噪声远小于 1 周,
        // 取 2 周在保证整周跳变捕获的同时抑制多普勒噪声引起的误报
        constexpr double kCycleSlipThresholdCycle = 2.0;
        if(std::abs(delta_cp_obs - delta_cp_pred - tdec_common) > kCycleSlipThresholdCycle)
        {
            mark_slip(9);
            ROS_INFO("Satellite %d: Delta CP Pred=%.3f cycles, Delta CP Obs=%.3f cycles", sat_id, delta_cp_pred, delta_cp_obs);
            return true;
        }

        return false;
    }

    // 求解结果发散防护: 正常情况下 FGO 与 SPP 初值偏差在数十米内。
    // 个别历元(如接收机钟跳步引发约束重构时)优化器会收敛到公里级错误盆地
    // 且 opt_time 飙升, 此时把最新历元状态回退为 SPP 初值, 防止异常输出与
    // 滑窗状态污染 (下一历元以正常约束重新求解, 观测上表现为单点恢复)。
    bool guardLatestStateAgainstDivergence(const Eigen::Vector3d &spp_init,
                                           double max_jump_m = 50.0)
    {
        if (measSize <= 0 || state_array.empty() || state_array[measSize - 1] == nullptr)
        {
            return false;
        }
        const int last = measSize - 1;
        const Eigen::Vector3d solved(state_array[last][0], state_array[last][1], state_array[last][2]);
        if (!std::isfinite(solved.norm()))
        {
            ROS_WARN("[FGO Guard] non-finite state, fallback to SPP init");
            state_array[last][0] = spp_init.x();
            state_array[last][1] = spp_init.y();
            state_array[last][2] = spp_init.z();
            return true;
        }
        const double jump = (solved - spp_init).norm();
        if (jump > max_jump_m)
        {
            ROS_WARN("[FGO Guard] state jumped %.0f m from SPP init, fallback (solver divergence)", jump);
            state_array[last][0] = spp_init.x();
            state_array[last][1] = spp_init.y();
            state_array[last][2] = spp_init.z();
            return true;
        }
        return false;
    }

    // 根据优化结果更新载波相位变化残差 (dop_cp, 单位: 米, 与 psr 残差量纲统一)
    // 残差 = λ × [Δcp + doppler_mean·Δt] —— 载波相位与该星多普勒积分的相干性 (TDEC):
    // 观测链路健康时为厘米级, 周跳表现为 λ 整数倍的跳变, 多径/跟踪退化时增大。
    // 卫星运动与接收机钟漂均已含在 doppler 中, 故不依赖星历几何,
    // 与 cycleSlipDetect 的 TDEC 周跳检验保持同一物理模型。
    void CPresidualsUpdate()
    {
        residuals.clear();
        auto curr_obs_it = gnss_raw_map.find(time_frame_now);
        if (curr_obs_it == gnss_raw_map.end())
        {
            return;
        }
        for (const auto &obs : curr_obs_it->second)
        {
            if (!obs) continue;
            int sat = static_cast<int>(obs->sat);
            residuals[sat] = std::map<std::string, double>({{"dop_cp", 0.0}});
        }

        auto prev_obs_it = gnss_raw_map.find(time_frame_last);
        if (prev_obs_it == gnss_raw_map.end())
        {
            return;
        }
        const auto curr_sv_it = sv_info_window_map.find(time_frame_now);
        if (curr_sv_it == sv_info_window_map.end())
        {
            return;
        }

        for (const auto &obs : curr_obs_it->second)
        {
            if (!obs) continue;
            const int sat = static_cast<int>(obs->sat);

            int l1_idx = -1;
            L1_freq(obs, &l1_idx);
            if (l1_idx < 0 ||
                l1_idx >= static_cast<int>(obs->cp.size()) ||
                l1_idx >= static_cast<int>(obs->dopp.size()))
            {
                continue;
            }
            const double curr_cp_cycle = obs->cp[l1_idx];
            const double dopp_curr = obs->dopp[l1_idx];
            const double curr_time_sec = obs->time.time + obs->time.sec;

            // 上一历元同卫星的 L1 载波相位与多普勒
            double prev_cp_cycle = 0.0, dopp_prev = 0.0, prev_time_sec = 0.0;
            bool prev_found = false;
            for (const auto &prev_obs : prev_obs_it->second)
            {
                if (prev_obs && static_cast<int>(prev_obs->sat) == sat)
                {
                    int prev_l1_idx = -1;
                    L1_freq(prev_obs, &prev_l1_idx);
                    if (prev_l1_idx >= 0 &&
                        prev_l1_idx < static_cast<int>(prev_obs->cp.size()) &&
                        prev_l1_idx < static_cast<int>(prev_obs->dopp.size()))
                    {
                        prev_cp_cycle = prev_obs->cp[prev_l1_idx];
                        dopp_prev = prev_obs->dopp[prev_l1_idx];
                        prev_time_sec = prev_obs->time.time + prev_obs->time.sec;
                        prev_found = true;
                    }
                    break;
                }
            }
            if (!prev_found || prev_cp_cycle == 0.0 || curr_cp_cycle == 0.0)
            {
                continue;
            }
            const double dt_sec = curr_time_sec - prev_time_sec;
            if (!std::isfinite(dt_sec) || dt_sec <= 0.0)
            {
                continue;
            }

            const auto lam_it = curr_sv_it->second.find(sat);
            if (lam_it == curr_sv_it->second.end())
            {
                continue;
            }
            const double lamda_l1 = lam_it->second.lamda_l1;
            if (!std::isfinite(lamda_l1) || lamda_l1 <= 0.0)
            {
                continue;
            }

            const double resid_cycles = (curr_cp_cycle - prev_cp_cycle) + 0.5 * (dopp_curr + dopp_prev) * dt_sec;
            const double dop_cp_residual_m = lamda_l1 * resid_cycles;
            if (std::isfinite(dop_cp_residual_m) && std::abs(dop_cp_residual_m) < 1e4)
            {
                residuals[sat]["dop_cp"] = dop_cp_residual_m;
            }
        }
    }
    // 输出ENU结果估计的协方差
    bool printLatestPosCovarianceENU(Eigen::Vector3d& result_enu) const
    {
        if (covMatrix.rows() < 3 || covMatrix.cols() < 3)
        {
            ROS_WARN("[TRDDCPFactorGraph::printLatestPosCovarianceENU] covariance matrix is unavailable");
            return false;
        }

        if (!std::isfinite(enu_ref_llh.norm()) || enu_ref_llh.norm() <= 0.0)
        {
            ROS_WARN("[TRDDCPFactorGraph::printLatestPosCovarianceENU] ENU reference is unavailable");
            return false;
        }

        const Eigen::Matrix3d cov_ecef = covMatrix.block<3, 3>(0, 0);
        const Eigen::Matrix3d R_enu_ecef = gnss_comm::geo2rotation(enu_ref_llh).transpose();
        const Eigen::Matrix3d cov_enu = R_enu_ecef * cov_ecef * R_enu_ecef.transpose();

        std::ostringstream oss;
        #if SHOW_COVMAX
        oss << std::scientific << std::setprecision(6)
            << "Latest ENU covariance [m^2]:\n"
            << cov_enu;
        ROS_INFO_STREAM(oss.str());
        #endif
        ROS_INFO("Latest ENU sigma [m]: E=%.4f N=%.4f U=%.4f",
                 std::sqrt(std::max(0.0, cov_enu(0, 0))),
                 std::sqrt(std::max(0.0, cov_enu(1, 1))),
                 std::sqrt(std::max(0.0, cov_enu(2, 2))));
        result_enu = Eigen::Vector3d(std::sqrt(std::max(0.0, cov_enu(0, 0))),
                                    std::sqrt(std::max(0.0, cov_enu(1, 1))),
                                    std::sqrt(std::max(0.0, cov_enu(2, 2))));
        return true;
    }

};
