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
#define ar_size 100
#define residual_gate_threshold 0.1

class SelfAdjustedFactorGraph: public FactorGraph{
    
public:
    std::vector<TRDDMeasurement> tr_measurments; // 当前窗口内的TR测量列表

    std::map<double, Eigen::Matrix<double, 7, 1>> Ps_full_state; // 保存完整状态（位置+钟差）

    std::map<int,int> sat_lock_count_l1,sat_lock_count_l2;  // 卫星连续观测计数器L1,L2频段分开统计
    std::map<int,bool> sat_lock_check_l1,sat_lock_check_l2; // 卫星连续观测检查标志
    std::map<int,double> sat_cp_const_l1,sat_cp_const_l2;   // 卫星TR锁定状态
    std::map<int,double> sat_const_residual_l1,sat_const_residual_l2; // 预拟合残差基线（按卫星维护）

    int MaxTRFactorNum = 100;   // TR因子数量上限，动态调整以控制优化时间
    int TRFactorCount;          // 统计优化规模

    bool SAME_RELIABLE = false;
    bool SAME_TIME_WEIGHT = false;

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

        //遍历窗口，挑选候选TR因子
        auto curr_iter = iter;
        auto gnss_data = curr_iter->second;
        for(const auto& obs : gnss_data)
        {
            if(!obs) // 连续观测不足2个历元的卫星不考虑构建TR因子
            {
                continue;
            }
            else
            {
                auto ref_iter = curr_iter;
                int i = 1;
                while (i <= windowSize && it != gnss_raw_map.begin())
                {
                    // ROS_INFO("FIND TR FACTOR: Satellite %d, Epoch Gap: %d", obs->sat, i);
                    --ref_iter; // move to previous epoch; i means epoch-gap (>=1)
                    //构建双差因子L1/L2：
                    for(int j=0;j<obs->freq.size();j++)
                    {
                        int l1_idx=-1,l2_idx=-1;
                        L1_freq(obs,&l1_idx);
                        L2_freq(obs,&l2_idx);

                        int freq_idx = -1;

                        // l1_idx/l2_idx 可能都是1
                        if(l1_idx >= 0 && j == l1_idx)
                        {
                            // L1频段的双差因子构建
                            freq_idx = 1;
                        }

                        else if(l2_idx >= 0 && j == l2_idx)
                        {
                            // L2频段的双差因子构建
                            freq_idx = 2;
                        }
                    }
                    TRDDMeasurement tr_meas;
                    tr_meas.freq_idx = freq_idx;
                    
                    tr_meas.u_master_SV = obs;

                    // 不满足观测连续性
                    if( ((freq_idx == 1) && (sat_lock_count_l1[obs->sat] <= i)) || 
                        ((freq_idx == 2) && (sat_lock_count_l2[obs->sat] <= i)) )
                    {
                        break;
                    }
                    
                    // 在连续锁定的卫星中寻找副卫星
                    std::map<int,int> sat_lock_count;
                    if(freq_idx == 1) sat_lock_count = sat_lock_count_l1;
                    else if (freq_idx == 2) sat_lock_count = sat_lock_count_l2;
                    for(auto pair: sat_lock_count)
                    {
                        // ROS_INFO("CHECKING :  Satellite %d lock count: %d", pair.first, pair.second);
                        if(pair.first == obs->sat)
                        {
                            continue; // 跳过主卫星
                        }

                        if (satsys(pair.first, nullptr) != satsys(obs->sat, nullptr))
                        {
                            continue; // 仅在同一星座/系统内构建TR双差，避免系统间钟差引入偏置 TODO:加入卫星枢纽，取消系统限制
                        }

                        if(pair.second > (i + 1)) // 连续观测至少i+1个历元（含当前），增加TR因子稳定性
                        {
                            // printf("selected: sat_id: %d,  lock_count: %d \n", pair.first, pair.second);
                            // 找到一个**同频段**满足条件的副卫星，构建TR因子

                            //TODO: 此时的副卫星没有任何特征筛选机制，
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
                                // ROS_INFO("skip reason 1");
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
                                // ROS_INFO("skip reason 2");
                                continue;
                            }

                            // 将预测的移动距离关联到TR测量中，后续时间折扣函数使用
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
                                // ROS_INFO("skip reason 3");
                                continue; // 已存在相同卫星组合的TR因子，跳过以避免冗余
                            }

                            // 因子权重自适应调整（TODO）
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

                            tr_meas.tr_score = pow(p_rel,1.0/6.0) * time_discount ; // 综合可靠性得分，作为TR因子的权重依据
                            tr_measurments.push_back(tr_meas);   //添加到候选双差因子列表


                            TRFactorCount++;
                            // 达到数量上限，退出
                            if(TRFactorCount >= MaxTRFactorNum)
                            {
                                break; 
                            }
                            
                        }
                    }
                    ++i;
                    if(TRFactorCount >= MaxTRFactorNum)
                    {
                        break; 
                    }
                }
            }
            if(TRFactorCount >= MaxTRFactorNum)
            {
               break;
            }
        }

        for(int i =0; i<tr_measurments.size(); i++)
        {
            // 双差伪距因子
            auto reference_sv_info = sv_info_window_map[tr_measurments[i].prev_time];
            auto current_sv_info = sv_info_window_map[tr_measurments[i].curr_time];
            int prev_epoch_index = tr_measurments[i].prev_epoch_index;
            int epoch_index = tr_measurments[i].curr_epoch_index;
            

            ceres::CostFunction* dd_pr_function =
            new ceres::AutoDiffCostFunction<DDPseudorangeFactor, 1, state_size, state_size>(
                new DDPseudorangeFactor(tr_measurments[i], current_sv_info, reference_sv_info, 1.0, tr_measurments[i].freq_idx));
            
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
        if(run_time < 10.0)
        {
            MaxTRFactorNum = 100;
            return; // 优化时间在可接受范围内，无需调整
        }
        MaxTRFactorNum = (MaxTRFactorNum / (run_time / max_running_time_ms)); // 如果优化时间过长，减少TR因子数量上限
        if (MaxTRFactorNum < 10)
        {
            MaxTRFactorNum = 10;
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

    bool ObservationCheck()
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
                    cycleSlipDetect(sat_id, 1)
                }
                if(l2_idx >= 0)
                {
                    cycleSlipDetect(sat_id, 2)
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

    // for single freq check
    bool cycleSlipDetect(int sat_id, int freq)
    {
        auto mark_slip = [this, sat_id, int freq]() {
            if(freq == 1)
            {
                sat_lock_count_l1[sat_id] = 0;
                sat_lock_check_l1[sat_id] = true;
                sat_cp_const_l1.erase(sat_id);
                sat_const_residual_l1.erase(sat_id);
            }
            else
            {
                sat_lock_count_l2[sat_id] = 0;
                sat_lock_check_l2[sat_id] = true;
                sat_cp_const_l2.erase(sat_id);
                sat_const_residual_l2.erase(sat_id);
            }
            
        };

        if (gnss_raw_map.size() < 2 || measSize < 2)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 1");
            return true;
        }

        // 只与上一个历元的观测对比
        const auto allobs_now = std::prev(gnss_raw_map.end());
        const auto allobs_prev = std::prev(allobs_now);

        gnss_comm::ObsPtr prev_obs;
        gnss_comm::ObsPtr curr_obs;
        findSatellitewithSameId(sat_id, allobs_prev->second, prev_obs);
        findSatellitewithSameId(sat_id, allobs_now->second, curr_obs);
        if (!prev_obs || !curr_obs)
        {
            // 没有该卫星的观测
            mark_slip();
            return true;
        }

        double prev_cp_cycle = 0.0;
        double curr_cp_cycle = 0.0;
        if (!getValidCarrierPhase(prev_obs, prev_cp_cycle,freq) ||
            !getValidCarrierPhase(curr_obs, curr_cp_cycle,freq))
        {
            // 载波相位、观测值异常
            mark_slip();
            return true;
        }

        const double prev_time_sec = prev_obs->time.time + prev_obs->time.sec;
        const double curr_time_sec = curr_obs->time.time + curr_obs->time.sec;
        const double dt = curr_time_sec - prev_time_sec;
        if (!std::isfinite(dt) || dt <= 0.0)
        {
            // 时间顺序错误
            mark_slip();
            return true;
        }
       
        const auto prev_sv_it = sv_info_window_map.find(allobs_prev->first);
        const auto curr_sv_it = sv_info_window_map.find(allobs_now->first);
        if (prev_sv_it == sv_info_window_map.end() || curr_sv_it == sv_info_window_map.end())
        {
            // 没有最新的卫星星历
            mark_slip();
            return true;
        }

        const auto prev_sat_it = prev_sv_it->second.find(sat_id);
        const auto curr_sat_it = curr_sv_it->second.find(sat_id);
        if (prev_sat_it == prev_sv_it->second.end() || curr_sat_it == curr_sv_it->second.end())
        {
            // 没有该卫星星历
            mark_slip();
            return true;
        }


        const sv_info& prev_sv = prev_sat_it->second;
        const sv_info& curr_sv = curr_sat_it->second;

        double curr_lambda,prev_lambda;
        if(freq == 1)
        {
            curr_lambda = (curr_sv.lamda_l1 > 0.0) ? curr_sv.lamda_l1 : 0;
            prev_lambda = (prev_sv.lamda_l1 > 0.0) ? prev_sv.lamda_l1 : 0;
        }
        else if(freq == 2)
        {
            curr_lambda = (curr_sv.lamda_l2 > 0.0) ? curr_sv.lamda_l2 : 0;
            prev_lambda = (prev_sv.lamda_l2 > 0.0) ? prev_sv.lamda_l2 : 0;
        }
        
        if (!std::isfinite(curr_lambda) || curr_lambda <= 0.0 || !std::isfinite(prev_lambda) || prev_lambda <= 0.0 || prev_lambda != curr_lambda)
        {
            // 卫星波长无效
            mark_slip();
            return true;
        }

        //更新sat_cp_const以适应卫星状态的变化，避免过时的常数导致误判周跳
        std::map<int,double> sat_cp_const;
        if(freq == 1) sat_cp_const = sat_cp_const_l1;
        else if(freq == 2) sat_cp_const = sat_cp_const_l2;
        sat_cp_const[sat_id] = curr_cp_cycle - prev_cp_cycle;
        
        const Eigen::Vector3d prev_vel{doppler_map[allobs_prev->first].twist.twist.linear.x,
                                    doppler_map[allobs_prev->first].twist.twist.linear.y,
                                    doppler_map[allobs_prev->first].twist.twist.linear.z};
        const Eigen::Vector3d curr_vel{doppler_map[allobs_now->first].twist.twist.linear.x,
                                    doppler_map[allobs_now->first].twist.twist.linear.y,
                                    doppler_map[allobs_now->first].twist.twist.linear.z};
                

        // 积分多普勒辅助的载波相位检查
        const Eigen::Vector3d prev_rcv(state_array[measSize - 2][0], state_array[measSize - 2][1], state_array[measSize - 2][2]);
        Eigen::Vector3d curr_rcv;
        // 时间戳被放大10倍保存，所以/10
        curr_rcv.x() = prev_rcv.x() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[allobs_now->first].twist.twist.linear.x;
        curr_rcv.y() = prev_rcv.y() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[allobs_now->first].twist.twist.linear.y;
        curr_rcv.z() = prev_rcv.z() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[allobs_now->first].twist.twist.linear.z;
        const Eigen::Vector3d prev_sat(prev_sv.pos[0], prev_sv.pos[1], prev_sv.pos[2]);
        const Eigen::Vector3d curr_sat(curr_sv.pos[0], curr_sv.pos[1], curr_sv.pos[2]);
        
        // LOS卫地距离
        const Eigen::Vector3d prev_los = prev_sat - prev_rcv;
        const Eigen::Vector3d curr_los = curr_sat - curr_rcv;
        if (!std::isfinite(prev_los.norm()) || !std::isfinite(curr_los.norm()) ||
            prev_los.norm() <= 0.0 || curr_los.norm() <= 0.0)
        {
            // 卫地距离异常
            mark_slip();
            return true;
        }

        Eigen::Vector3d los = prev_los + curr_los;
        if (!std::isfinite(los.norm()) || los.norm() <= 0.0)
        {
            los = curr_los;
        }
        else
        {
            los.normalize();
        }

        // 载波相位变化预测值
        const double delta_cp_pred = ((curr_vel + prev_vel).dot(los) * (curr_time_sec - prev_time_sec) / 10.0)/(2 * lambda);
        // 载波相位变化观测值

        double sat_cp_const = (freq==1) ? sat_cp_const_l1[sat_id] : ((freq==2) ? sat_cp_const_l2[sat_id] : 0.0)

        const double delta_cp_obs = (curr_cp_cycle - prev_cp_cycle) - sat_cp_const;
        // printf(" | obs delta cp: %f \n",(curr_cp_cycle - prev_cp_cycle));
        
        constexpr double kCycleSlipThresholdCycle = 1.0; // 预测和观测的载波相位变化超过1周期则判定为可能发生了周跳
        if(std::abs(delta_cp_obs - delta_cp_pred) > kCycleSlipThresholdCycle)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 10");
            ROS_INFO("Satellite %d: Delta CP Pred=%.3f cycles, Delta CP Obs=%.3f cycles", sat_id, delta_cp_pred, delta_cp_obs);
            return true;
        }

        return false;
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
