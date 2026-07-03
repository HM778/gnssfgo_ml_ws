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
#include "bucket_struct.hpp"
#include "../tools/checking.hpp"

#define state_size 7 // x,y,z, clk_gps, clk_glo, clk_gal, clk_bds
#define ar_size 100
#define residual_gate_threshold 0.1

struct BucketContributionRecord
{
    BucketKey key;
    double residual = 0.0;
    double gap_epoch = 0.0;
    double delta_pos = 0.0;
    bool success = false;
};

struct EpochBucketSnapshot
{
    int epoch_index = -1;
    std::vector<BucketContributionRecord> contributions;
};

class BinningFactorGraph: public FactorGraph{
    
public:
    BucketManager bucket_manager;
    BucketIndexer indexer;
    std::vector<TRfeature> tr_factor_features;
    std::vector<TRRTKMeasurement> tr_measurments; // 当前窗口内的TR测量列表

    std::map<double, Eigen::Matrix<double, 7, 1>> Ps_full_state; // 保存完整状态（位置+钟差）

    std::map<int,int> sat_lock_count; // 卫星连续观测计数器
    std::map<int,bool> sat_lock_check; // 卫星连续观测检查标志
    std::map<int,double> sat_cp_const; // 卫星TR锁定状态
    std::map<int,double> sat_const_residual; // 预拟合残差基线（按卫星维护）
    std::map<int, EpochBucketSnapshot> epoch_bucket_history;

    int MaxTRFactorNum = 100; // TR因子数量上限，动态调整以控制优化时间
    int TRFactorCount;
    int recent_history_epoch_window = 30; // 仅统计最近 N 个历元的桶统计窗口

    bool SAME_RELIABLE = false;
    bool SAME_TIME_WEIGHT = false;

public:
    BinningFactorGraph()
    {
        tr_factor_features.reserve(windowSize); // 预分配空间，减少动态扩容的开销
        loss_function = nullptr;
        problem_options.cost_function_ownership = ceres::TAKE_OWNERSHIP;
        problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
        problem_options.local_parameterization_ownership = ceres::TAKE_OWNERSHIP;
        problem.~Problem();
        new (&problem) ceres::Problem(problem_options);
    }


    ~BinningFactorGraph()
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

            for (const auto& saved_state : Ps)
            {
                if (saved_state.first == time)
                {
                    state_array[m][0] = saved_state.second.x();
                    state_array[m][1] = saved_state.second.y();
                    state_array[m][2] = saved_state.second.z();
                    if (m > 0)
                    {
                        state_array[m][3] = state_array[m - 1][3];
                        state_array[m][4] = state_array[m - 1][4];
                        state_array[m][5] = state_array[m - 1][5];
                        state_array[m][6] = state_array[m - 1][6];
                    }
                    else
                    {
                        state_array[m][3] = 0.0;
                        state_array[m][4] = 0.0;
                        state_array[m][5] = 0.0;
                        state_array[m][6] = 0.0;
                    }
                    break;
                }
            }
        }
        // Debug: compare this initialization with previously saved optimal states
        if (!Ps_full_state.empty())
        {
            int idx = 0;
            auto it = gnss_raw_map.begin();
            double global_max_diff = 0.0;
            int mismatch_count = 0;
            for (; idx < measSize && it != gnss_raw_map.end(); ++idx, ++it)
            {
                const double time = it->first;
                auto prev_it = Ps_full_state.find(time);
                if (prev_it != Ps_full_state.end())
                {
                    Eigen::Matrix<double, 7, 1> prev = prev_it->second;
                    Eigen::Matrix<double, 7, 1> cur;
                    cur << state_array[idx][0], state_array[idx][1], state_array[idx][2],
                            state_array[idx][3], state_array[idx][4], state_array[idx][5], state_array[idx][6];
                    Eigen::Matrix<double, 7, 1> diff = (cur - prev).cwiseAbs();
                    double local_max = diff.maxCoeff();
                    if (local_max > 1e-9)
                    {
                        ++mismatch_count;
                        // ROS_INFO("[TRDDCP INIT DEBUG] epoch_idx=%d time=%.6f max_abs_diff=%.12f",
                                //  idx, time, local_max);
                    }
                    if (local_max > global_max_diff)
                    {
                        global_max_diff = local_max;
                    }
                }
                else
                {
                    // ROS_INFO("[TRDDCP INIT DEBUG] epoch_idx=%d time=%.6f no previous saved optimal state", idx, it->first);
                }
            }
            // ROS_INFO("[TRDDCP INIT DEBUG] compared=%d mismatches=%d global_max_diff=%.12f", measSize, mismatch_count, global_max_diff);
        }
        else
        {
            // ROS_INFO("[TRDDCP INIT DEBUG] Ps_full_state empty: no previous optimal states to compare");
        }

        return true;
    }

    bool saveGraphStateToVector(bool allsave = false)
    {
        FactorGraph::saveGraphStateToVector(allsave);

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
                              state_array[idx][3], state_array[idx][4],
                              state_array[idx][5], state_array[idx][6];
                Ps_full_state[time] = full_state;
            }
            ++iter_pr;
            ++idx;
        }
        return true;
    }

    
    // 计算TR双差因子残差的逻辑，只保留历史时刻的信息，不涉及当前时刻的状态
    double computeTRResidual(const TRRTKMeasurement& dd_measurement, Eigen::Vector3d& prev_pos)
    {
        // 计算TR双差因子残差的逻辑，只保留历史时刻的信息，不涉及当前时刻的状态
        // 这里需要根据双差测量的定义，计算出理论值和观测值之间的差异
        // 具体的计算方法会依赖于双差测量的定义和使用的观测类型（如伪距、载波相位等）
        Eigen::Vector3d now_pos(state_array[measSize-1][0], state_array[measSize-1][1], state_array[measSize-1][2]);

        const auto* curr_master_sv_ptr = getSatelliteInfoPtr(sv_info_map, dd_measurement.u_master_SV);
        const auto* curr_i_sv_ptr = getSatelliteInfoPtr(sv_info_map, dd_measurement.u_iSV);
        const auto reference_sv_info_it = sv_info_window_map.find(dd_measurement.prev_time);
        if (!curr_master_sv_ptr || !curr_i_sv_ptr || reference_sv_info_it == sv_info_window_map.end())
        {
            return 100.0;
        }

        const auto& reference_sv_info_map = reference_sv_info_it->second;
        const auto* prev_master_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_master_SV);
        const auto* prev_i_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_iSV);
        if (!prev_master_sv_ptr || !prev_i_sv_ptr)
        {
            return 100.0;
        }

        const sv_info& curr_master_sv = *curr_master_sv_ptr;
        const sv_info& curr_i_sv = *curr_i_sv_ptr;
        const sv_info& prev_master_sv = *prev_master_sv_ptr;
        const sv_info& prev_i_sv = *prev_i_sv_ptr;


        Eigen::Vector3d curr_master_pos(curr_master_sv.pos[0], curr_master_sv.pos[1], curr_master_sv.pos[2]);
        Eigen::Vector3d curr_i_pos(curr_i_sv.pos[0], curr_i_sv.pos[1], curr_i_sv.pos[2]);
        Eigen::Vector3d prev_master_pos(prev_master_sv.pos[0], prev_master_sv.pos[1], prev_master_sv.pos[2]);
        Eigen::Vector3d prev_i_pos(prev_i_sv.pos[0], prev_i_sv.pos[1], prev_i_sv.pos[2]);

        double est_prev_master = sqrt((prev_pos.x() - prev_master_pos.x()) * (prev_pos.x() - prev_master_pos.x()) +
                                 (prev_pos.y() - prev_master_pos.y()) * (prev_pos.y() - prev_master_pos.y()) +
                                 (prev_pos.z() - prev_master_pos.z()) * (prev_pos.z() - prev_master_pos.z()));
        est_prev_master = est_prev_master + (OMGE_ / CLIGHT_) *
            (prev_master_pos.x() * prev_pos[1] - (prev_master_pos.y()) * prev_pos[0]);

        double est_prev_i = sqrt((prev_pos.x() - prev_i_pos.x()) * (prev_pos.x() - prev_i_pos.x()) +
                            (prev_pos.y() - prev_i_pos.y()) * (prev_pos.y() - prev_i_pos.y()) +
                            (prev_pos.z() - prev_i_pos.z()) * (prev_pos.z() - prev_i_pos.z()));
        est_prev_i = est_prev_i + (OMGE_ / CLIGHT_) *
            (prev_i_pos.x() * prev_pos[1] - prev_i_pos.y() * prev_pos[0]);

        double est_curr_master = sqrt((now_pos.x() - curr_master_pos.x()) * (now_pos.x() - curr_master_pos.x()) +
                                 (now_pos.y() - curr_master_pos.y()) * (now_pos.y() - curr_master_pos.y()) +
                                 (now_pos.z() - curr_master_pos.z()) * (now_pos.z() - curr_master_pos.z()));
        est_curr_master = est_curr_master + (OMGE_ / CLIGHT_) *
            (curr_master_pos.x() * now_pos[1] - curr_master_pos.y() * now_pos[0]);

        double est_curr_i = sqrt((now_pos.x() - curr_i_pos.x()) * (now_pos.x() - curr_i_pos.x()) +
                            (now_pos.y() - curr_i_pos.y()) * (now_pos.y() - curr_i_pos.y()) +
                            (now_pos.z() - curr_i_pos.z()) * (now_pos.z() - curr_i_pos.z()));
        est_curr_i = est_curr_i + (OMGE_ / CLIGHT_) *
            (curr_i_pos.x() * now_pos[1] - curr_i_pos.y() * now_pos[0]);

        // 双差载波相位-估计值：
        double est_dd_cp = (est_curr_i - est_prev_i) - (est_curr_master - est_prev_master);

        int l1_idx_prev_master = -1;
        int l1_idx_prev_i = -1;
        int l1_idx_curr_master = -1;
        int l1_idx_curr_i = -1;
        L1_freq(dd_measurement.r_master_SV, &l1_idx_prev_master);
        L1_freq(dd_measurement.r_iSV, &l1_idx_prev_i);
        L1_freq(dd_measurement.u_master_SV, &l1_idx_curr_master);
        L1_freq(dd_measurement.u_iSV, &l1_idx_curr_i);
        if(l1_idx_prev_master < 0 || l1_idx_prev_i < 0 || l1_idx_curr_master < 0 || l1_idx_curr_i < 0)
        {
            // 无法找到L1频点的观测数据，无法计算载波相位残差，返回一个较大的残差值以降低该因子的权重
            return 100.0; 
        }

        double prev_master_cp = dd_measurement.r_master_SV->cp[l1_idx_prev_master] * prev_master_sv.lamda;
        double prev_i_cp = dd_measurement.r_iSV->cp[l1_idx_prev_i] * prev_i_sv.lamda;
        double curr_master_cp = dd_measurement.u_master_SV->cp[l1_idx_curr_master] * curr_master_sv.lamda;
        double curr_i_cp = dd_measurement.u_iSV->cp[l1_idx_curr_i] * curr_i_sv.lamda;
        // 双差载波相位-观测值
        double dd_cp = (curr_i_cp - prev_i_cp) - (curr_master_cp - prev_master_cp);
        double residuals_raw = (est_dd_cp - dd_cp);
        return residuals_raw; 
    }

    void applyEpochSnapshot(const EpochBucketSnapshot& snapshot)
    {
        for (const auto& contribution : snapshot.contributions)
        {
            bucket_manager.update(contribution.key, contribution.residual, contribution.success);
            bucket_manager.updateTimeDiscount(snapshot.epoch_index,
                                              contribution.gap_epoch,
                                              contribution.residual,
                                              contribution.delta_pos);
        }
    }

    void rollbackEpochSnapshot(const EpochBucketSnapshot& snapshot)
    {
        for (const auto& contribution : snapshot.contributions)
        {
            bucket_manager.remove(contribution.key, contribution.residual, contribution.success);
        }
        bucket_manager.removeTimeDiscountEpoch(snapshot.epoch_index);
    }

    //update history info 一轮优化后使用
    bool updateHistoryInfo()
    {
        if (measSize <= 0)
        {
            return true;
        }

        const int latest_epoch_index = measSize - 1;
        const int earliest_epoch_index = (recent_history_epoch_window <= 0)
            ? 0
            : std::max(0, latest_epoch_index - std::max(1, recent_history_epoch_window) + 1);

        auto existing_snapshot_it = epoch_bucket_history.find(latest_epoch_index);
        if (existing_snapshot_it != epoch_bucket_history.end())
        {
            rollbackEpochSnapshot(existing_snapshot_it->second);
            epoch_bucket_history.erase(existing_snapshot_it);
        }

        EpochBucketSnapshot current_snapshot;
        current_snapshot.epoch_index = latest_epoch_index;

        for (auto& iter : tr_measurments)
        {
            if (iter.curr_epoch_index != latest_epoch_index)
            {
                continue;
            }

            const auto& tr_factor = iter.tr_feature;
            if (iter.curr_epoch_index < earliest_epoch_index || iter.prev_epoch_index < earliest_epoch_index)
            {
                continue; // 只统计最近 N 个历元内的卫星特征数据和残差
            }

            int gap = tr_factor.gap_epoch;
            if(gap >= recent_history_epoch_window)
            {
                continue; // 超出历史窗口范围，不更新
            }

            if(measSize-gap-1 >= 0)
            {
                Eigen::Vector3d prev_pos(state_array[measSize-gap-1][0], state_array[measSize-gap-1][1], state_array[measSize-gap-1][2]);
                Eigen::Vector3d curr_pos(state_array[measSize-1][0], state_array[measSize-1][1], state_array[measSize-1][2]);
                double delta_pos = (curr_pos - prev_pos).norm();
                double residual = computeTRResidual(iter,prev_pos);

                BucketKey key = indexer.makeKey(
                tr_factor.min_SNR,
                tr_factor.max_SNR,
                tr_factor.min_elev,
                tr_factor.max_elev,
                tr_factor.delta_azm);

                BucketContributionRecord contribution;
                contribution.key = key;
                contribution.residual = residual;
                contribution.gap_epoch = tr_factor.gap_epoch;
                contribution.delta_pos = delta_pos;
                contribution.success = residual < 0.05;
                current_snapshot.contributions.push_back(contribution);
            }
            
            
        }

        applyEpochSnapshot(current_snapshot);
        epoch_bucket_history.emplace(latest_epoch_index, std::move(current_snapshot));

        if (recent_history_epoch_window > 0)
        {
            while (!epoch_bucket_history.empty())
            {
                const auto oldest_it = epoch_bucket_history.begin();
                if (oldest_it->first >= earliest_epoch_index)
                {
                    break;
                }

                rollbackEpochSnapshot(oldest_it->second);
                epoch_bucket_history.erase(oldest_it);
            }
        }

        return true;
    }

    bool addTRDDCPFactors()
    {
        tr_factor_features.clear();
        tr_measurments.clear();
        TRFactorCount = 0;

        if(measSize < 2)
        {
            return false; // 需要至少两个历元的数据才能构建TR因子
        }

        auto iter = gnss_raw_map.end();
        --iter;
        // ROS_INFO("CURRENT TIME: %f  ---- LATEST ONS TIME: %f", iter->first, time_frame_now);
        // 重置卫星连续观测检查标志
        for (auto& pair : sat_lock_check) 
        {
            pair.second = false;
        }
        if(iter!=gnss_raw_map.begin())
        {
            // 维护信息表
            for(int i=0;i<iter->second.size();i++)
            {
                int sat_id = iter->second[i]->sat;
                if(sat_lock_check.find(sat_id) == sat_lock_check.end() )
                {
                    sat_lock_check[sat_id] = true;
                    sat_lock_count[sat_id] = 1;
                    continue;
                }
                
                if(cycleSlipDetect(sat_id))
                {
                    // ROS_WARN("Cycle slip detected for satellite %d at time %.3f", sat_id, iter->first);
                    sat_lock_check[sat_id] = true;
                    sat_lock_count[sat_id] = 1;
                }
                else
                {
                    sat_lock_count[sat_id]++;
                    sat_lock_check[sat_id] = true;
                }
                // ROS_INFO("Satellite %d lock count: %d", sat_id, sat_lock_count[sat_id]);
            }
            //失锁判定
            for(auto& pair : sat_lock_check)
            {
                if(pair.second == false)
                {
                    sat_lock_count[pair.first] = 0;
                    // ROS_INFO("Satellite %d lock LOSE: %d", pair.first, sat_lock_count[pair.first]);
                }
            }
        }

        //1. 遍历窗口，挑选候选TR因子
        auto curr_iter = iter;
        auto gnss_data = curr_iter->second;
        for(const auto& obs : gnss_data)
        {
            if(!obs) // 连续观测不足2个历元的卫星不考虑构建TR因子，增加稳定性
            {
                continue;
            }
            else
            {
                auto it = curr_iter;
                int i = 1;
                while (i <= windowSize && it != gnss_raw_map.begin())
                {
                    // ROS_INFO("FIND TR FACTOR: Satellite %d, Epoch Gap: %d", obs->sat, i);
                    --it; // move to previous epoch; i means epoch-gap (>=1)
                    //构建双差因子：
                    TRRTKMeasurement tr_meas;
                    tr_meas.u_master_SV = obs;

                    if(sat_lock_count[obs->sat] <= i)
                    {
                        break;
                    }
                    
                    // 在连续锁定的卫星中寻找副卫星
                    for(auto pair: sat_lock_count)
                    {
                        // ROS_INFO("CHECKING :  Satellite %d lock count: %d", pair.first, pair.second);
                        if(pair.first == obs->sat)
                        {
                            continue; // 跳过主卫星
                        }

                        if (satsys(pair.first, nullptr) != satsys(obs->sat, nullptr))
                        {
                            continue; // 仅在同一星座/系统内构建TR双差，避免系统间钟差引入偏置
                        }

                        if(pair.second > (i + 1)) // 连续观测至少i+1个历元（含当前），增加TR因子稳定性
                        {
                            // printf("selected: sat_id: %d,  lock_count: %d \n", pair.first, pair.second);
                            // 找到一个满足条件的副卫星，构建TR因子
                            findSatellitewithSameId(pair.first, gnss_data, tr_meas.u_iSV);
                            findSatellitewithSameId(pair.first, it->second, tr_meas.r_iSV);
                            findSatellitewithSameId(obs->sat, it->second, tr_meas.r_master_SV);

                            tr_meas.prev_epoch_index = measSize - 1 - i;
                            tr_meas.curr_epoch_index = measSize - 1;
                            tr_meas.prev_time = it->first;
                            tr_meas.curr_time = curr_iter->first;

                            if (tr_meas.prev_epoch_index < 0 || tr_meas.prev_epoch_index >= measSize ||
                                tr_meas.curr_epoch_index < 0 || tr_meas.curr_epoch_index >= measSize ||
                                tr_meas.prev_epoch_index == tr_meas.curr_epoch_index)
                            {
                                // ROS_INFO("skip reason 1");
                                continue;
                            }

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

                            Eigen::Vector3d prev_pos(state_array[tr_meas.prev_epoch_index][0], state_array[tr_meas.prev_epoch_index][1], state_array[tr_meas.prev_epoch_index][2]);
                            Eigen::Vector3d pred_pos(state_array[tr_meas.curr_epoch_index-1][0], state_array[tr_meas.curr_epoch_index-1][1], state_array[tr_meas.curr_epoch_index-1][2]);
                            Eigen::Vector3d dop_vel ;
                            dop_vel.x() = doppler_map[tr_meas.curr_time].twist.twist.linear.x;
                            dop_vel.y() = doppler_map[tr_meas.curr_time].twist.twist.linear.y;
                            dop_vel.z() = doppler_map[tr_meas.curr_time].twist.twist.linear.z;
                            pred_pos = pred_pos + dop_vel*((time_frame_now - time_frame_last)/10.0);
                            // ROS_INFO("Predicted Move Distance for TR factor: %.2f m (Doppler Velocity: [%.2f, %.2f, %.2f] m/s, Time Gap: %.1f s)", 
                            //          (pred_pos - prev_pos).norm(), dop_vel.x(), dop_vel.y(), dop_vel.z(), (time_frame_now - time_frame_last)/10.0);
                            double pred_move = (pred_pos - prev_pos).norm();
                            tr_meas.pred_move = pred_move; // 将预测的移动距离关联到TR测量中，方便后续时间折扣函数使用


                           if(checkRedundantPair(tr_meas))
                            {
                                // ROS_INFO("skip reason 3");
                                continue; // 已存在相同卫星组合的TR因子，跳过以避免冗余
                            }

                            if(tr_meas.u_iSV && tr_meas.r_master_SV && tr_meas.r_iSV)
                            {
                                TRfeature tr_feature_cap;
                                tr_feature_cap = buildTRFactor(tr_meas);
                                double time_discount = 1.0;
                                double p_rel = 1.0;
                                if(!SAME_RELIABLE)
                                {
                                    p_rel = bucket_manager.getReliability(indexer.makeKey(
                                                                tr_feature_cap.min_SNR,
                                                                tr_feature_cap.max_SNR,
                                                                tr_feature_cap.min_elev,
                                                                tr_feature_cap.max_elev,
                                                                tr_feature_cap.delta_azm));
                                }
                                
                                if(!SAME_TIME_WEIGHT)
                                {
                                    time_discount = bucket_manager.getTimeDiscount(tr_feature_cap.gap_epoch, tr_meas.pred_move); // 这里暂时不考虑delta_pos对时间折扣的影响
                                
                                } 

                                tr_meas.tr_feature = tr_feature_cap; // 将TR因子信息关联到TR测量中，方便后续更新历史信息时使用
                                
                                tr_meas.tr_score = pow(p_rel,1.0/6.0) * time_discount ; // 综合可靠性得分，作为TR因子的权重依据

                            
                                // ROS_INFO("Found TR factor: Rel_Prob=%.3f, |Gap: %d, Time_Discount=%.3f | TR_Score=%.3f", p_rel, time_discount,tr_meas.curr_epoch_index-tr_meas.prev_epoch_index, tr_meas.tr_score);

                                if(TRFactorCount >= MaxTRFactorNum)
                                {
                                    // ROS_INFO("Reached max TR factor limit (%d), skipping remaining candidates.", MaxTRFactorNum);
                                    break; // 达到TR因子数量上限，停止添加更多的TR因子
                                }
                                tr_factor_features.push_back(tr_feature_cap); //TR因子信息列表
                                tr_measurments.push_back(tr_meas);   //候选双差因子列表
                                TRFactorCount++;
                            }
                        }
                    }
                    ++i;
                    if(TRFactorCount >= MaxTRFactorNum)
                    {
                        // ROS_INFO("Reached max TR factor limit (%d), skipping remaining candidates.", MaxTRFactorNum);
                        break; // 达到TR因子数量上限，停止添加更多的TR因子
                    }
                }
            }
            if(TRFactorCount >= MaxTRFactorNum)
            {
                // ROS_INFO("Reached max TR factor limit (%d), skipping remaining candidates.", MaxTRFactorNum);
                break; // 达到TR因子数量上限，停止添加更多的TR因子
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
                new DDPseudorangeFactor(tr_measurments[i], current_sv_info, reference_sv_info, 1.0));
            
            problem.AddResidualBlock(dd_pr_function, loss_function, state_array[prev_epoch_index], state_array[epoch_index]);
            
            
            ceres::CostFunction* dd_cp_function =
                new ceres::AutoDiffCostFunction<TRDDCPFactor, 1, state_size, state_size>(
                    new TRDDCPFactor(tr_measurments[i], current_sv_info, reference_sv_info, tr_measurments[i].tr_score));
            problem.AddResidualBlock(dd_cp_function, loss_function, state_array[prev_epoch_index], state_array[epoch_index]);
        }
        int cnt_gps=0, cnt_glo=0, cnt_gal=0, cnt_bds=0;
        int cnt_l1=0, cnt_l2=0;
        for (const auto& tr_m : tr_measurments) {
            if (!tr_m.u_master_SV) continue;
            const int sys = satsys(tr_m.u_master_SV->sat, nullptr);
            if (sys==SYS_GPS) cnt_gps++; else if (sys==SYS_GLO) cnt_glo++;
            else if (sys==SYS_GAL) cnt_gal++; else if (sys==SYS_BDS) cnt_bds++;
            // Check freq availability for master SV
            int l1=-1,l2=-1; L1_freq(tr_m.u_master_SV,&l1); L2_freq(tr_m.u_master_SV,&l2);
            if (l1>=0) cnt_l1++; if (l2>=0) cnt_l2++;
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

    bool checkRedundantPair(TRRTKMeasurement tr_meas)
    {
        for(int i=0;i<tr_measurments.size();i++)
        {
            if(tr_meas.r_master_SV->sat == tr_measurments[i].r_master_SV->sat 
                && tr_meas.r_iSV->sat == tr_measurments[i].r_iSV->sat 
                && tr_meas.curr_epoch_index == tr_measurments[i].curr_epoch_index 
                && tr_meas.prev_epoch_index == tr_measurments[i].prev_epoch_index)
            {
                // ROS_INFO("Redundant pair found: Master SV %d, iSV %d, Epoch gap %d", tr_meas.r_master_SV->sat, tr_meas.r_iSV->sat, tr_meas.curr_epoch_index - tr_meas.prev_epoch_index);
                return true;
            }
            else if(tr_meas.r_master_SV->sat == tr_measurments[i].r_iSV->sat 
                && tr_meas.r_iSV->sat == tr_measurments[i].r_master_SV->sat
                && tr_meas.curr_epoch_index == tr_measurments[i].curr_epoch_index
                && tr_meas.prev_epoch_index == tr_measurments[i].prev_epoch_index)
            {
                // ROS_INFO("Redundant pair found (reversed): Master SV %d, iSV %d, Epoch gap %d", tr_meas.r_master_SV->sat, tr_meas.r_iSV->sat, tr_meas.curr_epoch_index - tr_meas.prev_epoch_index);
                return true;
            }
        }
        return false;
    }

    static bool hasValidCarrierPhaseForTR(const gnss_comm::ObsPtr& obs, int freq_sel = -1)
    {
        if (!obs)
        {
            return false;
        }

        int freq_idx = -1;
        double lamda_d = 0.0, lamda_l2_d = 0.0;
        getFreqIndex(obs, freq_sel, freq_idx, lamda_d, lamda_l2_d);
        if (freq_idx < 0 ||
            freq_idx >= static_cast<int>(obs->cp.size()) ||
            freq_idx >= static_cast<int>(obs->LLI.size()) ||
            freq_idx >= static_cast<int>(obs->status.size()))
        {
            return false;
        }

        if (obs->cp[freq_idx] < 10.0)
        {
            return false;
        }

        const bool carrier_phase_valid = (obs->status[freq_idx] & 0x2u) != 0u;
        const bool lost_lock = obs->LLI[freq_idx] != 0u;
        return carrier_phase_valid && !lost_lock;
    }

    // 复用WCP前置筛选思想：仅当载波相位预拟合残差在合理范围时才允许通过。
    bool hasReasonableWcpPrefitResidual(const Eigen::Matrix<double, 7, 1>& state,
                                        const gnss_comm::ObsPtr& obs,
                                        const sv_info& info,
                                        double max_abs_residual_m,
                                        int freq_sel = -1)
    {
        if (!obs || info.lamda <= 0.0 || !std::isfinite(max_abs_residual_m) || max_abs_residual_m <= 0.0)
        {
            return false;
        }

        int freq_idx = -1;
        double lamda_d = 0.0, lamda_l2_d = 0.0;
        getFreqIndex(obs, freq_sel, freq_idx, lamda_d, lamda_l2_d);
        if (freq_idx < 0 || freq_idx >= static_cast<int>(obs->cp.size()) || obs->cp[freq_idx] < 10.0)
        {
            return false;
        }

        const Eigen::Vector3d sat_pos(info.pos[0], info.pos[1], info.pos[2]);
        if (!std::isfinite(state(0)) || !std::isfinite(state(1)) || !std::isfinite(state(2)) ||
            !std::isfinite(state(3)) || !std::isfinite(state(4)) ||
            !std::isfinite(state(5)) || !std::isfinite(state(6)) ||
            !std::isfinite(sat_pos.x()) || !std::isfinite(sat_pos.y()) || !std::isfinite(sat_pos.z()) ||
            !std::isfinite(info.dt))
        {
            return false;
        }

        const double dx = state(0) - sat_pos.x();
        const double dy = state(1) - sat_pos.y();
        const double dz = state(2) - sat_pos.z();
        double rho = std::sqrt(dx * dx + dy * dy + dz * dz);
        rho += (OMGE_ / CLIGHT_) * (sat_pos.x() * state(1) - sat_pos.y() * state(0));

        const int sat_sys = satsys(static_cast<int>(obs->sat), nullptr);
        const double recv_clk = (sat_sys == SYS_GPS) ? state(3) :
                                 ((sat_sys == SYS_GLO) ? state(4) :
                                 ((sat_sys == SYS_GAL) ? state(5) :
                                 ((sat_sys == SYS_BDS) ? state(6) : 0.0)));
        const double lambda_m = getWavelengthForFreq(&info, freq_sel);
        const double model_m = rho + recv_clk - CLIGHT_ * info.dt;
        const double meas_m = obs->cp[freq_idx] * lambda_m;
        const double residual_m = model_m - meas_m;

        auto baseline_it = sat_const_residual.find(obs->sat);
        if (baseline_it == sat_const_residual.end() || std::abs(sat_const_residual[obs->sat]) > 100.0 ||
            (std::abs(residual_m) - std::abs(sat_const_residual[obs->sat])) > 1.0)
        {
            sat_const_residual[obs->sat] = residual_m;
            return true;
        }

        const double baseline_residual_m = baseline_it->second;
        const double residual_delta_m = std::abs(residual_m) - std::abs(baseline_residual_m);
        return std::isfinite(residual_delta_m) && std::abs(residual_m) <= max_abs_residual_m;
    }

    bool cycleSlipDetect(int sat_id)
    {
        auto mark_slip = [this, sat_id]() {
            sat_lock_count[sat_id] = 0;
            sat_lock_check[sat_id] = true;
            sat_cp_const.erase(sat_id);
            sat_const_residual.erase(sat_id);
        };

        if (gnss_raw_map.size() < 2 || measSize < 2)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 1");
            return true;
        }

        const auto iter_now = std::prev(gnss_raw_map.end());
        const auto iter_prev = std::prev(iter_now);

        gnss_comm::ObsPtr prev_obs;
        gnss_comm::ObsPtr curr_obs;
        findSatellitewithSameId(sat_id, iter_prev->second, prev_obs);
        findSatellitewithSameId(sat_id, iter_now->second, curr_obs);
        if (!prev_obs || !curr_obs)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 2");
            return true;
        }

        int prev_l1_idx = -1;
        int curr_l1_idx = -1;
        double prev_cp_cycle = 0.0;
        double curr_cp_cycle = 0.0;
        if (!getValidL1CarrierPhase(prev_obs, prev_l1_idx, prev_cp_cycle) ||
            !getValidL1CarrierPhase(curr_obs, curr_l1_idx, curr_cp_cycle) ||
            !hasValidCarrierPhaseForTR(prev_obs) ||
            !hasValidCarrierPhaseForTR(curr_obs))
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 3");
            return true;
        }

        const double prev_time_sec = prev_obs->time.time + prev_obs->time.sec;
        const double curr_time_sec = curr_obs->time.time + curr_obs->time.sec;
        const double dt = curr_time_sec - prev_time_sec;
        if (!std::isfinite(dt) || dt <= 0.0)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 5");
            return true;
        }

        // ROS_INFO("CURRENT TIME: %f  ---- PREV TIME: %f", iter_now->first, iter_prev->first);
        const auto prev_sv_it = sv_info_window_map.find(iter_prev->first);
        const auto curr_sv_it = sv_info_window_map.find(iter_now->first);
        if (prev_sv_it == sv_info_window_map.end() || curr_sv_it == sv_info_window_map.end())
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 6");
            return true;
        }

        const auto prev_sat_it = prev_sv_it->second.find(sat_id);
        const auto curr_sat_it = curr_sv_it->second.find(sat_id);
        if (prev_sat_it == prev_sv_it->second.end() || curr_sat_it == curr_sv_it->second.end())
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 7");
            return true;
        }


        const sv_info& prev_sv = prev_sat_it->second;
        const sv_info& curr_sv = curr_sat_it->second;


        const double lambda = (curr_sv.lamda > 0.0) ? curr_sv.lamda : prev_sv.lamda;
        const double prev_lambda = (prev_sv.lamda > 0.0) ? prev_sv.lamda : lambda;
        if (!std::isfinite(lambda) || lambda <= 0.0 || !std::isfinite(prev_lambda) || prev_lambda <= 0.0 || prev_lambda != lambda)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 8");
            return true;
        }

        //更新sat_cp_const以适应卫星状态的变化，避免过时的常数导致误判周跳
        sat_cp_const[sat_id] = curr_cp_cycle - prev_cp_cycle;
        // printf("sat_id %d: const: %f cycles. | ", sat_id, sat_cp_const[sat_id]);
        
        const Eigen::Vector3d prev_vel{doppler_map[iter_prev->first].twist.twist.linear.x,
                                    doppler_map[iter_prev->first].twist.twist.linear.y,
                                    doppler_map[iter_prev->first].twist.twist.linear.z};
        const Eigen::Vector3d curr_vel{doppler_map[iter_now->first].twist.twist.linear.x,
                                    doppler_map[iter_now->first].twist.twist.linear.y,
                                    doppler_map[iter_now->first].twist.twist.linear.z};
                

        const Eigen::Vector3d prev_rcv(state_array[measSize - 2][0], state_array[measSize - 2][1], state_array[measSize - 2][2]);
        Eigen::Vector3d curr_rcv;
        curr_rcv.x() = prev_rcv.x() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[iter_now->first].twist.twist.linear.x;
        curr_rcv.y() = prev_rcv.y() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[iter_now->first].twist.twist.linear.y;
        curr_rcv.z() = prev_rcv.z() + ((curr_time_sec - prev_time_sec) / 10.0) * doppler_map[iter_now->first].twist.twist.linear.z;
        const Eigen::Vector3d prev_sat(prev_sv.pos[0], prev_sv.pos[1], prev_sv.pos[2]);
        const Eigen::Vector3d curr_sat(curr_sv.pos[0], curr_sv.pos[1], curr_sv.pos[2]);
        // ROS_INFO("Prev RCV: [%.3f, %.3f, %.3f], Curr RCV: [%.3f, %.3f, %.3f]", prev_rcv.x(), prev_rcv.y(), prev_rcv.z(), curr_rcv.x(), curr_rcv.y(), curr_rcv.z());

        const Eigen::Vector3d prev_los = prev_sat - prev_rcv;
        const Eigen::Vector3d curr_los = curr_sat - curr_rcv;
        if (!std::isfinite(prev_los.norm()) || !std::isfinite(curr_los.norm()) ||
            prev_los.norm() <= 0.0 || curr_los.norm() <= 0.0)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 9");
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

        const double delta_cp_pred = ((curr_vel + prev_vel).dot(los) * (curr_time_sec - prev_time_sec) / 10.0)/(2 * lambda);
        const double delta_cp_obs = (curr_cp_cycle - prev_cp_cycle) - sat_cp_const[sat_id];
        // printf(" | obs delta cp: %f \n",(curr_cp_cycle - prev_cp_cycle));
        
        constexpr double kCycleSlipThresholdCycle = 1.0; // 预测和观测的载波相位变化超过0.5周期则判定为可能发生了周跳
        if(std::abs(delta_cp_obs - delta_cp_pred) > kCycleSlipThresholdCycle)
        {
            mark_slip();
            // ROS_INFO("SLIP REASON 10");
            ROS_INFO("Satellite %d: Delta CP Pred=%.3f cycles, Delta CP Obs=%.3f cycles", sat_id, delta_cp_pred, delta_cp_obs);
            return true;
        }

        return false;
    }


    TRfeature buildTRFactor(const TRRTKMeasurement& tr_meas)
    {
        // 计算TR因子的信噪比特征，可以根据主卫星和副卫星的信噪比进行计算
        int prev_epoch_index = tr_meas.prev_epoch_index;
        double prev_time = tr_meas.prev_time;
        auto prev_obs = gnss_raw_map[prev_time];
        const auto prev_sv_info_it = sv_info_window_map.find(prev_time);
        if (prev_sv_info_it == sv_info_window_map.end())
        {
            return TRfeature{};
        }
        const auto& prev_sv_info = prev_sv_info_it->second;

        double snr_m = 30.0, snr_i = 30.0;
        double elev_m = 15.0, elev_i = 15.0;
        Eigen::Vector3d sat_m_pos, sat_i_pos;
        
        for(const auto& obs : prev_obs)
        {
            if(obs->sat == tr_meas.r_master_SV->sat)
            {
                int l1_idx = -1;
                L1_freq(obs, &l1_idx);
                if(l1_idx >= 0)
                {
                    snr_m = obs->CN0[l1_idx];
                    const auto it_sv = prev_sv_info.find(obs->sat);
                    if (it_sv != prev_sv_info.end())
                    {
                        elev_m = it_sv->second.elevation * R2D;
                        sat_m_pos = it_sv->second.pos;
                    }
                }
            }
            else if(obs->sat == tr_meas.r_iSV->sat)
            {
                int l1_idx = -1;
                L1_freq(obs, &l1_idx);
                if(l1_idx >= 0)
                {
                    snr_i = obs->CN0[l1_idx];
                    const auto it_sv = prev_sv_info.find(obs->sat);
                    if (it_sv != prev_sv_info.end())
                    {
                        elev_i = it_sv->second.elevation * R2D;
                        sat_i_pos = it_sv->second.pos;
                    }
                }
            }
        }

        Eigen::Vector3d now_pos{state_array[measSize-2][0],state_array[measSize-2][1],state_array[measSize-2][2]};
        Eigen::Vector3d los_e_m = (sat_m_pos - now_pos).normalized();
        Eigen::Vector3d los_e_i = (sat_i_pos - now_pos).normalized();
        Eigen::Vector3d zenith = now_pos.normalized();
        Eigen::Vector3d horiz_m = los_e_m - (los_e_m.dot(zenith)) * zenith;
        Eigen::Vector3d horiz_i = los_e_i - (los_e_i.dot(zenith)) * zenith;
        double azm_m = (horiz_m-horiz_i).norm();


        TRfeature tr_factor;
        tr_factor.min_SNR = std::min(snr_m, snr_i); // 取主副卫星中较小的信噪比作为TR因子的特征
        tr_factor.max_SNR = std::max(snr_m, snr_i); // 主副卫星信噪比差值作为另一个特征
        tr_factor.min_elev = std::min(elev_m, elev_i); //
        tr_factor.max_elev = std::max(elev_m, elev_i); // 计算地理仰角特征
        tr_factor.delta_elev = std::abs(elev_m - elev_i);
        tr_factor.delta_azm = azm_m;
        tr_factor.gap_epoch = std::max(0, tr_meas.curr_epoch_index - tr_meas.prev_epoch_index); // 计算历元间隔特征
        

        return tr_factor; // 取主副卫星中较小的信噪比作为TR因子的特征
    }

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


