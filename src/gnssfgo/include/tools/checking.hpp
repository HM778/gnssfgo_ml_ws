#ifndef CHECKING_HPP
#define CHECKING_HPP

#include <cmath>
#include <limits>
#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_constant.hpp>
#include "../datatype.h"

// 检查观测数据中载波相位是否有效
bool getValidCarrierPhase(const gnss_comm::ObsPtr& obs,int& freq_idx, double& cp_cycle, int freq)
{
    if (!obs)
    {
        return false;
    }

    const auto is_valid_cp = [](double v) {
        // Keep sign (TDCP may be positive or negative), reject only invalid/near-zero placeholders.
        return std::isfinite(v) && std::abs(v) >= 1.0e-3;
    };

    freq_idx = -1;
    cp_cycle = 0.0;

    if(freq == 1)
    {
        L1_freq(obs, &freq_idx);
    }
    else if(freq == 2)
    {
        L2_freq(obs, &freq_idx);
    }
    
    if (freq_idx >= 0 && freq_idx < static_cast<int>(obs->cp.size()))
    {
        cp_cycle = obs->cp[freq_idx];
        if(obs->cp[freq_idx] < 10.0)
        {
            // 异常值
            return false;
        }

        if((obs->status[freq_idx] & 0x2u != 0u) || (obs->LLI[freq_idx] != 0u))
        {
            // 信号显式的异常指标检查
            return false;
        }
      
        if (is_valid_cp(cp_cycle))
        {
            // ROS_INFO("Found valid L2 carrier phase: satellite %d, cp=%.3f cycles", obs->sat, cp_cycle);
            return true;
        }
    }
    
    return false;
}

// 确认观测数据中是否包含对应卫星的信息，确保双差观测因子中所需的卫星信息完整，避免因缺失卫星信息导致的优化问题
bool hasSatelliteInfo(const std::map<int, sv_info>& local_sv_info_map, const gnss_comm::ObsPtr& obs)
{
    return obs && (local_sv_info_map.find(static_cast<int>(obs->sat)) != local_sv_info_map.end());
}

inline const sv_info* getSatelliteInfoPtr(const std::map<int, sv_info>& local_sv_info_map,
                                          const gnss_comm::ObsPtr& obs)
{
    if (!obs)
    {
        return nullptr;
    }

    const auto it = local_sv_info_map.find(static_cast<int>(obs->sat));
    return (it != local_sv_info_map.end()) ? &it->second : nullptr;
}

inline const sv_info* getSatelliteInfoPtr(const std::map<int, sv_info>& local_sv_info_map, int sat_id)
{
    const auto it = local_sv_info_map.find(sat_id);
    return (it != local_sv_info_map.end()) ? &it->second : nullptr;
}

// idx用于区分L1/L2，信号本身包含的失锁信息
bool hasLLISlip(const gnss_comm::ObsPtr& obs, int idx) 
{
    if (!obs || idx < 0 || idx >= static_cast<int>(obs->LLI.size()))
    {
        return false;
    }
    int lli = static_cast<int>(obs->LLI[idx]);
    
    // 低两位通常用于表示周跳/半周等异常状态，保守处理为“疑似失锁”
    return (lli & 0x03) != 0;
}

/**
 * @brief check the wheather all the consided satellite have carrier-phase
 * @param gnss double-difference measurements
 * @return true or false
*/
bool checkCarrierPhaseConsistency(DDMeasurement DD_measurement)
{
    int check_freq = DD_measurement.freq_idx;

    int idx_u_m = -1;
    int idx_u_i = -1;
    int idx_r_m = -1;
    int idx_r_i = -1;
    double cp_u_m = 0.0;
    double cp_u_i = 0.0;
    double cp_r_m = 0.0;
    double cp_r_i = 0.0;

    if (!getValidCarrierPhase(DD_measurement.u_master_SV, idx_u_m, cp_u_m,check_freq) ||
        !getValidCarrierPhase(DD_measurement.u_iSV, idx_u_i, cp_u_i,check_freq) ||
        !getValidCarrierPhase(DD_measurement.r_master_SV, idx_r_m, cp_r_m,check_freq) ||
        !getValidCarrierPhase(DD_measurement.r_iSV, idx_r_i, cp_r_i,check_freq))
    {
        return false;
    }

    // Reject obvious loss-of-lock flags.
    if (hasLLISlip(DD_measurement.u_master_SV, idx_u_m) ||
        hasLLISlip(DD_measurement.u_iSV, idx_u_i) ||
        hasLLISlip(DD_measurement.r_master_SV, idx_r_m) ||
        hasLLISlip(DD_measurement.r_iSV, idx_r_i))
    {
        return false;
    }

    return true;
}

#endif // CHECKING_HPP
