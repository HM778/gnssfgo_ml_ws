#ifndef CHECKING_HPP
#define CHECKING_HPP

#include <cmath>
#include <limits>
#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_constant.hpp>
#include "../datatype.h"

// 检查观测数据中L1载波相位是否有效
bool getValidL1CarrierPhase(const gnss_comm::ObsPtr& obs, int& l1_idx, double& cp_cycle)
    {
        l1_idx = -1;
        cp_cycle = 0.0;
        if (!obs)
        {
            return false;
        }
        L1_freq(obs, &l1_idx);
        const auto is_valid_cp = [](double v) {
            // Keep sign (TDCP may be positive or negative), reject only invalid/near-zero placeholders.
            return std::isfinite(v) && std::abs(v) >= 1.0e-3;
        };

        if (l1_idx >= 0 && l1_idx < static_cast<int>(obs->cp.size()))
        {
            cp_cycle = obs->cp[l1_idx];
            if (is_valid_cp(cp_cycle))
            {
                // ROS_INFO("Found valid L1 carrier phase: satellite %d, cp=%.3f cycles", obs->sat, cp_cycle);
                return true;
            }
        }

        // // Fallback: if L1 index is unavailable/invalid, try any valid carrier-phase entry.
        // for (int i = 0; i < static_cast<int>(obs->cp.size()); ++i)
        // {
        //     if (is_valid_cp(obs->cp[i]))
        //     {
        //         l1_idx = i;
        //         cp_cycle = obs->cp[i];
        //         return true;
        //     }
        // }

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

bool hasLLISlip(const gnss_comm::ObsPtr& obs, int l1_idx) 
    {
        if (!obs || l1_idx < 0 || l1_idx >= static_cast<int>(obs->LLI.size()))
        {
            return false;
        }
        const int lli = static_cast<int>(obs->LLI[l1_idx]);
        // 低两位通常用于表示周跳/半周等异常状态，保守处理为“疑似失锁”
        return (lli & 0x03) != 0;
    }

bool detectCycleSlipOrLoss(const gnss_comm::ObsPtr& prev_obs, const gnss_comm::ObsPtr& curr_obs) 
{
    int prev_l1 = -1;
    int curr_l1 = -1;
    double prev_cp = 0.0;
    double curr_cp = 0.0;
    if (!getValidL1CarrierPhase(prev_obs, prev_l1, prev_cp) || !getValidL1CarrierPhase(curr_obs, curr_l1, curr_cp))
    {
        return true;
    }

    if (hasLLISlip(prev_obs, prev_l1) || hasLLISlip(curr_obs, curr_l1))
    {
        return true;
    }

    // 极端跳变保护（阈值取大，避免误判正常动态）
    const double cp_jump = std::abs(curr_cp - prev_cp);
    if (!std::isfinite(cp_jump) || cp_jump > 1.0e5)
    {
        return true;
    }
    return false;
}

inline double obsTimeSeconds(const gnss_comm::ObsPtr& obs)
{
    if (!obs)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return obs->time.time + obs->time.sec;
}

inline bool isLikelyAdjacentEpochPair(const gnss_comm::ObsPtr& prev_obs,
                                      const gnss_comm::ObsPtr& curr_obs,
                                      double min_dt_sec = 0.02,
                                      double max_dt_sec = 2.0)
{
    const double t_prev = obsTimeSeconds(prev_obs);
    const double t_curr = obsTimeSeconds(curr_obs);
    if (!std::isfinite(t_prev) || !std::isfinite(t_curr))
    {
        return false;
    }
    const double dt = std::abs(t_curr - t_prev);
    return (dt >= min_dt_sec) && (dt <= max_dt_sec);
}

/**
 * @brief check the wheather all the consided satellite have carrier-phase
 * @param gnss double-difference measurements
 * @return true or false
*/
bool checkCarrierPhaseConsistency(DDMeasurement DD_measurement)
{
    int l1_idx_u_m = -1;
    int l1_idx_u_i = -1;
    int l1_idx_r_m = -1;
    int l1_idx_r_i = -1;
    double cp_u_m = 0.0;
    double cp_u_i = 0.0;
    double cp_r_m = 0.0;
    double cp_r_i = 0.0;

    if (!getValidL1CarrierPhase(DD_measurement.u_master_SV, l1_idx_u_m, cp_u_m) ||
        !getValidL1CarrierPhase(DD_measurement.u_iSV, l1_idx_u_i, cp_u_i) ||
        !getValidL1CarrierPhase(DD_measurement.r_master_SV, l1_idx_r_m, cp_r_m) ||
        !getValidL1CarrierPhase(DD_measurement.r_iSV, l1_idx_r_i, cp_r_i))
    {
        return false;
    }

    // Reject obvious loss-of-lock flags.
    if (hasLLISlip(DD_measurement.u_master_SV, l1_idx_u_m) ||
        hasLLISlip(DD_measurement.u_iSV, l1_idx_u_i) ||
        hasLLISlip(DD_measurement.r_master_SV, l1_idx_r_m) ||
        hasLLISlip(DD_measurement.r_iSV, l1_idx_r_i))
    {
        return false;
    }

    return true;
}

/**
 * @brief check the wheather all the consided satellite have carrier-phase 
 * or have cycle slip?
 * @param gnss double-difference measurements
 * @return true or false
*/
bool checkCarrierPhaseConsistencyCycleSlip(DDMeasurement DD_measurement)
{
    if (!checkCarrierPhaseConsistency(DD_measurement))
    {
        return false;
    }

    // For TDCP-style DD construction (same receiver, adjacent epochs),
    // additionally reject satellites with cycle slip/loss between the two epochs.
    // 相邻对检查
    const bool master_is_adjacent_epoch_pair =
        isLikelyAdjacentEpochPair(DD_measurement.r_master_SV, DD_measurement.u_master_SV);
    const bool i_is_adjacent_epoch_pair =
        isLikelyAdjacentEpochPair(DD_measurement.r_iSV, DD_measurement.u_iSV);

    if (master_is_adjacent_epoch_pair &&
        detectCycleSlipOrLoss(DD_measurement.r_master_SV, DD_measurement.u_master_SV))
    {
        return false;
    }
    if (i_is_adjacent_epoch_pair &&
        detectCycleSlipOrLoss(DD_measurement.r_iSV, DD_measurement.u_iSV))
    {
        return false;
    }

    return true;
}

// 载波相位失锁性检查
bool hasValidTRRTKCarrierPhase(const gnss_comm::ObsPtr& obs)
{
    int l1_idx = -1;
    double cp_cycle = 0.0;
    return getValidL1CarrierPhase(obs, l1_idx, cp_cycle);
}



	// 检查观测数据中L2载波相位是否有效
	inline bool getValidL2CarrierPhase(const gnss_comm::ObsPtr& obs, int& l2_idx, double& cp_cycle)
	{
	    l2_idx = -1;
	    cp_cycle = 0.0;
	    if (!obs)
	    {
	        return false;
	    }
	    gnss_comm::L2_freq(obs, &l2_idx);
	    const auto is_valid_cp = [](double v) {
	        return std::isfinite(v) && std::abs(v) >= 1.0e-3;
	    };

	    if (l2_idx >= 0 && l2_idx < static_cast<int>(obs->cp.size()))
	    {
	        cp_cycle = obs->cp[l2_idx];
	        if (is_valid_cp(cp_cycle))
	        {
	            return true;
	        }
	    }

	    return false;
	}

	// 检查某颗卫星是否同时具备L1和L2观测
	inline bool hasDualFreqObs(const gnss_comm::ObsPtr& obs)
	{
	    return gnss_comm::has_dual_freq(obs);
	}

	// freq_idx: -1 = L1 (default), -2 = L2, >= 0 = direct index
	inline double getFreqIndex(const gnss_comm::ObsPtr& obs, int freq_sel, int& idx, double& lamda, double& lamda_l2, const sv_info* info = nullptr)
	{
	    idx = -1;
	    lamda = 0.0;
	    lamda_l2 = 0.0;
	    double freq = -1.0;
	    if (freq_sel == -1)
	    {
	        freq = gnss_comm::L1_freq(obs, &idx);
	        if (info) lamda = info->lamda;
	    }
	    else if (freq_sel == -2)
	    {
	        freq = gnss_comm::L2_freq(obs, &idx);
	        if (info) lamda = info->lamda_l2;
	    }
	    else if (freq_sel >= 0)
	    {
	        idx = freq_sel;
	        if (obs && idx < static_cast<int>(obs->freqs.size()))
	            freq = obs->freqs[idx];
	        if (info)
	        {
	            lamda = (std::abs(info->freq_l1 - freq) < std::abs(info->freq_l2 - freq))
	                     ? info->lamda : info->lamda_l2;
	        }
	    }
	    return freq;
	}

	inline double getWavelengthForFreq(const sv_info* info, int freq_sel)
	{
	    if (!info) return 0.0;
	    return (freq_sel == -2 && info->lamda_l2 > 0.0) ? info->lamda_l2 : info->lamda;
	}

		// IF pseudorange: P_IF = (f1^2*P1 - f2^2*P2) / (f1^2 - f2^2). Returns 0 if dual-freq unavailable.
		inline double getIFPseudorange(const gnss_comm::ObsPtr& obs)
		{
		    if (!obs) return 0.0;
		    int l1 = -1, l2 = -1;
		    const double f1 = gnss_comm::L1_freq(obs, &l1);
		    const double f2 = gnss_comm::L2_freq(obs, &l2);
		    if (l1 < 0 || l2 < 0 || l1 >= (int)obs->psr.size() || l2 >= (int)obs->psr.size()) return 0.0;
		    if (!std::isfinite(obs->psr[l1]) || !std::isfinite(obs->psr[l2])) return 0.0;
		    if (obs->psr[l1] <= 0.0 || obs->psr[l2] <= 0.0) return 0.0;
		    return gnss_comm::psr_iono_free(f1, f2, obs->psr[l1], obs->psr[l2]);
		}

		// IF carrier phase (cycles): CP_IF = (f1^2*CP1 - f2^2*CP2) / (f1^2 - f2^2). Returns 0 if dual-freq unavailable.
		inline double getIFCarrierPhase(const gnss_comm::ObsPtr& obs)
		{
		    if (!obs) return 0.0;
		    int l1 = -1, l2 = -1;
		    const double f1 = gnss_comm::L1_freq(obs, &l1);
		    const double f2 = gnss_comm::L2_freq(obs, &l2);
		    if (l1 < 0 || l2 < 0 || l1 >= (int)obs->cp.size() || l2 >= (int)obs->cp.size()) return 0.0;
		    if (!std::isfinite(obs->cp[l1]) || !std::isfinite(obs->cp[l2])) return 0.0;
		    return gnss_comm::cp_iono_free(f1, f2, obs->cp[l1], obs->cp[l2]);
		}

		// IF pseudorange std: sigma_IF = noise_factor * sigma_L1
		inline double getIFPsrStd(const gnss_comm::ObsPtr& obs)
		{
		    if (!obs) return 3.0;
		    int l1 = -1, l2 = -1;
		    const double f1 = gnss_comm::L1_freq(obs, &l1);
		    const double f2 = gnss_comm::L2_freq(obs, &l2);
		    const double base_std = (l1 >= 0 && l1 < (int)obs->psr_std.size() &&
		                             std::isfinite(obs->psr_std[l1]) && obs->psr_std[l1] > 0.0)
		                             ? obs->psr_std[l1] : 3.0;
		    const double nf = (l2 >= 0 && f2 > 0.0 && f1 > 0.0 && std::abs(f1*f1 - f2*f2) > 1e-6)
		                       ? gnss_comm::if_noise_factor(f1, f2) : 1.0;
		    return base_std * nf;
		}

#endif // CHECKING_HPP
