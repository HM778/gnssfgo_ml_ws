/*******************************************************
 * Copyright (C) 2026, 

 *******************************************************/

#ifndef GNSS_COMM_EXTRA_HPP
#define GNSS_COMM_EXTRA_HPP

#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/GnssEphemMsg.h>
#include <gnss_comm/GnssGloEphemMsg.h>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_spp_extra.hpp>
#include <gnss_comm/gnss_constant.hpp>

#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include "../datatype.h"

using namespace Eigen;


#define use_fixed_cov_ar 1
namespace gnss_comm_extra{
    struct CorrectedPseudorangeMeasurement
    {
        int freq = -1;
        int sat = 0;
        std::string sat_sys;
        Eigen::Vector3d sat_pos = Eigen::Vector3d::Zero();
        double pseudorange = 0.0;
        double sigma = 0.0;
        double sv_dt_sec = 0.0;
        double tgd_sec = 0.0;
        double ion_delay_m = 0.0;
        double tro_delay_m = 0.0;
        double elevation_rad = 0.0;
        double psr_std = 0.0;
        bool valid = false;
    };

    inline void matchObsAndEphems(const std::vector<gnss_comm::ObsPtr> &obs,
                                  const std::vector<gnss_comm::EphemBasePtr> &ephems,
                                  std::vector<gnss_comm::ObsPtr> &matched_obs,
                                  std::vector<gnss_comm::EphemBasePtr> &matched_ephems)
    {
        matched_obs.clear();
        matched_ephems.clear();
        matched_obs.reserve(obs.size());
        matched_ephems.reserve(obs.size());

        for (const auto &single_obs : obs)
        {
            if (!single_obs)
            {
                continue;
            }
            for (const auto &single_ephem : ephems)
            {
                if (single_ephem && single_ephem->sat == single_obs->sat)
                {
                    matched_obs.push_back(single_obs);
                    matched_ephems.push_back(single_ephem);
                    break;
                }
            }
        }
    }

    inline double computePseudorangeWeight(const gnss_comm::ObsPtr &obs,
                                           const gnss_comm::EphemBasePtr &ephem_base,
                                           double /*elevation_rad*/,
                                           int freq)
    {
        double weight = 1.0;
        if(freq == 1)
        {
            int l1_idx;
            gnss_comm::L1_freq(obs,&l1_idx);
            if (obs && l1_idx >= 0 && l1_idx < static_cast<int>(obs->psr_std.size()) && obs->psr_std[l1_idx] > 0.0)
            {
                weight /= (obs->psr_std[l1_idx] / 0.16);
            }
        }
        else if(freq == 2)
        {
            int l2_idx;
            gnss_comm::L2_freq(obs,&l2_idx);
            if (obs && l2_idx >= 0 && l2_idx < static_cast<int>(obs->psr_std.size()) && obs->psr_std[l2_idx] > 0.0)
            {
                weight /= (obs->psr_std[l2_idx] / 0.16);
            }
        }
        const int obs_sys = obs ? gnss_comm::satsys(obs->sat, NULL) : SYS_NONE;
        if (gnss_comm::EphemPtr ephem = std::dynamic_pointer_cast<gnss_comm::Ephem>(ephem_base))
        {
            if (obs_sys == SYS_GPS || obs_sys == SYS_BDS)
            {
                weight /= std::max(1.0, ephem->ura - 1.0);
            }
            else if (obs_sys == SYS_GAL)
            {
                weight /= std::max(1.0, ephem->ura - 2.0);
            }
        }

        return std::isfinite(weight) && weight > 0.0 ? weight : 0.0;
    }

    inline std::vector<CorrectedPseudorangeMeasurement> buildCorrectedPseudorangeMeasurements(
        const std::vector<gnss_comm::ObsPtr> &obs,
        const std::vector<gnss_comm::EphemBasePtr> &ephems,
        const Eigen::Vector3d &state_guess,
        bool has_state_guess,
        const std::vector<double> &iono_params = std::vector<double>(8, 0.0),
        int freq_sel = -1)
    {
        std::vector<CorrectedPseudorangeMeasurement> corrected_measurements;

        std::vector<gnss_comm::ObsPtr> matched_obs;
        std::vector<gnss_comm::EphemBasePtr> matched_ephems;
        matchObsAndEphems(obs, ephems, matched_obs, matched_ephems);
        if (matched_obs.size() < 4 || matched_ephems.size() != matched_obs.size())
        {
            return corrected_measurements;
        }

        const std::vector<gnss_comm::SatStatePtr> sat_states = gnss_comm::sat_states(matched_obs, matched_ephems);
        if (sat_states.size() != matched_obs.size())
        {
            return corrected_measurements;
        }

        Eigen::Matrix<double, 7, 1> receiver_state = Eigen::Matrix<double, 7, 1>::Zero();
        if (has_state_guess)
        {
            receiver_state.head<3>() = state_guess;
        }

        Eigen::VectorXd residuals;
        Eigen::MatrixXd jacobian;
        std::vector<Eigen::Vector2d> atmos_delay;
        std::vector<Eigen::Vector2d> all_sv_azel;
        // saastamoninen model for tropospheric delay, 
        // klobuchar model for ionospheric delay, 
        gnss_comm::psr_res_extra(receiver_state, matched_obs, sat_states, iono_params, residuals, jacobian, atmos_delay, all_sv_azel);

        // 预留 L1/L2 频段
        corrected_measurements.reserve( 2* matched_obs.size());

        // 遍历同时具有观测和星历的卫星观测
        for (size_t sat_i = 0; sat_i < matched_obs.size(); ++sat_i)
        {
            CorrectedPseudorangeMeasurement measurement;
            const auto &single_obs = matched_obs[sat_i];
            const auto &sat_state = sat_states[sat_i];
            if (!single_obs || !sat_state)
            {
                continue;
            }

            // 双频段观测可用性检查
            int freq_idx_l1 = -1,freq_idx_l2=-1;
            bool valid_l1 = true, valid_l2 = true;

            gnss_comm::L1_freq(single_obs,&freq_idx_l1);
            gnss_comm::L2_freq(single_obs,&freq_idx_l2);

            if (freq_idx_l1 < 0 || freq_idx_l1 >= static_cast<int>(single_obs->psr.size()) || single_obs->psr[freq_idx_l1] <= 0.0)
            {
                valid_l1 = false;
            }
            if (freq_idx_l2 < 0 || freq_idx_l2 >= static_cast<int>(single_obs->psr.size()) || single_obs->psr[freq_idx_l2] <= 0.0)
            {
                valid_l2 = false;
            }
            if( !valid_l2 && !valid_l1)
            {
                continue;
            }

            if (!std::isfinite(sat_state->pos.x()) || !std::isfinite(sat_state->pos.y()) || !std::isfinite(sat_state->pos.z()))
            {
                continue;
            }

            // 这些数组按卫星索引，使用外层 sat_i
            const double elevation_rad = (sat_i < all_sv_azel.size()) ? all_sv_azel[sat_i](1) : M_PI / 2.0;
            const double ion_delay_m = (sat_i < atmos_delay.size()) ? atmos_delay[sat_i](0) : 0.0;
            const double tro_delay_m = (sat_i < atmos_delay.size()) ? atmos_delay[sat_i](1) : 0.0;
            const auto &matched_ephem = matched_ephems[sat_i];

            for (int freq_i = 0; freq_i < 2; ++freq_i)
            {
                int freq_idx = -1;
                if (freq_i == 0 && valid_l1)
                {
                    freq_idx = freq_idx_l1;
                    measurement.freq = 1;
                }
                else if (freq_i == 1 && valid_l2)
                {
                    freq_idx = freq_idx_l2;
                    measurement.freq = 2;
                }

                // 单频段无效时直接跳过，避免 freq_idx = -1 导致 psr[-1] UB
                if (freq_idx < 0)
                {
                    continue;
                }

                const double weight = computePseudorangeWeight(single_obs, matched_ephem, elevation_rad, freq_idx);
                if (weight <= 0.0)
                {
                    continue;
                }

                const int obs_sys = gnss_comm::satsys(single_obs->sat, NULL);
                measurement.sat = static_cast<int>(single_obs->sat);
                measurement.sat_sys = (obs_sys == SYS_GPS) ? "GPS" :
                                    ((obs_sys == SYS_GLO) ? "GLONASS" :
                                    ((obs_sys == SYS_GAL) ? "Galileo" :
                                    ((obs_sys == SYS_BDS) ? "BeiDou" : "Unknown")));
                measurement.sat_pos = sat_state->pos;
                measurement.pseudorange = single_obs->psr[freq_idx];
                measurement.sigma = std::sqrt(1.0 / weight);
                measurement.sv_dt_sec = sat_state->dt;
                measurement.tgd_sec = sat_state->tgd;
                measurement.ion_delay_m = ion_delay_m;
                measurement.tro_delay_m = tro_delay_m;
                measurement.elevation_rad = elevation_rad;
                measurement.valid = true;
                measurement.psr_std = single_obs->psr_std[freq_idx] <= 0.0 ? 3.0 : single_obs->psr_std[freq_idx];
                corrected_measurements.push_back(measurement);
            }
        }

        return corrected_measurements;
    }


    double getDistance(Eigen::Vector3d p1, Eigen::Vector3d p2)
    {
        double xMod = pow((p1.x() - p2.x()), 2);
        double yMod = pow((p1.y() - p2.y()), 2);
        double zMod = pow((p1.z() - p2.z()), 2);
        double mod = sqrt(xMod + yMod + zMod);
        return mod;
    }

    

    /* get variance for carrier-phase from a single satellite based on elevation/SNR */
    // 无论传入L1还是L2的观测数据，它都使用同一套经验模型和相同的固定参数来计算方差。
    /* TDCP 单差观测标准差模型 (米): 历元间载波相位差分的噪声
     * a: 天顶方向基线 (相位噪声 + 多普勒积分误差, ~1cm)
     * b/sin(el): 低仰角多径与大气残差放大 (15° 时约 4cm)
     * 返回值直接作为 TRDDCP 因子的 σ (米), 与伪距因子 (σ≈1.4m) 量纲统一 */
    double getVarofCp_ele_SNR(gnss_comm::ObsPtr single_sat_data,std::map<int,sv_info> sv_info_map)
    {
        constexpr double a = 0.01;
        constexpr double b = 0.01;
        const auto it = sv_info_map.find(int(single_sat_data->sat));
        const double elR = (it == sv_info_map.end())
            ? M_PI / 4.0
            : it->second.elevation; // radians
        const double sin_el2_d2 = std::isfinite(elR)
            ? std::max(std::sin(elR) * std::sin(elR), 1e-4)
            : 1.0;
        const double var_ele = a * a + b * b / sin_el2_d2;
        return std::sqrt(var_ele);
    }

};

#endif // GNSS_COMM_EXTRA_HPP
