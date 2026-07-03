/*******************************************************
 * Copyright (C) 2026
 * 
 * This file is reroduce version of gnssfgo Suzuki TR-RTK-GNSS.
 *
 * Author: Huimin ==> taro suzuki
 *******************************************************/

// #include <gnss_comm/gnss_constant.hpp>
#include <algorithm>
#include <ceres/ceres.h>
#include <Eigen/Eigen>
#include <map>
#include <utility>

#include "../datatype.h"
#include "../tools/checking.hpp"

struct pseudorangeFactor
{
    //伪距因子
    pseudorangeFactor(std::string sat_sys,
                    double s_g_x, double s_g_y, double s_g_z,
                    double pseudorange, double weight, double sv_dt_sec = 0.0, double tgd_sec = 0.0,
                    double ion_delay_m = 0.0, double tro_delay_m = 0.0)
                : sat_sys(std::move(sat_sys)),
                    s_g_x(s_g_x), s_g_y(s_g_y), s_g_z(s_g_z), pseudorange(pseudorange), weight(weight), sv_clk_m(sv_dt_sec * CLIGHT_), tgd_m(tgd_sec * CLIGHT_),
                    ion_delay_m(ion_delay_m), tro_delay_m(tro_delay_m) {}

    template <typename T>
    bool operator()(const T* state, T* residuals) const
    {
        residuals[0] = T(0);
        if (!std::isfinite(s_g_x) || !std::isfinite(s_g_y) || !std::isfinite(s_g_z) ||
            !std::isfinite(pseudorange) || !std::isfinite(weight) || weight <= 1e-6)
        {
            return true;
        }

        T est_pseudorange; 

        T delta_x = pow((state[0] - s_g_x),2);
        T delta_y = pow((state[1] - s_g_y),2);
        T delta_z = pow((state[2] - s_g_z),2);

        est_pseudorange = sqrt(delta_x+ delta_y + delta_z);

        est_pseudorange = est_pseudorange + OMGE_ * (s_g_x*state[1]-s_g_y*state[0]) / CLIGHT_;

        // Satellite clock bias correction (meters).
        // Observation model: rho = r + c*(dt_rcv - dt_sv) + ...
        // 加上估算的误差，只关注位置和钟差的增量，卫星钟差、大气延迟作为已知量进行校正
        est_pseudorange = est_pseudorange - T(sv_clk_m) + T(tgd_m) + T(ion_delay_m) + T(tro_delay_m);
        
        if(sat_sys == "GPS")
        {
            // Receiver clock bias in meters for GPS.
            est_pseudorange = est_pseudorange + state[3];
        }
        else if(sat_sys == "GLONASS")
        {
            // Receiver clock bias in meters for GLONASS.
            est_pseudorange = est_pseudorange + state[4];
        }
        else if(sat_sys == "Galileo")
        {
            // Receiver clock bias in meters for Galileo.
            est_pseudorange = est_pseudorange + state[5];
        }
        else if(sat_sys == "BeiDou")
        {
            // Receiver clock bias in meters for BeiDou.
            est_pseudorange = est_pseudorange + state[6];
        }

        const double confidence = std::max(1e-3, std::min(1.0, 1.0 / (weight * weight)));
        residuals[0] = (est_pseudorange - T(pseudorange)) * T(confidence);
        // printf("pseudorange factor residual|confidence: %f | %f \n", residuals[0], static_cast<double>(confidence));
        if (!ceres::IsFinite(residuals[0]))
        {
            residuals[0] = T(0);
            return true;
        }
        
        return true;
    }

    double s_g_x, s_g_y, s_g_z, pseudorange, weight;
    double sv_clk_m = 0.0;
    double tgd_m = 0.0;
    double ion_delay_m = 0.0;
    double tro_delay_m = 0.0;
    std::string sat_sys; // satellite system
};

struct DDPseudorangeFactor
{
    DDPseudorangeFactor(DDMeasurement dd_measurement,
                            std::map<int, sv_info> current_sv_info_map,
                            std::map<int, sv_info> reference_sv_info_map,
                            double out_conf = 1.0, int freq_idx = -1)
        : dd_measurement(dd_measurement),
          current_sv_info_map(std::move(current_sv_info_map)),
          reference_sv_info_map(std::move(reference_sv_info_map)),
          out_conf(out_conf), freq_idx(freq_idx) {}

    template <typename T>
    bool operator()(const T* prev_state, const T* curr_state, T* residuals) const
    {
        const auto* curr_master_sv_ptr = getSatelliteInfoPtr(current_sv_info_map, dd_measurement.u_master_SV);
        const auto* curr_i_sv_ptr = getSatelliteInfoPtr(current_sv_info_map, dd_measurement.u_iSV);
        const auto* prev_master_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_master_SV);
        const auto* prev_i_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_iSV);
        if (!curr_master_sv_ptr || !curr_i_sv_ptr || !prev_master_sv_ptr || !prev_i_sv_ptr)
        {
            residuals[0] = T(0);
            return true;
        }

        // 当前历元主卫星，次卫星
        const sv_info& curr_master_sv = *curr_master_sv_ptr;
        const sv_info& curr_i_sv = *curr_i_sv_ptr;
        // 参考历元主卫星，次卫星
        const sv_info& prev_master_sv = *prev_master_sv_ptr;
        const sv_info& prev_i_sv = *prev_i_sv_ptr;

        Eigen::Vector3d curr_master_pos(curr_master_sv.pos[0], curr_master_sv.pos[1], curr_master_sv.pos[2]);
        Eigen::Vector3d curr_i_pos(curr_i_sv.pos[0], curr_i_sv.pos[1], curr_i_sv.pos[2]);
        Eigen::Vector3d prev_master_pos(prev_master_sv.pos[0], prev_master_sv.pos[1], prev_master_sv.pos[2]);
        Eigen::Vector3d prev_i_pos(prev_i_sv.pos[0], prev_i_sv.pos[1], prev_i_sv.pos[2]);

        // 当前、参考历元主卫星，次卫星 sagnac项矫正
        T est_prev_master = sqrt((prev_state[0] - T(prev_master_pos.x())) * (prev_state[0] - T(prev_master_pos.x())) +
                                 (prev_state[1] - T(prev_master_pos.y())) * (prev_state[1] - T(prev_master_pos.y())) +
                                 (prev_state[2] - T(prev_master_pos.z())) * (prev_state[2] - T(prev_master_pos.z())));
        est_prev_master = est_prev_master + T(OMGE_ / CLIGHT_) *
            (T(prev_master_pos.x()) * prev_state[1] - T(prev_master_pos.y()) * prev_state[0]);

        T est_prev_i = sqrt((prev_state[0] - T(prev_i_pos.x())) * (prev_state[0] - T(prev_i_pos.x())) +
                            (prev_state[1] - T(prev_i_pos.y())) * (prev_state[1] - T(prev_i_pos.y())) +
                            (prev_state[2] - T(prev_i_pos.z())) * (prev_state[2] - T(prev_i_pos.z())));
        est_prev_i = est_prev_i + T(OMGE_ / CLIGHT_) *
            (T(prev_i_pos.x()) * prev_state[1] - T(prev_i_pos.y()) * prev_state[0]);

        T est_curr_master = sqrt((curr_state[0] - T(curr_master_pos.x())) * (curr_state[0] - T(curr_master_pos.x())) +
                                 (curr_state[1] - T(curr_master_pos.y())) * (curr_state[1] - T(curr_master_pos.y())) +
                                 (curr_state[2] - T(curr_master_pos.z())) * (curr_state[2] - T(curr_master_pos.z())));
        est_curr_master = est_curr_master + T(OMGE_ / CLIGHT_) *
            (T(curr_master_pos.x()) * curr_state[1] - T(curr_master_pos.y()) * curr_state[0]);

        T est_curr_i = sqrt((curr_state[0] - T(curr_i_pos.x())) * (curr_state[0] - T(curr_i_pos.x())) +
                            (curr_state[1] - T(curr_i_pos.y())) * (curr_state[1] - T(curr_i_pos.y())) +
                            (curr_state[2] - T(curr_i_pos.z())) * (curr_state[2] - T(curr_i_pos.z())));
        est_curr_i = est_curr_i + T(OMGE_ / CLIGHT_) *
            (T(curr_i_pos.x()) * curr_state[1] - T(curr_i_pos.y()) * curr_state[0]);

        // 双差伪距增量-估计值
        T est_dd_pr = (est_curr_i - est_prev_i) - (est_curr_master - est_prev_master);

        int l1_idx_prev_master = -1;
        int l1_idx_prev_i = -1;
        int l1_idx_curr_master = -1;
        int l1_idx_curr_i = -1;
        L1_freq(dd_measurement.r_master_SV, &l1_idx_prev_master);
        L1_freq(dd_measurement.r_iSV, &l1_idx_prev_i);
        L1_freq(dd_measurement.u_master_SV, &l1_idx_curr_master);
        L1_freq(dd_measurement.u_iSV, &l1_idx_curr_i);

        auto valid_psr_index = [](const gnss_comm::ObsPtr& obs, int idx) -> bool {
            return obs && idx >= 0 && idx < static_cast<int>(obs->psr.size());
        };
        if (!valid_psr_index(dd_measurement.r_master_SV, l1_idx_prev_master) ||
            !valid_psr_index(dd_measurement.r_iSV, l1_idx_prev_i) ||
            !valid_psr_index(dd_measurement.u_master_SV, l1_idx_curr_master) ||
            !valid_psr_index(dd_measurement.u_iSV, l1_idx_curr_i))
        {
            residuals[0] = T(0);
            return true;
        }

        const T prev_master_pr = T(dd_measurement.r_master_SV->psr[l1_idx_prev_master]);
        const T prev_i_pr = T(dd_measurement.r_iSV->psr[l1_idx_prev_i]);
        const T curr_master_pr = T(dd_measurement.u_master_SV->psr[l1_idx_curr_master]);
        const T curr_i_pr = T(dd_measurement.u_iSV->psr[l1_idx_curr_i]);
        // 双差伪距增量-测量值
        const T dd_pr = (curr_i_pr - prev_i_pr) - (curr_master_pr - prev_master_pr);

        // 观测值的方差估计，基于伪距观测的标准差进行传播，假设观测误差独立且服从正态分布
        auto obs_sigma = [](const gnss_comm::ObsPtr& obs, int l1_idx) -> double {
            if (!obs || l1_idx < 0 || l1_idx >= static_cast<int>(obs->psr_std.size()) || obs->psr_std[l1_idx] <= 0.0)
            {
                return 3.0;
            }
            return obs->psr_std[l1_idx];
        };

        const double sigma_prev_master = obs_sigma(dd_measurement.r_master_SV, l1_idx_prev_master);
        const double sigma_prev_i = obs_sigma(dd_measurement.r_iSV, l1_idx_prev_i);
        const double sigma_curr_master = obs_sigma(dd_measurement.u_master_SV, l1_idx_curr_master);
        const double sigma_curr_i = obs_sigma(dd_measurement.u_iSV, l1_idx_curr_i);
        const double sigma = std::max(1e-3,
            std::sqrt(sigma_prev_master * sigma_prev_master + sigma_prev_i * sigma_prev_i +
                      sigma_curr_master * sigma_curr_master + sigma_curr_i * sigma_curr_i));

        const double bounded_out_conf = std::max(1e-3, std::min(1.0, out_conf));
        const double sqrt_info = bounded_out_conf / sigma;
        residuals[0] = (est_dd_pr - dd_pr) * T(sqrt_info);
        // printf("sat_pair: %d & %d DDPR residuals: %f | confidence: %f |out_conf:%f \n",dd_measurement.u_master_SV->sat,dd_measurement.u_iSV->sat ,residuals[0], confidence,out_conf);
        return true;
    }

    DDMeasurement dd_measurement;
    std::map<int, sv_info> current_sv_info_map;
    std::map<int, sv_info> reference_sv_info_map;
    int freq_idx = -1;  // -1=L1, -2=L2
    double out_conf = 1.0; // confidence for outlier handling, currently not used
};
