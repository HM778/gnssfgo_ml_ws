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
#include <gnss_comm/gnss_constant.hpp>

#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include "../datatype.h"

using namespace Eigen;


#define use_fixed_cov_ar 1
namespace gnss_comm_extra{

    // Eigen::Matrix<double, 7, 1> psr_pos_extra(const std::vector<ObsPtr> &obs, 
    //     const std::vector<EphemBasePtr> &ephems, const std::vector<GloEphemPtr> &gephems, const std::vector<double> &iono_params)
    // {
        
    // }


    struct CorrectedPseudorangeMeasurement
    {
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
                                           int l1_idx)
    {
        double weight = 1.0;
        if (obs && l1_idx >= 0 && l1_idx < static_cast<int>(obs->psr_std.size()) && obs->psr_std[l1_idx] > 0.0)
        {
            weight /= (obs->psr_std[l1_idx] / 0.16);
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
        gnss_comm::psr_res(receiver_state, matched_obs, sat_states, iono_params, residuals, jacobian, atmos_delay, all_sv_azel);

        corrected_measurements.reserve(matched_obs.size());
        for (size_t i = 0; i < matched_obs.size(); ++i)
        {
            CorrectedPseudorangeMeasurement measurement;
            const auto &single_obs = matched_obs[i];
            const auto &sat_state = sat_states[i];
            if (!single_obs || !sat_state)
            {
                corrected_measurements.push_back(measurement);
                continue;
            }

            int freq_idx = -1;
            double lamda_dummy = 0.0, lamda_l2_dummy = 0.0;
            getFreqIndex(single_obs, freq_sel, freq_idx, lamda_dummy, lamda_l2_dummy);
            if (freq_idx < 0 || freq_idx >= static_cast<int>(single_obs->psr.size()) || single_obs->psr[freq_idx] <= 0.0)
            {
                corrected_measurements.push_back(measurement);
                continue;
            }

            if (!std::isfinite(sat_state->pos.x()) || !std::isfinite(sat_state->pos.y()) || !std::isfinite(sat_state->pos.z()))
            {
                corrected_measurements.push_back(measurement);
                continue;
            }

            const double elevation_rad = (i < all_sv_azel.size()) ? all_sv_azel[i](1) : M_PI / 2.0;
            const double weight = computePseudorangeWeight(single_obs, matched_ephems[i], elevation_rad, freq_idx);
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
            measurement.ion_delay_m = (i < atmos_delay.size()) ? atmos_delay[i](0) : 0.0;
            measurement.tro_delay_m = (i < atmos_delay.size()) ? atmos_delay[i](1) : 0.0;
            measurement.elevation_rad = elevation_rad;
            measurement.valid = true;
            measurement.psr_std = single_obs->psr_std[freq_idx] <= 0.0 ? 3.0 : single_obs->psr_std[freq_idx];
            corrected_measurements.push_back(measurement);
        }

        return corrected_measurements;
    }


    double getDistanceFrom2Points(Eigen::Vector3d p1, Eigen::Vector3d p2)
    {
        double xMod = pow((p1.x() - p2.x()), 2);
        double yMod = pow((p1.y() - p2.y()), 2);
        double zMod = pow((p1.z() - p2.z()), 2);
        double mod = sqrt(xMod + yMod + zMod);
        return mod;
    }

    /* get variance for carrier-phase from a single satellite based on elevation */
    double getVarofCp(double ele)
    {
        double a = 0.01; 
        double b = 0.01;
        double c = 0;
        double d = 0;
        double square_sigma = pow(a,2) + pow(b,2) / pow(sin(ele * D2R), 2) + pow(c,2) + pow(d,2);
        // return 0.004;
        return square_sigma;
    }

    
    /* get variance for pseudorange from a single satellite based on elevation */
    double getVarofPr(double ele)
    {
        double p_c_ratio = 1; // 10
        double a = 3 * p_c_ratio; 
        double b = 3 * p_c_ratio;
        double c = 0;
        double d = 0;
        double square_sigma = pow(a,2) + pow(b,2) / pow(sin(ele * D2R), 2) + pow(c,2) + pow(d,2);
        // return 0.4;
        // std::cout<<"var of pseudorange -> "<<var<<std::endl;
        return square_sigma;
    }

    /* get pseudorange double-differenced Jacobian matrix */
    void getPrDDJacobian(std::map<int,sv_info> sv_info_map,Eigen::Vector3d u_pose, Eigen::Vector3d base_pose, DDMeasurement dd_measurement, Eigen::MatrixXd& jacobian_matrix, int jac_row,Eigen::MatrixXd& weighting_matrix)
    {
        Eigen::Vector3d pose_r = base_pose;

        /* satellite position*/
        int u_m_sat = dd_measurement.u_master_SV->sat;
        Eigen::Vector3d u_pose_m(sv_info_map[u_m_sat].pos[0], sv_info_map[u_m_sat].pos[1], sv_info_map[u_m_sat].pos[2]);
        double var_u2m = getVarofCp(sv_info_map[u_m_sat].elevation);

        int u_i_sat = dd_measurement.u_iSV->sat;
        Eigen::Vector3d u_pose_i(sv_info_map[u_i_sat].pos[0], sv_info_map[u_i_sat].pos[1], sv_info_map[u_i_sat].pos[2]); 
        double var_u2i = getVarofCp(sv_info_map[u_i_sat].elevation);

        int r_m_sat = dd_measurement.r_master_SV->sat;
        Eigen::Vector3d r_pose_m(sv_info_map[r_m_sat].pos[0], sv_info_map[r_m_sat].pos[1], sv_info_map[r_m_sat].pos[2]);
        double var_r2m = getVarofCp(sv_info_map[r_m_sat].elevation);

        int r_i_sat = dd_measurement.r_iSV->sat;
        Eigen::Vector3d r_pose_i(sv_info_map[r_i_sat].pos[0], sv_info_map[r_i_sat].pos[1], sv_info_map[r_i_sat].pos[2]); 
        double var_r2i = getVarofCp(sv_info_map[r_i_sat].elevation);

        double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

        double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*u_pose(1)-u_pose_m(1)*u_pose(0))/CLIGHT_;
        

        double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*u_pose(1)-u_pose_i(1)*u_pose(0))/CLIGHT_;

        // jacobian_matrix(jac_row, 0) = (u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m;

        // jacobian_matrix(jac_row, 1) = (u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m;

        // jacobian_matrix(jac_row, 2) = (u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m;

        double factor = 1;
        jacobian_matrix(jac_row, 0) = factor * ((u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m);

        jacobian_matrix(jac_row, 1) = factor * ((u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m);

        jacobian_matrix(jac_row, 2) = factor * ((u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m);

        for(int i = 3; i < jacobian_matrix.cols(); i++)
        {
        jacobian_matrix(jac_row, i) = 0;
        }
        
        #if use_fixed_cov_ar
        weighting_matrix(jac_row,jac_row) = 1.0 / (pow(0.4, 2));
        #else
        weighting_matrix(jac_row,jac_row) = 1.0 / ((var_u2m + var_u2i + var_r2m + var_r2i)/4.0); 
        #endif 

    }

    /* get carrier-phase double-differenced Jacobian matrix */
    void getCpDDJacobian(std::map<int,sv_info> sv_info_map, Eigen::Vector3d u_pose, Eigen::Vector3d base_pose, DDMeasurement dd_measurement, Eigen::MatrixXd& jacobian_matrix, int jac_row, int carrier_phase_index, Eigen::MatrixXd& weighting_matrix)
    {
        Eigen::Vector3d pose_r = base_pose;

        /* satellite position*/
        int u_m_sat = dd_measurement.u_master_SV->sat;
        Eigen::Vector3d u_pose_m(sv_info_map[u_m_sat].pos[0], sv_info_map[u_m_sat].pos[1], sv_info_map[u_m_sat].pos[2]);
        double var_u2m = getVarofCp(sv_info_map[u_m_sat].elevation);

        int u_i_sat = dd_measurement.u_iSV->sat;
        Eigen::Vector3d u_pose_i(sv_info_map[u_i_sat].pos[0], sv_info_map[u_i_sat].pos[1], sv_info_map[u_i_sat].pos[2]); 
        double var_u2i = getVarofCp(sv_info_map[u_i_sat].elevation);

        int r_m_sat = dd_measurement.r_master_SV->sat;
        Eigen::Vector3d r_pose_m(sv_info_map[r_m_sat].pos[0], sv_info_map[r_m_sat].pos[1], sv_info_map[r_m_sat].pos[2]);
        double var_r2m = getVarofCp(sv_info_map[r_m_sat].elevation);

        int r_i_sat = dd_measurement.r_iSV->sat;
        Eigen::Vector3d r_pose_i(sv_info_map[r_i_sat].pos[0], sv_info_map[r_i_sat].pos[1], sv_info_map[r_i_sat].pos[2]); 
        double var_r2i = getVarofCp(sv_info_map[r_i_sat].elevation);
        // double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
        // double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);

        // double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
        // double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);

        double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

        double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*u_pose(1)-u_pose_m(1)*u_pose(0))/CLIGHT_;
        

        double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*u_pose(1)-u_pose_i(1)*u_pose(0))/CLIGHT_;

        double factor = 1;
        // LOG(INFO)<<"jacobian_matrix.cols()"<<jacobian_matrix.cols();
        // LOG(INFO)<<"jacobian_matrix.rows()"<<jacobian_matrix.rows();
        // LOG(INFO)<< "jac_row-> "<<jac_row;
        // LOG(INFO)<<"carrier_phase_index->" << carrier_phase_index;
        // LOG(INFO)<<"jacobian_matrix(jac_row, 0)->" << jacobian_matrix(jac_row, 0);
        jacobian_matrix(jac_row, 0) = factor * ((u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m);

        jacobian_matrix(jac_row, 1) = factor * ((u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m);

        jacobian_matrix(jac_row, 2) = factor * ((u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m);

        for(int i = 3; i < jacobian_matrix.cols(); i++)
        {
        jacobian_matrix(jac_row, i) = 0; 
        }
        
        // jacobian_matrix(jac_row, 3 + carrier_phase_index) = dd_measurement.r_master_SV.lamda;
        jacobian_matrix(jac_row, 3 + carrier_phase_index) = sv_info_map[r_m_sat].lamda;

        // weighting_matrix(jac_row,jac_row) = 1.0/(var_u2m + var_u2i + var_r2m + var_r2i);
        weighting_matrix(jac_row,jac_row) = 1.0/(pow(0.004, 2));

        #if use_fixed_cov_ar
        
        weighting_matrix(jac_row,jac_row) = 1.0/(pow(0.004, 2));
        #else
        weighting_matrix(jac_row,jac_row) = 1.0 / ((var_u2m + var_u2i + var_r2m + var_r2i)/4.0); 
        #endif 
    }

    /* get variance for carrier-phase from a single satellite based on elevation/SNR */
    double getVarofCp_ele_SNR(gnss_comm::ObsPtr single_sat_data,std::map<int,sv_info> sv_info_map)
    {
        Eigen::Matrix<double,4,1> parameters;
        parameters<<50.0, 30.0, 30.0, 10.0; // loosely coupled 
        // parameters<<50.0, 30.0, 20.0, 30.0; // loosely coupled 
        double snr_1 = parameters(0); // T = 50
        double snr_A = parameters(1); // A = 30
        double snr_a = parameters(2);// a = 30
        double snr_0 = parameters(3); // F = 10
        const double snr_R = (!single_sat_data || single_sat_data->CN0.empty())
            ? std::numeric_limits<double>::quiet_NaN()
            : single_sat_data->CN0[0];

        const auto it = sv_info_map.find(int(single_sat_data->sat));
        const double elR = (it == sv_info_map.end())
            ? std::numeric_limits<double>::quiet_NaN()
            : it->second.elevation; // radians

        // if(elR<15) elR = 30;
        if(std::isfinite(elR) && elR < (15.0 * D2R))
        {
        // LOG(INFO) << "satellite elevation -> " << elR;
        }

        const double sin_el = std::sin(elR);
        const double sin_el2 = sin_el * sin_el;
        const double q_R_1 = (sin_el2 > 1e-12) ? (1.0 / sin_el2) : 1e12;
        const double q_R_2 = std::pow(10.0, (-(snr_R - snr_1) / snr_a));
        const double denom = (std::pow(10.0, (-(snr_0 - snr_1) / snr_a)) - 1.0);
        const double q_R_3 = (((denom != 0.0) ? (snr_A / denom / (snr_0 - snr_1)) : 0.0) * (snr_R - snr_1) + 1.0);
        const double q_R = q_R_1 * (q_R_2 * q_R_3);
        const double var = (q_R > 1e-12) ? (1.0 / q_R) : 1e12; // larger -> larger uncertainty

        double a = 0.01; 
        double b = 0.01;
        double c = 0;
        double d = 0;
        const double sin_el2_d2 = sin_el2 > 1e-12 ? sin_el2 : 1e-12;
        const double var_ele = pow(a,2) + pow(b,2) / sin_el2_d2 + pow(c,2) + pow(d,2);
        // return 0.004;
        #if useEleVar
        return sqrt(var_ele);
        #else 
        return 0.0001 * std::sqrt(1.0/var);
        #endif 
    }

};

#endif // GNSS_COMM_EXTRA_HPP
