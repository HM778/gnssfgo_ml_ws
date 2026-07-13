/**
* This file is part of gnss_comm.
*
* Copyright (C) 2021 Aerial Robotics Group, Hong Kong University of Science and Technology
* Author: CAO Shaozu (shaozu.cao@gmail.com)
*
* gnss_comm is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* gnss_comm is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with gnss_comm. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GNSS_SPS_HPP_
#define GNSS_SPS_HPP_

#include <eigen3/Eigen/Dense>
#include "gnss_constant.hpp"
#include <mutex>

namespace gnss_comm
{
    /* filter observation, only keep L1-observed Obs -----------------------------
    * args   : std::vector<ObsPtr>&         obs                 I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems              I   GNSS ephemeris data
    * args   : std::vector<ObsPtr>&         L1_obs              I   GNSS L1 observation data
    *          std::vector<EphemBasePtr>&   L1_ephems           I   GNSS L1 ephemeris data
    * return : void
    *-----------------------------------------------------------------------------*/
    
    extern Eigen::Matrix<double,4,4> dopp_cov;
    extern std::mutex dopp_cov_mux;
    bool dopp_cov_set(Eigen::Matrix<double,4,4> cov);
    Eigen::Matrix<double,4,4> dopp_cov_get();

    void filter_L1(const std::vector<ObsPtr> &obs, const std::vector<EphemBasePtr> &ephems,
        std::vector<ObsPtr> &L1_obs, std::vector<EphemBasePtr> &L1_ephems);

    /* filter observation, only keep dual-frequency (L1+L2) Obs -------------------
    * args   : std::vector<ObsPtr>&         obs                 I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems              I   GNSS ephemeris data
    * args   : std::vector<ObsPtr>&         IF_obs              O   GNSS dual-freq observation data
    *          std::vector<EphemBasePtr>&   IF_ephems           O   GNSS dual-freq ephemeris data
    * return : void
    *-----------------------------------------------------------------------------*/
    void filter_dual_freq(const std::vector<ObsPtr> &obs, const std::vector<EphemBasePtr> &ephems,
        std::vector<ObsPtr> &IF_obs, std::vector<EphemBasePtr> &IF_ephems);

    /* calculate satellite states -----------------------------------------------
    * args   : std::vector<ObsPtr>&         obs                 I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems              I   GNSS ephemeris data
    * return : std::vector<SatStatePtr>     satellite states
    *-----------------------------------------------------------------------------*/
    std::vector<SatStatePtr> sat_states(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems);
    
    /* calculate pseudo-range residual and Jacobian -------------------------------
    * args   : Eigen::Matrix<double, 7, 1>&     rcv_state          I   receiver state
    *          std::vector<ObsPtr>&             obs                I   GNSS observations
    *          std::vector<SatStatePtr>&        all_sv_states      I   satellite states
    *          std::vector<double>&             iono_params        I   ionospheric parameters
    *          Eigen::VectorXd&                 res                O   pseudo-range residual
    *          Eigen::MatrixXd&                 J                  O   Jacobian
    *          std::vector<Eigen::Vector2d>&    atmos_delay        O   ion and tro delay
    *          std::vector<Eigen::Vector2d>&    all_sv_azel        O   satellite azimuth and elevation (radius)
    * return : void
    *-----------------------------------------------------------------------------*/
    void psr_res(const Eigen::Matrix<double, 7, 1> &rcv_state, const std::vector<ObsPtr> &obs, 
        const std::vector<SatStatePtr> &all_sv_states, const std::vector<double> &iono_params, 
        Eigen::VectorXd &res, Eigen::MatrixXd &J, std::vector<Eigen::Vector2d> &atmos_delay, 
        std::vector<Eigen::Vector2d> &all_sv_azel);

    /* positioning by pseudo-range localization ----------------------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          std::vector<double>&         iono_params I   ionosphere parameters
    * return : receiver position in ECEF and four clock bias for 4 constellations 
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 7, 1> psr_pos(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, const std::vector<double> &iono_params);

    /* positioning by pseudo-range localization (Ionosphere-Free) ------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data (L1 filtered)
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    * return : receiver position in ECEF and four clock bias for 4 constellations
    * notes  : uses L1+L2 ionosphere-free combination when dual-freq available,
    *          falls back to L1+Klobuchar for single-freq satellites
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 7, 1> psr_pos_IF(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems);

    /* calculate doppler residual and Jacobian -------------------------------
    * args   : Eigen::Matrix<double, 4, 1>&     rcv_state          I   receiver state
    *          Eigen::Vector3d&                 rcv_ecef           I   receiver ECEF position
    *          std::vector<ObsPtr>&             obs                I   GNSS observations
    *          std::vector<SatStatePtr>&        all_sv_states      I   satellite states
    *          Eigen::VectorXd&                 res                O   pseudo-range residual
    *          Eigen::MatrixXd&                 J                  O   Jacobian
    * return : void
    *-----------------------------------------------------------------------------*/
    void dopp_res(const Eigen::Matrix<double, 4, 1> &rcv_state, const Eigen::Vector3d &rcv_ecef,
                  const std::vector<ObsPtr> &obs, const std::vector<SatStatePtr> &all_sv_states, 
                  Eigen::VectorXd &res, Eigen::MatrixXd &J);

    /* calculate velocity by using Doppler measurement -------------------------------------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          Eigen::Vector3d&             ref_ecef    IO  reference ECEF position, (0,0,0) if unknown
    * return : receiver velocity in ECEF and clock bias changing rate + covariance(huimin)
    *----------------------------------------------------------------------------------------------------*/
    Eigen::Matrix<double, 4, 1> dopp_vel(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef);

    Eigen::Matrix<double, 4, 1> dopp_vel_GGLweight(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef);

    /* calculate velocity by using Ionosphere-Free Doppler --------------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data (L1 filtered)
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          Eigen::Vector3d&             ref_ecef    IO  reference ECEF position
    * return : receiver velocity in ECEF and clock bias changing rate
    * notes  : uses L1+L2 IF Doppler when dual-freq available
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 4, 1> dopp_vel_IF(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef);
}

#endif