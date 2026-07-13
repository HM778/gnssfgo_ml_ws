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
*
* gnss_spp_extra.hpp
* Extended SPP functions that use both L1 and L2 band data (prefer L1, fallback L2)
* for all four GNSS systems (GPS, GLONASS, Galileo, BeiDou).
* Compared to gnss_spp.hpp which uses L1-only, these functions can also leverage
* L2 observations when L1 is unavailable, increasing satellite availability.
*/

#ifndef GNSS_SPP_EXTRA_HPP_
#define GNSS_SPP_EXTRA_HPP_

#include <eigen3/Eigen/Dense>
#include "gnss_constant.hpp"

namespace gnss_comm
{
    /* filter observation, keep L1 or L2 Obs from four GNSS systems ------------
    * args   : std::vector<ObsPtr>&         obs                 I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems              I   GNSS ephemeris data
    * args   : std::vector<ObsPtr>&         L1L2_obs            O   filtered observation data
    *          std::vector<EphemBasePtr>&   L1L2_ephems         O   filtered ephemeris data
    * return : void
    *-----------------------------------------------------------------------------*/
    void filter_L1L2_extra(const std::vector<ObsPtr> &obs, const std::vector<EphemBasePtr> &ephems,
        std::vector<ObsPtr> &L1L2_obs, std::vector<EphemBasePtr> &L1L2_ephems);

    /* calculate satellite states -----------------------------------------------
    * args   : std::vector<ObsPtr>&         obs                 I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems              I   GNSS ephemeris data
    * return : std::vector<SatStatePtr>     satellite states
    *-----------------------------------------------------------------------------*/
    std::vector<SatStatePtr> sat_states_extra(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems);

    /* calculate pseudo-range residual and Jacobian (L1-prefer, L2-fallback) ----
    * Tries L1 frequency first; if unavailable, falls back to L2.
    * Ionospheric delay is scaled appropriately for the used frequency.
    * args   : Eigen::Matrix<double, 7, 1>&     rcv_state          I   receiver state
    *          std::vector<ObsPtr>&             obs                I   GNSS observations
    *          std::vector<SatStatePtr>&        all_sv_states      I   satellite states
    *          std::vector<double>&             iono_params        I   ionospheric parameters
    *          Eigen::VectorXd&                 res                O   pseudo-range residual
    *          Eigen::MatrixXd&                 J                  O   Jacobian
    *          std::vector<Eigen::Vector2d>&    atmos_delay        O   ion and tro delay
    *          std::vector<Eigen::Vector2d>&    all_sv_azel        O   satellite azimuth and elevation
    * return : void
    *-----------------------------------------------------------------------------*/
    void psr_res_extra(const Eigen::Matrix<double, 7, 1> &rcv_state, const std::vector<ObsPtr> &obs,
        const std::vector<SatStatePtr> &all_sv_states, const std::vector<double> &iono_params,
        Eigen::VectorXd &res, Eigen::MatrixXd &J, std::vector<Eigen::Vector2d> &atmos_delay,
        std::vector<Eigen::Vector2d> &all_sv_azel);

    /* positioning by pseudo-range localization (L1/L2) -------------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          std::vector<double>&         iono_params I   ionosphere parameters
    * return : receiver position in ECEF and four clock bias for 4 constellations
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 7, 1> psr_pos_extra(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, const std::vector<double> &iono_params);

    /* calculate doppler residual and Jacobian (L1-prefer, L2-fallback) ---------
    * args   : Eigen::Matrix<double, 4, 1>&     rcv_state          I   receiver state
    *          Eigen::Vector3d&                 rcv_ecef           I   receiver ECEF position
    *          std::vector<ObsPtr>&             obs                I   GNSS observations
    *          std::vector<SatStatePtr>&        all_sv_states      I   satellite states
    *          Eigen::VectorXd&                 res                O   doppler residual
    *          Eigen::MatrixXd&                 J                  O   Jacobian
    * return : void
    *-----------------------------------------------------------------------------*/
    void dopp_res_extra(const Eigen::Matrix<double, 4, 1> &rcv_state, const Eigen::Vector3d &rcv_ecef,
                      const std::vector<ObsPtr> &obs, const std::vector<SatStatePtr> &all_sv_states,
                      Eigen::VectorXd &res, Eigen::MatrixXd &J);

    /* calculate velocity by using Doppler measurement (L1/L2) ------------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          Eigen::Vector3d&             ref_ecef    IO  reference ECEF position, (0,0,0) if unknown
    * return : receiver velocity in ECEF and clock bias changing rate
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 4, 1> dopp_vel_extra(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef);

    /* calculate velocity by using Doppler with GGL weight (L1/L2) --------------
    * args   : std::vector<ObsPtr>&         obs         I   GNSS observation data
    *          std::vector<EphemBasePtr>&   ephems      I   GNSS ephemeris data
    *          Eigen::Vector3d&             ref_ecef    IO  reference ECEF position
    * return : receiver velocity in ECEF and clock bias changing rate
    *-----------------------------------------------------------------------------*/
    Eigen::Matrix<double, 4, 1> dopp_vel_GGLweight_extra(const std::vector<ObsPtr> &obs,
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef);
}

#endif
