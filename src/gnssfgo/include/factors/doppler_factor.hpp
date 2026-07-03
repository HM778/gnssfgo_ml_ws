/*******************************************************
 * Copyright (C) 2026
 * 
 * This file is reroduce version of gnssfgo Suzuki TR-RTK-GNSS.
 *
 * Author: Huimin ==> taro suzuki
 *******************************************************/

#include <string>
#include <cmath>
#include <Eigen/Eigen>
#include "../datatype.h"

/* doppler factor*/
struct dopplerFactor
{
    /** delta t: 参考历元到当前历元之间的时间
     *  v_x, v_y, v_z: 速度
     * state[i,j][x,y,z]:两个历元估计的相对于初始位置的位置增量
     * */

    dopplerFactor(double v_x, double v_y, double v_z, double delta_t, Eigen::Vector3d var_vector)
                :v_x(v_x),v_y(v_y), v_z(v_z), delta_t(delta_t),var_vector(var_vector){}

    template <typename T>
    bool operator()(const T* state_i, const T* state_j, T* residuals) const
    {
        T est_v_x = (state_j[0] - state_i[0])/ delta_t;
        T est_v_y = (state_j[1] - state_i[1])/ delta_t;
        T est_v_z = (state_j[2] - state_i[2])/ delta_t;

        // est_pseudorange = sqrt(delta_x+ delta_y + delta_z);

        // residuals[0] = (state_j[0] - state_i[0]) - v_x*delta_t;
        // residuals[1] = (state_j[1] - state_i[1]) - v_y*delta_t;
        // residuals[2] = (state_j[2] - state_i[2]) - v_z*delta_t;

        constexpr double kMinConfidence = 1e-3;
        double sx = var_vector(0);
        double sy = var_vector(1);
        double sz = var_vector(2);
        if (!std::isfinite(sx) || sx <= 0.0) sx = kMinConfidence;
        if (!std::isfinite(sy) || sy <= 0.0) sy = kMinConfidence;
        if (!std::isfinite(sz) || sz <= 0.0) sz = kMinConfidence;

        const double wx = std::max(kMinConfidence, std::min(1.0, sx));
        const double wy = std::max(kMinConfidence, std::min(1.0, sy));
        const double wz = std::max(kMinConfidence, std::min(1.0, sz));

        residuals[0] = (est_v_x - T(v_x)) * T(wx);
        residuals[1] = (est_v_y - T(v_y)) * T(wy);
        residuals[2] = (est_v_z - T(v_z)) * T(wz);
        
        // printf("doppler factor residuals|weights: [%f, %f, %f] | [%f, %f, %f]\n",residuals[0],residuals[1],residuals[2], wx, wy, wz);
        return true;
    }

    double v_x, v_y, v_z;
    
    double delta_t;
    Eigen::Vector3d var_vector;
    std::string sat_sys; // satellite system
};

