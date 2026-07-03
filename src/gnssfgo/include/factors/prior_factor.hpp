#include <string>
#include <cmath>
#include <Eigen/Eigen>
#include "../datatype.h"
#define state_size 7 // x,y,z, clk_gps, clk_glo, clk_gal, clk_bds

struct FirstStatePriorFactor
{
    FirstStatePriorFactor(const std::array<double, state_size> &mean, const Eigen::Matrix<double, state_size, state_size> &sqrt_info)
        : mean_(mean), sqrt_info_(sqrt_info)
    {
    }

    template <typename T>
    bool operator()(const T *const state, T *residual) const
    {
        T dx[state_size];
        for (int i = 0; i < state_size; ++i)
        {
            dx[i] = state[i] - T(mean_[i]);
        }

        for (int r = 0; r < state_size; ++r)
        {
            T value = T(0.0);
            for (int c = 0; c < state_size; ++c)
            {
                value += T(sqrt_info_(r, c)) * dx[c];
            }
            residual[r] = value;
        }
        return true;
    }

    std::array<double, state_size> mean_{};
    Eigen::Matrix<double, state_size, state_size> sqrt_info_ = Eigen::Matrix<double, state_size, state_size>::Identity();

};