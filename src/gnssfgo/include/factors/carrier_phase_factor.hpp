/*******************************************************

 *******************************************************/

#include <ceres/ceres.h>
#include <algorithm>
#include <cmath>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <gnss_comm/gnss_constant.hpp>
#include <map>
#include <utility>
#include <vector>
#include "../tools/gnss_comm_extra.h"
#include "../tools/checking.hpp"
#include "../datatype.h"
using namespace gnss_comm;


// tr_score∈[0,1] → 残差乘性权重: 0.1→≈19, 0.5→500, 1→≈993 (单位 1/米)
// 供因子内部(T=ceres::Jet)与因子图日志(T=double)共用, 保证两边数值一致
template <typename T>
inline T trddcpConfWeight(T score)
{
    const T x = (score - T(0.5)) * T(10.0);  // 0→-5, 0.5→0, 1→5
    return T(1.0) + (T(1.0) / (T(1.0) + exp(-x))) * T(999.0);
}

struct TRDDCPFactor
{
    TRDDCPFactor(TRDDMeasurement dd_measurement,
                std::map<int, sv_info> current_sv_info_map,
                std::map<int, sv_info> reference_sv_info_map,
                double sigma, int freq_idx = -1)
                : dd_measurement(dd_measurement),
                    current_sv_info_map(std::move(current_sv_info_map)),
                    reference_sv_info_map(std::move(reference_sv_info_map)),
                    sqrt_info(std::max(1e-3, sigma)), freq_idx(freq_idx) {}

    template <typename T>
    bool operator()(const T* prev_state, const T* curr_state, T* residuals) const
    {
        const auto* curr_master_sv_ptr = getSatelliteInfoPtr(current_sv_info_map, dd_measurement.u_master_SV);
        const auto* curr_i_sv_ptr = getSatelliteInfoPtr(current_sv_info_map, dd_measurement.u_iSV);
        const auto* prev_master_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_master_SV);
        const auto* prev_i_sv_ptr = getSatelliteInfoPtr(reference_sv_info_map, dd_measurement.r_iSV);

        const sv_info& curr_master_sv = *curr_master_sv_ptr;
        const sv_info& curr_i_sv = *curr_i_sv_ptr;
        const sv_info& prev_master_sv = *prev_master_sv_ptr;
        const sv_info& prev_i_sv = *prev_i_sv_ptr;

        Eigen::Vector3d curr_master_pos(curr_master_sv.pos[0], curr_master_sv.pos[1], curr_master_sv.pos[2]);
        Eigen::Vector3d curr_i_pos(curr_i_sv.pos[0], curr_i_sv.pos[1], curr_i_sv.pos[2]);
        Eigen::Vector3d prev_master_pos(prev_master_sv.pos[0], prev_master_sv.pos[1], prev_master_sv.pos[2]);
        Eigen::Vector3d prev_i_pos(prev_i_sv.pos[0], prev_i_sv.pos[1], prev_i_sv.pos[2]);

        // 参考历元：主卫星-接收机位置
        T est_prev_master = sqrt((prev_state[0] - T(prev_master_pos.x())) * (prev_state[0] - T(prev_master_pos.x())) +
                                 (prev_state[1] - T(prev_master_pos.y())) * (prev_state[1] - T(prev_master_pos.y())) +
                                 (prev_state[2] - T(prev_master_pos.z())) * (prev_state[2] - T(prev_master_pos.z())));
        est_prev_master = est_prev_master + T(OMGE_ / CLIGHT_) *
            (T(prev_master_pos.x()) * prev_state[1] - T(prev_master_pos.y()) * prev_state[0]);

        // 参考历元：副卫星-接收机位置
        T est_prev_i = sqrt((prev_state[0] - T(prev_i_pos.x())) * (prev_state[0] - T(prev_i_pos.x())) +
                            (prev_state[1] - T(prev_i_pos.y())) * (prev_state[1] - T(prev_i_pos.y())) +
                            (prev_state[2] - T(prev_i_pos.z())) * (prev_state[2] - T(prev_i_pos.z())));
        est_prev_i = est_prev_i + T(OMGE_ / CLIGHT_) *
            (T(prev_i_pos.x()) * prev_state[1] - T(prev_i_pos.y()) * prev_state[0]);

        // 当前历元：主卫星-接收机位置
        T est_curr_master = sqrt((curr_state[0] - T(curr_master_pos.x())) * (curr_state[0] - T(curr_master_pos.x())) +
                                 (curr_state[1] - T(curr_master_pos.y())) * (curr_state[1] - T(curr_master_pos.y())) +
                                 (curr_state[2] - T(curr_master_pos.z())) * (curr_state[2] - T(curr_master_pos.z())));
        est_curr_master = est_curr_master + T(OMGE_ / CLIGHT_) *
            (T(curr_master_pos.x()) * curr_state[1] - T(curr_master_pos.y()) * curr_state[0]);
        // 当前历元：副卫星-接收机位置
        T est_curr_i = sqrt((curr_state[0] - T(curr_i_pos.x())) * (curr_state[0] - T(curr_i_pos.x())) +
                            (curr_state[1] - T(curr_i_pos.y())) * (curr_state[1] - T(curr_i_pos.y())) +
                            (curr_state[2] - T(curr_i_pos.z())) * (curr_state[2] - T(curr_i_pos.z())));
        est_curr_i = est_curr_i + T(OMGE_ / CLIGHT_) *
            (T(curr_i_pos.x()) * curr_state[1] - T(curr_i_pos.y()) * curr_state[0]);

        // 双差载波相位-估计值：这里是双时差双差，常值整周模糊度会相互消除，
        // 因此不再单独引入 ambiguity state。
        T est_dd_cp = (est_curr_i - est_prev_i) - (est_curr_master - est_prev_master);

        // 参考历元主卫星载波相位，参考历元副卫星载波相位，当前历元主卫星载波相位，当前历元副卫星载波相位观测值
        T prev_master_cp(0), prev_i_cp(0), curr_master_cp(0), curr_i_cp(0);

        auto trddcp_cp_m = [this](const gnss_comm::ObsPtr& o, const sv_info& sv) -> double {
            int li = -1;
            if(freq_idx == 1)  
            {
                L1_freq(o, &li);
                return (li >= 0 && li < (int)o->cp.size()) ? (o->cp[li] * sv.lamda_l1) : 0.0;
            }
            else if(freq_idx == 2) 
            {
                L2_freq(o, &li);
                return (li >= 0 && li < (int)o->cp.size()) ? (o->cp[li] * sv.lamda_l2) : 0.0;
            }
            else
            {
                return 0.0;
            }
        };

        auto trddcp_sigma = [this](const gnss_comm::ObsPtr& o, const std::map<int, sv_info>& svmap) -> double {
            return gnss_comm_extra::getVarofCp_ele_SNR(o, svmap);
        };

        prev_master_cp = T(trddcp_cp_m(dd_measurement.r_master_SV, prev_master_sv));
        prev_i_cp      = T(trddcp_cp_m(dd_measurement.r_iSV, prev_i_sv));
        curr_master_cp = T(trddcp_cp_m(dd_measurement.u_master_SV, curr_master_sv));
        curr_i_cp      = T(trddcp_cp_m(dd_measurement.u_iSV, curr_i_sv));

        // 双差载波相位-观测值
        T dd_cp = (curr_i_cp - prev_i_cp) - (curr_master_cp - prev_master_cp);

        // 载波相位观测的标准差计算，用于求因子权重
        // TDCP 单差 σ 由高程模型给出 (天顶 ~1.4cm, 15° ~4cm);
        // 双差由 4 个独立单差观测平方传播 → ≈ 2× 单差, 与伪距因子同为米量纲
        double sigma_prev_master = 0.0, sigma_prev_i = 0.0, sigma_curr_master = 0.0, sigma_curr_i = 0.0;
        sigma_prev_master = trddcp_sigma(dd_measurement.r_master_SV, reference_sv_info_map);
        sigma_prev_i      = trddcp_sigma(dd_measurement.r_iSV, reference_sv_info_map);
        sigma_curr_master = trddcp_sigma(dd_measurement.u_master_SV, current_sv_info_map);
        sigma_curr_i      = trddcp_sigma(dd_measurement.u_iSV, current_sv_info_map);

        const double sigma_dd = std::sqrt(sigma_prev_master * sigma_prev_master +
                                          sigma_prev_i * sigma_prev_i +
                                          sigma_curr_master * sigma_curr_master +
                                          sigma_curr_i * sigma_curr_i);
     
        const double sigma = std::max(0.10, sigma_dd);

        if(0)
        {
            // 权重由 tr_score 经 trddcpConfWeight 映射而来, 量级见函数注释;
            // 注意: 这不是 1/σ 加权, 高程/SNR 的 sigma 模型仅保留在下方禁用分支中
            residuals[0] = T(est_dd_cp - dd_cp) * trddcpConfWeight(T(sqrt_info));
        }
        else
        {
            residuals[0] = (T(est_dd_cp - dd_cp) / sigma ) * T(sqrt_info);
        }
        
        return true;
    }

    TRDDMeasurement dd_measurement;
    std::map<int, sv_info> current_sv_info_map;
    std::map<int, sv_info> reference_sv_info_map;
    double sqrt_info;
    int freq_idx = -1;  // 1=L1, 2=L2
};


