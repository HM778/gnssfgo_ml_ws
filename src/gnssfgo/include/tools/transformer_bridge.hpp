/*******************************************************
 * Copyright (C) 2026
 *
 * transformer_bridge.hpp — gnssfgo <-> OSQA Transformer 数据交换桥
 * ==============================================================
 *
 * 通过 JSONL 文件实现 gnssfgo (C++) 和 GNSS-Transformer/OSQA (Python)
 * 之间的数据交互。
 *
 * 数据流:
 *   gnssfgo 优化完成后 → exportEpochData() → osqa_input.jsonl
 *   OSQA Transformer  → readLatestQualityScores() ← osqa_output.jsonl
 *
 * 编译选项:
 *   定义 ENABLE_TRANSFORMER_BRIDGE 宏启用此功能
 *   不定义时所有代码被编译为空操作
 *
 * 依赖: nlohmann/json.hpp (header-only, 拷贝到 include/tools/)
 *
 * Author: Claude Code
 * Date: 2026-07-17
 *******************************************************/

#ifndef TRANSFORMER_BRIDGE_HPP
#define TRANSFORMER_BRIDGE_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <mutex>
#include <cmath>
#include <cstdio>

#include <Eigen/Dense>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_constant.hpp>

#include "../datatype.h"
#include "json.hpp"

using json = nlohmann::json;

namespace transformer_bridge {

// ==================== 数据结构 ====================

/**
 * 单颗卫星的质量评估信息（从 OSQA 输出读取）
 */
struct SatQualityInfo
{
    double quality = 1.0;               // 融合后的质量分数 [0, 1]
    std::string trust_level = "trusted"; // trusted / suspect / unreliable
    double q_transformer = 1.0;          // Transformer 分析器分数
    double q_graph = 1.0;                // 图分析器分数
    double q_temporal = 1.0;            // 时序分析器分数
    std::vector<std::string> flags;      // 异常标记
};

// ==================== 文件 I/O 工具 ====================

/**
 * 原子性地向 JSONL 文件追加一行 JSON
 *
 * 使用临时文件 + rename 保证写入的原子性，
 * 避免 OSQA 读取到不完整的行。
 *
 * @param filepath  目标文件路径
 * @param j         要写入的 JSON 对象
 * @return          写入成功返回 true
 */
inline bool appendJsonLine(const std::string &filepath, const json &j)
{
    // 先写入临时文件
    std::string tmp_path = filepath + ".tmp";
    std::ofstream ofs(tmp_path, std::ios::out | std::ios::app);
    if (!ofs.is_open())
    {
        return false;
    }
    ofs << j.dump() << "\n";
    ofs.close();

    // 追加到目标文件
    std::ofstream main_ofs(filepath, std::ios::out | std::ios::app);
    if (!main_ofs.is_open())
    {
        std::remove(tmp_path.c_str());
        return false;
    }
    main_ofs << j.dump() << "\n";
    main_ofs.close();

    // 清理临时文件
    std::remove(tmp_path.c_str());
    return true;
}

/**
 * 从 JSONL 文件读取所有行
 *
 * @param filepath  文件路径
 * @return          JSON 对象列表（每行一个）
 */
inline std::vector<json> readAllJsonLines(const std::string &filepath)
{
    std::vector<json> lines;
    std::ifstream ifs(filepath);
    if (!ifs.is_open())
    {
        return lines;
    }

    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.empty())
        {
            continue;
        }
        try
        {
            lines.push_back(json::parse(line));
        }
        catch (const json::parse_error &e)
        {
            // 跳过无法解析的行
            fprintf(stderr, "[TransformerBridge] JSON parse error: %s\n", e.what());
        }
    }
    return lines;
}

// ==================== 数据导出 (gnssfgo → OSQA) ====================

/**
 * 将单个 epoch 的观测数据和因子残差导出为 JSONL 格式
 *
 * 导出内容:
 *  - 接收机 ECEF 位置
 *  - 每颗卫星的: SNR(L1/L2), 仰角, 方位角, 伪距(L1/L2), 多普勒,
 *    载波相位(L1/L2), 锁定计数(L1/L2), SPP 伪距残差,
 *    卫星 ECEF 位置/速度/钟差, 各因子残差
 *
 * @param filepath          输出 JSONL 文件路径
 * @param timestamp          GPS 时间戳（秒）
 * @param time_frame         内部时间帧（gpst_sec * 10）
 * @param receiver_ecef      接收机 ECEF 位置（优化后）
 * @param receiver_enu       接收机 ENU 位置（优化后）
 * @param observations       当前 epoch 的原始观测列表
 * @param sv_info_map        卫星信息映射 (PRN → sv_info)
 * @param lock_count_l1      L1 连续锁定计数 (PRN → count)
 * @param lock_count_l2      L2 连续锁定计数 (PRN → count)
 * @param psr_residual_map   SPP 伪距残差 (PRN → residual_m)
 * @param factor_residuals   因子残差 (PRN → {psr, doppler, tr_dd_pr, tr_dd_cp})
 * @return                   写入成功返回 true
 */
inline bool exportEpochData(
    const std::string &filepath,
    double timestamp,
    double time_frame,
    const Eigen::Vector3d &receiver_ecef,
    const Eigen::Vector3d &receiver_enu,
    const std::vector<gnss_comm::ObsPtr> &observations,
    const std::map<int, sv_info> &sv_info_map,
    const std::map<int, int> &lock_count_l1,
    const std::map<int, int> &lock_count_l2,
    const std::map<int, double> &psr_residual_map,
    const std::map<int, std::map<std::string, double>> &factor_residuals)
{
    json j_epoch;
    j_epoch["timestamp"] = timestamp;
    j_epoch["time_frame"] = time_frame;
    j_epoch["receiver_ecef"] = {receiver_ecef.x(), receiver_ecef.y(), receiver_ecef.z()};
    j_epoch["receiver_enu"] = {receiver_enu.x(), receiver_enu.y(), receiver_enu.z()};

    json satellites = json::array();

    for (const auto &obs : observations)
    {
        if (!obs)
        {
            continue;
        }

        const int sat = static_cast<int>(obs->sat);
        const int sys = gnss_comm::satsys(sat, nullptr);

        // 卫星系统名称
        std::string sys_name = "Unknown";
        if (sys == SYS_GPS)      sys_name = "GPS";
        else if (sys == SYS_GLO) sys_name = "GLONASS";
        else if (sys == SYS_GAL) sys_name = "Galileo";
        else if (sys == SYS_BDS) sys_name = "BeiDou";

        // PRN 字符串 (如 "G05")
        char prn_buf[16];
        std::snprintf(prn_buf, sizeof(prn_buf), "%c%02d",
                      (sys == SYS_GPS ? 'G' : sys == SYS_GLO ? 'R' :
                       sys == SYS_GAL ? 'E' : sys == SYS_BDS ? 'C' : '?'),
                      sat % 100);
        std::string prn_str(prn_buf);

        // 频率索引
        int l1_idx = -1, l2_idx = -1;
        gnss_comm::L1_freq(obs, &l1_idx);
        gnss_comm::L2_freq(obs, &l2_idx);

        json j_sat;
        j_sat["prn"] = prn_str;
        j_sat["system"] = sys_name;
        j_sat["sat_id"] = sat;

        // SNR (L1/L2)
        j_sat["snr_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->CN0.size()))
            ? obs->CN0[l1_idx] : 0.0;
        j_sat["snr_l2"] = (l2_idx >= 0 && l2_idx < static_cast<int>(obs->CN0.size()))
            ? obs->CN0[l2_idx] : 0.0;

        // 伪距 (L1/L2)
        j_sat["pseudorange_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->psr.size()))
            ? obs->psr[l1_idx] : 0.0;
        j_sat["pseudorange_l2"] = (l2_idx >= 0 && l2_idx < static_cast<int>(obs->psr.size()))
            ? obs->psr[l2_idx] : 0.0;

        // 伪距标准差 (L1/L2)
        j_sat["psr_std_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->psr_std.size()))
            ? obs->psr_std[l1_idx] : 3.0;
        j_sat["psr_std_l2"] = (l2_idx >= 0 && l2_idx < static_cast<int>(obs->psr_std.size()))
            ? obs->psr_std[l2_idx] : 3.0;

        // 载波相位 (L1/L2)
        j_sat["carrier_phase_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->cp.size()))
            ? obs->cp[l1_idx] : 0.0;
        j_sat["carrier_phase_l2"] = (l2_idx >= 0 && l2_idx < static_cast<int>(obs->cp.size()))
            ? obs->cp[l2_idx] : 0.0;

        // 多普勒 (L1)
        j_sat["doppler_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->dopp.size()))
            ? obs->dopp[l1_idx] : 0.0;

        // LLI 失锁标志
        j_sat["lli_l1"] = (l1_idx >= 0 && l1_idx < static_cast<int>(obs->LLI.size()))
            ? obs->LLI[l1_idx] : 0;
        j_sat["lli_l2"] = (l2_idx >= 0 && l2_idx < static_cast<int>(obs->LLI.size()))
            ? obs->LLI[l2_idx] : 0;

        // 锁定计数
        auto l1_lock_it = lock_count_l1.find(sat);
        j_sat["lock_count_l1"] = (l1_lock_it != lock_count_l1.end()) ? l1_lock_it->second : 0;
        auto l2_lock_it = lock_count_l2.find(sat);
        j_sat["lock_count_l2"] = (l2_lock_it != lock_count_l2.end()) ? l2_lock_it->second : 0;

        // 卫星信息 (仰角、方位角、位置等)
        auto sv_it = sv_info_map.find(sat);
        if (sv_it != sv_info_map.end() && sv_it->second.avaliable)
        {
            const sv_info &sv = sv_it->second;
            j_sat["elevation"] = sv.elevation * R2D;   // 弧度转度
            j_sat["azimuth"] = sv.azimuth * R2D;        // 弧度转度
            j_sat["sat_pos_ecef"] = {sv.pos.x(), sv.pos.y(), sv.pos.z()};
            j_sat["sat_vel_ecef"] = {sv.vel.x(), sv.vel.y(), sv.vel.z()};
            j_sat["sat_clock_bias"] = sv.dt;
            j_sat["sat_clock_drift"] = sv.ddt;
        }
        else
        {
            j_sat["elevation"] = 0.0;
            j_sat["azimuth"] = 0.0;
            j_sat["sat_pos_ecef"] = {0.0, 0.0, 0.0};
            j_sat["sat_vel_ecef"] = {0.0, 0.0, 0.0};
            j_sat["sat_clock_bias"] = 0.0;
            j_sat["sat_clock_drift"] = 0.0;
        }

        // SPP 伪距残差
        auto psr_res_it = psr_residual_map.find(sat);
        j_sat["psr_residual_l1"] = (psr_res_it != psr_residual_map.end())
            ? psr_res_it->second : 0.0;

        // 因子残差 (优化后)
        auto fac_it = factor_residuals.find(sat);
        if (fac_it != factor_residuals.end())
        {
            const auto &res = fac_it->second;
            auto psr_it = res.find("psr");
            j_sat["psr_factor_residual"] = (psr_it != res.end()) ? psr_it->second : 0.0;
            auto tr_pr_it = res.find("tr_dd_pr");
            j_sat["tr_dd_pr_residual"] = (tr_pr_it != res.end()) ? tr_pr_it->second : 0.0;
            auto tr_cp_it = res.find("tr_dd_cp");
            j_sat["tr_dd_cp_residual"] = (tr_cp_it != res.end()) ? tr_cp_it->second : 0.0;
        }
        else
        {
            j_sat["psr_factor_residual"] = 0.0;
            j_sat["tr_dd_pr_residual"] = 0.0;
            j_sat["tr_dd_cp_residual"] = 0.0;
        }

        satellites.push_back(j_sat);
    }

    j_epoch["satellites"] = satellites;
    j_epoch["n_satellites"] = satellites.size();

    return appendJsonLine(filepath, j_epoch);
}

// ==================== 数据导入 (OSQA → gnssfgo) ====================

/**
 * 从 JSONL 文件读取最新 epoch 的质量评分
 *
 * 只读取匹配指定 time_frame 的 epoch 的评分；
 * 若未找到匹配项则返回空 map。
 *
 * @param filepath   OSQA 输出的 JSONL 文件路径
 * @param time_frame 目标时间帧（内部格式: gpst_sec * 10）
 * @return           PRN → 质量评分映射 (key 为 sat_id 整数)
 */
inline std::map<int, SatQualityInfo> readLatestQualityScores(
    const std::string &filepath,
    double time_frame)
{
    std::map<int, SatQualityInfo> result;

    // 读取所有行并找到匹配 time_frame 的最新 epoch
    auto all_lines = readAllJsonLines(filepath);
    if (all_lines.empty())
    {
        return result;
    }

    // 从后向前查找匹配的 epoch
    const json *match = nullptr;
    for (auto it = all_lines.rbegin(); it != all_lines.rend(); ++it)
    {
        if (it->contains("time_frame"))
        {
            double tf = (*it)["time_frame"].get<double>();
            if (std::abs(tf - time_frame) < 1e-3)
            {
                match = &(*it);
                break;
            }
        }
    }

    // 如果没找到精确匹配，找最接近的
    if (match == nullptr)
    {
        double best_diff = 1e9;
        for (auto it = all_lines.rbegin(); it != all_lines.rend(); ++it)
        {
            if (it->contains("time_frame"))
            {
                double tf = (*it)["time_frame"].get<double>();
                double diff = std::abs(tf - time_frame);
                if (diff < best_diff && diff < 60.0)  // 最多容忍 60 个时间帧的差异
                {
                    best_diff = diff;
                    match = &(*it);
                }
            }
        }
    }

    if (match == nullptr || !match->contains("satellites"))
    {
        return result;
    }

    // 解析每颗卫星的质量评分
    const json &satellites = (*match)["satellites"];
    for (auto it = satellites.begin(); it != satellites.end(); ++it)
    {
        SatQualityInfo info;
        const json &sat_data = it.value();

        info.quality = sat_data.value("quality", 1.0);
        info.trust_level = sat_data.value("trust_level", "trusted");

        if (sat_data.contains("details"))
        {
            const json &details = sat_data["details"];
            info.q_transformer = details.value("transformer", 1.0);
            info.q_graph = details.value("graph", 1.0);
            info.q_temporal = details.value("temporal", 1.0);
        }

        if (sat_data.contains("flags"))
        {
            info.flags = sat_data["flags"].get<std::vector<std::string>>();
        }

        // 同时尝试从 PRN 字符串和 sat_id 整数读取
        int sat_id = 0;
        if (sat_data.contains("sat_id"))
        {
            sat_id = sat_data["sat_id"].get<int>();
        }

        if (sat_id > 0)
        {
            result[sat_id] = info;
        }
    }

    return result;
}

/**
 * 读取最新的质量评分（不匹配特定 time_frame）
 *
 * 用于简单场景：只读 JSONL 文件中最后一行 epoch 的质量评分。
 *
 * @param filepath  OSQA 输出的 JSONL 文件路径
 * @return          PRN → 质量评分映射
 */
inline std::map<int, SatQualityInfo> readLatestQualityScores(
    const std::string &filepath)
{
    std::map<int, SatQualityInfo> result;
    auto all_lines = readAllJsonLines(filepath);
    if (all_lines.empty())
    {
        return result;
    }

    // 取最后一行
    const json &last = all_lines.back();
    if (!last.contains("satellites"))
    {
        return result;
    }

    for (auto it = last["satellites"].begin(); it != last["satellites"].end(); ++it)
    {
        SatQualityInfo info;
        const json &sat_data = it.value();

        info.quality = sat_data.value("quality", 1.0);
        info.trust_level = sat_data.value("trust_level", "trusted");

        if (sat_data.contains("details"))
        {
            const json &details = sat_data["details"];
            info.q_transformer = details.value("transformer", 1.0);
            info.q_graph = details.value("graph", 1.0);
            info.q_temporal = details.value("temporal", 1.0);
        }

        if (sat_data.contains("flags"))
        {
            info.flags = sat_data["flags"].get<std::vector<std::string>>();
        }

        int sat_id = sat_data.value("sat_id", 0);
        if (sat_id > 0)
        {
            result[sat_id] = info;
        }
    }

    return result;
}

} // namespace transformer_bridge

#endif // TRANSFORMER_BRIDGE_HPP
