#ifndef _HOTSTART_HPP_
#define _HOTSTART_HPP_
#include <fstream>
#include <iostream>
#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/GnssEphemMsg.h>
#include <gnss_comm/GnssGloEphemMsg.h>
#include <gnss_comm/GnssTimeMsg.h>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <ros/ros.h>

bool loadEphemsMap(ros::NodeHandle& nh, const std::string& filename) 
{
    //模拟硬件热启动，加载已知的卫星数据，加快解算速度

    ros::Publisher load_ephem_pub;
    load_ephem_pub = nh.advertise<gnss_comm::GnssEphemMsg>("/ublox_driver/ephem",10);

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }
    
    std::string line; 
    bool reading_eph = false;
    int loaded_count = 0;
    gnss_comm::GnssEphemMsg eph_msg;  //加载到gnss_comm可用的星历数据库中，发布
    gnss_comm::Ephem current_eph;
    
    while (std::getline(file, line)) {
        // 跳过空行和注释
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // 检查星历开始标记
        if (line.find("BEGIN_EPHEMERIS SAT") == 0) {
            if (reading_eph) {
                // 保存前一个星历
                loaded_count++;
            }
            reading_eph = true;
            continue;
        }
        
        // 检查星历结束标记
        if (line == "---END_EPHEMERIS---") {
            if (reading_eph && eph_msg.sat != 0) {
                load_ephem_pub.publish(eph_msg);
                loaded_count++;
            }
            reading_eph = false;
            continue;
        }
        
        if (!reading_eph) continue;
        
        // 解析键值对
        size_t colon_pos = line.find(':');
        if (colon_pos == std::string::npos) {
            continue;
        }
        
        std::string key = line.substr(0, colon_pos);
        std::string value_str = line.substr(colon_pos + 1);
        
        // 去除首尾空格
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value_str.erase(0, value_str.find_first_not_of(" \t"));
        value_str.erase(value_str.find_last_not_of(" \t") + 1);
        
        // 根据关键字解析值
        try {
            if (key == "SAT") eph_msg.sat = std::stoi(value_str);
            else if (key == "IODE") eph_msg.iode = std::stoi(value_str);
            else if (key == "IODC") eph_msg.iodc = std::stoi(value_str);
            else if (key == "URA") eph_msg.ura = std::stod(value_str);
            else if (key == "SVH") eph_msg.health = std::stoi(value_str);
            else if (key == "WEEK") eph_msg.week = std::stoi(value_str);
            else if (key == "CODE") eph_msg.code = std::stoi(value_str);
            else if (key == "TOE_WEEK") eph_msg.toe.week = std::stod(value_str);
            else if (key == "TOE_TOW") eph_msg.toe.tow = std::stod(value_str);
            else if (key == "TOC_WEEK") eph_msg.toc.week = std::stod(value_str);
            else if (key == "TOC_TOW") eph_msg.toc.tow = std::stod(value_str);
            else if (key == "TTR_WEEK") eph_msg.ttr.week = std::stod(value_str);
            else if (key == "TTR_TOW") eph_msg.ttr.tow = std::stod(value_str);
            else if (key == "A") eph_msg.A = std::stod(value_str);
            else if (key == "E") eph_msg.e = std::stod(value_str);
            else if (key == "I0") eph_msg.i0 = std::stod(value_str);
            else if (key == "OMG0") eph_msg.OMG0 = std::stod(value_str);
            else if (key == "OMG") eph_msg.omg = std::stod(value_str);
            else if (key == "M0") eph_msg.M0 = std::stod(value_str);
            else if (key == "DELN") eph_msg.delta_n = std::stod(value_str);
            else if (key == "OMGD") eph_msg.OMG_dot = std::stod(value_str);
            else if (key == "IDOT") eph_msg.i_dot = std::stod(value_str);
            else if (key == "CRC") eph_msg.crc = std::stod(value_str);
            else if (key == "CRS") eph_msg.crs = std::stod(value_str);
            else if (key == "CUC") eph_msg.cuc = std::stod(value_str);
            else if (key == "CUS") eph_msg.cus = std::stod(value_str);
            else if (key == "CIC") eph_msg.cic = std::stod(value_str);
            else if (key == "CIS") eph_msg.cis = std::stod(value_str);
            else if (key == "F0") eph_msg.af0 = std::stod(value_str);
            else if (key == "F1") eph_msg.af1 = std::stod(value_str);
            else if (key == "F2") eph_msg.af2 = std::stod(value_str);
            else if (key == "ADOT") eph_msg.A_dot = std::stod(value_str);
            else if (key == "NDOT") eph_msg.n_dot = std::stod(value_str);
            else if (key == "TGD0") eph_msg.tgd0 = std::stod(value_str);
            else if (key == "TGD1") eph_msg.tgd1 = std::stod(value_str);
            
        } catch (const std::exception& e) {
            std::cerr << "解析错误 - 行: " << line << " - " << e.what() << std::endl;
        }
    }

    file.close();
    std::cout << "load " << loaded_count << " ephems from " << filename << std::endl;
    return loaded_count > 0;
}

bool saveEphemsMap(const std::map<int, gnss_comm::EphemPtr>& eph_map, const std::string& filename,bool append = false) 
{
    std::ios_base::openmode mode;
    if (append) {
        mode = std::ios::app;  // 追加模式，不清除内容
    } else {
        mode = std::ios::out | std::ios::trunc;  // 输出模式，清除内容
    }
    
    std::ofstream file(filename, mode);
    
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }
    
    // 设置高精度输出
    file << std::setprecision(15) << std::scientific;
    
    // 如果不是追加模式，写入文件头
    if (!append) {
        file << "# GPS Ephemeris Data Batch Export\n";
        file << "# Total satellites: " << eph_map.size() << "\n";
        file << "# Format: KEY VALUE\n";
        file << "# Separator: ---END_EPHEMERIS---\n\n";
    }
    
    int saved_count = 0;
    for (const auto& kv : eph_map) {
        const int sat = kv.first;
        const auto& eph = kv.second;
        file << "BEGIN_EPHEMERIS SAT" << "\n";
        // 整型数据直接保存，不需要小数和高精度
        file << "SAT: " << eph->sat << "\n";
        file << "IODE: " << eph->iode << "\n";
        file << "IODC: " << eph->iodc << "\n";
        file << "SVH: " << eph->health << "\n";
        file << "WEEK: " << eph->week << "\n";
        file << "CODE: " << eph->code << "\n";
        
        // 时间相关数据：保留3位小数足够
        file << std::fixed << std::setprecision(3);

        gnss_comm::GnssTimeMsg time_msg;
        time_msg.tow = gnss_comm::time2gpst(eph->toe,&time_msg.week);
        file << "TOE_WEEK: " << time_msg.week << "\n";
        file << "TOE_TOW: " << time_msg.tow << "\n";
        time_msg.tow = gnss_comm::time2gpst(eph->toc,&time_msg.week);
        file << "TOC_WEEK: " << time_msg.week << "\n";
        file << "TOC_TOW: " << time_msg.tow << "\n";
        time_msg.tow = gnss_comm::time2gpst(eph->ttr,&time_msg.week);
        file << "TTR_WEEK: " << time_msg.week << "\n";
        file << "TTR_TOW: " << time_msg.tow << "\n";
        
        
        // 轨道参数：根据数值大小设置合适的精度
        file << std::scientific << std::setprecision(6);
        file << "A: " << eph->A << "\n";
        file << "E: " << eph->e << "\n";
        file << "I0: " << eph->i0 << "\n";
        file << "OMG0: " << eph->OMG0 << "\n";
        file << "OMG: " << eph->omg << "\n";
        file << "M0: " << eph->M0 << "\n";
        
        // 小数值参数：科学计数法，6位精度足够
        file << "DELN: " << eph->delta_n << "\n";
        file << "OMGD: " << eph->OMG_dot << "\n";
        file << "IDOT: " << eph->i_dot << "\n";
        
        // 摄动参数：根据数值范围设置精度
        if (std::abs(eph->crc) > 1.0 || std::abs(eph->crs) > 1.0) {
            file << std::fixed << std::setprecision(3);
        } else {
            file << std::scientific << std::setprecision(6);
        }
        file << "CRC: " << eph->crc << "\n";
        file << "CRS: " << eph->crs << "\n";
        
        // 小摄动参数：科学计数法
        file << std::scientific << std::setprecision(6);
        file << "CUC: " << eph->cuc << "\n";
        file << "CUS: " << eph->cus << "\n";
        file << "CIC: " << eph->cic << "\n";
        file << "CIS: " << eph->cis << "\n";
        
        // 钟差参数：科学计数法
        file << "F0: " << eph->af0 << "\n";
        file << "F1: " << eph->af1 << "\n";
        file << "F2: " << eph->af2 << "\n";
        
        // 变化率参数
        file << "ADOT: " << eph->A_dot << "\n";
        file << "NDOT: " << eph->n_dot << "\n";
        
        // 群延迟参数：科学计数法
        file << "TGD0: " << eph->tgd[0] << "\n";
        file << "TGD1: " << eph->tgd[1] << "\n";
        
        // 精度参数：固定小数
        file << std::fixed << std::setprecision(1);
        file << "URA: " << eph->ura << "\n";
        
        file << "---END_EPHEMERIS---\n\n";
        saved_count++;
    }
    
    file.close();
    std::cout << "save " << saved_count << " ephems to " << filename 
            << " (mode: " << (append ? "append" : "overwrite") << ")" << std::endl;
    return true;
}

#endif // _HOTSTART_HPP_