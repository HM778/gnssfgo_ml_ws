#ifndef D2R
#define D2R         M_PI/180.0
#endif
#ifndef R2D
#define R2D         180.0/M_PI
#endif

#ifndef DATATYPE_H
#define DATATYPE_H
#define OMGE_       7.2921151467E-5
#define CLIGHT_     299792458.0

struct DDMeasurement
{
    /* the master satellite is found from user end*/
    gnss_comm::ObsPtr u_master_SV;
    gnss_comm::ObsPtr u_iSV;

    gnss_comm::ObsPtr r_master_SV;
    gnss_comm::ObsPtr r_iSV;

    int freq_idx;   //新增了L1/L2频段数据处理后，必须要区分频段

    double var_pr;
    double var_cp;
};

// 跨历元双差因子信息
struct TRDDMeasurement : DDMeasurement
{
    int prev_epoch_index;
    int curr_epoch_index;
    double prev_time;
    double curr_time;

    double pred_move; // 预测的移动距离（用于时间折扣函数）
    double tr_score;
};

/**
 * 卫星解析信息
 * 可用性
 * 卫星编码
 * 卫星系统
 * 位置、速度（ECEF）
 * 钟差、钟速
 * 方位角、仰角
 * 波长
 */
struct sv_info
{
    bool avaliable = false;
    int sat_prn = 0;
    int sys = 0;
    Eigen::Vector3d pos = Eigen::Vector3d::Zero(); // ECEF
    Eigen::Vector3d vel = Eigen::Vector3d::Zero(); // ECEF
    double dt = 0.0;
    double ddt = 0.0;
    double azimuth = 0.0;
    double elevation = 0.0;
    double lamda_l1 = 0.0;    // L1 wavelength (m)
    double lamda_l2 = 0.0;    // L2 wavelength (m)
    double freq_l1 = 0.0;     // L1 frequency (Hz)
    double freq_l2 = 0.0;     // L2 frequency (Hz)
};
#endif // DATATYPE_H