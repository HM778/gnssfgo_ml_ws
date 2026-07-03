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

    double var_pr;
    double var_cp;
};

struct TRfeature{
    //保存TR因子相关信息的结构体
    double min_SNR;    
    double max_SNR;   
    double min_elev;    
    double max_elev;    
    double delta_elev;  
    double delta_azm;      
    double gap_epoch; // 时间间隔

    double residual; // TR双差因子残差
    double pred_move; // 预测的移动距离（用于时间折扣函数）
};

struct TRRTKMeasurement : DDMeasurement
{
    /* the master satellite is found from user end*/
    int prev_epoch_index;
    int curr_epoch_index;
    double prev_time;
    double curr_time;
   
    TRfeature tr_feature; // 关联的TR因子信息
    double pred_move; // 预测的移动距离（用于时间折扣函数）
    double tr_score;
};

/**
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
    double lamda = 0.0;       // L1 wavelength (m)
    double lamda_l2 = 0.0;    // L2 wavelength (m)
    double freq_l1 = 0.0;     // L1 frequency (Hz)
    double freq_l2 = 0.0;     // L2 frequency (Hz)
};
#endif // DATATYPE_H