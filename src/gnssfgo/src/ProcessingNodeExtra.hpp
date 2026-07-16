#include <iostream>
#include <string>  
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>

#include <math.h>
#include <time.h>
#include <algorithm>

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <ros/ros.h>
#include <ros/package.h>

#include <ros/callback_queue.h>

#include <geometry_msgs/Point32.h>
#include <stdio.h>
#include <queue>
#include <map>
#include <queue>
#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/NavSatFix.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <stdarg.h>
#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_spp_extra.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_constant.hpp>

#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/GnssEphemMsg.h>
#include <gnss_comm/GnssGloEphemMsg.h>
#include <gnss_comm/gnss_ros.hpp>


#include "../include/datatype.h"
#include "../include/tools/hotstart.hpp"

class ProcessingNodeExtra
{
public:
    ros::NodeHandle nh;

    ros::CallbackQueue cbq_range_;
    ros::CallbackQueue cbq_ephem_;
    ros::NodeHandle nh_range_;
    ros::NodeHandle nh_ephem_;

    std::unique_ptr<ros::AsyncSpinner> range_spinner_;
    std::unique_ptr<ros::AsyncSpinner> ephem_spinner_;

    std::mutex m_gnss_raw_mux;

    ros::Subscriber sub_range_meas_;
    ros::Subscriber sub_ephem_;
    ros::Subscriber sub_glo_ephem_;
    ros::Subscriber sub_origin_node;

    std::vector<gnss_comm::ObsPtr> meas;
    std::vector<gnss_comm::ObsPtr> L1_meas;
    std::vector<gnss_comm::ObsPtr> L2_meas;
    nav_msgs::Odometry dop_vel_meas;
    std::vector<gnss_comm::EphemBasePtr> ephem_array;
    std::vector<gnss_comm::GloEphemPtr> gloephem_array;
    std::map<int, gnss_comm::EphemPtr> current_ephem_map;
    std::map<int, gnss_comm::GloEphemPtr> current_glo_ephem_map;
    std::map<int, sv_info> current_sv_info_map;

    std::atomic<bool> hasNewData{false};

    // ── CPU optimization: condition variable + pre-allocated buffers ──
    std::condition_variable cv_new_data_;
    std::mutex cv_wait_mutex_;

    // Pre-allocated buffers to avoid repeated malloc/free in hot callbacks
    std::vector<gnss_comm::ObsPtr> meas_buf_;
    std::vector<gnss_comm::ObsPtr> L1_meas_buf_;
    std::vector<gnss_comm::ObsPtr> L2_meas_buf_;

    // Backpressure: skip heavy SPP/Doppler computation when the optimizer
    // hasn't consumed the previous epoch yet, preventing callback pile-up.
    static constexpr int kMaxCallbackQueued = 2;
    std::atomic<int> pending_cb_count_{0};

    //  settings
    // 星历加载-热启动
    bool HOTSTART_ENABLE = false;
    bool load_ephem_flag = false;
    std::string filename = ros::package::getPath("gnssfgo") +  "/../ublox_driver_Ephems.txt";
    std::string output_filename = ros::package::getPath("gnssfgo") +  "/../outputfile.txt";
    //是否只处理GPS和BDS卫星 (default: false = all 4 systems)
    bool PART_SATSYS = false;
    // windowsize
    int windowSize = 0;
    // 是否启用边缘化
    bool marginal_enable = true;

    //保证时间同步
    gnss_comm::gtime_t current_sys_time;            //后处理时设置为星历时间
    double current_gpst_sec = -1.0;

    // 结果输出设置
    bool enu_ref_set = false;
    Eigen::Matrix<double, 3,1> enu_ref_llh, enu_ref_ecef;
    Eigen::Matrix<double, 3,1> latest_pos_ecef, latest_pos_llh, latest_pos_enu;  //最小二乘结果
    std::string WLS_LLH_TOPIC = "/gnssfgo/psr_llh_latest";
    std::string WLS_ENU_TOPIC = "/gnssfgo/psr_enu_latest";
    std::string FGO_LLH_TOPIC = "/gnssfgo/fgo_llh";
    std::string FGO_ENU_TOPIC = "/gnssfgo/fgo_enu";
    
    ros::Publisher pub_psr_llh_latest;
    ros::Publisher pub_psr_enu_latest;
    ros::Publisher pub_fgo_llh_latest;
    ros::Publisher pub_fgo_enu_latest;

    void StartSpinners()
    {
        range_spinner_ = std::make_unique<ros::AsyncSpinner>(1, &cbq_range_);
        ephem_spinner_ = std::make_unique<ros::AsyncSpinner>(1, &cbq_ephem_);
        range_spinner_->start();
        ephem_spinner_->start();
    }

    // ── CPU optimization: blocking wait to replace busy-poll ──
    // Call this instead of polling hasNewData in a tight loop.
    // Returns true if new data is available, false on timeout.
    // timeout_ms: how long to block waiting for data (default 50ms).
    //             Use -1 to block indefinitely until data arrives.
    inline bool waitForNewData(int timeout_ms = 50)
    {
        std::unique_lock<std::mutex> lk(cv_wait_mutex_);
        if (timeout_ms < 0)
        {
            cv_new_data_.wait(lk, [this] {
                return hasNewData.load(std::memory_order_acquire);
            });
            return true;
        }
        return cv_new_data_.wait_for(lk, std::chrono::milliseconds(timeout_ms),
                                     [this] {
                                         return hasNewData.load(std::memory_order_acquire);
                                     });
    }

    void InitialSubTopics()
    {
        nh_range_ = ros::NodeHandle(nh);
        nh_ephem_ = ros::NodeHandle(nh);
        nh_range_.setCallbackQueue(&cbq_range_);
        nh_ephem_.setCallbackQueue(&cbq_ephem_);
         /* subscriber of three topics -- UBLOX_DRIVER */
        sub_range_meas_ = nh_range_.subscribe<gnss_comm::GnssMeasMsg>("/ublox_driver/range_meas", 20, &ProcessingNodeExtra::range_meas_cb, this);
        sub_ephem_ = nh_ephem_.subscribe<gnss_comm::GnssEphemMsg>("/ublox_driver/ephem", 200, &ProcessingNodeExtra::ephem_cb, this);
        sub_glo_ephem_ = nh_ephem_.subscribe<gnss_comm::GnssGloEphemMsg>("/ublox_driver/glo_ephem", 200, &ProcessingNodeExtra::gephem_cb, this);
        sub_origin_node = nh_range_.subscribe<sensor_msgs::NavSatFix>("/ublox_driver/receiver_lla", 1000, &ProcessingNodeExtra::origin_node_cb, this);
    }

    void InitialPubTopics()
    {
        pub_psr_llh_latest = nh.advertise<sensor_msgs::NavSatFix>(WLS_LLH_TOPIC, 100); 
        pub_psr_enu_latest = nh.advertise<nav_msgs::Odometry>(WLS_ENU_TOPIC, 100); 
        
        pub_fgo_llh_latest = nh.advertise<sensor_msgs::NavSatFix>(FGO_LLH_TOPIC, 100);
        pub_fgo_enu_latest = nh.advertise<nav_msgs::Odometry>(FGO_ENU_TOPIC, 100); 
    }

    void loadParams()
    {
        nh = ros::NodeHandle("~"); 
        nh.param<bool>("hotstart_enable", HOTSTART_ENABLE, false);
        nh.param<bool>("GPS/BDS_only", PART_SATSYS, false);
        nh.param<std::string>("wls_llh_topic", WLS_LLH_TOPIC, "/gnssfgo/psr_llh_latest");
        nh.param<std::string>("wls_enu_topic", WLS_ENU_TOPIC, "/gnssfgo/psr_enu_latest");
        nh.param<std::string>("fgo_llh_topic", FGO_LLH_TOPIC, "/gnssfgo/fgo_llh");
        nh.param<std::string>("fgo_enu_topic", FGO_ENU_TOPIC, "/gnssfgo/fgo_enu");
        nh.param<std::string>("output_file", output_filename, "outputfile.txt");
        output_filename = ros::package::getPath("gnssfgo") + "/../" + output_filename;
        nh.param<int>("windowSize", windowSize, 2);
        nh.param<bool>("marginal_enable", marginal_enable, true);
        ROS_INFO("Parameters loaded: \n  hotstart_enable=%s, \n  GPS/BDS_only=%s, \n wls_llh_topic=%s, \n wls_enu_topic=%s, \n fgo_llh_topic=%s, \n fgo_enu_topic=%s \n windowSize:%d \n save to %s",
                 HOTSTART_ENABLE ? "true" : "false",
                  PART_SATSYS ? "true" : "false", 
                  WLS_LLH_TOPIC.c_str(), 
                  WLS_ENU_TOPIC.c_str(), 
                  FGO_LLH_TOPIC.c_str(), 
                  FGO_ENU_TOPIC.c_str(),
                  windowSize,
                  output_filename.c_str());
    }


    bool savetofile(Eigen::Vector3d lla, Eigen::Vector3d enu, Eigen::Vector3d cov_enu, double gpst_sec, double opt_time)
    {
        std::ofstream outfile(output_filename, std::ios::out | std::ios::app);
        if (!outfile.is_open())
        {
            ROS_ERROR("Failed to open output file: %s", output_filename.c_str());
            return false;
        }

        // Write header if file is new
        if (outfile.tellp() == 0)
        {
            outfile << "timestamp,latitude,longitude,altitude,enu_x,enu_y,enu_z,opt_time\n";
        }

        // Write data
        outfile << std::fixed << std::setprecision(6) << gpst_sec << ","
                << lla[0] << ","
                << lla[1] << ","
                << lla[2] << ","
                << enu[0] << ","
                << enu[1] << ","
                << enu[2] << ","
                << cov_enu[0] << "," 
                << cov_enu[1] << "," 
                << cov_enu[2] << ","
                << opt_time << "\n";

        outfile.close();
        return true;
    }

    void range_meas_cb(const gnss_comm::GnssMeasMsgConstPtr &meas_msg)
    {
        // ── Backpressure: if the optimizer is falling behind (>kMaxCallbackQueued
        //     epochs unconsumed), skip heavy SPP/Doppler computation in this callback
        //     to shed CPU load. Raw measurements are still stored so no data is lost.
        //     This prevents callback pile-up from saturating the CPU. ──
        const int queued = pending_cb_count_.fetch_add(1, std::memory_order_acq_rel);
        const bool skip_heavy = (queued >= kMaxCallbackQueued);

        hasNewData.store(false, std::memory_order_release);

        // ── Reuse pre-allocated buffers (swap + clear instead of re-allocate) ──
        const size_t n_obs = meas_msg->meas.size();
        meas_buf_.clear();
        L1_meas_buf_.clear();
        L2_meas_buf_.clear();
        meas_buf_.reserve(n_obs);
        L1_meas_buf_.reserve(n_obs);
        L2_meas_buf_.reserve(n_obs);

        for (size_t i = 0; i < n_obs; ++i)
        {
            const gnss_comm::GnssObsMsg& obs_msg = meas_msg->meas[i];
            gnss_comm::ObsPtr obs(new gnss_comm::Obs());
            obs->time       = gnss_comm::gpst2time(obs_msg.time.week, obs_msg.time.tow);
            obs->sat        = obs_msg.sat;
            obs->freqs      = obs_msg.freqs;
            obs->CN0        = obs_msg.CN0;
            obs->LLI        = obs_msg.LLI;
            obs->code       = obs_msg.code;
            obs->psr        = obs_msg.psr;
            obs->psr_std    = obs_msg.psr_std;
            obs->cp         = obs_msg.cp;
            obs->cp_std     = obs_msg.cp_std;
            obs->dopp       = obs_msg.dopp;
            obs->dopp_std   = obs_msg.dopp_std;
            obs->status     = obs_msg.status;

            if (PART_SATSYS)
            {
                const int sys = gnss_comm::satsys(obs->sat, nullptr);
                if (sys != SYS_BDS && sys != SYS_GPS)
                {
                    continue;
                }
            }

            meas_buf_.push_back(obs);

            // ── Combined L1+L2 filtering in a single pass (was two separate loops) ──
            int l1_idx = -1, l2_idx = -1;
            L1_freq(obs, &l1_idx);
            L2_freq(obs, &l2_idx);
            if (l1_idx >= 0) L1_meas_buf_.push_back(obs);
            if (l2_idx >= 0) L2_meas_buf_.push_back(obs);
        }

        if (meas_buf_.empty())
        {
            pending_cb_count_.fetch_sub(1, std::memory_order_acq_rel);
            return;
        }

        const gnss_comm::gtime_t msg_time = meas_buf_[0]->time;
        const double msg_gpst_sec = static_cast<double>(meas_msg->meas[0].time.week) * WEEK_SECONDS
                                    + meas_msg->meas[0].time.tow;

        // ── Swap buffers under lock (O(1) instead of O(n) copy) ──
        std::vector<gnss_comm::EphemBasePtr> local_ephem_array;
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            meas.swap(meas_buf_);
            L1_meas.swap(L1_meas_buf_);
            L2_meas.swap(L2_meas_buf_);
            current_sys_time = msg_time;
            current_gpst_sec = msg_gpst_sec;
            local_ephem_array = ephem_array;
        }

        if (HOTSTART_ENABLE && !load_ephem_flag)
        {
            loadEphemsMap(nh, filename);
            load_ephem_flag = true;
        }

        const size_t n_l1 = L1_meas.size();
        const size_t n_eph = local_ephem_array.size();
        const bool can_doppler = (n_l1 >= 4) && (n_eph >= 4);
        const bool can_psr     = (n_l1 >= 6) && (n_eph >= 6);

        // ── Throttle: skip CPU-heavy SPP & Doppler when the optimizer is behind ──
        if (skip_heavy)
        {
            // Only store raw data; skip SPP, Doppler, and publishing.
            // The optimizer will still get the raw measurements via meas/L1_meas/...,
            // and init values will come from the previous epoch's optimized state.
            ROS_DEBUG_THROTTLE(1.0, "CB overload — skipping SPP/Doppler (queued=%d)", queued);
            hasNewData.store(true, std::memory_order_release);
            cv_new_data_.notify_one();
            pending_cb_count_.fetch_sub(1, std::memory_order_acq_rel);
            return;
        }

        Eigen::Vector3d pos_ecef = Eigen::Vector3d::Zero();
        Eigen::Vector3d pos_llh  = Eigen::Vector3d::Zero();
        bool pos_ok = false;

        if (can_psr)
        {
            const Eigen::Matrix<double, 7, 1> psr_result =
                gnss_comm::psr_pos_extra(meas, local_ephem_array);
            pos_ecef = psr_result.head<3>();
            if (pos_ecef.norm() > 1e-3)
            {
                pos_llh = gnss_comm::ecef2geo(pos_ecef);
                pos_ok  = posValid(pos_llh);
            }
        }

        if (!enu_ref_set && pos_ok)
        {
            enu_ref_llh  = pos_llh;
            enu_ref_ecef = gnss_comm::geo2ecef(enu_ref_llh);
            enu_ref_set  = true;
        }

        Eigen::Vector3d pos_enu = Eigen::Vector3d::Zero();
        if (pos_ok)
        {
            pos_enu = gnss_comm::ecef2enu(enu_ref_llh, pos_ecef - enu_ref_ecef);
        }

        // ── Doppler velocity ──
        if (!enu_ref_set || !can_doppler)
        {
            hasNewData.store(true, std::memory_order_release);
            cv_new_data_.notify_one();
            pending_cb_count_.fetch_sub(1, std::memory_order_acq_rel);
            return;
        }

        const Eigen::Vector4d doppler_est =
            gnss_comm::dopp_vel_extra(meas, local_ephem_array, enu_ref_ecef);
        const Eigen::Vector3d vel_xyz(doppler_est[0], doppler_est[1], doppler_est[2]);

        nav_msgs::Odometry local_dop_meas;
        bool doppler_updated = false;
        if (velValid(vel_xyz))
        {
            const Eigen::MatrixXd covM = gnss_comm::dopp_cov_get();
            const bool cov_ok = (covM.rows() >= 3 && covM.cols() >= 3);

            local_dop_meas.header.stamp    = ros::Time::now();
            local_dop_meas.header.frame_id = "map";
            local_dop_meas.child_frame_id  = "map";
            local_dop_meas.twist.twist.linear.x = doppler_est[0];
            local_dop_meas.twist.twist.linear.y = doppler_est[1];
            local_dop_meas.twist.twist.linear.z = doppler_est[2];
            local_dop_meas.twist.covariance[0]  = cov_ok ? covM(0,0) : 1.0;
            local_dop_meas.twist.covariance[1]  = cov_ok ? covM(1,1) : 1.0;
            local_dop_meas.twist.covariance[2]  = cov_ok ? covM(2,2) : 1.0;
            doppler_updated = true;
        }
        else
        {
            ROS_WARN_THROTTLE(1.0, "Doppler invalid this epoch; keeping last valid velocity.");
        }

        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            if (pos_ok)
            {
                latest_pos_ecef = pos_ecef;
                latest_pos_llh  = pos_llh;
                latest_pos_enu  = pos_enu;
            }
            if (doppler_updated)
            {
                dop_vel_meas = local_dop_meas;
            }
        }

        if (pos_ok)
        {
            sensor_msgs::NavSatFix llh_msg;
            llh_msg.header.frame_id = "map";
            llh_msg.latitude  = pos_llh(0);
            llh_msg.longitude = pos_llh(1);
            llh_msg.altitude  = pos_llh(2);
            pub_psr_llh_latest.publish(llh_msg);

            nav_msgs::Odometry enu_msg;
            enu_msg.header.frame_id = "map";
            enu_msg.pose.pose.position.x = pos_enu(0);
            enu_msg.pose.pose.position.y = pos_enu(1);
            enu_msg.pose.pose.position.z = pos_enu(2);
            pub_psr_enu_latest.publish(enu_msg);

            printf("\033[1;32mCURRENT POS | [%.8f, %.8f, %.3f]| VEL:[%.2f, %.2f, %.2f]\033[0m \n",
                   pos_llh(0), pos_llh(1), pos_llh(2), vel_xyz(0), vel_xyz(1), vel_xyz(2));
        }
        else
        {
            ROS_WARN_THROTTLE(1.0, "SPP invalid this epoch; keeping last valid position.");
        }

        // ── Notify waiting optimization thread ──
        hasNewData.store(true, std::memory_order_release);
        cv_new_data_.notify_one();
        pending_cb_count_.fetch_sub(1, std::memory_order_acq_rel);
    }

    void ephem_cb(const gnss_comm::GnssEphemMsgConstPtr &ephem_msg)
    {
        gnss_comm::EphemPtr ud_ephem(new gnss_comm::Ephem());
        ud_ephem->sat = ephem_msg->sat;
        ud_ephem->ttr = gnss_comm::gpst2time(ephem_msg->ttr.week, ephem_msg->ttr.tow);
        ud_ephem->toe = gnss_comm::gpst2time(ephem_msg->toe.week, ephem_msg->toe.tow);
        ud_ephem->toc = gnss_comm::gpst2time(ephem_msg->toc.week, ephem_msg->toc.tow);
        ud_ephem->toe_tow = ephem_msg->toe_tow;
        ud_ephem->week = ephem_msg->week;
        ud_ephem->iode = ephem_msg->iode;
        ud_ephem->iodc = ephem_msg->iodc;
        ud_ephem->health = ephem_msg->health;
        ud_ephem->code = ephem_msg->code;
        ud_ephem->ura = ephem_msg->ura;
        ud_ephem->A = ephem_msg->A;
        ud_ephem->e = ephem_msg->e;
        ud_ephem->i0 = ephem_msg->i0;
        ud_ephem->omg = ephem_msg->omg;
        ud_ephem->OMG0 = ephem_msg->OMG0;
        ud_ephem->M0 = ephem_msg->M0;
        ud_ephem->delta_n = ephem_msg->delta_n;
        ud_ephem->OMG_dot = ephem_msg->OMG_dot;
        ud_ephem->i_dot = ephem_msg->i_dot;
        ud_ephem->cuc = ephem_msg->cuc;
        ud_ephem->cus = ephem_msg->cus;
        ud_ephem->crc = ephem_msg->crc;
        ud_ephem->crs = ephem_msg->crs;
        ud_ephem->cic = ephem_msg->cic;
        ud_ephem->cis = ephem_msg->cis;
        ud_ephem->af0 = ephem_msg->af0;
        ud_ephem->af1 = ephem_msg->af1;
        ud_ephem->af2 = ephem_msg->af2;
        ud_ephem->tgd[0] = ephem_msg->tgd0;
        ud_ephem->tgd[1] = ephem_msg->tgd1;
        ud_ephem->A_dot = ephem_msg->A_dot;
        ud_ephem->n_dot = ephem_msg->n_dot;

        if(PART_SATSYS)
        {
            int sys = gnss_comm::satsys(ud_ephem->sat, NULL);
            if(sys != SYS_BDS && sys != SYS_GPS)
            {
                return;
            }
        }

        gnss_comm::gtime_t local_time;
        double local_gpst_sec = -1.0;
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            current_ephem_map[ud_ephem->sat] = ud_ephem; //sat已经经过ublox_driver处理，不冲突
            local_time = current_sys_time;
            local_gpst_sec = current_gpst_sec;
        }
        
        // 系统时间未就绪
        if (!gpsTimeValid(local_gpst_sec) && !sysTimeValid(local_time))
        {
            return;
        }

        // std::cout << "ud_ephem->toe: " << ud_ephem->toe.time + ud_ephem->toe.sec << std::endl;
        double tk = time_diff(local_time, ud_ephem->toe);

        if (tk > WEEK_SECONDS/2)  tk -= WEEK_SECONDS;
        else if (tk < -WEEK_SECONDS/2)  tk += WEEK_SECONDS;

        if (std::abs(tk) > EPH_VALID_SECONDS)
        {
            LOG(WARNING) << "Ephemeris is not valid anymore : " << ud_ephem->sat << "->time out:" << std::abs(tk);
            return;
        }

        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            int x;
            for(x = 0; x < (int)ephem_array.size(); x++)
            {
                if(ephem_array[x]->sat == ud_ephem->sat)
                {
                    ephem_array[x] = ud_ephem;
                    break;
                }
            }
            if(x >= (int)ephem_array.size())
            {
                ephem_array.push_back(ud_ephem);
            }
        }
        return;
    }

    void gephem_cb(const gnss_comm::GnssGloEphemMsgConstPtr &gephem_msg)
    {
        if(PART_SATSYS)
        {
            return;
        }
        
        gnss_comm::GloEphemPtr ud_gephem = gnss_comm::msg2glo_ephem(gephem_msg);

        // --- get system time reference ---
        gnss_comm::gtime_t local_time;
        double local_gpst_sec = -1.0;
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            current_glo_ephem_map[ud_gephem->sat] = ud_gephem;
            local_time     = current_sys_time;
            local_gpst_sec = current_gpst_sec;
        }

        // System time not ready yet — defer processing
        if (!gpsTimeValid(local_gpst_sec) && !sysTimeValid(local_time))
        {
            return;
        }

        // --- time validity check ---
        // GLONASS ephemeris valid window: ~30 min (1800s). Use the same 7200s
        // window as GPS/BDS for consistency since toe is already relative.
        double tk = time_diff(local_time, ud_gephem->toe);
        if (tk > WEEK_SECONDS/2)  tk -= WEEK_SECONDS;
        else if (tk < -WEEK_SECONDS/2)  tk += WEEK_SECONDS;

        if (std::abs(tk) > EPH_VALID_SECONDS)
        {
            LOG(WARNING) << "GLONASS ephemeris expired: sat=" << ud_gephem->sat
                         << " tk=" << std::abs(tk) << "s (limit " << EPH_VALID_SECONDS << "s)";
            return;
        }

        // --- update both storage arrays ---
        // 1) gloephem_array: typed storage for direct GloEphem access
        // 2) ephem_array: base-class pointer vector used by psr_pos / dopp_vel_GGLweight
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            int x;

            // Update gloephem_array (typed)
            for (x = 0; x < (int)gloephem_array.size(); x++)
            {
                if (gloephem_array[x]->sat == ud_gephem->sat)
                {
                    gloephem_array[x] = ud_gephem;
                    break;
                }
            }
            if (x >= (int)gloephem_array.size())
            {
                gloephem_array.push_back(ud_gephem);
            }

            // Update ephem_array (base pointer — used by SPP/Doppler solvers)
            for (x = 0; x < (int)ephem_array.size(); x++)
            {
                if (ephem_array[x]->sat == ud_gephem->sat)
                {
                    ephem_array[x] = ud_gephem;
                    break;
                }
            }
            if (x >= (int)ephem_array.size())
            {
                ephem_array.push_back(ud_gephem);
            }
        }

        ROS_DEBUG("GLONASS ephemeris updated: sat=%u, freqo=%d, toe=%.2f",
                  ud_gephem->sat, ud_gephem->freqo,
                  ud_gephem->toe.time + ud_gephem->toe.sec);

        return;
    }

    // ENU original point set
    void origin_node_cb(const sensor_msgs::NavSatFixConstPtr &msg)
    {
        if (!msg)
            return;
        if (!std::isfinite(msg->latitude) || !std::isfinite(msg->longitude) || !std::isfinite(msg->altitude))
            return;
        if (std::abs(msg->latitude) > 90.0 || std::abs(msg->longitude) > 180.0)
            return;
        if (msg->status.status < 0)
            return;

        Eigen::Vector3d llh(msg->latitude, msg->longitude, msg->altitude);
        if (!posValid(llh))
            return;

        std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
        if (enu_ref_set)
            return;

        enu_ref_llh = llh;
        enu_ref_ecef = gnss_comm::geo2ecef(enu_ref_llh);
        latest_pos_llh = llh;
        latest_pos_ecef = enu_ref_ecef;
        latest_pos_enu = Eigen::Vector3d::Zero();
        enu_ref_set = true;
        ROS_INFO("[ProcessingNodeExtra] ENU origin set from /ublox_driver/receiver_lla: (%.9f, %.9f, %.3f)",
                 enu_ref_llh(0), enu_ref_llh(1), enu_ref_llh(2));
    }

    void update_fgo_enu(double e, double n, double h)
    {
        nav_msgs::Odometry fgo_enu_msg;
        fgo_enu_msg.header.frame_id = "map";
        fgo_enu_msg.pose.pose.position.x = e;
        fgo_enu_msg.pose.pose.position.y = n;
        fgo_enu_msg.pose.pose.position.z = h;
        pub_fgo_enu_latest.publish(fgo_enu_msg);
    }

    void update_fgo_llh(double lat, double lon, double alt)
    {
        sensor_msgs::NavSatFix fgo_llh_msg;
        fgo_llh_msg.header.frame_id = "map";
        fgo_llh_msg.latitude = lat;
        fgo_llh_msg.longitude = lon;
        fgo_llh_msg.altitude = alt;
        printf("\033[1;32m   FGO LLH  | [%.8f, %.8f, %.3f]\033[0m | ", lat, lon, alt);
        pub_fgo_llh_latest.publish(fgo_llh_msg);
    }

    bool gpsTimeValid(double time)
    {
        return std::isfinite(time) && time > 0.0;
    }

    bool sysTimeValid(gnss_comm::gtime_t time)
    {
        return (time.time > 0) || (time.time == 0 && time.sec >= 0.0);
    }

    bool can_solve()
    {
        Eigen::Vector3d pos_llh;
        Eigen::Vector3d vel_xyz;
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            pos_llh = latest_pos_llh;
            vel_xyz = {dop_vel_meas.twist.twist.linear.x, dop_vel_meas.twist.twist.linear.y, dop_vel_meas.twist.twist.linear.z};
        }

        return posValid(pos_llh) && velValid(vel_xyz);
    }

    bool posValid(Eigen::Vector3d pos)
    {
        if (!std::isfinite(pos(0)) || !std::isfinite(pos(1)) || !std::isfinite(pos(2)))
        {
            return false;
        }
        if( pos(0) < -90.0 || pos(0) > 90.0 ||
            pos(1) < -180.0 || pos(1) > 180.0 ||
            pos(2) < -450.0 || pos(2) > 100000.0)
        {
            return false;
        }
        else if((pos(0) == 0.0 && pos(1) == 0.0 && pos(2) == 0.0))
        {
            return false;
        }
        else
            return true;
    }
    
    bool velValid(Eigen::Vector3d vel)
    {
        if (!std::isfinite(vel(0)) || !std::isfinite(vel(1)) || !std::isfinite(vel(2)))
        {
            return false;
        }
        if(fabs(vel(0)) > 20.0 || fabs(vel(1)) > 20.0 || fabs(vel(2)) > 20.0)
        {
            return false;
        }
        else if((vel(0) == 0.0 && vel(1) == 0.0 && vel(2) == 0.0))
        {
            return false;
        }
        else
            return true;
    }

    // useless for now, as we directly update sv_info in ephem_cb. 
    // But may be useful in the future if we want to do more refined satellite position/velocity calculation using current receiver position.
    bool update_sat_info(const std::vector<int> &ephem_indices)
    {
        for(int idx : ephem_indices)
        {
            if(idx < 0 || idx >= (int)ephem_array.size())
            {
                continue;
            }

            gnss_comm::EphemPtr ephem = std::dynamic_pointer_cast<gnss_comm::Ephem>(ephem_array[idx]);
            if(!ephem)
            {
                continue;
            }

            const int sat_i = ephem->sat;
            auto &sv = current_sv_info_map[sat_i];
            sv.sat_prn = sat_i;
            sv.sys = gnss_comm::satsys(sat_i, NULL);
            sv.ddt = gnss_comm::eph2svdt(current_sys_time, ephem);
            sv.pos = gnss_comm::eph2pos(current_sys_time, ephem, &sv.dt);
            sv.vel = gnss_comm::eph2vel(current_sys_time, ephem, &sv.ddt);

            if(isnan(sv.pos[0]) || sv.pos[0] == 0)
            {
                sv.avaliable = false;
                continue;
            }
            sv.avaliable = true;

            if(posValid(latest_pos_llh) && enu_ref_set)
            {
                // O(1) refine using current receiver position
                double d_sr = sqrt(pow(latest_pos_ecef(0) - sv.pos[0], 2) +
                                   pow(latest_pos_ecef(1) - sv.pos[1], 2) +
                                   pow(latest_pos_ecef(2) - sv.pos[2], 2));
                double tau = d_sr / LIGHT_SPEED;
                gnss_comm::gtime_t transmit_time = gnss_comm::time_add(current_sys_time, -tau);
                sv.pos = eph2pos(transmit_time, current_ephem_map[sat_i], &sv.dt);

                Eigen::Vector3d sat_pos = {sv.pos[0], sv.pos[1], sv.pos[2]};
                Eigen::Vector3d rev2sat_ecef = (sat_pos - latest_pos_ecef).normalized();
                Eigen::Vector3d rev2sat_enu = gnss_comm::ecef2enu(enu_ref_llh, rev2sat_ecef);
                sv.azimuth = rev2sat_ecef.head<2>().norm() < 1e-12 ? 0.0 : atan2(rev2sat_enu.x(), rev2sat_enu.y());
                sv.azimuth += (sv.azimuth < 0 ? 2*M_PI : 0);
                sv.elevation = asin(rev2sat_enu.z());

                if(sv.elevation > 1.0 || sv.elevation < -1.0)
                {
                    sv.avaliable = false;
                    // std::cout<<" SAT " << sv.sat_prn << "low elevation, not avaliable "<< std::endl;
                    sv.pos = {0,0,0};
                    continue;
                }
            }
        }
        return true;
    }

    void freeSpinner()
    {
        if (range_spinner_) range_spinner_->stop();
        if (ephem_spinner_) ephem_spinner_->stop();
    }

    void saveEphems()
    {
        std::map<int, gnss_comm::EphemPtr> ephem_map_copy;
        {
            std::lock_guard<std::mutex> lk(m_gnss_raw_mux);
            ephem_map_copy = current_ephem_map;
        }
        saveEphemsMap(ephem_map_copy, filename, false);
    }
};
