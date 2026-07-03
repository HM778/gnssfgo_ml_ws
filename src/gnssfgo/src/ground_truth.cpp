#include <cmath>
#include <string>
#include <iostream>
#include <Eigen/Eigen>
#include <fstream>
#include <iomanip>
#include <mutex>

#include <ros/ros.h>
#include <ros/package.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/NavSatFix.h>

#include <gnss_comm/GnssMeasMsg.h>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/GnssPVTSolnMsg.h>

#include <GeographicLib/LocalCartesian.hpp>

class GroundTruthEnuNode
{
public:
	GroundTruthEnuNode(ros::NodeHandle nh, ros::NodeHandle pnh)
	: nh_(std::move(nh)),
		pnh_(std::move(pnh)),
		origin_set_(false),
		local_cart_()
	{
		pnh_.param<std::string>("input_topic", input_topic_, std::string("/ublox_driver/receiver_pvt"));
		pnh_.param<std::string>("output_topic", output_topic_, std::string("/gnssfgo/gt_enu"));
		pnh_.param<std::string>("frame_id", frame_id_, std::string("map"));
		pnh_.param<std::string>("child_frame_id", child_frame_id_, std::string("gt"));
		pnh_.param<std::string>("output_file", output_file_, std::string("gt_lla.csv"));
		pnh_.param<bool>("recording",recording,true);
		pnh_.param<bool>("use_gps_time", use_gps_time_, true);
		pnh_.param<std::string>("gps_time_topic", gps_time_topic_, std::string("/ublox_driver/range_meas"));
		if (!output_file_.empty() && output_file_.front() != '/')
		{
			output_file_ = ros::package::getPath("gnssfgo") + "/../" + output_file_;
		}

		pub_ = nh_.advertise<nav_msgs::Odometry>(output_topic_, 10);
		sub_ = nh_.subscribe(input_topic_, 50, &GroundTruthEnuNode::gnssPVTSolnCb, this);
		if (use_gps_time_)
		{
			gps_time_sub_ = nh_.subscribe(gps_time_topic_, 50, &GroundTruthEnuNode::gpsTimeCb, this);
			ROS_INFO("[ground_truth_node] GPS time source: %s", gps_time_topic_.c_str());
		}

		ROS_INFO("[ground_truth_node] Subscribing: %s", input_topic_.c_str());
		ROS_INFO("[ground_truth_node] Publishing:  %s", output_topic_.c_str());
		ROS_INFO("[ground_truth_node] Output file: %s", output_file_.c_str());
	}

private:
	static bool isValidFix(const sensor_msgs::NavSatFix &msg)
	{
		if (!std::isfinite(msg.latitude) || !std::isfinite(msg.longitude) || !std::isfinite(msg.altitude))
			return false;
		if (std::abs(msg.latitude) > 90.0 || std::abs(msg.longitude) > 180.0)
			return false;
		if (msg.status.status < 0)
			return false;
		return true;
	}

	void setOrigin(const sensor_msgs::NavSatFix &msg)
	{
		origin_lla_lat_ = msg.latitude;
		origin_lla_lon_ = msg.longitude;
		origin_lla_alt_ = msg.altitude;
		local_cart_.Reset(origin_lla_lat_, origin_lla_lon_, origin_lla_alt_);
		origin_set_ = true;
		ROS_INFO("[ground_truth_node] Origin set LLA=(%.9f, %.9f, %.3f)", origin_lla_lat_, origin_lla_lon_, origin_lla_alt_);
	}

	bool savetofile(Eigen::Vector3d lla, double gpst_sec, int fixtype)
    {
        std::ofstream outfile(output_file_, std::ios::out | std::ios::app);
        if (!outfile.is_open())
        {
            ROS_ERROR("Failed to open output file: %s", output_file_.c_str());
            return false;
        }

        // Write header if file is new
        if (outfile.tellp() == 0)
        {
			outfile << "timestamp,latitude,longitude,altitude,fixtype\n";
        }

        // Write data
        outfile << std::fixed << std::setprecision(6) << gpst_sec << ","
                << lla[0] << ","
                << lla[1] << ","
				<< lla[2] << ","
				<< fixtype << "\n";

        outfile.close();
        return true;
    }

	void gpsTimeCb(const gnss_comm::GnssMeasMsgConstPtr &meas_msg)
	{
		if (!meas_msg || meas_msg->meas.empty())
		{
			return;
		}
		const auto &time_msg = meas_msg->meas.front().time;
		double gpst_sec = static_cast<double>(time_msg.week) * WEEK_SECONDS + time_msg.tow;
		if (!std::isfinite(gpst_sec) || gpst_sec <= 0.0)
		{
			return;
		}
		std::lock_guard<std::mutex> lk(gpst_mutex_);
		latest_gpst_sec_ = gpst_sec;
		latest_gpst_valid_ = true;
	}

	bool getLatestGpsSec(double &gpst_sec)
	{
		std::lock_guard<std::mutex> lk(gpst_mutex_);
		if (!latest_gpst_valid_)
		{
			return false;
		}
		gpst_sec = latest_gpst_sec_;
		return true;
	}

	void gnssPVTSolnCb(const gnss_comm::GnssPVTSolnMsgConstPtr &msg)
	{
		if (!msg)
		{
			return;
		}
		if (!std::isfinite(msg->latitude) || !std::isfinite(msg->longitude) || !std::isfinite(msg->altitude))
		{
			ROS_WARN_THROTTLE(1.0, "[ground_truth_node] Invalid PVT (lat/lon/alt)");
			return;
		}
		if (std::abs(msg->latitude) > 90.0 || std::abs(msg->longitude) > 180.0)
		{
			ROS_WARN_THROTTLE(1.0, "[ground_truth_node] Invalid PVT (range)");
			return;
		}

		sensor_msgs::NavSatFix fix_msg;
		fix_msg.latitude = msg->latitude;
		fix_msg.longitude = msg->longitude;
		fix_msg.altitude = msg->altitude;
		fix_msg.status.status = msg->valid_fix ? 0 : -1;

		if (!origin_set_)
		{
			setOrigin(fix_msg);
		}

		double east = 0.0, north = 0.0, up = 0.0;
		local_cart_.Forward(msg->latitude, msg->longitude, msg->altitude, east, north, up);

		double time_sec = msg->time.week * WEEK_SECONDS + msg->time.tow;
		if (!std::isfinite(time_sec) || time_sec <= 0.0)
		{
			time_sec = ros::Time::now().toSec();
		}

		const int fixtype = (msg->valid_fix ? 3 : 0) + (msg->diff_soln ? 1 : 0) + msg->carr_soln;

		if (recording)
		{
			ROS_INFO("recording ground truth");
			savetofile(Eigen::Vector3d(msg->latitude, msg->longitude, msg->altitude), time_sec, fixtype);
		}

		nav_msgs::Odometry odom;
		odom.header.stamp = ros::Time().fromSec(time_sec);
		odom.header.frame_id = frame_id_;
		odom.child_frame_id = child_frame_id_;

		odom.pose.pose.position.x = east;
		odom.pose.pose.position.y = north;
		odom.pose.pose.position.z = up;
		odom.pose.pose.orientation.w = 1.0;
		odom.pose.pose.orientation.x = 0.0;
		odom.pose.pose.orientation.y = 0.0;
		odom.pose.pose.orientation.z = 0.0;

		pub_.publish(odom);
	}

	void navSatFixCb(const sensor_msgs::NavSatFixConstPtr &msg)
	{
		if (!msg)
			return;
		if (!isValidFix(*msg))
		{
			ROS_WARN_THROTTLE(1.0, "[ground_truth_node] Invalid NavSatFix (lat/lon/alt/status)");
			return;
		}

		if (!origin_set_)
			setOrigin(*msg);

		double east = 0.0, north = 0.0, up = 0.0;
		local_cart_.Forward(msg->latitude, msg->longitude, msg->altitude, east, north, up);
		
		double time_sec = msg->header.stamp.toSec();
		if (use_gps_time_)
		{
			double gpst_sec = -1.0;
			if (getLatestGpsSec(gpst_sec))
			{
				time_sec = gpst_sec;
			}
		}

		if(recording)
		{
			ROS_INFO("recording ground truth");
			savetofile(Eigen::Vector3d(msg->latitude, msg->longitude, msg->altitude), time_sec, 0);
		}

		nav_msgs::Odometry odom;
		odom.header.stamp = ros::Time().fromSec(time_sec);
		odom.header.frame_id = frame_id_;
		odom.child_frame_id = child_frame_id_;

		odom.pose.pose.position.x = east;
		odom.pose.pose.position.y = north;
		odom.pose.pose.position.z = up;
		odom.pose.pose.orientation.w = 1.0;
		odom.pose.pose.orientation.x = 0.0;
		odom.pose.pose.orientation.y = 0.0;
		odom.pose.pose.orientation.z = 0.0;

		// If the NavSatFix covariance is provided, map ENU (E,N,U) into pose covariance.
		// NavSatFix position_covariance is documented as ENU in meters^2 when known.
		if (msg->position_covariance_type != sensor_msgs::NavSatFix::COVARIANCE_TYPE_UNKNOWN)
		{
			for (double &c : odom.pose.covariance)
				c = 0.0;
			odom.pose.covariance[0] = msg->position_covariance[0];  // xx (E)
			odom.pose.covariance[1] = msg->position_covariance[1];
			odom.pose.covariance[2] = msg->position_covariance[2];
			odom.pose.covariance[6] = msg->position_covariance[3];
			odom.pose.covariance[7] = msg->position_covariance[4];  // yy (N)
			odom.pose.covariance[8] = msg->position_covariance[5];
			odom.pose.covariance[12] = msg->position_covariance[6];
			odom.pose.covariance[13] = msg->position_covariance[7];
			odom.pose.covariance[14] = msg->position_covariance[8]; // zz (U)
		}

		pub_.publish(odom);
	}

private:
	ros::NodeHandle nh_;
	ros::NodeHandle pnh_;
	ros::Subscriber sub_;
	ros::Subscriber gps_time_sub_;
	ros::Publisher pub_;

	std::string input_topic_;
	std::string output_topic_;
	std::string frame_id_;
	std::string child_frame_id_;
	std::string output_file_;
	bool recording;
	bool use_gps_time_{true};
	std::string gps_time_topic_;
	std::mutex gpst_mutex_;
	double latest_gpst_sec_{-1.0};
	bool latest_gpst_valid_{false};

	bool origin_set_;
	double origin_lla_lat_{};
	double origin_lla_lon_{};
	double origin_lla_alt_{};
	GeographicLib::LocalCartesian local_cart_;
};

int main(int argc, char **argv)
{
	ros::init(argc, argv, "ground_truth_node");
	ros::NodeHandle nh;
	ros::NodeHandle pnh("~");
	GroundTruthEnuNode node(nh, pnh);
	ros::spin();
	return 0;
}

