#ifndef DIFFERENTIAL_FUNCTIONS_HPP
#define DIFFERENTIAL_FUNCTIONS_HPP

#include <gnss_comm/gnss_spp.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_constant.hpp>
#include "../datatype.h"
#include "checking.hpp"

/**
 * @brief find master satellites from user end gnss data,
 * @param satellite id, gnss data array (one epoch), gnss data (one satellite)
 * @return none
 */
bool findMasterSatellite(int sat_id, std::vector<gnss_comm::ObsPtr> user_gnss_data, gnss_comm::ObsPtr& master_sv, gnss_comm::ObsPtr& i_sv,std::map<int,sv_info> sv_info_map)
{   
    //确定sat_id所属的卫星系统(GPS/BDS/GAL/GLO)
    /* get navigation satellite system */
    int sys=gnss_comm::satsys(sat_id,NULL);
    if(sys==SYS_GPS)
    {
        // LOG(INFO) << "GPS Satellite  ";
    }
    else if(sys==SYS_GLO)
    {
        // LOG(INFO) << "GLONASS Satellite  ";
    }
    else if(sys==SYS_GAL)
    {
        // LOG(INFO) << "Galileo Satellite  ";
    }
    else if(sys==SYS_BDS)
    {
        // LOG(INFO) << "BeiDou Satellite   ";
    }
    else
    {
        // LOG(INFO) << "Unknown Satellite   ";
    }
    
    int sv_cnt = user_gnss_data.size(); 
    int same_satsystem_cnt = 0;
    double master_sv_id = -1;
    double max_elevation = -1;
    for(int i = 0; i < sv_cnt; i++)
    {
        if(sys == (gnss_comm::satsys(user_gnss_data[i]->sat,NULL))) // same satellite system
        {
            same_satsystem_cnt++;
            //TODO ： 有遮挡情况下这里另当别议
            //确定同一系统卫星的最大高度角--具有最大高度角的为主卫星（信号最良好）
            int sat = user_gnss_data[i]->sat;
            if(sv_info_map[sat].elevation >= max_elevation)
            {
                max_elevation = sv_info_map[sat].elevation;
                master_sv_id = i;
            }

            // same satellite system and same satellite id
            // 具有相同卫星系统和相同卫星id的卫星为i卫星（双差的另一个卫星）
            if(sat_id ==int(user_gnss_data[i]->sat))
            {
                //更新信息
                i_sv = user_gnss_data[i];
                // i_sv->elevation = sv_info_map[sat].elevation;
                // if(i_sv.elevation==0)
                // {
                //     LOG(INFO) << "satellite with zero elevation angle---";
                // }
            }
        }
    }

    if(same_satsystem_cnt>=2) // 同一系统下2颗以上卫星 主卫星才有意义 at least 2 satellites with same sat system
    {
        master_sv = user_gnss_data[master_sv_id];
        // LOG(INFO) << "elevation of master satellite " << max_elevation;
        /*check if you have find the satellite with same ID in user end*/
        if(sat_id==master_sv->sat || (i_sv->psr[0]<10)) // the sat_id is same as master satellite (double-difference should not be done between master and master)
        {
            // LOG(INFO) << "Warning!!! master satellite is itself!!!!!";
            return false;
        }
        return true;
    }
    else 
        return false;

}

/**
 * @brief find closest gnss measurement from user end gnss data
 * @param time, gnss map from user end, gnss data (one epoch)
 * @return none
*/
void findClosestEpoch(double t, std::map<double, std::vector<gnss_comm::ObsPtr>> gnss_raw_map, std::vector<gnss_comm::ObsPtr>& cloest_epoch_gnss)
{
    std::map<double, std::vector<gnss_comm::ObsPtr>>::iterator gnss_iter;
    gnss_iter = gnss_raw_map.begin();
    int length = gnss_raw_map.size();
    double time_diff = 100000;
    for(int i = 0;  i < length; i++,gnss_iter++) // initialize
    {
        if((fabs(t - gnss_iter->first))<time_diff)
        {
            // with smallest time difference
            time_diff = fabs(t - gnss_iter->first);
            cloest_epoch_gnss = gnss_iter->second;
        }
        

    }
    // LOG(INFO) << "time_diff " << time_diff;
}

/**
 * @brief find satellite with same id
 * @param id, gnss data in one epoch, one satellite
 * @return none
*/
void findSatellitewithSameId(double id, std::vector<gnss_comm::ObsPtr> gnss_data, gnss_comm::ObsPtr& same_id_sv, int freq_idx)
{
    int length = gnss_data.size();
    for(int i = 0; i < length; i++)
    {
        if(gnss_data[i]->sat == id )
        {
            gnss_comm::L1_freq(gnss_data[i],&search_l1);
            gnss_comm::L2_freq(gnss_data[i],&search_l2);
            if((freq_idx == 1 && search_l1 >= 0) || (freq_idx == 2 && search_l2 >= 0))
            {
                same_id_sv = gnss_data[i];
                break;
            }
        }
    }
}



#endif // DIFFERENTIAL_FUNCTIONS_HPP
