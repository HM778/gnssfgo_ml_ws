cmake_minimum_required(VERSION 2.8.3)
project(gnssfgo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Debug)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

find_package(catkin REQUIRED COMPONENTS
  roslib
  roscpp
  rospy
  std_msgs
  nav_msgs
  sensor_msgs
  gnss_comm
)


find_package(Eigen3)
find_package(Ceres REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GEOGRAPHICLIB REQUIRED geographiclib)

catkin_package(
  CATKIN_DEPENDS 
    roslib
    roscpp 
    std_msgs 
    nav_msgs
    sensor_msgs
    gnss_comm
  DEPENDS 
    Eigen3 
    Ceres
)

include_directories(include)
# add_subdirectory(./ThirdParty/GeographicLib/)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
#   ./ThirdParty/GeographicLib/include/
)

add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(ground_truth_node src/ground_truth.cpp)

# trbin trrtk
add_executable(trbin_node src/trbin_main.cpp)
target_compile_definitions(trbin_node PRIVATE
  EIGEN_MALLOC_ALREADY_ALIGNED=0
  EIGEN_MAX_ALIGN_BYTES=16
)
if(ENABLE_ASAN)
  target_compile_options(trbin_node PRIVATE -fsanitize=address -fno-omit-frame-pointer -O1 -g)
  target_link_options(trbin_node PRIVATE -fsanitize=address)
endif()


# libGeographiccc
target_link_libraries(trbin_node ${catkin_LIBRARIES} ${CERES_LIBRARIES})

target_link_libraries(ground_truth_node ${catkin_LIBRARIES} ${GEOGRAPHICLIB_LIBRARIES})

install(PROGRAMS
  scripts/rosbag_play_delayed.sh
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
