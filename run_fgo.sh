#!/bin/bash

set -euo pipefail

WS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OSQA_INPUT_PATH="${WS_DIR}/osqa_input.jsonl"
OSQA_OUTPUT_PATH="${WS_DIR}/osqa_output.jsonl"
ENABLE_OSQA="false"

# 每次运行前清空交换文件，避免历史数据影响实时处理与性能
: > "${OSQA_INPUT_PATH}"
: > "${OSQA_OUTPUT_PATH}"

echo "OSQA/Transformer disabled. Running gnssfgo only."

# 在新终端窗口中启动 gnssfgo
echo "Opening new terminal for ROS launch..."
gnome-terminal -- bash -c "cd \"${WS_DIR}\" && source devel/setup.bash && roslaunch gnssfgo trbin.launch rviz:=1 play_bag:=1 osqa_enabled:=${ENABLE_OSQA}; exec bash"
