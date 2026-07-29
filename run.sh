#!/bin/bash

set -euo pipefail

WS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OSQA_INPUT_PATH="${WS_DIR}/osqa_input.jsonl"
OSQA_OUTPUT_PATH="${WS_DIR}/osqa_output.jsonl"
ENABLE_OSQA="true"

# 每次运行前清空交换文件，避免历史数据影响实时处理与性能
: > "${OSQA_INPUT_PATH}"
: > "${OSQA_OUTPUT_PATH}"

if [ "${ENABLE_OSQA}" = "true" ]; then
    # 可选：开启 Transformer/OSQA（传参 1）
    echo "Running GNSS Quality Analyzer with visualization..."
    gnome-terminal -- bash -c "cd \"${WS_DIR}/GNSS-Transformer/gnss_quality_analyzer\" && python3 run_analyzer.py --gnssfgo-input \"${OSQA_INPUT_PATH}\" --gnssfgo-output \"${OSQA_OUTPUT_PATH}\" --urban --vis; exec bash"
    sleep 3s
else
    echo "OSQA/Transformer disabled. Running gnssfgo only."
fi

# 第二个命令在新终端窗口中运行
echo "Opening new terminal for ROS launch..."
gnome-terminal -- bash -c "cd \"${WS_DIR}\" && source devel/setup.bash && roslaunch gnssfgo trbin.launch rviz:=1 play_bag:=1 osqa_enabled:=${ENABLE_OSQA}; exec bash"
