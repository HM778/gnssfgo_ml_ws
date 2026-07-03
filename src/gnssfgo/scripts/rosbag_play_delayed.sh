#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <delay_sec> <rate> <start_time> <use_clock:true|false> <bag_path> [duration_sec] [extra rosbag args...]" >&2
  exit 2
fi

delay_sec="$1"
rate="$2"
start_time="$3"
use_clock="$4"
bag_path="$5"
duration_sec="${6:-}"
shift 5

if [[ -n "${duration_sec}" && "${duration_sec}" != "0" && "${duration_sec}" != "0.0" ]]; then
  rosbag_args+=(--duration "${duration_sec}")
  shift 1
fi

sleep "$delay_sec"

rosbag_args=(play -r "$rate" -s "$start_time")

if [[ "$use_clock" == "true" || "$use_clock" == "1" ]]; then
  rosbag_args+=(--clock)
fi

rosbag_args+=("$bag_path")

exec rosbag "${rosbag_args[@]}" "$@"
