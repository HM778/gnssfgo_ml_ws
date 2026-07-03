#!/usr/bin/env python3
import csv
import math
import os
from bisect import bisect_right

import rospy
import yaml

WGS84_A = 6378137.0
WGS84_E2 = 6.69437999014e-3


def lla_to_ecef(lat_deg, lon_deg, alt):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    x = (N + alt) * cos_lat * cos_lon
    y = (N + alt) * cos_lat * sin_lon
    z = (N * (1.0 - WGS84_E2) + alt) * sin_lat
    return (x, y, z)


def ecef_to_enu(x, y, z, ref_lat_deg, ref_lon_deg, ref_alt):
    ref_x, ref_y, ref_z = lla_to_ecef(ref_lat_deg, ref_lon_deg, ref_alt)
    dx = x - ref_x
    dy = y - ref_y
    dz = z - ref_z

    lat = math.radians(ref_lat_deg)
    lon = math.radians(ref_lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    t = -sin_lon * dx + cos_lon * dy
    e = t
    t = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    n = t
    t = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    u = t
    return (e, n, u)


def read_csv(path, expected_fields):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row:
                continue
            if not all(k in row for k in expected_fields):
                continue
            try:
                rows.append({k: float(row[k]) for k in expected_fields})
            except ValueError:
                continue
    return rows


def detect_time_scale(gt_times, sol_times):
    if not gt_times or not sol_times:
        return 1.0
    gt_max = max(gt_times)
    sol_max = max(sol_times)
    if gt_max <= 0 or sol_max <= 0:
        return 1.0
    ratio = sol_max / gt_max
    if ratio > 5.0 and ratio < 15.0:
        return 0.1
    if ratio > 50.0 and ratio < 150.0:
        return 0.01
    if ratio < 0.2 and ratio > 0.06:
        return 10.0
    return 1.0


def interpolate_linear(times, values, t, max_gap):
    idx = bisect_right(times, t)
    if idx == 0 or idx >= len(times):
        return None
    t0 = times[idx - 1]
    t1 = times[idx]
    if t1 - t0 <= 0:
        return None
    if (t1 - t0) > max_gap:
        return None
    v0 = values[idx - 1]
    v1 = values[idx]
    ratio = (t - t0) / (t1 - t0)
    return (
        v0[0] + ratio * (v1[0] - v0[0]),
        v0[1] + ratio * (v1[1] - v0[1]),
        v0[2] + ratio * (v1[2] - v0[2]),
    )


def compute_stats(errors):
    if not errors:
        return (math.nan, math.nan, math.nan)
    abs_vals = [abs(e) for e in errors]
    mean_abs = sum(abs_vals) / len(abs_vals)
    max_abs = max(abs_vals)
    mean = sum(errors) / len(errors)
    var = sum((e - mean) ** 2 for e in errors) / len(errors)
    std = math.sqrt(var)
    return (max_abs, mean_abs, std)


def compute_stats_mag(magnitudes):
    if not magnitudes:
        return (math.nan, math.nan, math.nan)
    mean = sum(magnitudes) / len(magnitudes)
    var = sum((m - mean) ** 2 for m in magnitudes) / len(magnitudes)
    return (max(magnitudes), mean, math.sqrt(var))


def percentile(values, q):
    if not values:
        return math.nan
    q = min(max(q, 0.0), 1.0)
    s = sorted(values)
    if len(s) == 1:
        return s[0]
    pos = (len(s) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    weight = pos - lo
    return s[lo] * (1.0 - weight) + s[hi] * weight


def main():
    rospy.init_node("compare_errors", anonymous=True)

    gt_path = rospy.get_param("~gt_path", "/home/nvidia/Work/gnssfgo_ws/src/gt_lla.csv")
    output_path = rospy.get_param("~output_path", "/home/nvidia/Work/gnssfgo_ws/src/result.csv")
    max_gap = float(rospy.get_param("~max_gap", 2.0))
    std_poor_quantile = float(rospy.get_param("~std_poor_quantile", 0.8))

    scheme_paths = rospy.get_param("~scheme_paths", {
        "std": "/home/nvidia/Work/gnssfgo_ws/src/standard_output.csv",
        "tdcp": "/home/nvidia/Work/gnssfgo_ws/src/tdcp_output.csv",
        "tddcp": "/home/nvidia/Work/gnssfgo_ws/src/tddcp_output.csv",
        "trrtk": "/home/nvidia/Work/gnssfgo_ws/src/trrtk_output.csv",
        "wcp": "/home/nvidia/Work/gnssfgo_ws/src/wcp_output.csv",
        "bin": "/home/nvidia/Work/gnssfgo_ws/src/trbin_output.csv",
    })

    if isinstance(scheme_paths, str):
        try:
            scheme_paths = yaml.safe_load(scheme_paths)
        except yaml.YAMLError:
            rospy.logerr("Invalid scheme_paths format; expected YAML map string.")
            return

    if not isinstance(scheme_paths, dict):
        rospy.logerr("scheme_paths must be a dict, got: %s", type(scheme_paths))
        return

    if not os.path.exists(gt_path):
        rospy.logerr("GT file not found: %s", gt_path)
        return

    gt_rows = read_csv(gt_path, ["timestamp", "latitude", "longitude", "altitude"])
    if len(gt_rows) <= 1:
        rospy.logerr("GT file too short: %s", gt_path)
        return

    gt_rows = gt_rows[1:]
    gt_rows = gt_rows[:10000]
    gt_rows.sort(key=lambda r: r["timestamp"])
    ref = gt_rows[0]
    ref_lat, ref_lon, ref_alt = ref["latitude"], ref["longitude"], ref["altitude"]

    gt_times = []
    gt_enu = []
    for r in gt_rows:
        x, y, z = lla_to_ecef(r["latitude"], r["longitude"], r["altitude"])
        enu = ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt)
        gt_times.append(r["timestamp"])
        gt_enu.append(enu)

    results = []
    series_by_scheme = {}
    for scheme, path in scheme_paths.items():
        if not os.path.exists(path):
            rospy.logwarn("Missing scheme file: %s", path)
            continue
        rows = read_csv(path, ["timestamp", "latitude", "longitude", "altitude"])
        if not rows:
            rospy.logwarn("Empty scheme file: %s", path)
            continue
        rows.sort(key=lambda r: r["timestamp"])
        sol_times_raw = [r["timestamp"] for r in rows]
        time_scale = detect_time_scale(gt_times, sol_times_raw)

        sol_times = []
        sol_enu = []
        for r in rows:
            t = r["timestamp"] * time_scale
            x, y, z = lla_to_ecef(r["latitude"], r["longitude"], r["altitude"])
            enu = ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt)
            sol_times.append(t)
            sol_enu.append(enu)

        e_err = []
        n_err = []
        u_err = []
        h_err = []
        v_err = []
        count = 0

        e_series = []
        n_series = []
        u_series = []
        h_series = []
        v_series = []
        valid_series = []

        for t, gt in zip(gt_times, gt_enu):
            sol_interp = interpolate_linear(sol_times, sol_enu, t, max_gap)
            if sol_interp is None:
                e_series.append(None)
                n_series.append(None)
                u_series.append(None)
                h_series.append(None)
                v_series.append(None)
                valid_series.append(False)
                continue
            de = sol_interp[0] - gt[0]
            dn = sol_interp[1] - gt[1]
            du = sol_interp[2] - gt[2]
            eh = math.hypot(de, dn)
            ev = abs(du)

            e_err.append(de)
            n_err.append(dn)
            u_err.append(du)
            h_err.append(eh)
            v_err.append(ev)
            count += 1

            e_series.append(de)
            n_series.append(dn)
            u_series.append(du)
            h_series.append(eh)
            v_series.append(ev)
            valid_series.append(True)

        max_e, mean_e, std_e = compute_stats(e_err)
        max_n, mean_n, std_n = compute_stats(n_err)
        max_u, mean_u, std_u = compute_stats(u_err)
        max_h, mean_h, std_h = compute_stats_mag(h_err)
        max_v, mean_v, std_v = compute_stats_mag(v_err)

        series_by_scheme[scheme] = {
            "valid": valid_series,
            "h": h_series,
            "v": v_series,
        }

        results.append({
            "scheme": scheme,
            "count": count,
            "availability": count / len(gt_times) if gt_times else math.nan,
            "time_scale": time_scale,
            "max_abs_e": max_e,
            "mean_abs_e": mean_e,
            "std_e": std_e,
            "max_abs_n": max_n,
            "mean_abs_n": mean_n,
            "std_n": std_n,
            "max_abs_u": max_u,
            "mean_abs_u": mean_u,
            "std_u": std_u,
            "max_h": max_h,
            "mean_h": mean_h,
            "std_h": std_h,
            "max_v": max_v,
            "mean_v": mean_v,
            "std_v": std_v,
        })

    std_series = series_by_scheme.get("std")
    std_h_poor_threshold = math.nan
    if std_series:
        std_h_vals = [h for h, v in zip(std_series["h"], std_series["valid"]) if v and h is not None]
        std_h_poor_threshold = percentile(std_h_vals, std_poor_quantile)

        for r in results:
            scheme_series = series_by_scheme.get(r["scheme"])
            if not scheme_series or not std_series or not std_h_vals:
                r.update({
                    "std_h_poor_threshold": std_h_poor_threshold,
                    "count_poor_std": 0,
                    "availability_poor_std": math.nan,
                    "mean_h_poor_std": math.nan,
                    "mean_v_poor_std": math.nan,
                    "p95_h_poor_std": math.nan,
                    "p95_v_poor_std": math.nan,
                })
                continue

            poor_h = []
            poor_v = []
            poor_total = 0
            poor_with_scheme = 0
            for std_valid, std_h, sch_valid, sch_h, sch_v in zip(
                std_series["valid"], std_series["h"],
                scheme_series["valid"], scheme_series["h"], scheme_series["v"]):
                if not std_valid or std_h is None:
                    continue
                if std_h <= std_h_poor_threshold:
                    continue
                poor_total += 1
                if sch_valid and sch_h is not None and sch_v is not None:
                    poor_with_scheme += 1
                    poor_h.append(sch_h)
                    poor_v.append(sch_v)

            r.update({
                "std_h_poor_threshold": std_h_poor_threshold,
                "count_poor_std": poor_with_scheme,
                "availability_poor_std": (poor_with_scheme / poor_total) if poor_total > 0 else math.nan,
                "mean_h_poor_std": (sum(poor_h) / len(poor_h)) if poor_h else math.nan,
                "mean_v_poor_std": (sum(poor_v) / len(poor_v)) if poor_v else math.nan,
                "p95_h_poor_std": percentile(poor_h, 0.95),
                "p95_v_poor_std": percentile(poor_v, 0.95),
            })

    def fmt(value):
        if isinstance(value, float):
            if math.isnan(value) or math.isinf(value):
                return "nan"
            return f"{value:.5f}"
        return value

    with open(output_path, "w", newline="") as f:
        fieldnames = [
            "scheme",
            "count",
            "availability",
            "time_scale",
            "max_abs_e",
            "mean_abs_e",
            "std_e",
            "max_abs_n",
            "mean_abs_n",
            "std_n",
            "max_abs_u",
            "mean_abs_u",
            "std_u",
            "max_h",
            "mean_h",
            "std_h",
            "max_v",
            "mean_v",
            "std_v",
            "std_h_poor_threshold",
            "count_poor_std",
            "availability_poor_std",
            "mean_h_poor_std",
            "mean_v_poor_std",
            "p95_h_poor_std",
            "p95_v_poor_std",
        ]
        f.write("\t".join(fieldnames) + "\n")
        for r in results:
            row = [fmt(r.get(k, "")) for k in fieldnames]
            f.write("\t".join(map(str, row)) + "\n")

    rospy.loginfo("Saved results to %s", output_path)


if __name__ == "__main__":
    main()
