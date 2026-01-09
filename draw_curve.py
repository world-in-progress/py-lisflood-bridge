import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# 配置
# -----------------------------
REAL_DATA_DIR = "real_data"
STATION_NAME = "Tsim Bei Tsui"
ANUGA_CSV = "Tsim_Bei_Tsui_stage.csv"
UVH_FILE = "uvh_ei_36409.txt"

UVH_START_TIME = pd.Timestamp("2023-09-07 17:00:00")
UVH_DT_SECONDS = 300

# -----------------------------
# 读取实测潮位数据
# -----------------------------
all_files = glob.glob(os.path.join(REAL_DATA_DIR, "*.csv"))
obs_list = []

for file in all_files:
    df = pd.read_csv(file)
    df_station = df[df["Tide Station"] == STATION_NAME].copy()

    if not df_station.empty:
        df_station["datetime"] = pd.to_datetime(
            df_station["Date"] + " " + df_station["Time"]
        )
        df_station["Height"] = df_station["Height(m)"].astype(float)
        obs_list.append(df_station[["datetime", "Height"]])

if obs_list:
    obs_df = pd.concat(obs_list, ignore_index=True)
    obs_df = obs_df.sort_values("datetime")
else:
    raise RuntimeError(f"No observation data found for {STATION_NAME}")

# -----------------------------
# 读取 ANUGA 模拟结果
# -----------------------------
sim_df = pd.read_csv(ANUGA_CSV, parse_dates=["datetime"])
sim_df = sim_df[["datetime", "depth_m"]].copy()

# -----------------------------
# 读取 UVH 数据 + 构造时间轴
# -----------------------------
uvh_depth = pd.read_csv(UVH_FILE, header=None, names=["depth_uvh"])

uvh_time = UVH_START_TIME + pd.to_timedelta(
    np.arange(len(uvh_depth)) * UVH_DT_SECONDS, unit="s"
)

uvh_df = pd.DataFrame({
    "datetime": uvh_time,
    "depth_uvh": uvh_depth["depth_uvh"].values
})

# ============================================================
# 第一张大图：三类数据【公共时间段】
# ============================================================
start_time = max(
    obs_df["datetime"].min(),
    sim_df["datetime"].min(),
    uvh_df["datetime"].min()
)
end_time = min(
    obs_df["datetime"].max(),
    sim_df["datetime"].max(),
    uvh_df["datetime"].max()
)

obs_df_overlap = obs_df.query("@start_time <= datetime <= @end_time")
sim_df_overlap = sim_df.query("@start_time <= datetime <= @end_time")
uvh_df_overlap = uvh_df.query("@start_time <= datetime <= @end_time")

# 插值到观测时间
sim_interp = (
    sim_df_overlap
    .set_index("datetime")
    .reindex(obs_df_overlap["datetime"])
    .interpolate(method="time")
)

uvh_interp = (
    uvh_df_overlap
    .set_index("datetime")
    .reindex(obs_df_overlap["datetime"])
    .interpolate(method="time")
)

plt.figure(figsize=(12, 6))
plt.plot(
    obs_df_overlap["datetime"],
    obs_df_overlap["Height"],
    label="Observed Height (m)",
    linewidth=2
)
plt.plot(
    obs_df_overlap["datetime"],
    sim_interp["depth_m"],
    label="ANUGA Depth (m)",
    linestyle="--"
)
plt.plot(
    obs_df_overlap["datetime"],
    uvh_interp["depth_uvh"],
    label="UVH Depth (m)",
    linestyle=":"
)

plt.title(f"{STATION_NAME} – Common Time Range Comparison")
plt.xlabel("Time")
plt.ylabel("Height / Depth (m)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ============================================================
# 第二张大图：三类数据【完整时间序列】
# ============================================================
plt.figure(figsize=(12, 6))

plt.plot(
    obs_df["datetime"],
    obs_df["Height"],
    label="Observed Height (m)",
    linewidth=2
)

plt.plot(
    sim_df["datetime"],
    sim_df["depth_m"],
    label="ANUGA Depth (m)",
    linestyle="--"
)

plt.plot(
    uvh_df["datetime"],
    uvh_df["depth_uvh"],
    label="UVH Depth (m)",
    linestyle=":"
)

plt.title(f"{STATION_NAME} – Full Time Series")
plt.xlabel("Time")
plt.ylabel("Height / Depth (m)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
