import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.spatial import cKDTree
import matplotlib.tri as mtri
from datetime import datetime, timedelta

# ===============================
# 文件路径
# ===============================
SWW_FILE = "domain.sww"
STATION_CSV = "coord_2326.csv"

# 最大允许匹配距离 (m)，可根据 DEM 网格间距设置
MAX_DIST = 50  # 超过 50 m 就警告

# 初始时间
START_TIME = datetime(2023, 9, 7, 17, 0, 0)

# ===============================
# 读取潮位站坐标（2326）
# ===============================
stations = pd.read_csv(STATION_CSV)
# CSV 至少包含：Station, Easting_2326_m, Northing_2326_m

# ===============================
# 打开 SWW 文件
# ===============================
ds = Dataset(SWW_FILE, "r")

x = ds.variables["x"][:]  # 节点 X
y = ds.variables["y"][:]  # 节点 Y
time = ds.variables["time"][:]  # 时间序列 (秒)
stage = ds.variables["stage"][:]  # (nt, nn)
elevation = ds.variables["elevation"][:]  # (nn,)
tris = ds.variables["triangles"][:] if "triangles" in ds.variables else None  # 可选

print(f"X range: {x.min()} ~ {x.max()}")
print(f"Y range: {y.min()} ~ {y.max()}")
print(f"✔ SWW loaded: nodes={len(x)}, timesteps={len(time)}")

# ===============================
# KDTree + Triangulation
# ===============================
pts = np.column_stack((x, y))
kdtree = cKDTree(pts)

tri_obj = None
trifinder = None
if tris is not None:
    tris = np.array(tris, dtype=int)
    tri_obj = mtri.Triangulation(x, y, triangles=tris)
    trifinder = tri_obj.get_trifinder()
    print(f"✔ Triangulation loaded: triangles={len(tris)}")
else:
    print("⚠ SWW 文件没有 triangles 信息，无法找所在三角形。")

# ===============================
# 工具：找最近节点
# ===============================
def find_nearest_node(px, py, kdtree_obj):
    """
    返回最近节点索引和距离 (m)
    """
    dist, idx = kdtree_obj.query([px, py], k=1)
    return idx, dist

# ===============================
# 对每个潮位站提取 stage + 三角形信息
# ===============================
for _, row in stations.iterrows():
    name = row["Station"].replace(" ", "_")
    px = row["Easting_2326_m"]
    py = row["Northing_2326_m"]

    # 最近节点
    node_id, dist = find_nearest_node(px, py, kdtree)
    if dist > MAX_DIST:
        print(f"⚠ WARNING: {name} matched node {node_id} is far away! dist={dist:.2f} m")
    dx = x[node_id] - px
    dy = y[node_id] - py
    print(f"▶ {name}: nearest node {node_id}, dx={dx:.2f} m, dy={dy:.2f} m, dist={dist:.2f} m")

    # 所在三角形
    if trifinder is not None:
        tri_id = trifinder(px, py)
        if tri_id != -1:
            tri_verts_idx = tris[tri_id]
            tri_coords = [(x[i], y[i]) for i in tri_verts_idx]
            print(f"  所在三角形索引: {tri_id}")
            for i, (vx, vy) in enumerate(tri_coords):
                print(f"    顶点{i}: ({vx:.3f}, {vy:.3f})")
        else:
            print(f"  ⚠ {name} 坐标不在任何三角形内")

    # 提取时间序列 stage
    station_stage = stage[:, node_id]
    elevation_stage = elevation[node_id]

    # 计算水深，保证不为负
    depth = np.maximum(station_stage - elevation_stage, 0)

    # 组成 DataFrame：时间加上初始时间
    df_out = pd.DataFrame({
        "datetime": [START_TIME + timedelta(seconds=t) for t in time],
        "time_seconds": time,
        "stage_m": station_stage,
        "depth_m": depth
    })

    out_csv = f"{name}_stage.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"  ✔ written: {out_csv}")

ds.close()
print("✅ 所有潮位站提取完成")
