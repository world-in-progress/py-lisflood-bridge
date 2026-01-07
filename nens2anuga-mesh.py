#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
功能：
- 读取 NS/NE 文件
- 生成带悬挂点处理的三角网格（你已有的算法）
- 导出 PLY（保留）
- 将网格转为 ANUGA Domain 并进行简单演化（5s）
"""

import math
from typing import List, Dict, Tuple
import numpy as np
import anuga
import csv
from datetime import datetime, timedelta

from infiltration_operator import TestOperator

# ---------------------------
# （保留你原来的类/解析函数）
# ---------------------------
class HydroSide:
    def __init__(self, index:int, direction:int, bounds:Tuple[float,float,float,float],
                 grid_index_a:int, grid_index_b:int, altitude:float=0.0, type_:int=0):
        self.index = index
        self.direction = direction
        self.bounds = bounds
        self.grid_index_a = grid_index_a
        self.grid_index_b = grid_index_b
        self.altitude = altitude
        self.type = type_

def parse_ns_line(row: list) -> HydroSide:
    index = int(row[0])
    direction_ns = int(row[1])
    left = int(row[2])
    right = int(row[3])
    bottom = int(row[4])
    top = int(row[5])
    length = float(row[6])
    cx, cy, cz = float(row[7]), float(row[8]), float(row[9])
    type_ = int(row[10])

    direction = 1 if direction_ns == 1 else 0

    if direction == 1:
        min_x = cx - length / 2.0
        max_x = cx + length / 2.0
        min_y = cy
        max_y = cy
        grid_index_a = top
        grid_index_b = bottom
    else:
        min_x = cx
        max_x = cx
        min_y = cy - length / 2.0
        max_y = cy + length / 2.0
        grid_index_a = left
        grid_index_b = right

    bounds = (min_x, min_y, max_x, max_y)
    return HydroSide(index, direction, bounds, grid_index_a, grid_index_b, cz, type_)

def load_ns_file(path: str) -> Dict[int, HydroSide]:
    side_map = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            row = [float(p) if '.' in p else int(p) for p in parts]
            side = parse_ns_line(row)
            side_map[side.index] = side
    return side_map

class HydroElement:
    def __init__(self, index:int, bounds:Tuple[float,float,float,float],
                 left_edges:List[int], right_edges:List[int],
                 bottom_edges:List[int], top_edges:List[int],
                 altitude:float, type_:int):
        self.index = index
        self.bounds = bounds
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.bottom_edges = bottom_edges
        self.top_edges = top_edges
        self.altitude = altitude
        self.type = type_

    @property
    def center(self):
        min_x, min_y, max_x, max_y = self.bounds
        return ((min_x + max_x)/2.0, (min_y + max_y)/2.0, self.altitude)

def parse_ne_line(row: list, side_map: Dict[int, HydroSide]) -> HydroElement:
    idx = 0
    index = int(row[idx]); idx += 1

    Ln = int(row[idx]); idx += 1
    Rn = int(row[idx]); idx += 1
    Bn = int(row[idx]); idx += 1
    Tn = int(row[idx]); idx += 1

    left_edges = [int(row[idx + i]) for i in range(Ln)]; idx += Ln
    right_edges = [int(row[idx + i]) for i in range(Rn)]; idx += Rn
    bottom_edges = [int(row[idx + i]) for i in range(Bn)]; idx += Bn
    top_edges = [int(row[idx + i]) for i in range(Tn)]; idx += Tn

    cx = float(row[idx]); idx += 1
    cy = float(row[idx]); idx += 1
    cz = float(row[idx]); idx += 1
    type_ = int(row[idx]); idx += 1

    if left_edges:
        min_x = side_map[left_edges[0]].bounds[0]
    else:
        min_x = cx

    if right_edges:
        max_x = side_map[right_edges[0]].bounds[2]
    else:
        max_x = cx

    if bottom_edges:
        min_y = side_map[bottom_edges[0]].bounds[1]
    else:
        min_y = cy

    if top_edges:
        max_y = side_map[top_edges[0]].bounds[3]
    else:
        max_y = cy

    bounds = (min_x, min_y, max_x, max_y)

    return HydroElement(index, bounds, left_edges, right_edges, bottom_edges, top_edges, cz, type_)

def load_ne_file(path: str, side_map: Dict[int, HydroSide]) -> List[HydroElement]:
    elements = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            row = [float(p) if '.' in p else int(p) for p in parts]
            elements.append(parse_ne_line(row, side_map))
    return elements

# -------------------------
# 你的三角化函数（保留）
# -------------------------
def tri_area(a, b, c) -> float:
    (x1,y1) = a; (x2,y2) = b; (x3,y3) = c
    return abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))/2.0

# -------------------------
# 你的三角化函数（修改后）
# -------------------------
def elements_to_triangles(elements: List[HydroElement], side_map: Dict[int, HydroSide]):
    points: List[Tuple[float, float]] = []
    elevations: List[float] = []
    triangles: List[Tuple[int,int,int]] = []
    point_index_map = {}

    def add_point(x, y, z):
        key = (round(float(x), 9), round(float(y), 9))
        if key in point_index_map:
            return point_index_map[key]
        idx = len(points)
        points.append((x, y))
        elevations.append(z)
        point_index_map[key] = idx
        return idx

    def side_endpoints(side: HydroSide):
        min_x, min_y, max_x, max_y = side.bounds
        return [(min_x, min_y), (max_x, max_y)]

    def extract_hanging(edge_ids, is_vertical):
        if len(edge_ids) <= 1:
            return []
        flat = []
        for sid in edge_ids:
            s = side_map[sid]
            a, b = side_endpoints(s)
            flat.append(a); flat.append(b)
        flat = list(set(flat))
        if is_vertical:
            flat = sorted(flat, key=lambda p: p[1])
        else:
            flat = sorted(flat, key=lambda p: p[0])
        return flat[1:-1]

    for elem in elements:
        min_x, min_y, max_x, max_y = elem.bounds
        z = elem.altitude

        left_hang   = extract_hanging(elem.left_edges,   True)
        right_hang  = extract_hanging(elem.right_edges,  True)
        bottom_hang = extract_hanging(elem.bottom_edges, False)
        top_hang    = extract_hanging(elem.top_edges,    False)

        has_hanging = (
            len(left_hang) > 0 or
            len(right_hang) > 0 or
            len(bottom_hang) > 0 or
            len(top_hang) > 0
        )

        if has_hanging:
            cx, cy, cz = elem.center
            center_idx = add_point(cx, cy, cz)

            boundary = []
            boundary.append((min_x, min_y))
            boundary.extend(bottom_hang)
            boundary.append((max_x, min_y))
            boundary.extend(right_hang)
            boundary.append((max_x, max_y))
            boundary.extend(reversed(top_hang))
            boundary.append((min_x, max_y))
            boundary.extend(reversed(left_hang))

            b_idx = [add_point(px, py, z) for (px, py) in boundary]
            for i in range(len(b_idx)):
                a, b, c = center_idx, b_idx[i], b_idx[(i+1) % len(b_idx)]
                # 确保三个顶点不重复
                if a != b and b != c and c != a and tri_area(points[a], points[b], points[c]) > 1e-12:
                    triangles.append((a, b, c))
            continue

        # no hanging: simple split
        p0 = add_point(min_x, min_y, z)
        p1 = add_point(max_x, min_y, z)
        p2 = add_point(max_x, max_y, z)
        p3 = add_point(min_x, max_y, z)
        if p0 != p1 and p1 != p2 and p2 != p0 and tri_area(points[p0], points[p1], points[p2]) > 1e-12:
            triangles.append((p0, p1, p2))
        if p0 != p2 and p2 != p3 and p3 != p0 and tri_area(points[p0], points[p2], points[p3]) > 1e-12:
            triangles.append((p0, p2, p3))

    return points, triangles, elevations


# -------------------------
# 按 (x,y) 去重 + 过滤退化三角形（修改后）
# -------------------------
def remove_duplicate_vertices_by_xy(points, triangles, elevations, round_digits=9):
    pts = np.asarray(points, dtype=float)    # (P,2)
    elevs = np.asarray(elevations, dtype=float)  # (P,)
    tris = np.asarray(triangles, dtype=np.int64) # (T,3)

    keys = np.round(pts, round_digits)
    structured = np.core.records.fromarrays([keys[:,0], keys[:,1]], names='x,y')
    unique_vals, unique_idx, inverse_idx = np.unique(structured, return_index=True, return_inverse=True)

    coords_unique = pts[unique_idx]
    N = coords_unique.shape[0]

    # 合并高度
    elevs_acc = np.zeros(N, dtype=float)
    counts = np.zeros(N, dtype=int)
    for orig_i, map_i in enumerate(inverse_idx):
        elevs_acc[map_i] += elevs[orig_i]
        counts[map_i] += 1
    counts[counts==0] = 1
    elevs_unique = elevs_acc / counts

    # 重建三角形索引
    try:
        tris_mapped = inverse_idx[tris]
    except Exception as e:
        valid_tris = []
        P = pts.shape[0]
        for tri in tris:
            if np.any(tri < 0) or np.any(tri >= P):
                continue
            valid_tris.append(tuple(inverse_idx[list(tri)]))
        tris_mapped = np.array(valid_tris, dtype=np.int64)

    # 删除退化三角形（面积接近 0 或重复顶点）
    valid_tris = []
    degenerate_count = 0
    for tri in tris_mapped:
        i, j, k = tri
        if i == j or j == k or k == i:
            degenerate_count += 1
            continue
        a, b, c = coords_unique[i], coords_unique[j], coords_unique[k]
        area = abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])) / 2.0
        if area < 1e-12:
            degenerate_count += 1
            continue
        valid_tris.append(tri)
    print(f"过滤了 {degenerate_count} 个退化三角形")
    tris_mapped = np.array(valid_tris, dtype=np.int64)

    return coords_unique.astype(np.float64), tris_mapped.astype(np.int64), elevs_unique.astype(np.float64)


# -------------------------
# PLY 导出（保留）
# -------------------------
def export_ply(points, triangles, elevations, path="mesh_output.ply"):
    with open(path, "w") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write(f"element vertex {len(points)}\n")
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write(f"element face {len(triangles)}\n")
        f.write("property list uchar int vertex_indices\n")
        f.write("end_header\n")
        for i, (x, y) in enumerate(points):
            f.write(f"{x} {y} {elevations[i]}\n")
        for (i,j,k) in triangles:
            f.write(f"3 {i} {j} {k}\n")
    print(f"PLY 已导出: {path}")



# --------- 替换：按 (x,y) 去重并重建三角形索引 ----------
def remove_duplicate_vertices_by_xy(points, triangles, elevations, round_digits=9):
    """
    按 (x,y) 去重：
    - points: list of (x,y)
    - triangles: list of (i,j,k)
    - elevations: list of z (one per original point)
    返回: coords_unique (N,2), tris_unique (M,3), elevs_unique (N,)
    高程合并策略：对重复 (x,y) 的点取平均高度
    """
    pts = np.asarray(points, dtype=float)    # (P,2)
    elevs = np.asarray(elevations, dtype=float)  # (P,)
    tris = np.asarray(triangles, dtype=np.int64) # (T,3)

    if pts.ndim != 2 or pts.shape[1] != 2:
        raise RuntimeError("points 应为 (P,2) 形状")

    # 建立键（四舍五入，避免浮点抖动导致不匹配）
    keys = np.round(pts, round_digits)
    # 将 (x,y) 转为可 unique 的 structured array
    structured = np.core.records.fromarrays([keys[:,0], keys[:,1]], names='x,y')
    unique_vals, unique_idx, inverse_idx = np.unique(structured, return_index=True, return_inverse=True)

    coords_unique = pts[unique_idx]  # (N,2)

    # 合并高度（平均）
    N = coords_unique.shape[0]
    elevs_acc = np.zeros(N, dtype=float)
    counts = np.zeros(N, dtype=int)
    for orig_i, map_i in enumerate(inverse_idx):
        elevs_acc[map_i] += elevs[orig_i]
        counts[map_i] += 1
    # 防止除零(理论上 counts 都>0)
    counts[counts==0] = 1
    elevs_unique = elevs_acc / counts

    # 重建三角形索引
    # inverse_idx: 原顶点 i -> 新顶点 index
    try:
        tris_mapped = inverse_idx[tris]  # (T,3)  通过 numpy 花式索引
    except Exception as e:
        # 如果 tris 中有越界或负值，会报错；在此做保护性处理
        print("⚠ tris 中可能包含非法索引，尝试逐项映射并跳过非法面...")
        valid_tris = []
        P = pts.shape[0]
        for tri in tris:
            if np.any(tri < 0) or np.any(tri >= P):
                continue
            valid_tris.append(tuple(inverse_idx[list(tri)]))
        tris_mapped = np.array(valid_tris, dtype=np.int64)

    return coords_unique.astype(np.float64), tris_mapped.astype(np.int64), elevs_unique.astype(np.float64)



def read_rain_csv(path):
    """
    返回时间序列降雨量列表 [(datetime, rain_m_per_s), ...]
    假设 CSV 单位是 mm/step，转换为 m/s
    """
    rain_data = []
    with open(path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            dt = datetime.strptime(row['DateTime'], "%Y/%m/%d %H:%M")
            rain_mm = float(row['rainfall'])
            rain_m = rain_mm / 1000.0   # mm -> m
            rain_data.append((dt, rain_m))
    return rain_data

def read_tide_csv(path):
    """
    返回 [(datetime, stage_m), ...]
    自动处理 BOM
    """
    data = []
    with open(path, newline='', encoding='utf-8-sig') as f:  # 注意 'utf-8-sig'
        reader = csv.DictReader(f)
        print("CSV 表头:", reader.fieldnames)
        for row in reader:
            dt = datetime.strptime(
                row['date'].strip() + " " + row['time'].strip(),
                "%m/%d/%Y %H:%M:%S"
            )
            stage = float(row['chaowei'])
            data.append((dt, stage))
    data.sort(key=lambda x: x[0])
    return data


# ======================================================
# 工具：按时间取“最近 ≤ 当前”的值
# ======================================================
def lookup_time_series(series, current_dt):
    """
    series: [(dt, value), ...] 已排序
    """
    valid = [v for dt, v in series if dt <= current_dt]
    return valid[-1] if valid else series[0][1]

def read_tide_csv(path):
    """
    返回 [(datetime, stage_m), ...]
    自动处理 BOM
    期望表头包含 'date','time','chaowei'（如果有 BOM 则用 utf-8-sig）
    """
    data = []
    with open(path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        print("Tide CSV headers:", reader.fieldnames)
        for row in reader:
            dt = datetime.strptime(row['date'].strip() + " " + row['time'].strip(), "%m/%d/%Y %H:%M:%S")
            stage = float(row['chaowei'])
            data.append((dt, stage))
    data.sort(key=lambda x: x[0])
    return data


# -------------------------
# 从 tris 里提取所有外边界（并保持方向与 tris 出现的一致）
# 返回 set of (i,j) with directed edges that are boundary edges
# -------------------------
def get_boundary_edges_oriented(tris):
    """
    tris: iterable of (i,j,k) with integers (assumed to be a numpy array or list)
    返回 set((i,j),...) 的有向外边集合（方向按 tri 中出现的一致）
    """
    tris = np.asarray(tris, dtype=np.int64)
    # 无向计数
    from collections import defaultdict
    undirected_count = defaultdict(int)
    # 同时记录有向边集合（出现过的方向）
    oriented_set = set()

    for tri in tris:
        i, j, k = int(tri[0]), int(tri[1]), int(tri[2])
        for a, b in ((i, j), (j, k), (k, i)):
            key = tuple(sorted((a, b)))
            undirected_count[key] += 1
            oriented_set.add((a, b))

    boundary_edges_oriented = set()
    for key, cnt in undirected_count.items():
        if cnt == 1:
            a, b = key
            # 优先保持 tri 中的方向（有向集合里哪个存在就取哪个）
            if (a, b) in oriented_set:
                boundary_edges_oriented.add((a, b))
            elif (b, a) in oriented_set:
                boundary_edges_oriented.add((b, a))
            else:
                # 理论上不应该发生
                pass

    return boundary_edges_oriented



# -------------------------
# main: 生成 + ANUGA（修正 Domain 初始化）
# -------------------------
def main():
    NS_PATH = "ns.txt"
    NE_PATH = "ne.txt"

    print("读取 NS 文件...")
    side_map = load_ns_file(NS_PATH)

    print("读取 NE 文件...")
    elements = load_ne_file(NE_PATH, side_map)

    print("生成三角网格...")
    points, triangles, elevations = elements_to_triangles(elements, side_map)
    print(f"mesh 顶点={len(points)}, 面={len(triangles)}")

    print("导出 PLY...")
    export_ply(points, triangles, elevations, "mesh_output.ply")

    # --- 去重并重建索引 ---
    coords, tris, elevs = remove_duplicate_vertices_by_xy(points, triangles, elevations)
    print("去重后 coords shape:", coords.shape)
    print("去重后 tris  shape:", tris.shape)
    print("去重后 elevs shape:", elevs.shape)

    # --- 准备 ANUGA Domain ---
    coords_list = coords.tolist()
    tris_list = tris.tolist()

    domain = anuga.Domain(coords.tolist(), tris.tolist())

    # 打印可用的 boundary tags（ANUGA API）
    print("Available boundary tags:", domain.get_boundary_tags())

    # 打印每个标签下的边（基于我们传入的 boundary_dict）
    # from collections import defaultdict
    # tag_to_edges = defaultdict(list)
    # for edge, tag in boundary_dict.items():
    #     tag_to_edges[tag].append(edge)

    # print("=== boundary tag summary (from mesh_builder's boundary_dict) ===")
    # for tag, edges in tag_to_edges.items():
    #     print(f"tag='{tag}', num_edges={len(edges)}, sample_edges={edges[:5]}")

    # ------------------
    # 设置数量场
    # ------------------
    domain.set_quantity("elevation", elevs.tolist())
    domain.set_quantity("stage", expression="elevation")
    domain.set_quantity("friction", 0.03)

    operator = TestOperator(domain, 3.0, 0.1, 2.0)

    # ------------------
    # 读取潮位（Flather）
    # ------------------
    tide_series = read_tide_csv("test_tide.csv")
    tide_series.sort(key=lambda x: x[0])
    tide_t0 = tide_series[0][0]

    def tide_stage_func(t):
        current_dt = tide_t0 + timedelta(seconds=t)
        return lookup_time_series(tide_series, current_dt)

    # 使用 Flather 外部海面边界（返回外海水位，外侧速度设为零）
    tide_bc = anuga.Flather_external_stage_zero_velocity_boundary(domain, tide_stage_func)

    # ------------------
    # 构造 boundary_map：只对 'ocean' 施加潮汐，其它标签用 Transmissive_boundary
    # ------------------
    boundary_map = {}
    for tag in domain.get_boundary_tags():
        if tag == 'ocean':
            boundary_map[tag] = tide_bc
        else:
            # 不设置潮汐的其他边使用透传边界（保持内部值）
            boundary_map[tag] = anuga.Transmissive_boundary(domain)

    domain.set_boundary(boundary_map)

    # ------------------
    # 读取降雨（保留你原有逻辑）
    # ------------------
    rain_series = read_rain_csv("test_rain.csv")
    rain_series.sort(key=lambda x: x[0])
    rain_t0 = rain_series[0][0]
 # 仿真时间步长（秒），与你的 CSV 采样步长应一致或者你需根据 CSV 调整
    yieldstep = 300  # 5 min

    # 起始绝对时间（用降雨的第一条时间作为仿真起点）
    rain_t0 = rain_series[0][0] if rain_series else datetime.now()


    # 降雨速率函数：返回 m/s（ANUGA 需要单位 m/s）
    def rain_rate_func(t):
        current_dt = rain_t0 + timedelta(seconds=t)
        depth_m = lookup_time_series(rain_series, current_dt)  # 深度（m）对应当前步
        # 如果你的 CSV 每个记录本身就是“该步累积深度（m）”，则 rate = depth / duration
        rate_m_per_s = depth_m / yieldstep
        return rate_m_per_s

    # 把降雨 forcing 添加到 domain
    rain_forcing = anuga.Rainfall(domain, rate=rain_rate_func, default_rate=0.0)
    domain.forcing_terms.append(rain_forcing)

    # ------------------
    # 演化：每一步打印详尽的降雨信息
    # ------------------
    finaltime = yieldstep * 100  # 仅示例
    print("开始演化（只含降雨 forcing），每步会打印降雨强度及相关信息...")
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        # 仿真相对时间 t (秒) -> 绝对时间
        cur_dt = rain_t0 + timedelta(seconds=t)
        # 当前步的深度（m）与强度（m/s）
        depth_m = lookup_time_series(rain_series, cur_dt) if rain_series else 0.0
        rate_m_s = depth_m / yieldstep if yieldstep > 0 else 0.0
        # 额外以更熟悉的单位打印： mm / mm/hr
        depth_mm = depth_m * 1000.0
        rate_mm_per_hr = rate_m_s * 1000.0 * 3600.0

        print(
            f"t={t:8.1f}s | time={cur_dt} | "
            f"depth_this_step={depth_m:.6e} m ({depth_mm:.3f} mm) | "
            f"rate={rate_m_s:.6e} m/s ({rate_mm_per_hr:.3f} mm/hr)"
        )

        # 继续输出 domain 的时间步统计（可选）
        print(domain.timestepping_statistics())

    print(">>> 模拟完成")



if __name__ == "__main__":
    main()
