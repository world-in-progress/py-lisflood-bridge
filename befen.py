
"""
功能：
- 读取 NS / NE 文件
- 生成支持悬挂点的三角网格
- 在【三角化阶段】直接生成 boundary（基于 NS 拓扑）
- 顶点去重时同步 remap boundary
- 输出可直接用于 anuga.Domain 的 points / triangles / boundary
"""

import numpy as np
from typing import List, Dict, Tuple

# ============================================================
# NS 相关
# ============================================================

class HydroSide:
    def __init__(
        self,
        index: int,
        direction: int,
        bounds: Tuple[float, float, float, float],
        left_grid_index: int,
        right_grid_index: int,
        altitude: float = 0.0,
        type_: int = 0
    ):
        self.index = index
        self.direction = direction
        self.bounds = bounds
        self.left_grid_index = left_grid_index
        self.right_grid_index = right_grid_index
        self.altitude = altitude
        self.type = type_
        # 是否为边界边
        self.is_boundary = (left_grid_index == 0 or right_grid_index == 0)


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

    if direction == 1:  # 水平边
        min_x = cx - length / 2
        max_x = cx + length / 2
        min_y = max_y = cy
    else:               # 垂直边
        min_y = cy - length / 2
        max_y = cy + length / 2
        min_x = max_x = cx

    bounds = (min_x, min_y, max_x, max_y)
    return HydroSide(index, direction, bounds, left, right, cz, type_)


def load_ns_file(path: str) -> Dict[int, HydroSide]:
    side_map = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            row = [float(p) if "." in p else int(p) for p in parts]
            s = parse_ns_line(row)
            side_map[s.index] = s
    return side_map


# ============================================================
# NE 相关
# ============================================================

class HydroElement:
    def __init__(
        self,
        index: int,
        bounds: Tuple[float, float, float, float],
        left_edges: List[int],
        right_edges: List[int],
        bottom_edges: List[int],
        top_edges: List[int],
        altitude: float,
        type_: int
    ):
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
        return ((min_x + max_x) / 2, (min_y + max_y) / 2, self.altitude)


def parse_ne_line(row: list, side_map: Dict[int, HydroSide]) -> HydroElement:
    idx = 0
    index = int(row[idx]); idx += 1
    Ln = int(row[idx]); idx += 1
    Rn = int(row[idx]); idx += 1
    Bn = int(row[idx]); idx += 1
    Tn = int(row[idx]); idx += 1

    left_edges = [int(row[idx+i]) for i in range(Ln)]; idx += Ln
    right_edges = [int(row[idx+i]) for i in range(Rn)]; idx += Rn
    bottom_edges = [int(row[idx+i]) for i in range(Bn)]; idx += Bn
    top_edges = [int(row[idx+i]) for i in range(Tn)]; idx += Tn

    cx = float(row[idx]); idx += 1
    cy = float(row[idx]); idx += 1
    cz = float(row[idx]); idx += 1
    type_ = int(row[idx])

    min_x = side_map[left_edges[0]].bounds[0] if left_edges else cx
    max_x = side_map[right_edges[0]].bounds[2] if right_edges else cx
    min_y = side_map[bottom_edges[0]].bounds[1] if bottom_edges else cy
    max_y = side_map[top_edges[0]].bounds[3] if top_edges else cy

    return HydroElement(
        index,
        (min_x, min_y, max_x, max_y),
        left_edges, right_edges, bottom_edges, top_edges,
        cz, type_
    )


def load_ne_file(path: str, side_map: Dict[int, HydroSide]) -> List[HydroElement]:
    elems = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(",")]
            row = [float(p) if "." in p else int(p) for p in parts]
            elems.append(parse_ne_line(row, side_map))
    return elems


# ============================================================
# 三角化 + boundary 生成（核心）
# ============================================================

def tri_area(a, b, c):
    return abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])) / 2


def elements_to_triangles(elements, side_map):
    points = []
    elevations = []
    triangles = []
    tri_boundary_flags = []   # 每个 triangle 一个 dict: edge_id -> tag
    point_map = {}

    def add_point(x, y, z):
        key = (round(x, 9), round(y, 9))
        if key in point_map:
            return point_map[key]
        idx = len(points)
        points.append((x, y))
        elevations.append(z)
        point_map[key] = idx
        return idx

    def is_boundary(edge_ids):
        return any(side_map[sid].is_boundary for sid in edge_ids)

    for elem in elements:
        min_x, min_y, max_x, max_y = elem.bounds
        z = elem.altitude

        p0 = add_point(min_x, min_y, z)
        p1 = add_point(max_x, min_y, z)
        p2 = add_point(max_x, max_y, z)
        p3 = add_point(min_x, max_y, z)

        # triangle 1
        if tri_area(points[p0], points[p1], points[p2]) > 1e-12:
            triangles.append((p0, p1, p2))
            tri_boundary_flags.append({
                2: 'ocean' if is_boundary(elem.bottom_edges) else None,
                1: 'ocean' if is_boundary(elem.right_edges) else None
            })

        # triangle 2
        if tri_area(points[p0], points[p2], points[p3]) > 1e-12:
            triangles.append((p0, p2, p3))
            tri_boundary_flags.append({
                0: 'ocean' if is_boundary(elem.top_edges) else None,
                2: 'ocean' if is_boundary(elem.left_edges) else None
            })

    return points, triangles, elevations, tri_boundary_flags


# ============================================================
# 去重 + boundary remap
# ============================================================

def remove_duplicate_vertices(points, triangles, elevations, tri_boundary_flags):
    pts = np.asarray(points)
    elevs = np.asarray(elevations)
    tris = np.asarray(triangles, dtype=int)

    keys = np.round(pts, 9)
    _, unique_idx, inverse = np.unique(
        keys.view([('x', float), ('y', float)]),
        return_index=True, return_inverse=True
    )

    new_pts = pts[unique_idx]
    new_elevs = np.zeros(len(unique_idx))
    counts = np.zeros(len(unique_idx))

    for i, ni in enumerate(inverse):
        new_elevs[ni] += elevs[i]
        counts[ni] += 1
    new_elevs /= counts

    new_tris = []
    new_boundary = {}

    for old_tid, tri in enumerate(tris):
        nt = tuple(inverse[i] for i in tri)
        tid = len(new_tris)
        new_tris.append(nt)
        for edge_id, tag in tri_boundary_flags[old_tid].items():
            if tag:
                new_boundary[(tid, edge_id)] = tag

    return new_pts, np.array(new_tris), new_elevs, new_boundary


# ============================================================
# 对外接口（你直接用这个）
# ============================================================

def build_mesh_from_ns_ne(ns_path, ne_path):
    side_map = load_ns_file(ns_path)
    elements = load_ne_file(ne_path, side_map)

    pts, tris, elevs, tri_bflags = elements_to_triangles(elements, side_map)
    pts, tris, elevs, boundary = remove_duplicate_vertices(
        pts, tris, elevs, tri_bflags
    )

    print(f"mesh: points={len(pts)}, triangles={len(tris)}, boundary edges={len(boundary)}")
    return pts, tris, elevs, boundary