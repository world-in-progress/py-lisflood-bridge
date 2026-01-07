import os
os.environ['OMP_NUM_THREADS'] = '16'
from mesh_build import build_mesh_from_ns_ne
import anuga
from datetime import datetime, timedelta
from anuga import Domain, Operator
from rain_process import read_rain_csv, process_rain_data
from tide_process import read_tide_csv, prepare_tide_series, make_tide_stage_func 
from infiltration_operator import HortonInfiltrationByTypeOperator

# ======================================================
# 工具：按时间取“最近 ≤ 当前”的值
# ======================================================
def lookup_time_series(series, current_dt):
    """
    series: [(dt, value), ...] 已排序
    """
    valid = [v for dt, v in series if dt <= current_dt]
    return valid[-1] if valid else series[0][1]


def main():

    ns_path = "ns.txt"
    ne_path = "ne.txt"
    rain_path = "test_rain.csv"
    tide_path = "test_tide.csv"
    # ------------------
    # 构造网格（Mesh）
    # ------------------

    coords, tris, elevs ,boundary, tri_types= build_mesh_from_ns_ne(ns_path, ne_path)

    domain = anuga.Domain(coords.tolist(), tris.tolist(),boundary=boundary)
    
    # ------------------
    # 构建模拟区（Domain）
    # ------------------

    domain.set_quantity("elevation", elevs.tolist())
    domain.set_quantity("stage", expression="elevation")
    domain.set_quantity("friction", 0.03)

    # ------------------
    # 读取潮位（Flather）
    # ------------------
    tts, tls = prepare_tide_series(tide_path)
    tide_stage_func = make_tide_stage_func(tts, tls)

    # 使用 Flather 外部海面边界（返回外海水位，外侧速度设为零）
    tide_bc = anuga.Flather_external_stage_zero_velocity_boundary(domain, tide_stage_func)

    boundary_map = {}
    for tag in domain.get_boundary_tags():
        if tag == 'ocean':
            boundary_map[tag] = tide_bc
        else:
            # 不设置潮汐的其他边使用透传边界（保持内部值）
            boundary_map[tag] = anuga.Transmissive_boundary(domain)

    domain.set_boundary(boundary_map)

    # ------------------
    # 处理降雨：每一步打印详尽的降雨信息
    # ------------------
    rain_rate_func, rain_t0,rts = process_rain_data(rain_path)

    # 把降雨 forcing 添加到 domain
    rain_forcing = anuga.Rainfall(domain, rate=rain_rate_func, default_rate=0.0)#anuga接受的降雨率变量单位是mm/s
    domain.forcing_terms.append(rain_forcing)

    # ------------------
    # 处理下渗：利用降雨和type信息
    # ------------------

    horton_op = HortonInfiltrationByTypeOperator(domain, tri_types)
    domain.forcing_terms.append(horton_op)

    # ------------------
    # 演化：每一步打印详尽的降雨信息
    # ------------------
    yieldstep = 300   # 秒 —— 你要求的仿真步长
    finaltime = rts[-1]  # 仿真终点（秒）

    print("开始演化（只含降雨 forcing），每步会打印降雨强度及累积降雨深度与潮位...")

    cumulative_rain_mm = 0.0  # 初始化累积降雨（单位 mm）

    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        cur_dt = rain_t0 + timedelta(seconds=t)

        # ⭐ 当前步雨强（mm/s） → ANUGA 内部会自动转 m/s
        rate_mm_s = rain_rate_func(t)
        rate_mm_hr = rate_mm_s * 3600.0         # mm/hr
        depth_step_mm = rate_mm_s * yieldstep  # 本步降雨深度 mm

        # 累积降雨深度
        cumulative_rain_mm += depth_step_mm

        # ⭐ 当前潮位（m）
        tide_m = tide_stage_func(t)

        print(
            f"t={t:8.1f}s | time={cur_dt} | "
            f"rain_rate={rate_mm_s:.3e} mm/s ({rate_mm_hr:.3f} mm/hr) | "
            f"step_depth={depth_step_mm:.3f} mm | "
            f"cumulative_depth={cumulative_rain_mm:.3f} mm | "
            f"tide={tide_m:.3f} m"
        )

        # 打印 ANUGA 内部水深/速度统计
        stats = domain.timestepping_statistics()
        print(stats)


if __name__ == "__main__":
    main()

