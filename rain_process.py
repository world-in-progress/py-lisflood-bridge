from datetime import datetime, timedelta
import csv

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
            rain_m = rain_mm  
            rain_data.append((dt, rain_m))
    return rain_data

# ------------------
# 构造与参考代码行为一致的 rain_rate_func（返回 m/s）
# ------------------
def make_rain_rate_func(rts, rqs):
    """
    返回一个函数 rain_rate_func(t)，t 为仿真开始后经过的秒数，返回单位 m/s。
    使用闭包维护 current index，行为等价于你给的参考片段（分段常值插值）。
    """
    current_idx = 0
    r_num_local = len(rts)

    def rain_rate_func(t):
        nonlocal current_idx

        # 如果在降雨记录之前或之后：返回 0
        if t < rts[0] or t > rts[-1]:
            return 0.0

        # 若 t 超过下一个时间节点，则推进索引（与参考代码的 current_rain_idx += 1 等价）
        # 用 while 以防仿真步长跳跃超过多个节点
        while current_idx < r_num_local - 2 and t > rts[current_idx + 1]:
            current_idx += 1

        # 计算该区间内的均匀降雨强度（m / s）
        dt = rts[current_idx + 1] - rts[current_idx]
        if dt <= 0:
            return 0.0

        # 区间总雨量用 rqs[current_idx + 1] （已是 m）
        rate_m_per_s = rqs[current_idx + 1] / dt
        return rate_m_per_s

    return rain_rate_func

def process_rain_data(rain_path):
    rain_series = read_rain_csv(rain_path)
    rain_series.sort(key=lambda x: x[0])

    if not rain_series:
        raise ValueError("Rain series is empty")
    
    # 起始绝对时间（用降雨的第一条时间作为仿真起点）
    rain_t0 = rain_series[0][0]

    # 把时间点变为相对秒（相对于 rain_t0）
    rain_times = [rt for rt, _ in rain_series]
    rain_values = [rv for _, rv in rain_series]  # CSV 中每条表示该时间段内的累计降雨量（你说的是 5 min 的累计）

    # rts: 秒（从 0 开始）
    rts = [(rt - rain_t0).total_seconds() for rt in rain_times]
    r_num = len(rts)

    # rqs: 对应每个时间节点的降雨量（单位：m）
    # 假设 CSV 的 rain_values 单位是 mm（常见），所以转换为 m：
    rqs = [v for v in rain_values] 

    rain_rate_func = make_rain_rate_func(rts, rqs)
    return rain_rate_func, rain_t0,rts
