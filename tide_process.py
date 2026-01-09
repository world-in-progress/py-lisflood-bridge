from datetime import datetime, timedelta
import csv
import numpy as np


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
            dt = datetime.strptime(row['date'].strip() + " " + row['time'].strip(), "%d/%m/%Y %H:%M:%S")
            stage = float(row['chaowei'])
            data.append((dt, stage))
    data.sort(key=lambda x: x[0])
    return data

# -------------------------------
# 线性插值函数
# -------------------------------
def lerp(y0, y1, alpha):
    """线性插值"""
    return y0 + alpha * (y1 - y0)


# -------------------------------
# 读取潮位 CSV 后处理成时间序列
# -------------------------------
def prepare_tide_series(tide_csv_path):
    tide_series = read_tide_csv(tide_csv_path)  # 返回 [(datetime, value), ...]
    tide_series.sort(key=lambda x: x[0])

    # 拆分成两个列表，便于插值
    tts = [dt for dt, v in tide_series]  # 时间列表
    tls = [v for dt, v in tide_series]   # 潮位列表

    return tts, tls


# -------------------------------
# 构造 tide_stage_func（线性插值）
# -------------------------------
def make_tide_stage_func(tts, tls):
    def tide_stage_func(t):
        current_time = tts[0] + timedelta(seconds=t)

        # 找到当前时间所在区间索引
        current_tide_idx = np.searchsorted(tts, current_time) - 1
        if current_tide_idx < 0:
            return tls[0]
        if current_tide_idx >= len(tls) - 1:
            return tls[-1]

        t0 = tts[current_tide_idx]
        t1 = tts[current_tide_idx + 1]
        y0 = tls[current_tide_idx]
        y1 = tls[current_tide_idx + 1]

        alpha = (current_time - t0) / (t1 - t0)
        return lerp(y0, y1, alpha)

    return tide_stage_func

