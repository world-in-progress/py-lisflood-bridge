import csv
import re
from pyproj import Transformer


# ---------------------------------------
# 坐标转换器：WGS84 -> HK1980 (EPSG:2326)
# ---------------------------------------
transformer = Transformer.from_crs(
    "EPSG:4326",
    "EPSG:2326",
    always_xy=True
)


def parse_deg(value: str) -> float:
    """
    解析形如 '22.291111°N' 或 '114.213333°E' 的字符串
    返回十进制度数（float）
    """
    value = value.strip()

    match = re.match(r"([0-9.]+)°([NSEW])", value)
    if not match:
        raise ValueError(f"无法解析坐标值: {value}")

    deg = float(match.group(1))
    hemi = match.group(2)

    if hemi in ("S", "W"):
        deg = -deg

    return deg


def main():
    input_csv = "coord.csv"
    output_csv = "coord_2326.csv"

    results = []

    with open(input_csv, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 3:
                continue

            name = row[0].strip()
            lat_str = row[1]
            lon_str = row[2]

            lat = parse_deg(lat_str)
            lon = parse_deg(lon_str)

            # 转换坐标
            E, N = transformer.transform(lon, lat)

            results.append([
                name,
                lon,
                lat,
                round(E, 3),
                round(N, 3)
            ])

    # 写出新 CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Station",
            "Longitude_WGS84",
            "Latitude_WGS84",
            "Easting_2326_m",
            "Northing_2326_m"
        ])
        writer.writerows(results)

    print(f"✅ 转换完成，结果已写入: {output_csv}")


if __name__ == "__main__":
    main()
