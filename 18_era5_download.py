#!/usr/bin/env python3
"""
18_era5_download.py — 从 CDS 下载 ERA5-Land 日尺度辐射(可选风速)

为什么用 ERA5-Land 而不是 CHELSA:
  CHELSA 是 1 km 逐日全球 COG，但没有区域子集接口，每天必须传约 71 MB 的
  中国窗口，实测 59 秒/天 → 2000-2022 年 5-9 月共 3519 天需 58 小时。
  ERA5-Land 支持服务端区域裁剪与日统计聚合，只传中国范围的日值，
  数据量小两个数量级。代价是分辨率 9 km（CHELSA 为 1 km）。
  对辐射而言这个取舍合理: 太阳辐射的空间变异主要由纬度、地形和云场驱动，
  且 1000 m 缓冲区本来就只覆盖 CHELSA 的 3-4 个像元，没用上 1 km 的精细结构。

前置: 需要 CDS 账号与 API key，写入 ~/.cdsapirc:
    url: https://cds.climate.copernicus.eu/api
    key: <你的 Personal Access Token>

用法:
    python3 18_era5_download.py rsds          # 仅辐射
    python3 18_era5_download.py rsds wind     # 辐射 + 风速
"""
import os, sys, time

import cdsapi

PROJ = "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
OUT = os.path.join(PROJ, "data_raw/era5_land")
YEARS = [str(y) for y in range(2000, 2023)]
MONTHS = ["05", "06", "07", "08", "09"]          # 与 SIF 的 5-9 月窗口一致
AREA = [54, 73, 17, 136]                          # N, W, S, E — 覆盖全部 924 站

# ssrd 是累积量: ERA5-Land 从每日 00 UTC 起累积，故 23:00 的值即当日总量(J m-2)，
# 只需取 1 个时次而非 24 个。风速是瞬时量，若要日均需另取多时次，这里同样取 23:00
# 作为代表时次(仅作区域背景控制变量)。
VARS = {
    "rsds":   "surface_solar_radiation_downwards",
    "wind_u": "10m_u_component_of_wind",
    "wind_v": "10m_v_component_of_wind",
}


def fetch(client, key, year):
    var = VARS[key]
    fp = os.path.join(OUT, f"era5land_{key}_{year}.nc")
    if os.path.exists(fp) and os.path.getsize(fp) > 10000:
        print(f"  {key} {year} 已存在，跳过", flush=True)
        return True
    req = {
        "variable": [var],
        "year": year,
        "month": MONTHS,
        "day": [f"{d:02d}" for d in range(1, 32)],
        "time": ["23:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": AREA,
    }
    for attempt in range(1, 4):
        try:
            client.retrieve("reanalysis-era5-land", req).download(fp)
            print(f"  {key} {year} 完成 ({os.path.getsize(fp)/1e6:.1f} MB)", flush=True)
            return True
        except Exception as e:                     # CDS 队列拥塞时常见，退避重试
            print(f"  {key} {year} 第{attempt}次失败: {str(e)[:110]}", flush=True)
            time.sleep(60 * attempt)
    return False


def main():
    which = sys.argv[1:] or ["rsds"]
    keys = ["rsds"] if which == ["rsds"] else (
        ["rsds", "wind_u", "wind_v"] if "wind" in which else which)

    if not os.path.exists(os.path.expanduser("~/.cdsapirc")):
        sys.exit("缺少 ~/.cdsapirc —— 请先注册 CDS 并写入 url 与 key")

    os.makedirs(OUT, exist_ok=True)
    c = cdsapi.Client()
    bad = []
    for k in keys:
        print(f"===== {k} =====", flush=True)
        for y in YEARS:
            if not fetch(c, k, y):
                bad.append(f"{k}/{y}")
    print("\n=== 结束 ===")
    print("失败:", ", ".join(bad) if bad else "无")


if __name__ == "__main__":
    main()
