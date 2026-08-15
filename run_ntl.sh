#!/usr/bin/env bash
# 夜间灯光流式处理: 下载(81MB) -> LZMA解压(11.7GB) -> 提取 -> 立即删除
# 注意: zip 成员用 LZMA(方法14)压缩且解压后 >4GB,
#   - GDAL 的 /vsizip/ 不支持方法14
#   - macOS 自带 unzip 处理不了 >4GB 成员, unzip -t 会假报损坏
# 因此用 Python zipfile 解压并自行校验。
set -uo pipefail
PROJ="/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
TMP="$PROJ/data_raw/ntl/tmp"; mkdir -p "$TMP"
URLS="${1:-$PROJ/data_raw/ntl/ntl_urls.txt}"
while read -r url name; do
  y=$(echo "$name" | grep -oE '[0-9]{4}' | tail -1)
  [ -f "$PROJ/data_raw/covariates_1km/ntl_station_${y}.csv" ] && { echo "== $y 已完成，跳过"; continue; }
  free=$(df -g "$PROJ" | awk 'NR==2{print $4}')
  [ "$free" -lt 14 ] && { echo "!! 磁盘剩余 ${free}G 不足 14G，中止"; break; }
  echo "===== $y (磁盘 ${free}G) ====="
  z="$TMP/$name"; ok=0
  for k in 1 2 3; do
    curl -sL --retry 5 --retry-delay 5 -o "$z" "$url" </dev/null || { sleep 20; continue; }
    python3 -c "import zipfile,sys; z=zipfile.ZipFile('$z'); sys.exit(0 if z.infolist() else 1)" </dev/null 2>/dev/null && { ok=1; break; }
    echo "  第${k}次下载/校验失败"; rm -f "$z"; sleep 20
  done
  [ "$ok" -eq 0 ] && { echo "!! $y 下载失败，跳过"; continue; }
  t="$TMP/ntl_$y.tif"
  python3 -c "
import zipfile,shutil
z=zipfile.ZipFile('$z'); i=[m for m in z.infolist() if m.filename.lower().endswith('.tif')][0]
with z.open(i) as f, open('$t','wb') as o: shutil.copyfileobj(f,o,1<<24)
" </dev/null || { echo "!! $y 解压失败"; rm -f "$z" "$t"; continue; }
  Rscript "$PROJ/16_ntl_1km.R" "$y" "$t" </dev/null > "$TMP/ext.log" 2>&1 || echo "  提取失败"
  grep -E "^\[|!!" "$TMP/ext.log" || true
  rm -f "$z" "$t"
done < "$URLS"
echo "=== 夜间灯光全部结束 ==="
