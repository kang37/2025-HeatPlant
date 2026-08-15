#!/usr/bin/env bash
# GLC_FCS30D 流式处理: 列瓦片 -> 只取覆盖站点的瓦片 -> 提取 -> 立即删除
#
# 为什么不下整包: Zenodo 上按 5° 经度带打包，12 个包约 50 GB、60+ 小时。
# 每个包内含 86 个 5°x5° 瓦片，站点只覆盖其中 46 个，按字节范围单独取回
# 约 9 GB。zip 是 ZIP64 且成员为 DEFLATE，由 17_glc_fetch_tile.py 处理。
set -uo pipefail
PROJ="/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
TMP="$PROJ/data_raw/glc/tmp"; PARTS="$PROJ/data_raw/glc/parts"; LIST="$PROJ/data_raw/glc/lists"
mkdir -p "$TMP" "$PARTS" "$LIST"
cd "$PROJ"

BANDS="70 80 90 100 110 120 130"   # 包按 10° 命名，每包含两列 5° 瓦片；覆盖经度 70-135

# --- 1. 列出各包的瓦片清单(缓存) ---
for e in $BANDS; do
  f="$LIST/E${e}.txt"
  [ -s "$f" ] && continue
  echo "== 列清单 E${e}"
  python3 17_glc_fetch_tile.py list "E${e}" > "$f" 2>"$LIST/E${e}.err" || { echo "  失败"; rm -f "$f"; }
  sleep 5
done

# --- 2. 站点需要哪些瓦片 ---
Rscript 17_glc_match.R || { echo "!! 瓦片匹配失败"; exit 1; }

# --- 3. 逐瓦片取回并提取 ---
# 两轮: zip 字节范围偶尔会取回坏数据(zlib "invalid block type")，
# 单纯跳过会永久缺该瓦片，故第一轮结束后对仍缺的瓦片再试一次。
#
# 注意: 这里所有子进程都要 </dev/null 并把输出写文件后再读。
#   曾三次踩坑: 子进程继承 while-read 的 stdin 会吞掉瓦片清单;
#   而 `cmd | tail -2` 在 pipefail 下会因 SIGPIPE 把成功误判为失败，删掉好文件。
for pass in 1 2; do
  echo "########## 第 $pass 轮 ##########"
  tail -n +2 "$PROJ/data_raw/glc/needed_tiles.csv" | while IFS=, read -r band name csize nst tid done_flag; do
    band=$(echo "$band" | tr -d '"'); name=$(echo "$name" | tr -d '"'); tid=$(echo "$tid" | tr -d '"')
    [ -s "$PARTS/glc_long_${tid}.csv" ] && { [ "$pass" = 1 ] && echo "== $tid 已完成，跳过"; continue; }
    free=$(df -g "$PROJ" | awk 'NR==2{print $4}')
    [ "$free" -lt 12 ] && { echo "!! 磁盘剩余 ${free}G 不足，中止"; break; }
    echo "===== $tid  ($(echo "scale=0;$csize/1048576"|bc) MB压缩, ${nst}站, 磁盘${free}G)"
    t="$TMP/${tid}.tif"
    if ! python3 17_glc_fetch_tile.py fetch "$band" "$name" "$t" </dev/null > "$TMP/fetch.log" 2>&1; then
      echo "  取回失败:"; tail -2 "$TMP/fetch.log"; rm -f "$t" "$t.part"; continue
    fi
    tail -1 "$TMP/fetch.log"
    [ -s "$t" ] || { echo "  文件为空，跳过"; continue; }
    Rscript 17_glc_extract.R "$t" "$tid" "$PARTS" </dev/null > "$TMP/ext.log" 2>&1 || echo "  提取失败"
    grep -E "^\[|注意" "$TMP/ext.log" || true
    rm -f "$t"
  done
done
echo "=== GLC 全部瓦片结束 ==="
Rscript 17_glc_match.R </dev/null > "$TMP/final.log" 2>&1
head -3 "$TMP/final.log"
ls -1 "$PARTS"/glc_long_*.csv 2>/dev/null | wc -l | xargs echo "分片文件数:"
