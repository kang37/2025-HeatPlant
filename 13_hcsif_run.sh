#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 13_hcsif_run.sh — HCSIF 流式下载 + 站点提取 + 立即清理
#
# 为什么要流式: 本机可用磁盘约 32 GiB，而 5-9 月全量约 44 GB。
# 逐年下载 -> 提取 924 个站点 -> 立刻删除栅格，峰值占用约 2 GB。
#
# 用法:
#   ./13_hcsif_run.sh <year|all> [manifest.txt]
#
# manifest 每行形如:
#   https://download.scidb.cn/download?fileId=<hash>&path=/V3/2000/2000196.tif&fileName=2000196.tif
# 注意文件名在查询串里，不能用 basename(URL) —— 那样只会得到 "download"。
#
# 只下 .tif: 实测 .tif 自带地理参考(EPSG:4326)，.tfw / .aux.xml 是冗余的，
# 且清单里只有 679/678 个(少于 1010 个 .tif)，强求会缺件。
# ---------------------------------------------------------------------------

set -uo pipefail

PROJ="/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
RAW_DIR="$PROJ/data_raw/hcsif/tmp"
OUT_DIR="${HCSIF_OUT_DIR:-$PROJ/data_raw/hcsif/station}"   # 可用环境变量覆盖
MANIFEST_DEFAULT="$PROJ/data_raw/hcsif/sif_download_link.txt"
MIN_FREE_GB=6          # 低于此值就停，避免把盘写满
DOY_MIN=121            # 5 月 1 日
DOY_MAX=273            # 9 月 30 日
# 服务器对单连接限速约 140 KB/s，但总吞吐随连接数线性增长(实测 3 并发 = 470 KB/s)。
# 取 10 并发约 1.4 MB/s，既能把 44 GB 压到约 9 小时，又不至于把公共服务打爆。
PARALLEL=10            # 并发下载数(约 1 GB 在途)

TARGET="${1:?用法: ./13_hcsif_run.sh <year|all> [manifest.txt]}"
MANIFEST="${2:-$MANIFEST_DEFAULT}"
[ -s "$MANIFEST" ] || { echo "找不到清单: $MANIFEST" >&2; exit 1; }

mkdir -p "$RAW_DIR" "$OUT_DIR"
free_gb() { df -g "$PROJ" | awk 'NR==2 {print $4}'; }

# 从清单里筛出某年 5-9 月的 .tif，输出 "url 空格 filename"
# (URL 和 HCSIF 文件名都不含空格，所以空格分隔是安全的)
select_year() {
  awk -F'fileName=' -v y="$1" -v lo="$DOY_MIN" -v hi="$DOY_MAX" '
    {
      fn = $2
      sub(/[&#].*$/, "", fn)                 # fileName 后面可能还有别的参数
      if (fn !~ /^[0-9][0-9][0-9][0-9][0-9][0-9][0-9]\.tif$/) next
      yr  = substr(fn, 1, 4)
      doy = substr(fn, 5, 3) + 0
      if (yr == y && doy >= lo && doy <= hi) print $1 "fileName=" fn " " fn
    }
  ' "$MANIFEST"
}

run_year() {
  local YEAR="$1"
  local SEL; SEL="$(mktemp)"
  select_year "$YEAR" > "$SEL"

  local N; N=$(wc -l < "$SEL" | tr -d ' ')
  echo "==> ${YEAR}: 命中 ${N} 个时相, 磁盘剩余 $(free_gb)G"
  if [ "$N" -eq 0 ]; then echo "   跳过(清单无匹配)"; rm -f "$SEL"; return 0; fi
  if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
    echo "   剩余磁盘不足 ${MIN_FREE_GB}G，中止。" >&2; rm -f "$SEL"; return 1
  fi

  # 并发下载。xargs -n 2 把每行的 url / filename 作为 $0 / $1 传给 bash -c。
  # 已完整存在的文件直接跳过，便于中断后重跑。
  #
  # 多轮重试: 实测服务端在持续大流量后会限流断连(2002 年曾 15/19 失败),
  # 而 curl 自身的 --retry 兜不住这类连接级失败。每轮之间退避，逐轮只补缺失的。
  local GOT=0 PASS
  for PASS in 1 2 3; do
    RAW_DIR="$RAW_DIR" MIN_BYTES=10000000 xargs -P "$PARALLEL" -n 2 bash -c '
      url="$0"; fn="$1"; dest="$RAW_DIR/$fn"
      if [ -s "$dest" ]; then exit 0; fi
      if curl -fsS --connect-timeout 30 --retry 5 --retry-delay 5 \
              --retry-all-errors --retry-max-time 180 -o "$dest.part" "$url"; then
        sz=$(wc -c < "$dest.part")
        if [ "$sz" -lt "$MIN_BYTES" ]; then
          rm -f "$dest.part"; echo "   BAD  $fn (仅 $sz 字节)" >&2
        else
          mv "$dest.part" "$dest"; echo "   ok   $fn ($(du -h "$dest" | cut -f1))"
        fi
      else
        rm -f "$dest.part"; echo "   FAIL $fn" >&2
      fi
    ' < "$SEL"

    GOT=$(ls -1 "$RAW_DIR"/${YEAR}*.tif 2>/dev/null | wc -l | tr -d ' ')
    echo "==> ${YEAR}: 第 ${PASS} 轮后到位 ${GOT}/${N}"
    [ "$GOT" -ge "$N" ] && break
    if [ "$PASS" -lt 3 ]; then
      echo "   缺 $((N-GOT)) 个，退避 60 秒后重试"; sleep 60
    fi
  done

  # 完整性闸门: 不完整就保留已下载的文件并跳过提取，
  # 否则会像 2002 年那样用 4/19 个时相生成一个看似正常、实则残缺的年度文件。
  if [ "$GOT" -lt "$N" ]; then
    echo "!! ${YEAR}: 仍缺 $((N-GOT)) 个时相，跳过提取并保留已下载文件(重跑本年可续)" >&2
    rm -f "$SEL"; return 1
  fi

  Rscript "$(dirname "$0")/13_hcsif_extract.R" "$YEAR" "$RAW_DIR" "$OUT_DIR" || {
    echo "   提取失败，保留原始文件以便排查" >&2; rm -f "$SEL"; return 1; }

  rm -f "$RAW_DIR"/${YEAR}*.tif "$RAW_DIR"/${YEAR}*.tfw "$RAW_DIR"/${YEAR}*.aux.xml
  echo "==> ${YEAR}: 完成并清理, 磁盘剩余 $(free_gb)G"
  rm -f "$SEL"
}

if [ "$TARGET" = "all" ]; then
  for y in $(seq 2000 2022); do
    run_year "$y" || echo "!! ${y} 未完成，继续下一年"
  done
  echo "=== 全部年份处理结束 ==="
  ls -1 "$OUT_DIR"/hcsif_station_*.csv 2>/dev/null | wc -l | xargs echo "已生成年份文件数:"
else
  run_year "$TARGET"
fi
