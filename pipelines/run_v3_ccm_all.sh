#!/usr/bin/env bash
# 等 station_v3 下载+提取完成，再依次跑 1/2/3 km CCM(N_SURR=0)，最后出对照与统计。
set -uo pipefail
PROJ="/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
cd "$PROJ"
LOG="$PROJ/data_proc/run_v3_ccm_all.log"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

say "等待 station_v3 提取完成…"
while pgrep -f "13_hcsif_run.sh" >/dev/null 2>&1; do sleep 300; done
NY=$(ls "$PROJ"/data_raw/hcsif/station_v3/hcsif_station_[0-9]*.csv 2>/dev/null | wc -l | tr -d ' ')
say "下载进程结束；station_v3 年份文件 = $NY"
if [ "$NY" -lt 20 ]; then
  say "!! 年份文件不足($NY<20)，可能下载未完整，仍继续但请核查。"
fi

run_ccm(){   # $1=BUF_R  $2=OUT_DIR
  local r="$1" out="$2"
  say ">>> buf${r} CCM 开始 → $out"
  HCSIF_BUF_R="$r" HCSIF_STAT_DIR="$PROJ/data_raw/hcsif/station_v3" \
  HCSIF_CCM_OUT="$out" HCSIF_N_SURR=0 \
    Rscript "$PROJ/pipelines/hcsif_buf${r}/01_ccm_buf${r}.R" >> "$LOG" 2>&1
  say "<<< buf${r} CCM 完成"
}
run_ccm 1000 "$PROJ/data_proc/ccm_hcsif_buf1000_v3"
run_ccm 2000 "$PROJ/data_proc/ccm_hcsif_buf2000"
run_ccm 3000 "$PROJ/data_proc/ccm_hcsif_buf3000"

say ">>> 汇总统计 + 1km 复现对照"
Rscript - <<'RS' >> "$LOG" 2>&1
suppressPackageStartupMessages(library(data.table))
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
latest <- function(dir, pat) { f <- list.files(dir, pattern=pat, full.names=TRUE)
  f <- f[!grepl("_raw", f)]; if(!length(f)) return(NA_character_); tail(sort(f),1) }
sets <- list(
  buf1000_v3 = latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_v3"), "ccm_hcsif_buf1000_vpd_.*csv$"),
  buf2000    = latest(file.path(PROJ,"data_proc/ccm_hcsif_buf2000"),    "ccm_hcsif_buf2000_vpd_.*csv$"),
  buf3000    = latest(file.path(PROJ,"data_proc/ccm_hcsif_buf3000"),    "ccm_hcsif_buf3000_vpd_.*csv$"))
cat("\n===== 三尺度 CCM 交付统计 =====\n")
for (nm in names(sets)) {
  p <- sets[[nm]]; if (is.na(p)) { cat(sprintf("%-11s : 无输出\n", nm)); next }
  d <- fread(p)
  n9 <- d[, .N, by=meteo_stat][N==9, .N]
  cat(sprintf("%-11s\n  路径: %s\n  站点数=%d | 完整9滞后站点=%d | rho中位=%.4f | |mean_coef|中位=%.4f\n",
      nm, p, uniqueN(d$meteo_stat), n9,
      median(d$rho, na.rm=TRUE), median(abs(d$mean_coef), na.rm=TRUE)))
}
# 复现对照: buf1000_v3 vs 既有 buf1000(v2)
old <- latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000"), "ccm_hcsif_buf1000_vpd_.*csv$")
if (!is.na(old) && !is.na(sets$buf1000_v3)) {
  a <- fread(sets$buf1000_v3)[, .(meteo_stat, tp, rho_new=rho, mc_new=mean_coef)]
  b <- fread(old)[, .(meteo_stat, tp, rho_old=rho, mc_old=mean_coef)]
  m <- merge(a, b, by=c("meteo_stat","tp"))
  cat(sprintf("\n===== 1km 复现对照 (v3 vs v2, 共 %d 行) =====\n", nrow(m)))
  cat(sprintf("  mean_coef: 最大绝对差=%.3e | 相关=%.6f\n",
      max(abs(m$mc_new-m$mc_old),na.rm=TRUE), cor(m$mc_new,m$mc_old,use="complete.obs")))
  cat(sprintf("  rho      : 最大绝对差=%.4f | 相关=%.6f\n",
      max(abs(m$rho_new-m$rho_old),na.rm=TRUE), cor(m$rho_new,m$rho_old,use="complete.obs")))
}
RS
say ">>> 全部完成"
touch "$PROJ/data_proc/.v3_ccm_done"
