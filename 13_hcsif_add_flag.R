#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_hcsif_add_flag.R
#
# 给站点提取结果追加 heavy_masked 标记列。
#
# 定义: 某站点在超过 THRESH_PCT 比例的时相上「750 m 缓冲区内全部像元都被掩膜」
#       (n_valid == 0)，即点值与邻域均值同时缺失，则标记为 TRUE。
#
# 为什么需要这一列: 这类站点集中在大城市(北京 4 个站、石家庄、济南、青岛、
# 南宁等)，而不是干旱区 —— 实测西部(经度<105)重度站占 4.9%、东部占 3.9%，
# 基本无差异。也就是说缺失最严重的恰好是城市化程度最高的站点，
# 对城市热环境研究构成选择性偏倚，必须能被单独识别和剔除。
#
# 用法:
#   Rscript 13_hcsif_add_flag.R [station_dir]
#   默认 data_raw/hcsif/station
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

PROJ       <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
THRESH_PCT <- 50

args     <- commandArgs(trailingOnly = TRUE)
STAT_DIR <- if (length(args) >= 1) args[1] else file.path(PROJ, "data_raw/hcsif/station")

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有找到提取结果: ", STAT_DIR)
log_msg("目录: ", STAT_DIR, "  (", length(fs), " 个年度文件)")

# --- 1. 汇总每站的全掩膜比例 ------------------------------------------------
all <- rbindlist(lapply(fs, function(f) fread(f, select = c("meteo_stat", "n_valid"))))
st  <- all[, .(n_tot = .N, n_allmask = sum(n_valid == 0)), by = meteo_stat]
st[, pct_allmask := round(100 * n_allmask / n_tot, 2)]
st[, heavy_masked := pct_allmask > THRESH_PCT]

log_msg("站点 ", nrow(st), " 个; heavy_masked = TRUE: ", sum(st$heavy_masked),
        sprintf(" (%.1f%%)", 100 * mean(st$heavy_masked)))

flag <- st[, .(meteo_stat, heavy_masked)]

# --- 2. 逐年追加列(先写临时文件再替换，避免中途失败损坏原文件) ---------------
for (f in fs) {
  d <- fread(f)
  if ("heavy_masked" %in% names(d)) d[, heavy_masked := NULL]
  d <- merge(d, flag, by = "meteo_stat", all.x = TRUE, sort = FALSE)
  d[is.na(heavy_masked), heavy_masked := FALSE]
  setorder(d, date, meteo_stat)

  tmp <- paste0(f, ".tmp")
  fwrite(d, tmp)
  if (file.info(tmp)$size > 0) file.rename(tmp, f) else stop("写入失败: ", f)
}
log_msg("已为 ", length(fs), " 个年度文件追加 heavy_masked 列")

# --- 3. 站点级汇总另存一份 --------------------------------------------------
fp <- file.path(STAT_DIR, "station_mask_flag.csv")
fwrite(st[order(-pct_allmask)], fp)
log_msg("写出站点级汇总: ", fp)

cat("\n=== heavy_masked 站点(按全掩膜比例排序) ===\n")
print(st[heavy_masked == TRUE][order(-pct_allmask)][, .(meteo_stat, n_allmask, n_tot, pct_allmask)])
