#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00_landcover_break_buf1000.R — CCM 前置筛选：识别缓冲区森林占比出现结构性
# 断点(疑似砍树/种树)的站点
#
# 背景: 之前(01_ccm_vpd.R / 03_run_full.R) 用 strucchange 直接在周尺度插值
# SIF 上做突变检测,效果不佳(679/898 站被判有突变——把干旱年/正常年际波动
# 也当成了断点,阈值形同虚设)。改法: 断点检测放到"地表覆盖"这个更干净、
# 有物理意义的信号上——GLC_FCS30D 2000-2022 逐年 1km 缓冲区森林类占比
# (data_raw/covariates_1km/glc_station_<年>.csv, LC04-LC13 十个森林细类
#  求和), 年际噪声远小于 SIF, 断点=真实地表覆盖切换的可能性更高。
#
# 方法: 每站 23 个年度森林占比点, 用 Chow 断点扫描(单变点最小二乘, 候选点
# 留 3 年边界)找最优断点 k*, 置换检验(打乱年份顺序 999 次)给 p 值——与本
# 项目 59_changepoint_test.R 同一套思路,不依赖 strucchange 包(未安装)。
# 判据: p<0.01 且断点前后均值差 > MAG_THRESH(默认 0.08,即 8 个百分点)。
#
# 输出: data_proc/output_hcsif_buf1000/landcover_break_stations.csv
#       每站一行: has_break, break_year, magnitude, pval, forest_frac 起止
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")

GLC_DIR     <- "data_raw/covariates_1km"
OUT_DIR     <- "data_proc/output_hcsif_buf1000"
YEARS       <- 2000:2022
FOREST_COLS <- sprintf("LC%02d", 4:13)   # 51,52,61,62,71,72,81,82,91,92 十个森林细类
MAG_THRESH  <- as.numeric(Sys.getenv("LC_MAG_THRESH", "0.08"))
P_THRESH    <- as.numeric(Sys.getenv("LC_P_THRESH", "0.01"))
NPERM       <- 999L
MARGIN      <- 3L   # 候选断点两端各留 3 年,保证段内均值可估
set.seed(42)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# --- 1. 读取逐年 GLC 站点占比,拼成 (stat_id, year, forest_frac) ------------
log_msg("读取 GLC_FCS30D 逐年站点森林占比 (", length(YEARS), " 年)...")
glc <- rbindlist(lapply(YEARS, function(y) {
  f <- file.path(GLC_DIR, sprintf("glc_station_%d.csv", y))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = c("stat_id", "year", FOREST_COLS))
  d[, forest_frac := rowSums(.SD), .SDcols = FOREST_COLS]
  d[, c("stat_id", "year", "forest_frac"), with = FALSE]
}))
log_msg("读到 ", uniqueN(glc$stat_id), " 站 x ", uniqueN(glc$year), " 年")

# --- 2. 断点扫描函数(Chow 检验风格,最小二乘) -------------------------------
chow_scan <- function(x) {
  n <- length(x)
  cand <- (MARGIN + 1):(n - MARGIN)
  if (length(cand) < 1) return(c(k = NA, mag = 0, f = -Inf))
  best <- list(k = NA_integer_, mag = 0, f = -Inf)
  for (k in cand) {
    m1 <- mean(x[1:k]); m2 <- mean(x[(k+1):n])
    rss1 <- sum((x[1:k] - m1)^2) + sum((x[(k+1):n] - m2)^2)
    rss0 <- sum((x - mean(x))^2)
    f <- if (rss1 > 0) ((rss0 - rss1) / 1) / (rss1 / (n - 2)) else Inf
    if (f > best$f) best <- list(k = k, mag = abs(m2 - m1), f = f)
  }
  c(k = best$k, mag = best$mag, f = best$f)
}

# --- 3. 逐站扫描 + 置换检验 -------------------------------------------------
stations <- sort(unique(glc$stat_id))
log_msg("逐站扫描 (", length(stations), " 站, 每站 ", NPERM, " 次置换)...")
res <- rbindlist(lapply(stations, function(sid) {
  d <- glc[stat_id == sid][order(year)]
  x <- d$forest_frac
  yrs <- d$year
  if (length(x) < 2 * MARGIN + 2 || anyNA(x))
    return(data.table(stat_id = sid, n_year = length(x), has_break = NA,
                       break_year = NA_integer_, magnitude = NA_real_,
                       pval = NA_real_, forest_frac_2000 = NA_real_,
                       forest_frac_2022 = NA_real_))
  obs <- chow_scan(x)
  fperm <- replicate(NPERM, chow_scan(sample(x))["f"])
  pval <- (1 + sum(fperm >= obs["f"] - 1e-9)) / (NPERM + 1)
  data.table(stat_id = sid, n_year = length(x),
             has_break = !is.na(pval) && pval < P_THRESH && obs["mag"] > MAG_THRESH,
             break_year = if (!is.na(obs["k"])) yrs[obs["k"]] else NA_integer_,
             magnitude = obs["mag"], pval = pval,
             forest_frac_2000 = x[1], forest_frac_2022 = x[length(x)])
}))

n_break <- sum(res$has_break, na.rm = TRUE)
log_msg("疑似地表覆盖断点站: ", n_break, " / ", nrow(res),
        sprintf(" (%.1f%%)", 100 * n_break / nrow(res)),
        "  [判据: p<", P_THRESH, " 且断点前后森林占比差>", MAG_THRESH, "]")
log_msg("magnitude 分位: ", paste(sprintf("%s=%.3f", c("P50","P75","P90","P95"),
        quantile(res$magnitude, c(.5,.75,.9,.95), na.rm = TRUE)), collapse = "  "))

fwrite(res, file.path(OUT_DIR, "landcover_break_stations.csv"))
log_msg("已写出 landcover_break_stations.csv")

# --- 4. 抽样可视化: 断点幅度最大的 12 站 ------------------------------------
top <- res[has_break == TRUE][order(-magnitude)][1:min(12, n_break)]
if (nrow(top) > 0) {
  pd <- glc[stat_id %in% top$stat_id]
  pd <- merge(pd, top[, .(stat_id, break_year)], by = "stat_id")
  p <- ggplot(pd, aes(year, forest_frac)) +
    geom_line() + geom_point(size = 0.8) +
    geom_vline(aes(xintercept = break_year), color = "red", linetype = "dashed") +
    facet_wrap(~stat_id, scales = "free_y") +
    labs(title = "疑似地表覆盖断点站点(森林占比年际序列, 红线=检测到的断点年)",
         x = "年份", y = "缓冲区森林占比") +
    theme_minimal(base_size = 11)
  ggsave(file.path(OUT_DIR, "landcover_break_examples.png"), p, width = 12, height = 9, dpi = 150)
  log_msg("已写出 landcover_break_examples.png (magnitude 最大的 ", nrow(top), " 站)")
}
log_msg("完成")
