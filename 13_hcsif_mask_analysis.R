#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_hcsif_mask_analysis.R
#
# 回答一个问题: HCSIF 的 0 值掩膜是"静态"的还是"随时间变化"的?
#
#   静态  = 一次性土地覆盖分类的产物 -> 被掩膜的站点永久不可用(点值尺度),
#           只能靠邻域均值补救，且缺失与城市化程度系统相关(选择性偏倚)。
#   时变  = 云/质控导致 -> 不同时相缺的站不同，多数站能凑出较完整的序列,
#           缺失更接近随机，对因果分析的威胁小得多。
#
# 用法: Rscript 13_hcsif_mask_analysis.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station")

fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有找到提取结果")

d <- rbindlist(lapply(fs, fread), fill = TRUE)
d[, masked := as.logical(pt_masked)]
setorder(d, date, meteo_stat)

n_t  <- uniqueN(d$date)
n_st <- uniqueN(d$meteo_stat)
cat("=== 数据规模 ===\n")
cat("年份    :", paste(sort(unique(d$year)), collapse = ", "), "\n")
cat("时相数  :", n_t, "\n")
cat("站点数  :", n_st, "\n")
cat("记录数  :", nrow(d), "\n\n")

# --- 1. 每站被掩膜的时相数 -------------------------------------------------
per_st <- d[, .(n_obs = .N, n_masked = sum(masked), n_valid = sum(!is.na(SIF))),
            by = meteo_stat]
per_st[, frac := n_masked / n_obs]

cat("=== 每站掩膜情况 ===\n")
cls <- per_st[, .(
  never     = sum(n_masked == 0),
  always    = sum(n_masked == n_obs),
  sometimes = sum(n_masked > 0 & n_masked < n_obs)
)]
cat("从不被掩膜  :", cls$never,     sprintf("(%.1f%%)\n", 100 * cls$never / n_st))
cat("始终被掩膜  :", cls$always,    sprintf("(%.1f%%)\n", 100 * cls$always / n_st))
cat("有时被掩膜  :", cls$sometimes, sprintf("(%.1f%%)\n", 100 * cls$sometimes / n_st))

cat("\n判定: ")
if (cls$sometimes == 0) {
  cat("完全静态 —— 每个站要么全程有值，要么全程被掩膜\n")
} else if (cls$sometimes / n_st < 0.05) {
  cat("以静态为主 —— 仅", cls$sometimes, "个站存在时相间变化\n")
} else {
  cat("时变 —— 有", cls$sometimes, "个站在不同时相间切换状态\n")
}

# --- 2. 时相之间掩膜集合的一致性(Jaccard) ----------------------------------
sets <- split(d[masked == TRUE]$meteo_stat, d[masked == TRUE]$date)
if (length(sets) >= 2) {
  ks <- names(sets)
  jac <- sapply(seq_len(length(ks) - 1), function(i) {
    a <- sets[[ks[i]]]; b <- sets[[ks[i + 1]]]
    if (length(union(a, b)) == 0) return(NA_real_)
    length(intersect(a, b)) / length(union(a, b))
  })
  cat("\n=== 相邻时相掩膜集合的 Jaccard 相似度 ===\n")
  cat(sprintf("均值 %.4f | 最小 %.4f | 最大 %.4f\n",
              mean(jac, na.rm = TRUE), min(jac, na.rm = TRUE), max(jac, na.rm = TRUE)))
  cat("(1.000 = 每个时相掩膜的站点完全相同, 即静态掩膜)\n")
}

# --- 3. 每个时相的掩膜站点数 -----------------------------------------------
per_t <- d[, .(n_masked = sum(masked), pct = round(100 * mean(masked), 2)), by = .(year, doy, date)]
setorder(per_t, date)
cat("\n=== 各时相掩膜站点数 ===\n")
print(per_t[, .(date, doy, n_masked, pct)])
cat(sprintf("\n跨时相极差: %d ~ %d 站 (相差 %d)\n",
            min(per_t$n_masked), max(per_t$n_masked),
            max(per_t$n_masked) - min(per_t$n_masked)))

# --- 4. 缓冲均值能补救多少 -------------------------------------------------
cat("\n=== 邻域均值的补救能力 ===\n")
cat("点值有效     :", sum(!is.na(d$SIF)),
    sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$SIF))))
cat("缓冲均值有效 :", sum(!is.na(d$SIF_buf750)),
    sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$SIF_buf750))))
rescued <- d[is.na(SIF) & !is.na(SIF_buf750), .N]
cat("被点值丢弃但缓冲可救回:", rescued,
    sprintf("(%.1f%% 的记录)\n", 100 * rescued / nrow(d)))
cat("两者都无值             :", d[is.na(SIF) & is.na(SIF_buf750), .N], "\n")

# --- 5. 与建成度的关系(选择性偏倚检验) -------------------------------------
bf_path <- file.path(PROJ, "data_raw/China_stations_buildings_with_coords30_local_rerun_filled.csv")
if (file.exists(bf_path)) {
  b <- fread(bf_path)
  b2 <- unique(b[, .(meteo_stat = as.integer(AirQualityStation),
                     bf = as.numeric(buildingFootprint))], by = "meteo_stat")[!is.na(meteo_stat)]
  m <- merge(per_st, b2, by = "meteo_stat")[!is.na(bf)]
  m[, grp := fifelse(n_masked == 0, "never", fifelse(n_masked == n_obs, "always", "sometimes"))]
  cat("\n=== 掩膜与建成度(30m 缓冲建筑占地, m2) ===\n")
  print(m[, .(n = .N, median_bf = round(median(bf), 1), mean_bf = round(mean(bf), 1)), by = grp][order(grp)])
  if (uniqueN(m$grp) >= 2 && all(c("never", "always") %in% m$grp)) {
    w <- wilcox.test(bf ~ grp, data = m[grp %in% c("never", "always")])
    cat("never vs always 的 Wilcoxon p =", format.pval(w$p.value, digits = 3), "\n")
  }
  cat("\n若 always 组建成度显著更高, 说明缺失是系统性的、与城市化正相关,\n")
  cat("直接用点值会把最城市化的站点选择性剔除 —— 这正是本研究关注的对象。\n")
}

fwrite(per_st, file.path(STAT_DIR, "mask_summary_by_station.csv"))
fwrite(per_t,  file.path(STAT_DIR, "mask_summary_by_date.csv"))
cat("\n写出: mask_summary_by_station.csv / mask_summary_by_date.csv\n")
