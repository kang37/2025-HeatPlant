#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_hcsif_qc.R — HCSIF 站点提取结果的质检
#
# 三条主线:
#   1. 缺失值:  点值 / 邻域均值的有效率, 缺失是否随机
#   2. 与原有 0.05° 产品的一致性: 相关、偏差、以及偏差的空间/城市化结构
#   3. 数值合理性: 量级、季节曲线、异常值
#
# 关键预期(而非"必须相符"):
#   HCSIF 与 0.05° 产品的差异应当在**城市化程度高的站点最大**。
#   若如此, 说明 500 m 确实分辨出了粗分辨率被稀释掉的城市信号 —— 这正是换数据的理由。
#   若差异与城市化无关且整体相关很低, 则更像是产品不一致或提取错误, 需要排查。
#
# 用法: Rscript 13_hcsif_qc.R [out_dir]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station")
MATCH_TOL <- 4L      # 与旧产品匹配的最大日期差(天); 两套 8 天网格固定相差 3 天

args    <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args)) args[1] else file.path(PROJ, "data_proc/ccm_hcsif_500m/qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

rule <- function(t) cat("\n", strrep("=", 66), "\n", t, "\n", strrep("=", 66), "\n", sep = "")

fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有提取结果")
new <- rbindlist(lapply(fs, fread), fill = TRUE)
new[, masked := as.logical(pt_masked)]
new[, date_d := as.Date(as.character(date), format = "%Y%m%d")]

rule("1. 数据规模")
cat("年份      :", paste(sort(unique(new$year)), collapse = ", "), "\n")
cat("时相数    :", uniqueN(new$date), "\n")
cat("站点数    :", uniqueN(new$meteo_stat), "\n")
cat("记录数    :", nrow(new), "\n")
cat("每站时相数:", paste(range(new[, .N, by = meteo_stat]$N), collapse = " ~ "), "\n")

# --------------------------------------------------------------------------
rule("2. 缺失值")
n <- nrow(new)
cat(sprintf("点值有效      : %6d (%.1f%%)\n", sum(!is.na(new$SIF)), 100*mean(!is.na(new$SIF))))
cat(sprintf("点值缺失      : %6d (%.1f%%)\n", sum(is.na(new$SIF)),  100*mean(is.na(new$SIF))))
cat(sprintf("  其中0值掩膜 : %6d\n", sum(new$masked, na.rm = TRUE)))
cat(sprintf("邻域均值有效  : %6d (%.1f%%)\n", sum(!is.na(new$SIF_buf750)), 100*mean(!is.na(new$SIF_buf750))))
cat(sprintf("两者皆缺      : %6d (%.1f%%)\n",
            new[is.na(SIF) & is.na(SIF_buf750), .N],
            100*new[is.na(SIF) & is.na(SIF_buf750), .N]/n))

per_st <- new[, .(n_obs = .N, n_miss = sum(is.na(SIF)), n_masked = sum(masked)), by = meteo_stat]
cat("\n按站点:\n")
cat("  全程有值  :", per_st[n_miss == 0, .N], "站\n")
cat("  部分缺失  :", per_st[n_miss > 0 & n_miss < n_obs, .N], "站\n")
cat("  全程缺失  :", per_st[n_miss == n_obs, .N], "站\n")

per_t <- new[, .(n_miss = sum(is.na(SIF)), pct = round(100*mean(is.na(SIF)), 2)), by = .(date, doy)]
setorder(per_t, date)
cat(sprintf("\n各时相缺失率: %.1f%% ~ %.1f%% (极差 %.1f 个百分点)\n",
            min(per_t$pct), max(per_t$pct), max(per_t$pct) - min(per_t$pct)))
cat("缺失率随 DOY 的变化(检验是否有季节性 -> 提示云/物候而非静态掩膜):\n")
print(per_t[, .(mean_pct = round(mean(pct), 2)), by = doy][order(doy)])

# --------------------------------------------------------------------------
rule("3. 数值合理性")
v <- new[!is.na(SIF)]$SIF
cat("点值分位数 (mW m-2 nm-1 sr-1):\n"); print(round(quantile(v, c(0,.01,.25,.5,.75,.99,1)), 4))
cat("\n负值比例:", sprintf("%.3f%%", 100*mean(v < 0)), "\n")
cat("超过 1.0 的比例:", sprintf("%.3f%%", 100*mean(v > 1)), "\n")
cat("\n季节曲线(按 DOY 的中位数):\n")
print(new[!is.na(SIF), .(median_sif = round(median(SIF), 4), n = .N), by = doy][order(doy)])

# --------------------------------------------------------------------------
rule("4. 与原有 0.05° 产品对比")
old <- fread(file.path(PROJ, "data_raw/meteo_stat_SIF_data.csv"))
old <- old[!is.na(SIF), .(meteo_stat, date_old = date, SIF_old = SIF)]
old[, date_d := as.Date(as.character(date_old), format = "%Y%m%d")]
old <- old[date_d >= min(new$date_d) - MATCH_TOL & date_d <= max(new$date_d) + MATCH_TOL]

# 按站点+最近日期匹配(rolling join)
setkey(old, meteo_stat, date_d)
nn <- new[!is.na(SIF), .(meteo_stat, date_d, doy, SIF, SIF_buf750)]
setkey(nn, meteo_stat, date_d)
mm <- old[nn, roll = "nearest"]
mm[, dgap := abs(as.integer(date_d - as.Date(as.character(date_old), format = "%Y%m%d")))]
mm <- mm[!is.na(SIF_old) & dgap <= MATCH_TOL]

cat("可比记录:", nrow(mm), " (日期差中位数", median(mm$dgap), "天)\n")
cat(sprintf("\nHCSIF   均值 %.4f  中位 %.4f  sd %.4f\n",
            mean(mm$SIF), median(mm$SIF), sd(mm$SIF)))
cat(sprintf("旧产品  均值 %.4f  中位 %.4f  sd %.4f\n",
            mean(mm$SIF_old), median(mm$SIF_old), sd(mm$SIF_old)))
cat(sprintf("比值(HCSIF/旧) 中位数 %.3f\n", median(mm$SIF / mm$SIF_old, na.rm = TRUE)))

cat(sprintf("\n全样本 Pearson r = %.3f | Spearman = %.3f\n",
            cor(mm$SIF, mm$SIF_old), cor(mm$SIF, mm$SIF_old, method = "spearman")))
cat(sprintf("邻域均值 vs 旧产品 r = %.3f\n",
            cor(mm$SIF_buf750, mm$SIF_old, use = "complete.obs")))

# 站点级相关(每站至少 10 对)
bs <- mm[, .(n = .N, r = if (.N >= 10) cor(SIF, SIF_old) else NA_real_,
             bias = mean(SIF - SIF_old)), by = meteo_stat][!is.na(r)]
cat("\n站点级相关系数分布 (n =", nrow(bs), "站):\n"); print(round(quantile(bs$r, c(0,.1,.25,.5,.75,.9,1)), 3))
cat("站点级偏差(HCSIF - 旧)分布:\n"); print(round(quantile(bs$bias, c(0,.25,.5,.75,1)), 4))

# --------------------------------------------------------------------------
rule("5. 差异是否集中在城市化站点(关键判据)")
bf_path <- file.path(PROJ, "data_raw/China_stations_buildings_with_coords30_local_rerun_filled.csv")
if (file.exists(bf_path)) {
  b  <- fread(bf_path)
  b2 <- unique(b[, .(meteo_stat = as.integer(AirQualityStation),
                     bf = as.numeric(buildingFootprint))], by = "meteo_stat")[!is.na(meteo_stat)]
  mb <- merge(bs, b2, by = "meteo_stat")[!is.na(bf)]
  # 用秩次分四组: 大量站点建筑占地为 0，直接按分位数切会因分位点重复而报错
  mb[, bf_grp := cut(frank(bf, ties.method = "first"),
                     breaks = 4, labels = c("Q1最低", "Q2", "Q3", "Q4最高"))]
  cat("按建成度四分位:\n")
  print(mb[, .(n = .N, median_r = round(median(r), 3),
               median_bias = round(median(bias), 4)), by = bf_grp][order(bf_grp)])
  ct <- suppressWarnings(cor.test(mb$bf, mb$r, method = "spearman"))
  cat(sprintf("\n建成度 vs 站点级相关: Spearman rho = %.3f, p = %s\n",
              ct$estimate, format.pval(ct$p.value, digits = 3)))
  cat("负相关说明: 越城市化的站点, 两套产品越不一致 —— 符合 500 m 分辨出城市信号的预期。\n")
} else cat("缺建筑数据, 跳过\n")

fwrite(bs, file.path(OUT_DIR, "qc_station_agreement.csv"))
fwrite(per_t, file.path(OUT_DIR, "qc_missing_by_date.csv"))
cat("\n写出:", file.path(OUT_DIR, "qc_station_agreement.csv"), "\n")
