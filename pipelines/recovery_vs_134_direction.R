#!/usr/bin/env Rscript
# =============================================================================
# recovery_vs_134_direction.R
#   Schwalm 式事件-恢复构建(同 recovery_schwalm_style_demo.R 的方法), 限定在
#   CCM统一因果确认判据(12_ccm_causal_confirmation.R)筛出的站上
#   (data_proc/smap_bivar_134/classify_134.csv, 脚本名沿用历史134命名, 实际
#   站数以该文件为准), 用该支线自己的双变量S-map方向标签(direction_2v:
#   Promote/Inhibit/Ambiguous)做分组生存分析。
#   替代已删除的 recovery_vs_ccm_direction.R(那个用的是44-60的grp2, 单变量
#   S-map、已知有偏, 见 docs/03_classification.md 2026-09-05 决定)。
#   已知局限(必须在结果解读时带上): direction_2v 本身已被
#   07_smap_robustness_theta_multivar.R 证明对状态空间维度不稳健(Inhibit组
#   加一维温度后近半反转), 见 docs/00_cross_module_issues.md 第1条。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(survival)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 站方向标签(站数以 classify_134.csv 实际行数为准) -------------------
stn <- fread("data_proc/smap_bivar_134/classify_134.csv")[, .(meteo_stat, direction_2v, optimal_tp)]
stn[, direction := factor(direction_2v, levels = c("Promote","Ambiguous","Inhibit"))]
cat("=== 站点数(按 direction_2v) ===\n"); print(stn[, .N, by = direction])
sel <- stn$meteo_stat

# ---- 2. HCSIF + VPD 距平构造(同前) -----------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)

deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE)
  a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]

# ---- 3. 逐站逐年: 事件(当季VPD峰值窗口) -> 基线(前2窗) -> 恢复/删失 ---------
events <- d[, {
  yy <- .SD[order(doy)]
  ev <- yy[which.max(vpd_z)]
  i0 <- which(yy$doy == ev$doy); n <- nrow(yy)
  base_idx <- if (i0 > 1) max(1, i0-2):(i0-1) else integer(0)
  baseline <- if (length(base_idx)) mean(yy$sif_anom[base_idx], na.rm = TRUE) else NA_real_
  post <- if (i0 < n) (i0+1):n else integer(0)
  rec_i <- if (is.na(baseline) || !length(post)) integer(0) else
    post[which(yy$sif_anom[post] >= baseline)[1]]
  list(baseline = baseline, n_post_avail = length(post),
       Trec_windows = if (length(rec_i)) rec_i - i0 else NA_integer_,
       censored = length(rec_i) == 0)
}, by = .(meteo_stat, year)]
events <- events[!is.na(baseline)]
events[, time  := ifelse(censored, n_post_avail, Trec_windows)]
events[, event := as.integer(!censored)]
events <- events[!is.na(time)]
events <- merge(events, stn, by = "meteo_stat")
fwrite(events, file.path(OUT, "recovery_vs_134_direction_events.csv"))
cat("\n事件总数:", nrow(events), " | 站数:", uniqueN(events$meteo_stat), "\n")
print(events[, .(n_events=.N, pct_censored=round(100*mean(censored),1)), by=direction])

# ---- 4. 生存分析: Promote 作参照组 ------------------------------------------
sd_  <- survdiff(Surv(time, event) ~ direction, data = events)
p_lr <- 1 - pchisq(sd_$chisq, length(sd_$n) - 1)
cx   <- coxph(Surv(time, event) ~ direction + cluster(meteo_stat), data = events)
cat("\n=== log-rank p =", signif(p_lr,3), "===\n")
print(summary(cx)$coefficients)

# ---- 5. 图: Kaplan-Meier, 三组 ----------------------------------------------
sf  <- survfit(Surv(time, event) ~ direction, data = events)
sfd <- as.data.table(summary(sf)[c("time","surv","strata")])
sfd[, direction := sub("direction=", "", as.character(strata))]

p <- ggplot(sfd, aes(time, 1-surv, color = direction)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c(Promote="#2166AC", Ambiguous="grey55", Inhibit="#C1121F")) +
  scale_x_continuous(breaks = 0:8) + scale_y_continuous(labels = scales::percent) +
  labs(title = "SIF recovery after peak-VPD event, by 134-station bivariate S-map direction",
       subtitle = sprintf("n=%d station-years, %d stations | log-rank p=%.3g | Cox ref=Promote, cluster-robust SE\nKnown caveat: direction_2v itself is unstable to state-space dimension (see docs/00_cross_module_issues.md #1)",
                           nrow(events), uniqueN(events$meteo_stat), p_lr),
       x = "t (8-day windows since peak-VPD event)", y = "Cumulative % recovered", color = NULL) +
  theme_bw(base_size=12) + theme(legend.position="top", plot.subtitle=element_text(size=8,color="grey35"))
ggsave(file.path(OUT, "recovery_vs_134_direction_KM.png"), p, width=9, height=6.5, dpi=300)
cat("\n-> recovery_vs_134_direction_KM.png\n")
