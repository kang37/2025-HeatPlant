#!/usr/bin/env Rscript
# =============================================================================
# recovery_ar1_pure.R
#   Forzieri et al. 2022 (Nature) 式"纯AR1"恢复力指标 vs De Keersmaecker 2015
#   式 ARX phi(控制了VPD) 的对比。都在年内滞后(不跨冬歇期)、去季节z-score距平
#   上估计，站点范围=统一CCM因果确认判据(pipelines/hcsif_buf1000/
#   12_ccm_causal_confirmation.R 的产出) —— 不涉及promote/inhibit方向标签，
#   不受方向不稳定问题影响。
#     模型1(纯AR1, Forzieri式):   sif_anom(t) = phi1*sif_anom(t-1) + e
#     模型2(ARX, De Keersmaecker式): sif_anom(t) = phi2*sif_anom(t-1) + beta*vpd_z(t) + e
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. CCM 确认因果耦合的站点清单(统一判据, 不涉及方向标签) ----------------
sel <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
cat("CCM 确认因果耦合的站数:", length(sel), "\n")

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
d[, lag1     := shift(sif_anom, 1), by = .(meteo_stat, year)]   # 年内滞后, 不跨越冬歇期

# ---- 3. 逐站拟合两个模型 ----------------------------------------------------
fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom) & !is.na(vpd_z)]
  if (nrow(sub) < 30) return(NULL)
  m1 <- tryCatch(lm(sif_anom ~ lag1, data = sub), error = function(e) NULL)
  m2 <- tryCatch(lm(sif_anom ~ lag1 + vpd_z, data = sub), error = function(e) NULL)
  if (is.null(m1) || is.null(m2)) return(NULL)
  co1 <- summary(m1)$coefficients; co2 <- summary(m2)$coefficients
  if (!("lag1" %in% rownames(co1)) || !("lag1" %in% rownames(co2))) return(NULL)
  list(phi1_pure = co1["lag1","Estimate"], phi1_se = co1["lag1","Std. Error"],
       phi2_arx  = co2["lag1","Estimate"], phi2_se = co2["lag1","Std. Error"],
       beta_arx  = co2["vpd_z","Estimate"], n = nrow(sub))
}
res <- d[, fit_one(.SD), by = meteo_stat]
fwrite(res, file.path(OUT, "recovery_ar1_pure_vs_arx.csv"))
cat("成功拟合的站数:", nrow(res), "\n")

# ---- 4. 关键结论 ------------------------------------------------------------
cat("\n=== phi1(纯AR1) vs phi2(ARX控制VPD) ===\n")
cat("Pearson r =", round(cor(res$phi1_pure, res$phi2_arx), 3), "\n")
cat("phi1 均值/中位数:", round(mean(res$phi1_pure),3), "/", round(median(res$phi1_pure),3), "\n")
cat("phi2 均值/中位数:", round(mean(res$phi2_arx),3),  "/", round(median(res$phi2_arx),3), "\n")
d_diff <- res$phi1_pure - res$phi2_arx
cat("phi1-phi2 差值 均值:", round(mean(d_diff),4), " | |差值|>0.05 的站数占比:",
    round(100*mean(abs(d_diff) > 0.05), 1), "%\n")
wt <- wilcox.test(res$phi1_pure, res$phi2_arx, paired = TRUE)
cat("配对 Wilcoxon (phi1 vs phi2): p =", signif(wt$p.value, 3), "\n")

# ---- 5. 图: phi1 vs phi2 一致性散点 -----------------------------------------
rng <- range(c(res$phi1_pure, res$phi2_arx))
p <- ggplot(res, aes(phi2_arx, phi1_pure)) +
  geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") +
  geom_point(alpha=.5, size=1.6, color="#2166AC") +
  geom_smooth(method="lm", color="black", linewidth=.6, se=FALSE) +
  coord_equal(xlim=rng, ylim=rng) +
  labs(title="Pure AR1 (Forzieri-style) vs VPD-controlled ARX phi (De Keersmaecker-style)",
       subtitle=sprintf("n=%d CCM-confirmed stations | Pearson r=%.2f | dashed=1:1 line\nphi1(pure)=sif_anom(t)~sif_anom(t-1); phi2(ARX)=sif_anom(t)~sif_anom(t-1)+vpd_z(t)",
                         nrow(res), cor(res$phi1_pure,res$phi2_arx)),
       x="phi2: ARX resilience (VPD controlled)", y="phi1: pure AR1 resilience (Forzieri-style)") +
  theme_bw(base_size=11)
ggsave(file.path(OUT, "recovery_ar1_pure_vs_arx.png"), p, width=7.5, height=7, dpi=300)
cat("\n-> recovery_ar1_pure_vs_arx.png\n")
