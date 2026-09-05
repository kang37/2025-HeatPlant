#!/usr/bin/env Rscript
# =============================================================================
# recovery_arx_index.R
#   De Keersmaecker et al. 2015 (GEB) 式 ARX 恢复力指数:
#   sif_anom(t) = phi*sif_anom(t-1) + beta*vpd_z(t) + e(t)   (年内滞后, 不跨冬歇期)
#   phi = resilience(自身异常的持续性), beta = resistance(VPD即时敏感性)
#   仅在 CCM 已确认因果耦合的站上估计(判据同 recovery_vs_ccm_direction.R:
#   ccm_hcsif_buf1000_surr 去趋势文件, p_surr<.1 & drho>0, 至少1个tp)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 站点清单: 方向标签 x CCM显著 --------------------------------------
L <- readRDS("data_proc/output_hcsif_buf1000/class_model_data.rds")$D[, .(stat_id, grp2)]
b <- fread("data_proc/ccm_hcsif_buf1000_surr/ccm_hcsif_buf1000_vpd_20260824_1601.csv")
b[, sig_loose := p_surr < 0.10 & (rho - rho_min) > 0]
sg    <- b[, .(any_loose = any(sig_loose)), by = .(stat_id = meteo_stat)]
coef0 <- b[tp == 0, .(stat_id = meteo_stat, coef_tp0 = mean_coef)]
stn   <- merge(merge(L, sg, by = "stat_id"), coef0, by = "stat_id")
stn[, direction := factor(ifelse(grp2 == 1, "Inhibit", "Promote"), levels = c("Promote","Inhibit"))]
sel <- stn[any_loose == TRUE, stat_id]
cat("CCM 确认因果耦合的站数:", length(sel), "\n")

# ---- 2. HCSIF + VPD 距平构造(同前几个脚本) --------------------------------
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

# ---- 3. 逐站 ARX 拟合 ------------------------------------------------------
fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom) & !is.na(vpd_z)]
  if (nrow(sub) < 30) return(NULL)
  m <- tryCatch(lm(sif_anom ~ lag1 + vpd_z, data = sub), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients
  if (!all(c("lag1","vpd_z") %in% rownames(co))) return(NULL)
  list(phi = co["lag1","Estimate"], phi_se = co["lag1","Std. Error"],
       beta = co["vpd_z","Estimate"], beta_se = co["vpd_z","Std. Error"],
       n = nrow(sub), r2 = summary(m)$r.squared)
}
arx <- d[, fit_one(.SD), by = meteo_stat]
arx <- merge(arx, stn, by.x = "meteo_stat", by.y = "stat_id")
fwrite(arx, file.path(OUT, "recovery_arx_index.csv"))
cat("成功拟合 ARX 的站数:", nrow(arx), "\n")

# ---- 4. 关键结论 ------------------------------------------------------------
cat("\n=== beta(ARX) vs mean_coef(tp=0, S-map) 一致性校验 ===\n")
cat("Pearson r =", round(cor(arx$beta, arx$coef_tp0), 3),
    " | Spearman rho =", round(cor(arx$beta, arx$coef_tp0, method="spearman"), 3), "\n")

cat("\n=== phi(恢复力) 按方向分组 ===\n")
print(arx[, .(n=.N, mean_phi=round(mean(phi),3), median_phi=round(median(phi),3),
              mean_absphi=round(mean(abs(phi)),3)), by=direction])
wt <- wilcox.test(phi ~ direction, data = arx)
wt_abs <- wilcox.test(abs(phi) ~ direction, data = arx)
cat("Wilcoxon phi ~ direction: W =", wt$statistic, " p =", signif(wt$p.value,3), "\n")
cat("Wilcoxon |phi| ~ direction: W =", wt_abs$statistic, " p =", signif(wt_abs$p.value,3), "\n")

# ---- 5. 图: 左=一致性校验散点, 右=phi按方向分组箱线图 ----------------------
p1 <- ggplot(arx, aes(coef_tp0, beta)) +
  geom_hline(yintercept=0,color="grey80",linewidth=.3) + geom_vline(xintercept=0,color="grey80",linewidth=.3) +
  geom_point(aes(color=direction), alpha=.6, size=1.6) +
  geom_smooth(method="lm", color="black", linewidth=.6, se=FALSE) +
  scale_color_manual(values=c(Promote="#2166AC", Inhibit="#C1121F")) +
  labs(title="(A) ARX beta vs S-map coef(tp=0): consistency check",
       subtitle=sprintf("Pearson r=%.2f, n=%d stations", cor(arx$beta,arx$coef_tp0), nrow(arx)),
       x="S-map mean_coef at tp=0", y="ARX beta (VPD sensitivity)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="top")

p2 <- ggplot(arx, aes(direction, phi, fill=direction)) +
  geom_boxplot(width=.5, outlier.alpha=.3) +
  scale_fill_manual(values=c(Promote="#2166AC", Inhibit="#C1121F")) +
  labs(title="(B) ARX resilience (phi) by CCM+S-map direction",
       subtitle=sprintf("Wilcoxon p=%.3g | higher phi = anomaly persists longer = lower resilience", wt$p.value),
       x=NULL, y="phi (AR(1) persistence of SIF anomaly)") +
  theme_bw(base_size=11) + theme(legend.position="none")

pg <- gridExtra::grid.arrange(p1, p2, ncol=2)
ggsave(file.path(OUT, "recovery_arx_index.png"), pg, width=11, height=5, dpi=300)
cat("\n-> recovery_arx_index.png\n")
