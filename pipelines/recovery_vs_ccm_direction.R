#!/usr/bin/env Rscript
# =============================================================================
# recovery_vs_ccm_direction.R
#   独立交叉验证: 只在 CCM 已确认存在因果耦合的站上, 用 Schwalm 式事件恢复构建
#   (recovery_schwalm_style_demo.R 的方法, 扩到全量站点), 检验"促进组恢复更快/
#   删失更少, 抑制组更慢/更易删失"是否成立 -- 作为 44-60 系列分类结果(grp2:
#   促进/抑制, 定义见 44_class_data.R)的一个独立外部佐证。
#
#   方向标签: class_model_data.rds (grp2, strict=样本A严格四分类, 全部=样本B)
#   因果确认: ccm_hcsif_buf1000_surr(去趋势, 与分类标签同一套系数, tp0系数相关=1)
#             的 p_surr, 站级"至少1个tp显著" —— 宽松 p_surr<.1&drho>0 / 严格 p_surr<.05&drho>.05
#   恢复构建: 每站每年生长季最大VPD窗口=事件, 事件前2窗均值=基线,
#             事件后首个 anom>=baseline 的窗口=恢复, 季末未恢复=右删失(time=可观测窗口数)
#   统计: Kaplan-Meier + log-rank + Cox(cluster=站点稳健SE), 促进=参照组
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(survival)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. 方向标签 x 因果确认 -> 站点清单
# ---------------------------------------------------------------------------
L <- readRDS("data_proc/output_hcsif_buf1000/class_model_data.rds")$D[, .(stat_id, grp2, strict)]

b <- fread("data_proc/ccm_hcsif_buf1000_surr/ccm_hcsif_buf1000_vpd_20260824_1601.csv")
b[, sig_loose  := p_surr < 0.10 & (rho - rho_min) > 0]
b[, sig_strict := p_surr < 0.05 & (rho - rho_min) > 0.05]
sg <- b[, .(any_loose = any(sig_loose), any_strict = any(sig_strict)), by = .(stat_id = meteo_stat)]

stn_all <- merge(L, sg, by = "stat_id")
stn_all[, direction := factor(ifelse(grp2 == 1, "Inhibit", "Promote"), levels = c("Promote","Inhibit"))]
cat("=== 样本量矩阵 (方向标签 x 因果确认) ===\n")
print(stn_all[, .N, by = .(sample = ifelse(strict,"A_strict4class","B_all"),
                           ccm_sig = ifelse(any_strict,"strict(p<.05)",
                                     ifelse(any_loose,"loose(p<.1)","not_sig")),
                           direction)][order(sample, ccm_sig, direction)])

sel_stations <- stn_all[any_loose == TRUE, stat_id]   # 最大需要集合, 之后按需再子集
cat("\n用于事件构建的站点数(样本B∩loose,覆盖后续所有子集需要):", length(sel_stations), "\n")

# ---------------------------------------------------------------------------
# 2. HCSIF + VPD 8天序列, 全量站点版(照搬 recovery_schwalm_style_demo.R 的构建)
# ---------------------------------------------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[meteo_stat %in% sel_stations, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel_stations]
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

# ---------------------------------------------------------------------------
# 3. 逐站逐年: 事件 -> 基线 -> 恢复/删失 (与 demo 脚本同一套规则)
# ---------------------------------------------------------------------------
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
# 生存时间: 恢复了就是恢复所需窗口数, 删失了就是"观测到的最后一个窗口数"(右删失时间)
events[, time  := ifelse(censored, n_post_avail, Trec_windows)]
events[, event := as.integer(!censored)]
events <- events[!is.na(time)]

events <- merge(events, stn_all, by.x = "meteo_stat", by.y = "stat_id")
fwrite(events, file.path(OUT, "recovery_vs_ccm_direction_events.csv"))
cat("\n事件总数(样本B∩loose):", nrow(events), " | 站数:", uniqueN(events$meteo_stat), "\n")

# ---------------------------------------------------------------------------
# 4. 生存分析: 4 个子集(样本A/B x loose/strict)都跑一遍, 汇总 HR/log-rank p
# ---------------------------------------------------------------------------
run_one <- function(ev, tag) {
  if (uniqueN(ev$direction) < 2 || nrow(ev) < 20) return(NULL)
  sd_ <- survdiff(Surv(time, event) ~ direction, data = ev)
  p_lr <- 1 - pchisq(sd_$chisq, length(sd_$n) - 1)
  cx <- coxph(Surv(time, event) ~ direction + cluster(meteo_stat), data = ev)
  co <- summary(cx)$coefficients
  data.table(subset = tag, n_events = nrow(ev), n_stations = uniqueN(ev$meteo_stat),
             n_inhibit_stn = uniqueN(ev[direction=="Inhibit"]$meteo_stat),
             n_promote_stn = uniqueN(ev[direction=="Promote"]$meteo_stat),
             pct_censored_promote = round(100*mean(ev[direction=="Promote"]$censored),1),
             pct_censored_inhibit = round(100*mean(ev[direction=="Inhibit"]$censored),1),
             HR_inhibit_vs_promote = round(exp(co[1,"coef"]), 3),
             HR_p_robust = signif(co[1,"Pr(>|z|)"], 3),
             logrank_p = signif(p_lr, 3))
}

subsets <- list(
  A_loose  = events[strict==TRUE  & any_loose==TRUE],
  B_loose  = events[                any_loose==TRUE],
  A_strict = events[strict==TRUE  & any_strict==TRUE],
  B_strict = events[                any_strict==TRUE]
)
res <- rbindlist(lapply(names(subsets), function(nm) run_one(subsets[[nm]], nm)), fill = TRUE)
cat("\n=== 生存分析汇总(HR>1 = 抑制组恢复更快[hazard更高]; HR<1 = 抑制组恢复更慢) ===\n")
print(res)
fwrite(res, file.path(OUT, "recovery_vs_ccm_direction_summary.csv"))

# ---------------------------------------------------------------------------
# 5. 图: 主分析(样本A∩loose) Kaplan-Meier 曲线
# ---------------------------------------------------------------------------
main_ev <- subsets$A_loose
sf <- survfit(Surv(time, event) ~ direction, data = main_ev)
sfd <- as.data.table(summary(sf)[c("time","surv","strata","n.risk")])
sfd[, direction := sub("direction=", "", as.character(strata))]
main_r <- res[subset == "A_loose"]

p <- ggplot(sfd, aes(time, 1 - surv, color = direction)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c(Promote = "#2166AC", Inhibit = "#C1121F")) +
  scale_x_continuous(breaks = 0:8) + scale_y_continuous(labels = scales::percent) +
  labs(title = "Cumulative probability of SIF recovery after the yearly peak-VPD event",
       subtitle = sprintf("Sample A (strict 4-class stations) x CCM-confirmed causal stations (loose p_surr<.1) | n=%d station-years, %d stations\nlog-rank p=%s | Cox HR(Inhibit vs Promote, station-clustered robust SE)=%.2f, p=%s | curve stays below 100%% = right-censored (season ended first)",
                          main_r$n_events, main_r$n_stations, formatC(main_r$logrank_p,format="e",digits=2),
                          main_r$HR_inhibit_vs_promote, formatC(main_r$HR_p_robust,format="e",digits=2)),
       x = "t (8-day windows since peak-VPD event)", y = "Cumulative % recovered", color = "CCM+S-map direction") +
  theme_bw(base_size = 12) + theme(legend.position = "top", plot.subtitle = element_text(size=8.5,color="grey35"))
ggsave(file.path(OUT, "recovery_vs_ccm_direction_KM.png"), p, width = 9, height = 6.5, dpi = 300)
cat("\n-> recovery_vs_ccm_direction_KM.png\n")
