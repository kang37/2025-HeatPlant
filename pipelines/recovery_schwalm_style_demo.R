#!/usr/bin/env Rscript
# =============================================================================
# recovery_schwalm_style_demo.R
#   可行性检验: 能否在本项目 HCSIF 8天合成数据上, 照搬 Schwalm et al. 2017 (Nature)
#   的"drought recovery time"构建方式——即直接在观测的 SIF 异常序列上定位一次
#   真实发生的 VPD 极端事件, 找异常最低点, 数到异常回正/回到扰动前基线所需的时长。
#   这和 recovery_indicator_demo.R/recovery_within_N.R 是两个不同的"recovery"概念:
#   那两个脚本是在已拟合好的 CCM/S-map 滞后响应曲线 beta(tp) 上算衰减, 依赖 CCM 结果;
#   这个脚本直接在原始 SIF/VPD 观测序列上找事件、算恢复, 不依赖 CCM/S-map 拟合结果,
#   是"能不能绕开 CCM 先做 resilience 分析"这个问题里 Schwalm 一路的真实还原。
#   10 个站点, 全英文图注(比照 recovery_indicator_demo.R 的约定)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. 读 HCSIF 站点序列 + VPD 8天缓存, 和 01_ccm_buf1000.R 完全一致的合并方式
# ---------------------------------------------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
d[, idx8 := as.integer(round(as.numeric(
      as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]

# ---------------------------------------------------------------------------
# 2. 选 10 站: 4 个沿用 recovery_indicator_demo.R 里已经定过型的站(强信号典型样例)
#    + 6 个按纬度均匀撒开、覆盖不同气候带, 检验方法对"普通站"是否也适用
# ---------------------------------------------------------------------------
curated <- c(`A_51814`=51814, `B_51567`=51567, `C_53487`=53487, `D_54943`=54943)
stn <- as.data.table(fread("data_raw/hcsif/stations_924.csv"))
stn <- stn[meteo_stat %in% unique(d$meteo_stat)][order(latitude)]
pick_lat <- stn[round(seq(1, .N, length.out = 8))]$meteo_stat
extra <- setdiff(pick_lat, curated)[1:6]
names(extra) <- paste0("E_", extra)
SEL <- c(curated, extra)
cat("=== 选中的 10 站 ===\n"); print(data.table(label=names(SEL), meteo_stat=SEL,
      stn[match(SEL, meteo_stat), .(latitude, longitude)]))

sub <- d[meteo_stat %in% SEL]
sub[, label := factor(names(SEL)[match(meteo_stat, SEL)], levels = names(SEL))]

# ---------------------------------------------------------------------------
# 3. 逐站: 去季节距平(z-score), 定位每个生长季里最极端的 VPD 窗口作为"事件",
#    往前2个窗口算基线, 往后数到异常回正(或回到>=基线)所需窗口数=recovery time,
#    若生长季在 9 月结束前未回正 -> censored(右删失, 这是本数据结构特有的问题)
# ---------------------------------------------------------------------------
deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE)
  a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
sub[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
sub[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]

events <- sub[, {
  yy <- .SD[order(doy)]
  ev <- yy[which.max(vpd_z)]
  i0 <- which(yy$doy == ev$doy)
  n  <- nrow(yy)
  base_idx <- max(1, i0 - 2):max(1, i0 - 1)
  baseline <- if (i0 > 1) mean(yy$sif_anom[base_idx], na.rm = TRUE) else NA_real_
  post <- if (i0 < n) (i0+1):n else integer(0)
  # Schwalm 原文口径: 异常回到扰动前基线水平(不是回到 0), 缺基线时不可判定
  rec_i <- if (is.na(baseline)) integer(0) else
    post[which(yy$sif_anom[post] >= baseline)[1]]
  list(event_doy = ev$doy, event_vpdz = ev$vpd_z, event_anom = ev$sif_anom,
       baseline = baseline, n_windows_season = n, n_post_avail = length(post),
       Trec_windows = if (length(rec_i)) rec_i - i0 else NA_integer_,
       censored = length(rec_i) == 0)
}, by = .(meteo_stat, label, year)]
events <- events[!is.na(event_doy)]

feas <- events[, .(n_years = .N,
                    avg_windows_per_season = round(mean(n_windows_season),1),
                    pct_censored = round(100*mean(censored),1),
                    median_Trec_8d = as.numeric(median(Trec_windows, na.rm=TRUE))),
                by = .(meteo_stat, label)]
cat("\n=== 逐站可行性汇总 ===\n"); print(feas)
fwrite(events, file.path(OUT, "recovery_schwalm_style_events.csv"))
fwrite(feas,   file.path(OUT, "recovery_schwalm_style_feasibility.csv"))

# ---------------------------------------------------------------------------
# 4. 图A: 10 站全时间序列总览(idx8 为真实时间轴, 冬歇期空档如实呈现),
#    每年事件窗口用竖虚线标出, 直观看"数据本身能不能支撑这种构建方式"
# ---------------------------------------------------------------------------
ev_mark <- merge(events[, .(meteo_stat, label, year, event_doy, censored)],
                  sub[, .(meteo_stat, year, doy, idx8)],
                  by.x=c("meteo_stat","year","event_doy"), by.y=c("meteo_stat","year","doy"))

p1 <- ggplot(sub, aes(idx8, sif_anom)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = .3) +
  geom_line(color = "grey55", linewidth = .3) +
  geom_point(size = .6, color = "grey40") +
  geom_vline(data = ev_mark, aes(xintercept = idx8, color = censored),
             linetype = "dashed", linewidth = .4, alpha=.8) +
  scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="#2166AC"),
                      labels = c(`TRUE`="event censored (season ends before recovery)",
                                 `FALSE`="event recovered within season"), name = NULL) +
  facet_wrap(~label, ncol = 2, scales = "free_x") +
  labs(title = "Full observed SIF anomaly series with yearly peak-VPD events marked (10 stations)",
       subtitle = "x-axis is real elapsed time (8-day steps since 2000-01-01) -- gaps ARE the Oct-Apr dormant season with no HCSIF data.\nEach dashed line = the most extreme VPD window of that growing season (candidate 'drought' event, Schwalm-style).",
       x = "idx8 (8-day step index since 2000-01-01)", y = "SIF anomaly (deseasoned z-score)") +
  theme_bw(base_size = 11) + theme(legend.position = "top",
    strip.text = element_text(face="bold"), plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_schwalm_style_overview.png"), p1, width = 13, height = 11, dpi = 300)

# ---------------------------------------------------------------------------
# 5. 图B: 事件复合(event-composite)轨迹 -- 把每年的事件对齐到 t=0,
#    画 t=-2..+6 的 SIF 异常轨迹, 叠加所有年份, 右删失(季末截断)的年份
#    用虚线+末端空心点标出、不假装它"没恢复", 只是"看不到后面"
# ---------------------------------------------------------------------------
traj <- sub[, {
  yy <- .SD[order(doy)]
  ev <- events[meteo_stat==.BY$meteo_stat & year==.BY$year]
  if (!nrow(ev) || is.na(ev$event_doy)) return(NULL)
  i0 <- which(yy$doy == ev$event_doy); n <- nrow(yy)
  rel <- (max(1,i0-2)):n
  data.table(t_rel = rel - i0, sif_anom = yy$sif_anom[rel],
             is_last = rel == n, censored = ev$censored)
}, by = .(meteo_stat, label, year)]
traj <- traj[t_rel <= 6]
last_pt <- traj[is_last == TRUE]

p2 <- ggplot(traj, aes(t_rel, sif_anom, group = year)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
  geom_line(alpha = .35, color = "grey40", linewidth = .4) +
  geom_point(data = last_pt[censored==TRUE], shape = 1, size = 2, color = "#C1121F") +
  stat_summary(aes(group=1), fun = mean, geom = "line", color = "#2166AC", linewidth = 1) +
  facet_wrap(~label, ncol = 2) +
  scale_x_continuous(breaks = -2:6) +
  labs(title = "Event-composite recovery trajectories: SIF anomaly aligned to each year's peak-VPD window (t=0)",
       subtitle = "Thin grey lines = individual years (2000-2022); thick blue = across-year mean; open red circles = last available\npoint before the growing season ends with recovery still undetermined (right-censored, not 'no recovery').",
       x = "t (8-day steps relative to peak-VPD window)", y = "SIF anomaly (deseasoned z-score)") +
  theme_bw(base_size = 11) + theme(strip.text = element_text(face="bold"),
    plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_schwalm_style_composite.png"), p2, width = 13, height = 11, dpi = 300)

cat("\n-> recovery_schwalm_style_overview.png\n-> recovery_schwalm_style_composite.png\n")
