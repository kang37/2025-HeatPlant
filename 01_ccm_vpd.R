# =============================================================================
# 01_ccm_vpd.R  (更新版)
# 修改点：
#   1. CCM 之前过滤 SIF 结构突变站点，并可视化正常 vs 异常站点的 SIF 序列
#   2. tp 延伸到 0..8
#   3. X 变量：vpd_mean（周均VPD）& heat_over_sum（周累积VPD超量）
# =============================================================================

pacman::p_load(dplyr, purrr, tidyr, ggplot2, targets, rEDM, tibble,
               stringr, strucchange, showtext, sysfonts)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUTPUT  <- "ccm_vpd_results.rds"   # 保存至根目录
OUT_FIG <- "data_proc"

# 字体
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 16
theme_cn <- function(bs = BS) {
  theme_minimal(base_size = bs) +
    theme(text = element_text(family = "heiti"),
          plot.title    = element_text(face = "bold", hjust = 0.5, size = bs + 4),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = bs),
          strip.text    = element_text(face = "bold", size = bs),
          axis.text     = element_text(size = bs - 2),
          axis.title    = element_text(size = bs))
}
theme_set(theme_cn())

# =============================================================================
# 1. 加载数据
# =============================================================================
tar_load(data_heat_sif_weekly)

safe_detrend <- function(x, t) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch(
    residuals(lm(x ~ t, na.action = na.exclude)),
    error = function(e) rep(NA_real_, length(x))
  )
}

cat("去趋势 vpd_mean 和 heat_over_sum...\n")
df_weekly <- data_heat_sif_weekly %>%
  group_by(meteo_stat_id) %>%
  arrange(year, week) %>%
  mutate(
    time_idx     = row_number(),
    vpd_mean_dt  = safe_detrend(vpd_mean,     time_idx),
    heat_over_dt = safe_detrend(heat_over_sum, time_idx)
  ) %>%
  ungroup()

# =============================================================================
# 2. SIF 结构突变检测（在 CCM 之前）
# =============================================================================
cat("\n检测 SIF 结构突变（在 CCM 前过滤）...\n")

detect_sif_break <- function(sid, df) {
  d <- df %>%
    filter(meteo_stat_id == sid, week %in% 20:39) %>%
    arrange(year, week) %>%
    filter(!is.na(sif_interp))
  if (nrow(d) < 30)
    return(tibble(meteo_stat_id = sid, has_break = NA_integer_,
                  break_pval = NA_real_, break_magnitude = NA_real_,
                  break_year = NA_real_))
  tryCatch({
    fs  <- strucchange::Fstats(sif_interp ~ 1, data = d, from = 0.15, to = 0.85)
    pv  <- strucchange::sctest(fs)$p.value
    bp  <- strucchange::breakpoints(sif_interp ~ 1, data = d, breaks = 1)
    idx <- bp$breakpoints
    if (!is.na(idx) && idx > 1 && idx < nrow(d)) {
      mu1 <- mean(d$sif_interp[1:idx],            na.rm = TRUE)
      mu2 <- mean(d$sif_interp[(idx+1):nrow(d)],  na.rm = TRUE)
      mag <- abs(mu2 - mu1) / sd(d$sif_interp, na.rm = TRUE)
      yr  <- d$year[idx]
    } else { mag <- 0; yr <- NA_real_ }
    tibble(meteo_stat_id  = sid,
           has_break      = as.integer(pv < 0.01 & mag > 1.0),
           break_pval     = pv,
           break_magnitude = mag,
           break_year     = yr)
  }, error = function(e)
    tibble(meteo_stat_id = sid, has_break = 0L,
           break_pval = NA_real_, break_magnitude = NA_real_, break_year = NA_real_))
}

all_sids <- df_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(sif_interp)) >= 30) %>%
  summarise(.groups = "drop") %>%
  pull(meteo_stat_id)

break_res <- map_dfr(all_sids, ~detect_sif_break(.x, df_weekly))
saveRDS(break_res, "data_proc/sif_break_results.rds")

n_break   <- sum(break_res$has_break == 1L, na.rm = TRUE)
n_clean   <- sum(break_res$has_break == 0L, na.rm = TRUE)
cat(sprintf("  突变站点: %d (%.1f%%)  |  正常站点: %d (%.1f%%)\n",
            n_break, n_break/length(all_sids)*100,
            n_clean, n_clean/length(all_sids)*100))

# =============================================================================
# 3. 可视化：正常站点 vs 突变站点的 SIF 时间序列
# =============================================================================
set.seed(7)
n_break_avail <- sum(break_res$has_break == 1L, na.rm = TRUE)
n_clean_avail <- sum(break_res$has_break == 0L, na.rm = TRUE)
sids_break <- break_res %>% filter(has_break == 1L) %>%
  slice_sample(n = min(4L, n_break_avail)) %>% pull(meteo_stat_id)
sids_clean <- break_res %>% filter(has_break == 0L) %>%
  slice_sample(n = min(4L, n_clean_avail)) %>% pull(meteo_stat_id)

plot_sids  <- c(sids_break, sids_clean)
plot_label <- tibble(
  meteo_stat_id = plot_sids,
  group_label   = c(rep("SIF突变站点（已排除）", length(sids_break)),
                    rep("正常站点（纳入CCM）",   length(sids_clean)))
)

sif_ts <- df_weekly %>%
  filter(meteo_stat_id %in% plot_sids, week %in% 20:39) %>%
  arrange(meteo_stat_id, year, week) %>%
  group_by(meteo_stat_id) %>%
  mutate(t_seq = row_number()) %>%
  ungroup() %>%
  left_join(plot_label, by = "meteo_stat_id") %>%
  left_join(break_res %>% dplyr::select(meteo_stat_id, break_year),
            by = "meteo_stat_id")

p_sif <- ggplot(sif_ts, aes(x = year + (week - 20) / 20, y = sif_interp)) +
  geom_line(color = "steelblue", linewidth = 0.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  geom_vline(aes(xintercept = break_year), color = "#D73027",
             linetype = "solid", linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ meteo_stat_id + group_label, scales = "free_y", ncol = 4) +
  labs(
    title    = "SIF 时间序列：突变站点 vs 正常站点",
    subtitle = "红色竖线 = 检测到的突变年份（p<0.01 且幅度>1SD）",
    x = "年份", y = "SIF（插值后）"
  ) +
  theme_cn() +
  theme(strip.text = element_text(size = BS - 4))

ggsave(file.path(OUT_FIG, "sif_break_vs_clean.png"),
       p_sif, width = 16, height = 8, dpi = 300)
cat("-> data_proc/sif_break_vs_clean.png\n")

# =============================================================================
# 4. 过滤突变站点，获得用于 CCM 的干净站点
# =============================================================================
stations_clean <- break_res %>%
  filter(has_break == 0L) %>%
  pull(meteo_stat_id)

stations_ok <- df_weekly %>%
  filter(week %in% 20:39, meteo_stat_id %in% stations_clean) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(vpd_mean_dt)) >= 30,
         sum(!is.na(sif_detrended)) >= 30) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(meteo_stat_id)

cat("突变过滤后候选站点:", length(stations_ok), "\n")

# =============================================================================
# 5. 分层抽样（测试版）
# =============================================================================
meta_all <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  distinct(meteo_stat_id, koppen_group, longitude, latitude, koppen_class) %>%
  filter(meteo_stat_id %in% stations_ok)

N_SAMPLE <- 25    # 25≈30min | 100≈2h | 全量≈14h

set.seed(42)
stations_sample <- meta_all %>%
  group_by(koppen_group) %>%
  slice_sample(prop = N_SAMPLE / nrow(meta_all), replace = FALSE) %>%
  ungroup() %>%
  { if (nrow(.) < N_SAMPLE)
      bind_rows(., anti_join(meta_all, ., by = "meteo_stat_id") %>%
                  slice_sample(n = N_SAMPLE - nrow(.)))
    else slice_sample(., n = N_SAMPLE) } %>%
  pull(meteo_stat_id)

cat("抽样站点数:", length(stations_sample), "\n")
cat("气候组分布:\n")
print(meta_all %>% filter(meteo_stat_id %in% stations_sample) %>% count(koppen_group))

stations_ok <- stations_sample

# =============================================================================
# 6. CCM 函数（tp 参数化，支持任意最大 lag）
# =============================================================================
perform_ccm_vpd <- function(station_id, df, x_col, tp_x = 0, min_pts = 40) {
  d <- df %>%
    filter(meteo_stat_id == station_id, week %in% 20:39) %>%
    arrange(year, week) %>%
    mutate(x_var = dplyr::lag(!!sym(x_col), n = tp_x)) %>%
    dplyr::select(sif = sif_detrended, heat_index = x_var) %>%
    filter(!is.na(sif), !is.na(heat_index)) %>%
    mutate(time = row_number(), .before = 1) %>%
    as.data.frame()

  if (nrow(d) < min_pts) return(NULL)
  n <- nrow(d)

  tryCatch({
    best_E <- rEDM::EmbedDimension(
      dataFrame = d, columns = "sif", target = "sif",
      lib = paste("1", n), pred = paste("1", n),
      maxE = min(8, floor(n / 10)), showPlot = FALSE
    ) %>% { .$E[which.max(.$rho)] }

    ccm_res <- rEDM::CCM(
      dataFrame = d, E = best_E, Tp = 0,
      columns = "heat_index", target = "sif",
      libSizes = paste(best_E + 2, n - best_E,
                       max(2, floor((n - 2*best_E - 2) / 10))),
      sample = 50, random = TRUE, showPlot = FALSE
    )
    ccm_sum <- ccm_res %>%
      group_by(LibSize) %>%
      summarise(rho_m = mean(`heat_index:sif`, na.rm = TRUE), .groups = "drop")

    final_rho <- ccm_sum$rho_m[nrow(ccm_sum)]
    trend_rho <- cor(ccm_sum$LibSize, ccm_sum$rho_m)

    sm <- rEDM::SMap(
      dataFrame = d, E = best_E, theta = 2,
      lib = paste("1", n), pred = paste("1", n),
      columns = "heat_index", target = "sif",
      embedded = FALSE
    )
    coef_col  <- which(grepl("heat", colnames(sm$coefficients),
                              ignore.case = TRUE))[1]
    if (is.na(coef_col)) coef_col <- 2
    mean_coef <- mean(sm$coefficients[, coef_col], na.rm = TRUE)

    tibble(meteo_stat_id = station_id, x_var = x_col, tp = tp_x,
           rho = final_rho, trend = trend_rho,
           mean_coef = mean_coef, n_obs = n)
  }, error = function(e) NULL)
}

# =============================================================================
# 7. 运行 CCM（两个 X 变量 × tp=0..8）
# =============================================================================
x_vars <- c("vpd_mean_dt", "heat_over_dt")
tp_seq <- 0:8    # ← 从0..5 改为 0..8

total_calls <- length(x_vars) * length(tp_seq) * length(stations_ok)
cat(sprintf("\n开始 CCM（%d站 × %d变量 × %d tp = %d次）...\n",
            length(stations_ok), length(x_vars), length(tp_seq), total_calls))

results_all <- map_dfr(x_vars, function(xv) {
  cat("\n>>> X:", xv, "\n")
  map_dfr(tp_seq, function(tp_now) {
    cat("  tp =", tp_now, "...\n")
    map_dfr(stations_ok, ~perform_ccm_vpd(.x, df_weekly, xv, tp_now))
  })
})

# 附加元数据
meta_coords <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  distinct(meteo_stat_id, longitude, latitude, koppen_class, koppen_group)

results_all <- results_all %>%
  left_join(meta_coords, by = "meteo_stat_id")

saveRDS(results_all, OUTPUT)
cat("\n结果已保存:", OUTPUT, "\n")
cat("总行数:", nrow(results_all), "\n")
print(results_all %>% count(x_var, tp))
