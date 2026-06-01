# =============================================================================
# 01_ccm_vpd.R
# CCM with raw VPD variables (vpd_mean, heat_over_sum)
# Tests two X variables: vpd_mean (weekly mean VPD in kPa)
#                        heat_over_sum (weekly cumulative VPD excess > 2.0 kPa)
# Output: ccm_vpd_results.rds  (columns same as results_weekly_0_5.rds + x_var)
# =============================================================================

pacman::p_load(dplyr, purrr, targets, rEDM, tibble)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUTPUT <- "analysis_vpd2/ccm_vpd_results.rds"

# ---- 1. 加载数据 ------------------------------------------------------------
tar_load(data_heat_sif_weekly)   # 需要 targets pipeline 已跑完

# 为 vpd_mean 和 heat_over_sum 添加去趋势版本
safe_detrend <- function(x, t) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  # na.exclude: 残差按原始行位置补NA，长度与输入一致
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
    time_idx      = row_number(),
    vpd_mean_dt   = safe_detrend(vpd_mean,     time_idx),
    heat_over_dt  = safe_detrend(heat_over_sum, time_idx)
  ) %>%
  ungroup()

cat("总周记录数:", nrow(df_weekly), "\n")

# ---- 2. CCM 函数 ------------------------------------------------------------
perform_ccm_vpd <- function(station_id, df, x_col, tp_x = 0, min_pts = 30) {
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
    # 最优 E
    best_E <- rEDM::EmbedDimension(
      dataFrame = d, columns = "sif", target = "sif",
      lib = paste("1", n), pred = paste("1", n),
      maxE = min(8, floor(n / 10)), showPlot = FALSE
    ) %>% { .$E[which.max(.$rho)] }

    # CCM
    ccm_res <- rEDM::CCM(
      dataFrame = d, E = best_E, Tp = 0,
      columns = "heat_index", target = "sif",
      libSizes = paste(best_E + 2, n - best_E, max(2, floor((n - best_E - best_E - 2) / 10))),
      sample = 50, random = TRUE, showPlot = FALSE
    )
    ccm_sum <- ccm_res %>%
      group_by(LibSize) %>%
      summarise(rho_m = mean(`heat_index:sif`, na.rm = TRUE), .groups = "drop")

    final_rho <- ccm_sum$rho_m[nrow(ccm_sum)]
    trend_rho  <- cor(ccm_sum$LibSize, ccm_sum$rho_m)

    # S-map
    sm <- rEDM::SMap(
      dataFrame = d, E = best_E, theta = 2,
      lib = paste("1", n), pred = paste("1", n),
      columns = "heat_index", target = "sif",
      embedded = FALSE
    )
    coef_col <- which(grepl("heat", colnames(sm$coefficients), ignore.case = TRUE))[1]
    if (is.na(coef_col)) coef_col <- 2
    coef_vec  <- sm$coefficients[, coef_col]
    mean_coef <- mean(coef_vec, na.rm = TRUE)

    tibble(
      meteo_stat_id = station_id,
      x_var         = x_col,
      tp            = tp_x,
      rho           = final_rho,
      trend         = trend_rho,
      mean_coef     = mean_coef,
      n_obs         = n
    )
  }, error = function(e) NULL)
}

# ---- 3. 获取候选站点 --------------------------------------------------------
stations_ok <- df_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(vpd_mean_dt)) >= 30, sum(!is.na(sif_detrended)) >= 30) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(meteo_stat_id)

cat("候选站点总数:", length(stations_ok), "\n")

# ---- 3b. 分层抽样（约30分钟版本）------------------------------------------
# 全量约14小时（747站 × 2变量 × 6tp）；按比例抽样到约25站，预计30分钟内完成
# 按 Köppen 气候组分层，保证各组有代表性
meta <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  distinct(meteo_stat_id, koppen_group) %>%
  filter(meteo_stat_id %in% stations_ok)

N_SAMPLE <- 25   # 调整此值控制运行时间：25≈30min，50≈1h，747≈14h

set.seed(42)
stations_sample <- meta %>%
  group_by(koppen_group) %>%
  slice_sample(prop = N_SAMPLE / nrow(meta), replace = FALSE) %>%
  ungroup() %>%
  # 确保至少取到 N_SAMPLE 个（因四舍五入可能少一点）
  { if (nrow(.) < N_SAMPLE) bind_rows(.,
      anti_join(meta, ., by="meteo_stat_id") %>% slice_sample(n = N_SAMPLE - nrow(.)))
    else slice_sample(., n = N_SAMPLE) } %>%
  pull(meteo_stat_id)

cat("抽样站点数:", length(stations_sample), "\n")
cat("气候组分布:\n")
print(meta %>% filter(meteo_stat_id %in% stations_sample) %>% count(koppen_group))

stations_ok <- stations_sample   # 用抽样结果覆盖

# ---- 4. 运行 CCM（两个 X 变量 × tp=0..5）-----------------------------------
x_vars <- c("vpd_mean_dt", "heat_over_dt")
tp_seq <- 0:5

cat("开始 CCM 分析（共", length(x_vars) * length(tp_seq) * length(stations_ok),
    "次调用）...\n")

results_all <- map_dfr(x_vars, function(xv) {
  cat("\n>>> X variable:", xv, "\n")
  map_dfr(tp_seq, function(tp_now) {
    cat("  tp =", tp_now, "...\n")
    map_dfr(stations_ok, ~perform_ccm_vpd(.x, df_weekly, xv, tp_now))
  })
})

# ---- 5. 附加元数据（坐标、气候）--------------------------------------------
# 从原始 results_weekly_0_5.rds 读取元数据
meta <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  distinct(meteo_stat_id, longitude, latitude, koppen_class, koppen_group)

results_all <- results_all %>%
  left_join(meta, by = "meteo_stat_id")

saveRDS(results_all, OUTPUT)
cat("\n结果已保存:", OUTPUT, "\n")
cat("总行数:", nrow(results_all), "\n")
print(results_all %>% count(x_var, tp))
