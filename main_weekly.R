# Weekly Sequence Analysis Script - Fixed Window (Weeks 20-39)
# This script analyzes data completeness for a fixed window of 20 weeks (May-Sep).

pacman::p_load(dplyr, ggplot2, tidyr, targets, lubridate, purrr)

# Load data
tar_load(data_heat_sif_weekly)

cat("【周度数据连续性分析 - 固定窗口 (第20-39周)】\n")

# 1. 筛选特定窗口 (5月-9月大约对应第20-39周)
target_weeks <- 20:39
n_target <- length(target_weeks)

weekly_window_data <- data_heat_sif_weekly %>%
  filter(week %in% target_weeks)

# 2. 分析每个站点每年在窗口内的完整性
station_year_completeness <- weekly_window_data %>%
  group_by(meteo_stat_id, year) %>%
  summarise(
    n_weeks = n(),
    is_complete = (n_weeks == n_target),
    .groups = "drop"
  )

# 3. 统计：哪些站点在所有年份都保证了这20周的完整性？
station_overall_stats <- station_year_completeness %>%
  group_by(meteo_stat_id) %>%
  summarise(
    total_years = n(),
    complete_years = sum(is_complete),
    perfect_record = (total_years == complete_years),
    completion_rate = complete_years / total_years,
    .groups = "drop"
  )

# 4. 留存率分析
n_total_stations <- nrow(station_overall_stats)
n_perfect_stations <- sum(station_overall_stats$perfect_record)
retention_rate <- (n_perfect_stations / n_total_stations) * 100

cat("\n--- 统计摘要 (窗口 20-39周) ---\n")
cat("总站点数:", n_total_stations, "\n")
cat("完美站点数 (所有年份均满20周):", n_perfect_stations, "\n")
cat("完美站点留存率:", round(retention_rate, 2), "%\n")

# 5. 高质量站点分析 (例如 80% 以上年份完整的站点)
n_high_quality <- sum(station_overall_stats$completion_rate >= 0.8)
cat("高质量站点数 (>=80% 年份完整):", n_high_quality, " (", round(n_high_quality/n_total_stations*100, 2), "%)\n")

# 6. 可视化：站点完整性分布
p_dist_city <- ggplot(station_overall_stats, aes(x = completion_rate)) +
  geom_histogram(bins = 20, fill = "darkorange", color = "white") +
  labs(title = "站点年度完整性比例分布", 
       subtitle = "反映了各站点有多少比例的年份达到了20周完整要求",
       x = "年份完整率 (1.0 代表完美)", y = "站点数量") +
  theme_minimal()

# 7. 可视化：每年达到完整的站点比例随时间变化
year_summary <- station_year_completeness %>%
  group_by(year) %>%
  summarise(
    pct_complete = mean(is_complete) * 100,
    .groups = "drop"
  )

p_year_trend <- ggplot(year_summary, aes(x = year, y = pct_complete)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "各年份达到20周完整要求的站点比例",
       x = "年份", y = "完整站点占比 (%)") +
  theme_minimal()

if(!dir.exists("data_proc")) dir.create("data_proc")
ggsave("data_proc/weekly_station_completion_dist.png", p_dist_city, width = 8, height = 6)
ggsave("data_proc/weekly_year_completion_trend.png", p_year_trend, width = 8, height = 6)

# 输出完美站点清单供后续分析参考
perfect_stations <- station_overall_stats %>% 
  filter(perfect_record) %>% 
  pull(meteo_stat_id)
cat("\n完美站点前10名:", head(perfect_stations, 10), "...\n")

# 选取前 20 个完美站点进行测试
test_stations <- head(perfect_stations, 20)

# 7. 周度热胁迫与 SIF 相关性分析 (20 站点散点图) ----
cat("\n【周度热胁迫与 SIF 相关性分析 (20 站点)】\n")

# 准备绘图数据：20 个测试站点的第 20-39 周数据
plot_data_20 <- data_heat_sif_weekly %>%
  filter(meteo_stat_id %in% test_stations, week %in% 20:39)

# 7.1 线性尺度 (去趋势)
p_corr_grid <- ggplot(plot_data_20, aes(x = heat_index_composite_detrended, y = sif_detrended)) +
  geom_point(alpha = 0.4, size = 1, color = "darkgreen") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ meteo_stat_id, scales = "free", ncol = 5) +
  labs(title = "周度热胁迫与 SIF 相关性 (线性尺度)",
       subtitle = "变量：去趋势后的残差 | 数据范围：第 20-39 周",
       x = "热胁迫指数 (去趋势)", y = "SIF (去趋势)") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("data_proc/weekly_heat_sif_correlation_20_linear.png", p_corr_grid, width = 15, height = 12, dpi = 300)

# 7.2 Log-Log 尺度 (原始物理量)
# 注意：log 空间要求数据为正，这里使用 heat_over_sum 和 sif_interp，并微调 0 值
plot_data_20_log <- plot_data_20 %>%
  filter(sif_interp > 0) %>%
  mutate(heat_val = heat_over_sum + 0.01) # 避开 0 以便取 log

p_corr_log_grid <- ggplot(plot_data_20_log, aes(x = heat_val, y = sif_interp)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkorange", linetype = "dashed", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ meteo_stat_id, scales = "free", ncol = 5) +
  labs(title = "周度热胁迫与 SIF 相关性 (Log-Log 尺度)",
       subtitle = "变量：累积 VPD 强度 (x) vs 原始 SIF (y) | 双对数坐标",
       x = "log10(累积 VPD 强度 + 0.01)", y = "log10(SIF)") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("data_proc/weekly_heat_sif_correlation_20_loglog.png", p_corr_log_grid, width = 15, height = 12, dpi = 300)

# 8. 周度 CCM 因果分析 (针对 20 个高质量站点测试) ----
cat("\n【周度 CCM 因果分析测试 (20 站点)】\n")

# 定义 CCM 函数 (逻辑同步自 _targets.R，但适配周度数据)
perform_ccm_weekly <- function(station_id, data, tp_x = 0) {
  # 准备数据：仅使用第 20-39 周
  df_ccm <- data %>%
    filter(meteo_stat_id == station_id, week %in% 20:39) %>%
    arrange(year, week) %>%
    # 滞后处理
    mutate(heat_index = dplyr::lag(heat_index_composite_detrended, n = tp_x)) %>% 
    select(sif = sif_detrended, heat_index) %>%
    filter(!is.na(sif), !is.na(heat_index)) %>%
    mutate(time = row_number(), .before = 1) %>%
    as.data.frame()
  
  n_data <- nrow(df_ccm)
  if (n_data < 30) return(NULL) # 即使是周数据，CCM 也需要足够点数
  
  tryCatch({
    # 1. 确定最优 E (基于 SIF)
    embed_sif <- EmbedDimension(dataFrame = df_ccm, columns = "sif", target = "sif", 
                                lib = paste("1", n_data), pred = paste("1", n_data),
                                maxE = 8, showPlot = FALSE)
    best_E <- embed_sif$E[which.max(embed_sif$rho)]
    
    # 2. 运行 CCM: Heat -> SIF
    ccm_res <- CCM(dataFrame = df_ccm, E = best_E, Tp = 0, 
                   columns = "heat_index", target = "sif",
                   libSizes = paste(best_E+2, n_data-best_E, 10), 
                   sample = 50, random = TRUE, showPlot = FALSE)
    # Bug: 这一步处理有必要吗？
    ccm_summary <- ccm_res %>%
      group_by(LibSize) %>%
      summarise(rho_mean = mean(`heat_index:sif`, na.rm = TRUE), .groups = "drop")
    
    final_rho <- ccm_summary$rho_mean[nrow(ccm_summary)]
    trend <- cor(ccm_summary$LibSize, ccm_summary$rho_mean)
    
    # 3. S-map 确定因果性质
    smap_res <- SMap(dataFrame = df_ccm, 
                     lib = paste("1", n_data),
                     pred = paste("1", n_data),
                     E = best_E, theta = 2,
                     columns = "heat_index", target = "sif", 
                     embedded = FALSE)
    
    coeffs <- smap_res$coefficients
    # 提取热胁迫系数 (通常在第2或第3列)
    coef_col <- which(grepl("heat", colnames(coeffs), ignore.case = TRUE))[1]
    if(is.na(coef_col)) coef_col <- 2
    
    mean_coef <- mean(coeffs[, coef_col], na.rm = TRUE)
    
    # 4. 判断逻辑
    # Bug：这个判定标准是否合理？
    is_causal <- (final_rho > 0.1 & trend > 0)
    effect_type <- case_when(
      !is_causal ~ "无因果",
      mean_coef > 0 ~ "促进",
      mean_coef < 0 ~ "抑制",
      # Bug: 其实是mean_coef = 0的部分。
      TRUE ~ "未知"
    )
    
    return(tibble(
      meteo_stat_id = station_id,
      tp = tp_x,
      rho = final_rho,
      trend = trend,
      effect_type = effect_type,
      mean_coef = mean_coef,
      n_obs = n_data
    ))
  }, error = function(e) return(NULL))
}

# 选取前 20 个完美站点进行测试
# Bug：要不要测试？
# test_stations <- head(perfect_stations, 20)
test_stations <- perfect_stations

results_weekly <- map_dfr(c(0, 1, 2), function(l) {
  map_dfr(test_stations, ~perform_ccm_weekly(.x, data_heat_sif_weekly, tp_x = l))
})

# 7. 可视化测试结果
if (nrow(results_weekly) > 0) {
  p_res <- ggplot(results_weekly, aes(x = factor(tp), fill = effect_type)) +
    geom_bar(position = "dodge") +
    labs(title = "周度 CCM 测试结果 (20个高质量站点)", 
         subtitle = "不同时间滞后 (Tp) 下的因果性质分布",
         x = "时间滞后 (周)", y = "站点数量", fill = "因果性质") +
    theme_minimal()
  
  ggsave("data_proc/weekly_ccm_test_results.png", p_res, width = 8, height = 6)
  
  # 保存结果以便后续分析
  saveRDS(results_weekly, "data_proc/results_weekly_test_20.rds")
  
  cat("\nCCM 测试完成。结果已保存至 data_proc/results_weekly_test_20.rds\n")
  print(results_weekly %>% group_by(tp, effect_type) %>% summarise(n = n(), .groups = "drop"))
}

