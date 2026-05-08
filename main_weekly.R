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
perfect_stations <- station_overall_stats %>% filter(perfect_record) %>% pull(meteo_stat_id)
cat("\n完美站点前10名:", head(perfect_stations, 10), "...\n")
