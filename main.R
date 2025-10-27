pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, tidyr, showtext)
showtext_auto()

# 读取数据。
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  # 计算SIF月均值。
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop") %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 画出各站点SIF每年各月份的变化。
ggplot(head(meteo_sif_data, 10000)) +
  geom_line(aes(x = month, y = sif, color = year, group = year)) +
  theme(legend.position = "none") + 
  facet_wrap(.~meteo_stat_id) 

# 读取气象站其他数据。
# 读取第一批气象站点数据。
# 函数：读取单个气象文件。
read_one_meteo <- function(path) {
  # 从文件名提取站点编号。
  station_id <- str_extract(basename(path), "\\d+")
  
  # 读取数据正文：从第2行开始。
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    # 把超大缺测码（如999990/999999等）统一置NA。
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1) %>%
    return()
}

# 读入目标站点并合并。
# 获取所有txt文件的文件名列表。
meteo_file_lin <- list.files(
  path = "data_raw/meteo_data_1961-2023",
  pattern = "\\.txt$", # 匹配所有以.txt结尾的文件。
  full.names = TRUE
) %>%
  .[!grepl("sta_lonlat_china.txt", .)]
# 气象变量。
var_meteo_lin <- c("tmax", "tmin", "tavg", "rh", "precip")
# 目标变量列表。
meteo_data_lin <- map(
  # 目标文件。
  meteo_file_lin[
    str_extract(basename(meteo_file_lin), "\\d+") %in%
      as.character(meteo_sif_data$meteo_stat_id)
  ],
  read_one_meteo
) %>%
  list_rbind() %>%
  rename_with(~ tolower(.x)) %>% 
  # 计算各气象变量每月均值。
  mutate(year = year(date), month = month(date)) %>% 
  group_by(meteo_stat_id, year, month) %>% 
  summarise(
    across(all_of(var_meteo_lin), ~mean(.x, na.rm = TRUE)), .groups = "drop"
  )

# 合并气象数据和SIF数据。
meteo_data <- meteo_sif_data %>%
  inner_join(meteo_data_lin, by = c("meteo_stat_id", "year", "month")) %>% 
  # 移除含有NA的行。
  drop_na() %>% 
  arrange(meteo_stat_id, year, month)

# ============================================================================
# 方案A：只分析6-8月数据，计算温度异常
# ============================================================================

# 1. 筛选6-8月数据并计算温度异常
meteo_data_summer <- meteo_data %>%
  filter(month %in% 6:8) %>%
  group_by(meteo_stat_id, month) %>%
  mutate(
    # 计算该站点该月份的90%分位数（历史阈值）
    temp_q90 = quantile(tmax, 0.9, na.rm = TRUE),
    # 计算温度异常（所有月份都计算）
    temp_anomaly = tmax - temp_q90,
    # 标记是否为热月
    is_hot_month = tmax > temp_q90
  ) %>%
  ungroup() %>%
  arrange(meteo_stat_id, year, month)

# 查看热月的分布
meteo_data_summer %>%
  group_by(meteo_stat_id) %>%
  summarise(
    n_total = n(),
    n_hot = sum(is_hot_month),
    pct_hot = round(n_hot / n_total * 100, 1),
    mean_anomaly_hot = mean(temp_anomaly[is_hot_month], na.rm = TRUE)
  ) %>%
  print(n = 20)

# ============================================================================
# 2. 可视化检查数据特征
# ============================================================================

# 2.1 温度异常的时间序列
meteo_data_summer %>%
  filter(meteo_stat_id %in% unique(meteo_stat_id)[1:9]) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>%
  ggplot(aes(x = date, y = temp_anomaly)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(color = is_hot_month), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  facet_wrap(~ meteo_stat_id, scales = "free") +
  labs(
    title = "Temperature Anomaly (6-8月)",
    subtitle = "红点 = 热月（超过90%分位数）",
    x = "Date", y = "Temperature Anomaly (°C)"
  ) +
  theme_minimal()

# 2.2 SIF的时间序列
meteo_data_summer %>%
  filter(meteo_stat_id %in% unique(meteo_stat_id)[1:9]) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>%
  ggplot(aes(x = date, y = sif)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  facet_wrap(~ meteo_stat_id, scales = "free") +
  labs(
    title = "SIF Time Series (6-8月)",
    subtitle = "红线 = 线性趋势",
    x = "Date", y = "SIF"
  ) +
  theme_minimal()

# ============================================================================
# 3. 检查长期趋势（添加错误处理）
# ============================================================================

trend_analysis <- meteo_data_summer %>%
  mutate(year_index = year - min(year)) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    # 数据点数
    n = n(),
    n_complete = sum(complete.cases(year_index, temp_anomaly, sif)),
    
    # 相关系数
    cor_temp_anomaly = cor(year_index, temp_anomaly, use = "complete.obs"),
    cor_sif = cor(year_index, sif, use = "complete.obs"),
    
    # 线性趋势显著性（添加错误处理）
    p_temp = tryCatch({
      if (sum(complete.cases(year_index, temp_anomaly)) >= 3) {
        cor.test(year_index, temp_anomaly)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    
    p_sif = tryCatch({
      if (sum(complete.cases(year_index, sif)) >= 3) {
        cor.test(year_index, sif)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
  ) %>%
  mutate(
    # 判断是否需要去趋势（相关系数>0.3 且 p<0.05）
    need_detrend_temp = abs(cor_temp_anomaly) > 0.3 & p_temp < 0.05 & !is.na(p_temp),
    need_detrend_sif = abs(cor_sif) > 0.3 & p_sif < 0.05 & !is.na(p_sif)
  )

# 查看趋势分析结果
print(trend_analysis, n = 20)

# 检查有问题的站点
trend_analysis %>%
  filter(is.na(p_temp) | is.na(p_sif) | n_complete < 10) %>%
  arrange(n_complete) %>%
  print()

# 统计需要去趋势的站点数量
trend_summary <- trend_analysis %>%
  summarise(
    n_stations_total = n(),
    n_stations_valid = sum(n_complete >= 10, na.rm = TRUE),
    n_need_detrend_temp = sum(need_detrend_temp, na.rm = TRUE),
    n_need_detrend_sif = sum(need_detrend_sif, na.rm = TRUE),
    n_need_detrend_both = sum(need_detrend_temp & need_detrend_sif, na.rm = TRUE),
    pct_need_detrend_temp = round(n_need_detrend_temp / n_stations_valid * 100, 1),
    pct_need_detrend_sif = round(n_need_detrend_sif / n_stations_valid * 100, 1)
  )

print(trend_summary)

# 可视化趋势分布（只包含有效站点）
trend_analysis %>%
  filter(n_complete >= 10) %>%  # 过滤掉数据点太少的站点
  ggplot(aes(x = cor_temp_anomaly, y = cor_sif)) +
  geom_point(aes(color = need_detrend_temp | need_detrend_sif), 
             size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                     labels = c("无需去趋势", "需要去趋势")) +
  labs(
    title = "站点趋势分析",
    subtitle = paste0("只显示数据点≥10的站点 (n=", 
                      sum(trend_analysis$n_complete >= 10), ")"),
    x = "Temperature Anomaly - Year Correlation",
    y = "SIF - Year Correlation",
    color = ""
  ) +
  theme_minimal()

# 可视化趋势分布的直方图
library(gridExtra)

p1 <- trend_analysis %>%
  filter(n_complete >= 10) %>%
  ggplot(aes(x = cor_temp_anomaly)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
  labs(title = "Temperature Anomaly Trend Distribution",
       x = "Correlation with Year", y = "Count") +
  theme_minimal()

p2 <- trend_analysis %>%
  filter(n_complete >= 10) %>%
  ggplot(aes(x = cor_sif)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
  labs(title = "SIF Trend Distribution",
       x = "Correlation with Year", y = "Count") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)


# ============================================================================
# 4. 去趋势处理
# ============================================================================

# 方法1：对所有站点统一去趋势
meteo_data_detrended_all <- meteo_data_summer %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    temp_anomaly_detrended = residuals(lm(temp_anomaly ~ time_idx)),
    sif_detrended = residuals(lm(sif ~ time_idx))
  ) %>%
  ungroup()


# ============================================================================
# 5. 对比去趋势前后的效果
# ============================================================================

# 选择几个典型站点可视化
example_stations <- unique(meteo_data_summer$meteo_stat_id)[1:6]

# 温度异常对比
meteo_data_detrended_all %>%
  filter(meteo_stat_id %in% example_stations) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>%
  pivot_longer(
    cols = c(temp_anomaly, temp_anomaly_detrended),
    names_to = "type",
    values_to = "value"
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(
    values = c("temp_anomaly" = "blue", "temp_anomaly_detrended" = "red"),
    labels = c("原始异常", "去趋势后")
  ) +
  facet_wrap(~ meteo_stat_id, scales = "free") +
  labs(
    title = "温度异常：去趋势前后对比",
    x = "Date", y = "Temperature Anomaly (°C)",
    color = ""
  ) +
  theme_minimal()

# SIF对比
meteo_data_detrended_all %>%
  filter(meteo_stat_id %in% example_stations) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>%
  pivot_longer(
    cols = c(sif, sif_detrended),
    names_to = "type",
    values_to = "value"
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(
    values = c("sif" = "darkgreen", "sif_detrended" = "orange"),
    labels = c("原始SIF", "去趋势后")
  ) +
  facet_wrap(~ meteo_stat_id, scales = "free") +
  labs(
    title = "SIF：去趋势前后对比",
    x = "Date", y = "SIF",
    color = ""
  ) +
  theme_minimal()

# ============================================================================
# 6. 导出处理后的数据
# ============================================================================

# 根据你的选择，导出相应的数据集
# 使用去趋势后的数据（推荐）
data_for_ccm <- meteo_data_detrended_all %>%
  select(meteo_stat_id, year, month, 
         temp_anomaly_original = temp_anomaly,
         temp_anomaly = temp_anomaly_detrended,
         sif_original = sif,
         sif = sif_detrended,
         is_hot_month)
