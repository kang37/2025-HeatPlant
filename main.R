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

# ============================================================================
# CCM分析准备（修正版 - 解决库长度问题）
# ============================================================================

library(rEDM)

# Example station ----
# ============================================================================
# 1. 单个站点的CCM分析示例
# ============================================================================

# ============================================================================
# 完整的滞后因果分析：t月温度异常 → t+1月、t+2月SIF
# ============================================================================

# ============================================================================
# 1. 数据准备
# ============================================================================

# 选择一个数据点充足的站点
example_station <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  slice(1) %>%
  pull(meteo_stat_id)

# 提取单个站点数据（去趋势后的数据）
single_station_data <- data_for_ccm %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  select(temp_anomaly, sif) %>%
  as.data.frame()

# 添加时间索引列
single_station_data <- data.frame(
  time = 1:nrow(single_station_data),
  temp_anomaly = single_station_data$temp_anomaly,
  sif = single_station_data$sif
)

cat("=== 站点信息 ===\n")
cat("站点ID:", example_station, "\n")
cat("数据点数:", nrow(single_station_data), "\n")
cat("时间跨度: 6-8月的数据\n\n")

n_data <- nrow(single_station_data)

# ============================================================================
# 2. 创建滞后变量（用于相关性分析）
# ============================================================================

lag_data <- single_station_data %>%
  mutate(
    # t+1月的SIF（lead函数：向前看1个位置）
    sif_lag1 = lead(sif, 1),
    # t+2月的SIF
    sif_lag2 = lead(sif, 2),
    # 同时也保留当月SIF（t+0）
    sif_lag0 = sif
  )

cat("=== 滞后变量创建完成 ===\n")
cat("可用于分析的数据点:\n")
cat("- Tp=0 (同期):", sum(!is.na(lag_data$sif_lag0)), "个\n")
cat("- Tp=1 (滞后1月):", sum(!is.na(lag_data$sif_lag1)), "个\n")
cat("- Tp=2 (滞后2月):", sum(!is.na(lag_data$sif_lag2)), "个\n\n")

# ============================================================================
# 3. 确定最优嵌入维度
# ============================================================================

lib_pred_size <- min(100, n_data)

embed_temp <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  maxE = 10,
  showPlot = FALSE
)

embed_sif <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "sif",
  target = "sif",
  maxE = 10,
  showPlot = FALSE
)

best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
E_for_ccm <- max(best_E_temp, best_E_sif)

cat("=== 嵌入维度 ===\n")
cat("Temperature Anomaly 最优 E:", best_E_temp, "\n")
cat("SIF 最优 E:", best_E_sif, "\n")
cat("CCM使用 E:", E_for_ccm, "\n\n")

# ============================================================================
# 4. CCM分析：测试不同滞后期的因果关系
# ============================================================================

# 计算库长度参数
tau <- 1
embedding_loss <- (E_for_ccm - 1) * tau
max_available_lib <- n_data - embedding_loss
lib_start <- max(E_for_ccm + 2, 10)
lib_end <- max_available_lib
lib_step <- max(2, floor((lib_end - lib_start) / 15))
libSizes_str <- paste(lib_start, lib_end, lib_step)

cat("=== 开始CCM分析 ===\n")
cat("测试滞后期: Tp = 0, 1, 2\n")
cat("嵌入维度 E:", E_for_ccm, "\n")
cat("库长度范围:", libSizes_str, "\n\n")

# 存储CCM结果
ccm_results_list <- list()
ccm_summary_list <- list()

for (tp in 0:2) {
  cat("正在分析 Tp =", tp, "...\n")
  
  # 调整最大库长度（因为滞后会损失数据点）
  adjusted_max_lib <- max_available_lib - tp
  adjusted_lib_end <- min(lib_end, adjusted_max_lib)
  adjusted_libSizes_str <- paste(lib_start, adjusted_lib_end, lib_step)
  
  tryCatch({
    ccm_result <- CCM(
      dataFrame = single_station_data,
      E = E_for_ccm,
      Tp = tp,
      columns = "temp_anomaly",
      target = "sif",
      libSizes = adjusted_libSizes_str,
      sample = 100,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 保存原始结果
    ccm_results_list[[paste0("Tp", tp)]] <- ccm_result
    
    # 提取汇总
    ccm_summary <- ccm_result %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Tp = tp)
    
    ccm_summary_list[[paste0("Tp", tp)]] <- ccm_summary
    
  }, error = function(e) {
    cat("  ✗ Tp =", tp, "分析失败:", e$message, "\n")
  })
}

# 合并所有Tp的汇总数据
ccm_all_summary <- bind_rows(ccm_summary_list)

cat("\n✓ CCM分析完成\n\n")

# ============================================================================
# 5. 相关性分析：测试不同滞后期的相关性
# ============================================================================

cat("=== 相关性分析 ===\n\n")

correlation_results <- list()

# Tp = 0 (同期)
cor_tp0 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag0)
correlation_results[["Tp0"]] <- tibble(
  Tp = 0,
  r = cor_tp0$estimate,
  p_value = cor_tp0$p.value,
  ci_lower = cor_tp0$conf.int[1],
  ci_upper = cor_tp0$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag0))
)

cat("Tp = 0 (同期: t月温度 vs t月SIF)\n")
cat("  Pearson r =", round(cor_tp0$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp0$p.value), "\n")
cat("  95% CI: [", round(cor_tp0$conf.int[1], 3), ",", 
    round(cor_tp0$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp0$estimate < 0, 
                 "✓ 负相关：温度↑ SIF↓", 
                 "✗ 正相关：温度↑ SIF↑"), "\n\n")

# Tp = 1 (滞后1月)
cor_tp1 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag1)
correlation_results[["Tp1"]] <- tibble(
  Tp = 1,
  r = cor_tp1$estimate,
  p_value = cor_tp1$p.value,
  ci_lower = cor_tp1$conf.int[1],
  ci_upper = cor_tp1$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag1))
)

cat("Tp = 1 (滞后1月: t月温度 vs t+1月SIF)\n")
cat("  Pearson r =", round(cor_tp1$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp1$p.value), "\n")
cat("  95% CI: [", round(cor_tp1$conf.int[1], 3), ",", 
    round(cor_tp1$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp1$estimate < 0, 
                 "✓ 负相关：t月温度↑ → t+1月SIF↓", 
                 "✗ 正相关：t月温度↑ → t+1月SIF↑"), "\n\n")

# Tp = 2 (滞后2月)
cor_tp2 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag2)
correlation_results[["Tp2"]] <- tibble(
  Tp = 2,
  r = cor_tp2$estimate,
  p_value = cor_tp2$p.value,
  ci_lower = cor_tp2$conf.int[1],
  ci_upper = cor_tp2$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag2))
)

cat("Tp = 2 (滞后2月: t月温度 vs t+2月SIF)\n")
cat("  Pearson r =", round(cor_tp2$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp2$p.value), "\n")
cat("  95% CI: [", round(cor_tp2$conf.int[1], 3), ",", 
    round(cor_tp2$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp2$estimate < 0, 
                 "✓ 负相关：t月温度↑ → t+2月SIF↓", 
                 "✗ 正相关：t月温度↑ → t+2月SIF↑"), "\n\n")

# 合并相关性结果
correlation_summary <- bind_rows(correlation_results)

# ============================================================================
# 6. 提取CCM的最终ρ值（用于综合报告）
# ============================================================================

ccm_final_rho <- ccm_all_summary %>%
  group_by(Tp) %>%
  filter(lib_size == max(lib_size)) %>%
  summarise(
    rho_final = mean(rho_mean),
    .groups = "drop"
  )

# 计算CCM趋势（收敛性）
ccm_trends <- ccm_all_summary %>%
  group_by(Tp) %>%
  summarise(
    trend = cor(lib_size, rho_mean),
    .groups = "drop"
  )

ccm_final_summary <- ccm_final_rho %>%
  left_join(ccm_trends, by = "Tp")

# ============================================================================
# 7. 综合结果表格
# ============================================================================

comprehensive_results <- correlation_summary %>%
  left_join(ccm_final_summary, by = "Tp") %>%
  mutate(
    causality = case_when(
      rho_final > 0.15 & trend > 0 ~ "有因果",
      TRUE ~ "无因果"
    ),
    effect_direction = case_when(
      r < 0 ~ "负向(温度↑ SIF↓)",
      r > 0 ~ "正向(温度↑ SIF↑)",
      TRUE ~ "无关"
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("=== 综合结果汇总 ===\n\n")
print(comprehensive_results %>%
        select(Tp, r, p_value, significance, rho_final, trend, 
               causality, effect_direction))

# ============================================================================
# 8. 可视化：相关性散点图
# ============================================================================

# 准备绘图数据
plot_data_tp0 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag0) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=0 (同期)")

plot_data_tp1 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag1) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=1 (滞后1月)")

plot_data_tp2 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag2) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=2 (滞后2月)")

plot_data_all <- bind_rows(plot_data_tp0, plot_data_tp1, plot_data_tp2) %>%
  mutate(Tp = factor(Tp, levels = c("Tp=0 (同期)", 
                                    "Tp=1 (滞后1月)", 
                                    "Tp=2 (滞后2月)")))

# 创建散点图
p_scatter <- ggplot(plot_data_all, aes(x = temp_anomaly, y = sif)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  facet_wrap(~ Tp, ncol = 3) +
  labs(
    title = paste("站点", example_station, "- 相关性分析"),
    subtitle = "温度异常 vs SIF（不同滞后期）",
    x = "温度异常 (°C)",
    y = "SIF",
    caption = "蓝线：线性拟合；阴影：95%置信区间"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11)
  )

# 添加相关系数标注
correlation_labels <- comprehensive_results %>%
  mutate(
    label = paste0("r = ", round(r, 3), significance, 
                   "\np = ", format.pval(p_value, digits = 2)),
    Tp = factor(paste0("Tp=", Tp, " (", 
                       c("同期", "滞后1月", "滞后2月")[Tp + 1], ")"),
                levels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)"))
  )

p_scatter_final <- p_scatter +
  geom_text(data = correlation_labels,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, 
            color = "red", fontface = "bold")

print(p_scatter_final)

# ============================================================================
# 9. 可视化：CCM收敛图
# ============================================================================

# CCM图：显示收敛性
p_ccm <- ggplot(ccm_all_summary, aes(x = lib_size, y = rho_mean, 
                                     color = factor(Tp))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = rho_mean - rho_sd, 
                  ymax = rho_mean + rho_sd,
                  fill = factor(Tp)),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.15, linetype = "dashed", 
             color = "gray", alpha = 0.7) +
  scale_color_manual(
    values = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A"),
    labels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)")
  ) +
  scale_fill_manual(
    values = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A"),
    labels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)")
  ) +
  labs(
    title = paste("站点", example_station, "- CCM因果分析"),
    subtitle = "温度异常 → SIF（不同滞后期的因果强度）",
    x = "Library Size",
    y = "Cross Map Skill (ρ)",
    color = "滞后期",
    fill = "滞后期",
    caption = "ρ > 0.15 且上升趋势 → 有因果关系"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_ccm)


# ============================================================================
# 10. 可视化：相关系数和CCM ρ随滞后期变化
# ============================================================================

# 准备综合对比数据
comparison_data <- comprehensive_results %>%
  select(Tp, r, rho_final) %>%
  pivot_longer(cols = c(r, rho_final),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    "r" = "相关系数 (Pearson r)",
                    "rho_final" = "CCM因果强度 (ρ)")
  )

p_comparison <- ggplot(comparison_data, aes(x = Tp, y = value, 
                                            color = metric, 
                                            group = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0.15, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = 0:2, 
                     labels = c("0\n(同期)", "1\n(滞后1月)", "2\n(滞后2月)")) +
  scale_color_manual(values = c("相关系数 (Pearson r)" = "#E41A1C",
                                "CCM因果强度 (ρ)" = "#377EB8")) +
  labs(
    title = paste("站点", example_station, "- 滞后效应综合对比"),
    subtitle = "相关性 vs 因果性随滞后期的变化",
    x = "滞后期 (月)",
    y = "强度",
    color = "指标",
    caption = "红色虚线(0.15)：CCM因果显著性阈值"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_comparison)


# ============================================================================
# 11. 创建组合图
# ============================================================================

combined_plot <- (p_scatter_final) / (p_ccm + p_comparison) +
  plot_annotation(
    title = paste("站点", example_station, 
                  "- 温度异常对SIF的滞后因果效应完整分析"),
    subtitle = "上：相关性分析（因果性质） | 左下：CCM收敛性（因果方向） | 右下：综合对比",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)


# ============================================================================
# 12. 最终结论报告
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("                   最终因果关系结论\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("站点:", example_station, "\n")
cat("数据: 去趋势后的6-8月温度异常和SIF\n\n")

for (tp in 0:2) {
  result <- comprehensive_results %>% filter(Tp == tp)
  
  cat("【Tp =", tp, "】", 
      c("同期效应", "滞后1月效应", "滞后2月效应")[tp + 1], "\n")
  
  # 因果方向（CCM）
  if (result$causality == "有因果") {
    cat("  ✓ 因果关系: 温度异常 → SIF (ρ =", 
        round(result$rho_final, 3), ")\n")
  } else {
    cat("  ✗ 无显著因果 (ρ =", round(result$rho_final, 3), ")\n")
  }
  
  # 因果性质（相关性）
  if (result$p_value < 0.05) {
    if (result$r < 0) {
      cat("  ✓ 效应方向: 负向影响 - 温度异常↑ → SIF↓\n")
      cat("    (r =", round(result$r, 3), result$significance, ")\n")
      cat("  → 符合研究预期！\n")
    } else {
      cat("  ✗ 效应方向: 正向影响 - 温度异常↑ → SIF↑\n")
      cat("    (r =", round(result$r, 3), result$significance, ")\n")
      cat("  → 不符合研究预期\n")
    }
  } else {
    cat("  - 相关性不显著 (r =", round(result$r, 3), ", p > 0.05)\n")
  }
  
  cat("\n")
}

# 找出最强的滞后效应
best_lag <- comprehensive_results %>%
  filter(causality == "有因果", r < 0, p_value < 0.05) %>%
  arrange(desc(abs(r))) %>%
  slice(1)

if (nrow(best_lag) > 0) {
  cat("*** 最佳滞后期 ***\n")
  cat("Tp =", best_lag$Tp, "月滞后效应最强\n")
  cat("- CCM因果强度: ρ =", round(best_lag$rho_final, 3), "\n")
  cat("- 负向影响强度: r =", round(best_lag$r, 3), 
      best_lag$significance, "\n")
  cat("→ t月温度异常每升高1°C，t+", best_lag$Tp, 
      "月SIF预期降低约", round(abs(best_lag$r), 3), "个标准差\n")
} else {
  cat("*** 警告 ***\n")
  cat("未发现显著的负向滞后因果效应！\n")
  cat("建议：\n")
  cat("1. 检查是否只分析了真正的热月数据\n")
  cat("2. 尝试更长的滞后期（Tp = 3, 4, 5...）\n")
  cat("3. 检查数据质量和站点选择\n")
}

cat("\n", rep("=", 70), "\n", sep = "")

# ============================================================================
# 批量分析所有站点：Tp=1滞后因果分析（增强版）
# ============================================================================

library(rEDM)
library(tidyverse)
library(patchwork)

# ============================================================================
# 1. 增强的CCM分析函数（包含相关性检验）
# ============================================================================

perform_ccm_with_correlation <- function(station_id, data, min_points = 30) {
  
  # 提取该站点数据
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif) %>%
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  # 检查数据点数
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    # 创建滞后变量
    station_data_with_lag <- station_data %>%
      mutate(
        sif_lag1 = lead(sif, 1)  # t+1月的SIF
      )
    
    # 添加时间列（用于CCM）
    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif
    )
    
    lib_pred_size <- min(100, n_data)
    
    # ===== 确定最优E =====
    embed_temp <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "temp_anomaly",
      target = "temp_anomaly",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    E_ccm <- max(best_E_temp, best_E_sif)
    
    # ===== 计算最大库长度 =====
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    tp <- 1  # 滞后1月
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) {
      return(NULL)
    }
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # ===== CCM分析（Tp=1）=====
    ccm_result <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = 1,  # 滞后1月
      columns = "temp_anomaly",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 提取汇总结果
    ccm_summary <- ccm_result %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_temp_to_sif_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_temp_to_sif_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sif_to_temp_mean = mean(`sif:temp_anomaly`, na.rm = TRUE),
        rho_sif_to_temp_sd = sd(`sif:temp_anomaly`, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 获取最终rho值
    final_rho_temp_to_sif <- ccm_summary$rho_temp_to_sif_mean[nrow(ccm_summary)]
    final_rho_sif_to_temp <- ccm_summary$rho_sif_to_temp_mean[nrow(ccm_summary)]
    
    # 计算趋势（收敛性）
    trend_temp_to_sif <- cor(ccm_summary$lib_size, 
                             ccm_summary$rho_temp_to_sif_mean)
    trend_sif_to_temp <- cor(ccm_summary$lib_size, 
                             ccm_summary$rho_sif_to_temp_mean)
    
    # ===== 相关性分析 =====
    cor_result <- cor.test(station_data_with_lag$temp_anomaly, 
                           station_data_with_lag$sif_lag1)
    
    correlation_r <- cor_result$estimate
    correlation_p <- cor_result$p.value
    
    # ===== 判断因果关系 =====
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    temp_causes_sif <- (final_rho_temp_to_sif > ccm_threshold_rho & 
                          trend_temp_to_sif > ccm_threshold_trend)
    sif_causes_temp <- (final_rho_sif_to_temp > ccm_threshold_rho & 
                          trend_sif_to_temp > ccm_threshold_trend)
    
    causality_type <- case_when(
      temp_causes_sif & !sif_causes_temp ~ "Temp → SIF",
      !temp_causes_sif & sif_causes_temp ~ "SIF → Temp",
      temp_causes_sif & sif_causes_temp ~ "双向因果",
      TRUE ~ "无显著因果"
    )
    
    # 判断效应方向
    effect_direction <- case_when(
      correlation_p >= 0.05 ~ "不显著",
      correlation_r < 0 ~ "负向(温度↑ SIF↓)",
      correlation_r > 0 ~ "正向(温度↑ SIF↑)",
      TRUE ~ "无关"
    )
    
    # ===== 创建CCM图表 =====
    ccm_plot <- ggplot(ccm_summary) +
      geom_line(aes(x = lib_size, y = rho_temp_to_sif_mean, 
                    color = "Temp → SIF"), linewidth = 1) +
      geom_point(aes(x = lib_size, y = rho_temp_to_sif_mean, 
                     color = "Temp → SIF"), size = 2) +
      geom_ribbon(aes(x = lib_size, 
                      ymin = rho_temp_to_sif_mean - rho_temp_to_sif_sd,
                      ymax = rho_temp_to_sif_mean + rho_temp_to_sif_sd,
                      fill = "Temp → SIF"),
                  alpha = 0.2) +
      geom_line(aes(x = lib_size, y = rho_sif_to_temp_mean, 
                    color = "SIF → Temp"), linewidth = 1) +
      geom_point(aes(x = lib_size, y = rho_sif_to_temp_mean, 
                     color = "SIF → Temp"), size = 2) +
      geom_ribbon(aes(x = lib_size, 
                      ymin = rho_sif_to_temp_mean - rho_sif_to_temp_sd,
                      ymax = rho_sif_to_temp_mean + rho_sif_to_temp_sd,
                      fill = "SIF → Temp"),
                  alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = ccm_threshold_rho, linetype = "dashed", 
                 color = "red", alpha = 0.3) +
      scale_color_manual(values = c("Temp → SIF" = "#377EB8", 
                                    "SIF → Temp" = "#E41A1C")) +
      scale_fill_manual(values = c("Temp → SIF" = "#377EB8", 
                                   "SIF → Temp" = "#E41A1C")) +
      labs(
        title = paste("站点", station_id, "- CCM分析 (Tp=1)"),
        subtitle = paste0(
          "因果: ", causality_type, " | ",
          "ρ(T→S)=", round(final_rho_temp_to_sif, 2), " | ",
          "ρ(S→T)=", round(final_rho_sif_to_temp, 2)
        ),
        x = "Library Size",
        y = "Cross Map Skill (ρ)",
        color = "方向",
        fill = "方向"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # ===== 创建相关性散点图 =====
    cor_plot <- ggplot(station_data_with_lag, 
                       aes(x = temp_anomaly, y = sif_lag1)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
      labs(
        title = paste("站点", station_id, "- 相关性分析 (Tp=1)"),
        subtitle = paste0(
          "效应: ", effect_direction, " | ",
          "r=", round(correlation_r, 3), " | ",
          "p=", format.pval(correlation_p, digits = 2)
        ),
        x = "t月温度异常 (°C)",
        y = "t+1月SIF"
      ) +
      theme_minimal()
    
    # ===== 返回结果 =====
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      max_lib_size = max_available_lib,
      E = E_ccm,
      best_E_temp = best_E_temp,
      best_E_sif = best_E_sif,
      rho_temp_to_sif = final_rho_temp_to_sif,
      trend_temp_to_sif = trend_temp_to_sif,
      rho_sif_to_temp = final_rho_sif_to_temp,
      trend_sif_to_temp = trend_sif_to_temp,
      temp_causes_sif = temp_causes_sif,
      sif_causes_temp = sif_causes_temp,
      causality_type = causality_type,
      correlation_r = correlation_r,
      correlation_p = correlation_p,
      effect_direction = effect_direction,
      ccm_summary_data = list(ccm_summary),
      ccm_raw_data = list(ccm_result),
      correlation_data = list(station_data_with_lag),
      ccm_plot = list(ccm_plot),
      cor_plot = list(cor_plot)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# ============================================================================
# 2. 批量分析所有站点
# ============================================================================

all_stations <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("\n准备分析", length(all_stations), "个站点...\n")
cat("分析内容: Tp=1 (t月温度 → t+1月SIF)\n\n")

# 简单进度显示（不使用progress包）
ccm_results_all <- map_dfr(seq_along(all_stations), function(i) {
  if (i %% 50 == 0) {
    cat("已完成:", i, "/", length(all_stations), "\n")
  }
  perform_ccm_with_correlation(all_stations[i], data_for_ccm, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_all), "/", length(all_stations), "\n")

# ============================================================================
# 3. 结果统计和分类
# ============================================================================

cat("\n=== CCM分析结果统计 ===\n\n")

causality_stats <- table(ccm_results_all$causality_type)
print(causality_stats)

cat("\n各因果类型占比:\n")
print(round(prop.table(causality_stats) * 100, 1))

effect_stats <- ccm_results_all %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  count(effect_direction) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\n温度→SIF的效应方向统计:\n")
print(effect_stats)

detailed_stats <- ccm_results_all %>%
  mutate(
    has_temp_to_sif = causality_type %in% c("Temp → SIF", "双向因果"),
    is_negative = correlation_r < 0 & correlation_p < 0.05
  ) %>%
  summarise(
    n_total = n(),
    n_temp_to_sif = sum(has_temp_to_sif),
    n_negative_effect = sum(has_temp_to_sif & is_negative),
    n_positive_effect = sum(has_temp_to_sif & !is_negative & correlation_p < 0.05),
    pct_temp_to_sif = round(n_temp_to_sif / n_total * 100, 1),
    pct_negative = round(n_negative_effect / n_temp_to_sif * 100, 1)
  )

cat("\n=== 综合统计 ===\n")
cat("总站点数:", detailed_stats$n_total, "\n")
cat("有温度→SIF因果的站点:", detailed_stats$n_temp_to_sif, 
    "(", detailed_stats$pct_temp_to_sif, "%)\n")
cat("  其中负向影响:", detailed_stats$n_negative_effect,
    "(", detailed_stats$pct_negative, "%)\n")
cat("  其中正向影响:", detailed_stats$n_positive_effect, "\n\n")

# ============================================================================
# 4. 因果方向 × 相关性质的交叉统计
# ============================================================================

cat("\n=== 因果方向 × 相关性质 交叉统计 ===\n\n")

detailed_classification <- ccm_results_all %>%
  mutate(
    effect_type = case_when(
      correlation_p >= 0.05 ~ "无显著相关",
      correlation_r < -0.1 ~ "显著负相关",
      correlation_r > 0.1 ~ "显著正相关",
      TRUE ~ "弱相关"
    ),
    causality_type = factor(causality_type, 
                            levels = c("Temp → SIF", "SIF → Temp", 
                                       "双向因果", "无显著因果"))
  )

cross_table <- detailed_classification %>%
  count(causality_type, effect_type) %>%
  pivot_wider(names_from = effect_type, 
              values_from = n, 
              values_fill = 0) %>%
  arrange(causality_type)

cat("【交叉统计表】\n")
print(cross_table)

cross_table_pct <- detailed_classification %>%
  count(causality_type, effect_type) %>%
  group_by(causality_type) %>%
  mutate(
    pct = round(n / sum(n) * 100, 1),
    label = paste0(n, " (", pct, "%)")
  ) %>%
  ungroup() %>%
  select(-n, -pct) %>%
  pivot_wider(names_from = effect_type, 
              values_from = label,
              values_fill = "0 (0%)")

cat("\n【交叉统计表（含百分比）】\n")
print(cross_table_pct)

cat("\n【详细统计】\n\n")

for (cause_type in c("Temp → SIF", "SIF → Temp", "双向因果", "无显著因果")) {
  
  subset_data <- detailed_classification %>%
    filter(causality_type == cause_type)
  
  if (nrow(subset_data) == 0) next
  
  cat("★", cause_type, "★\n")
  cat("总数:", nrow(subset_data), "个站点\n")
  
  stats <- subset_data %>%
    summarise(
      n_neg_sig = sum(correlation_r < 0 & correlation_p < 0.05),
      n_pos_sig = sum(correlation_r > 0 & correlation_p < 0.05),
      n_no_sig = sum(correlation_p >= 0.05),
      mean_r_neg = mean(correlation_r[correlation_r < 0 & correlation_p < 0.05], 
                        na.rm = TRUE),
      mean_r_pos = mean(correlation_r[correlation_r > 0 & correlation_p < 0.05], 
                        na.rm = TRUE)
    )
  
  cat("  • 显著负相关:", stats$n_neg_sig, 
      "(", round(stats$n_neg_sig/nrow(subset_data)*100, 1), "%)",
      ifelse(is.na(stats$mean_r_neg), "", 
             paste0(" - 平均r=", round(stats$mean_r_neg, 3))), "\n")
  
  cat("  • 显著正相关:", stats$n_pos_sig, 
      "(", round(stats$n_pos_sig/nrow(subset_data)*100, 1), "%)",
      ifelse(is.na(stats$mean_r_pos), "", 
             paste0(" - 平均r=", round(stats$mean_r_pos, 3))), "\n")
  
  cat("  • 无显著相关:", stats$n_no_sig, 
      "(", round(stats$n_no_sig/nrow(subset_data)*100, 1), "%)\n")
  
  cat("\n")
}

# ============================================================================
# 5. 温度→SIF因果关系的空间分布可视化
# ============================================================================

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

cat("=== 准备空间分析 ===\n")

# 读取站点坐标
station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

cat("读取到", nrow(station_coords), "个站点的坐标\n")

# 合并CCM结果和坐标
temp_to_sif_stations <- ccm_results_all %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    effect_category = case_when(
      correlation_r < 0 & correlation_p < 0.05 ~ "显著负相关",
      correlation_r > 0 & correlation_p < 0.05 ~ "显著正相关",
      TRUE ~ "无显著相关"
    ),
    effect_simple = case_when(
      correlation_r < 0 & correlation_p < 0.05 ~ "负向影响\n(温度↑ SIF↓)",
      correlation_r > 0 & correlation_p < 0.05 ~ "正向影响\n(温度↑ SIF↑)",
      TRUE ~ "无显著相关"
    ),
    effect_strength = abs(correlation_r)
  )

cat("\n有温度→SIF因果关系且有坐标的站点:", nrow(temp_to_sif_stations), "\n")
cat("  • 负向影响:", sum(temp_to_sif_stations$effect_category == "显著负相关"), "\n")
cat("  • 正向影响:", sum(temp_to_sif_stations$effect_category == "显著正相关"), "\n")
cat("  • 无显著相关:", sum(temp_to_sif_stations$effect_category == "无显著相关"), "\n\n")

library(rnaturalearth)
# 获取中国地图
china_map <- ne_countries(country = "china", scale = "medium", 
                          returnclass = "sf")

# 地图范围
bbox <- st_bbox(c(
  xmin = min(temp_to_sif_stations$lon) - 2,
  xmax = max(temp_to_sif_stations$lon) + 2,
  ymin = min(temp_to_sif_stations$lat) - 2,
  ymax = max(temp_to_sif_stations$lat) + 2
))

# 主地图
ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = temp_to_sif_stations,
             aes(x = longitude, y = latitude, 
                 color = effect_simple,
                 size = effect_strength),
             alpha = 0.7) +
  scale_color_manual(
    values = c(
      "负向影响\n(温度↑ SIF↓)" = "#d73027",
      "正向影响\n(温度↑ SIF↑)" = "#4575b4",
      "无显著相关" = "#999999"
    ),
    name = "效应方向"
  ) +
  scale_size_continuous(
    range = c(1, 3),
    name = "相关强度\n|r|",
    breaks = c(0.2, 0.4, 0.6)
  ) +
  labs(
    title = "热月温度异常对SIF的因果效应空间分布",
    subtitle = paste0("有因果关系的站点 (n=", nrow(temp_to_sif_stations), 
                      ") | Tp=1 (滞后1月)"),
    x = "经度",
    y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    panel.grid = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

