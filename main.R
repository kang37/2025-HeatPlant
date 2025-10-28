pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, tidyr)

# 读取数据
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop") %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 读取气象数据
read_one_meteo <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1) %>%
    return()
}

meteo_file_lin <- list.files(
  path = "data_raw/meteo_data_1961-2023",
  pattern = "\\.txt$",
  full.names = TRUE
) %>%
  .[!grepl("sta_lonlat_china.txt", .)]

var_meteo_lin <- c("tmax", "tmin", "tavg", "rh", "precip")

meteo_data_lin <- map(
  meteo_file_lin[
    str_extract(basename(meteo_file_lin), "\\d+") %in%
      as.character(meteo_sif_data$meteo_stat_id)
  ],
  read_one_meteo
) %>%
  list_rbind() %>%
  rename_with(~ tolower(.x)) %>% 
  mutate(year = year(date), month = month(date)) %>% 
  group_by(meteo_stat_id, year, month) %>% 
  summarise(
    across(all_of(var_meteo_lin), ~mean(.x, na.rm = TRUE)), .groups = "drop"
  )

# 合并数据
meteo_data <- meteo_sif_data %>%
  inner_join(meteo_data_lin, by = c("meteo_stat_id", "year", "month")) %>% 
  drop_na() %>% 
  arrange(meteo_stat_id, year, month)

# 筛选6-8月数据
meteo_data_summer <- meteo_data %>%
  filter(month %in% 6:8) %>%
  group_by(meteo_stat_id, month) %>%
  mutate(
    temp_q90 = quantile(tavg, 0.9, na.rm = TRUE),
    temp_anomaly = tavg - temp_q90,
    is_hot_month = tavg > temp_q90
  ) %>%
  ungroup() %>%
  arrange(meteo_stat_id, year, month)

# 检查趋势
trend_analysis <- meteo_data_summer %>%
  mutate(year_index = year - min(year)) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    n = n(),
    n_complete = sum(complete.cases(year_index, temp_anomaly, sif)),
    cor_temp_anomaly = cor(year_index, temp_anomaly, use = "complete.obs"),
    cor_sif = cor(year_index, sif, use = "complete.obs"),
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
    need_detrend_temp = abs(cor_temp_anomaly) > 0.3 & p_temp < 0.05 & !is.na(p_temp),
    need_detrend_sif = abs(cor_sif) > 0.3 & p_sif < 0.05 & !is.na(p_sif)
  )

# 线性去趋势
meteo_data_detrended_linear <- meteo_data_summer %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    temp_anomaly_detrended = residuals(lm(temp_anomaly ~ time_idx)),
    sif_detrended = residuals(lm(sif ~ time_idx))
  ) %>%
  ungroup()

data_for_ccm_linear <- meteo_data_detrended_linear %>%
  select(meteo_stat_id, year, month, 
         temp_anomaly_original = temp_anomaly,
         temp_anomaly = temp_anomaly_detrended,
         sif_original = sif,
         sif = sif_detrended,
         is_hot_month)

# ============================================================================
# 差分去趋势完整分析
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("               差分去趋势 vs 线性去趋势 完整对比分析\n")
cat(rep("=", 70), "\n\n", sep = "")

# 创建差分去趋势数据
meteo_data_detrended_diff <- meteo_data_summer %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    temp_anomaly_diff = temp_anomaly - lag(temp_anomaly),
    sif_diff = sif - lag(sif)
  ) %>%
  ungroup() %>%
  filter(!is.na(temp_anomaly_diff), !is.na(sif_diff))

data_for_ccm_diff <- meteo_data_detrended_diff %>%
  select(meteo_stat_id, year, month, 
         temp_anomaly = temp_anomaly_diff,
         sif = sif_diff,
         is_hot_month)

cat("线性去趋势数据点数:", nrow(data_for_ccm_linear), "\n")
cat("差分去趋势数据点数:", nrow(data_for_ccm_diff), 
    "(损失", nrow(data_for_ccm_linear) - nrow(data_for_ccm_diff), "个数据点)\n\n")

# ============================================================================
# 1. 单站点详细对比
# ============================================================================

example_station <- data_for_ccm_linear %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  slice(1) %>%
  pull(meteo_stat_id)

cat("=== 示例站点:", example_station, "===\n\n")

# 提取数据
station_linear <- meteo_data_detrended_linear %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-")))

station_diff <- meteo_data_detrended_diff %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-")))

station_original <- meteo_data_summer %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  mutate(date = as.Date(paste(year, month, 1, sep = "-")))

# 可视化对比
cat("【1. 温度异常的不同去趋势方法对比】\n")

p_temp_compare <- ggplot() +
  geom_line(data = station_original,
            aes(x = date, y = temp_anomaly, color = "原始数据"),
            alpha = 0.5, linewidth = 0.8) +
  geom_line(data = station_linear,
            aes(x = date, y = temp_anomaly_detrended, color = "线性去趋势"),
            linewidth = 0.8) +
  scale_color_manual(values = c(
    "原始数据" = "gray50",
    "线性去趋势" = "blue"
  )) +
  labs(
    title = paste("站点", example_station, "- 温度异常去趋势对比"),
    subtitle = "线性去趋势 vs 原始数据",
    x = "时间",
    y = "温度异常 (°C)",
    color = "方法"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_temp_compare)

p_temp_diff <- ggplot(station_diff, aes(x = date, y = temp_anomaly_diff)) +
  geom_line(color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = paste("站点", example_station, "- 差分去趋势结果"),
    subtitle = "逐月变化量（Δ温度异常）",
    x = "时间",
    y = "温度异常变化量 (°C)",
    caption = "差分后表示的是变化量"
  ) +
  theme_minimal()

print(p_temp_diff)

cat("\n【2. SIF的不同去趋势方法对比】\n")

p_sif_compare <- ggplot() +
  geom_line(data = station_original,
            aes(x = date, y = sif, color = "原始数据"),
            alpha = 0.5, linewidth = 0.8) +
  geom_line(data = station_linear,
            aes(x = date, y = sif_detrended, color = "线性去趋势"),
            linewidth = 0.8) +
  scale_color_manual(values = c(
    "原始数据" = "gray50",
    "线性去趋势" = "darkgreen"
  )) +
  labs(
    title = paste("站点", example_station, "- SIF去趋势对比"),
    subtitle = "线性去趋势 vs 原始数据",
    x = "时间",
    y = "SIF",
    color = "方法"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_sif_compare)

p_sif_diff <- ggplot(station_diff, aes(x = date, y = sif_diff)) +
  geom_line(color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = paste("站点", example_station, "- 差分去趋势结果"),
    subtitle = "逐月变化量（ΔSIF）",
    x = "时间",
    y = "SIF变化量"
  ) +
  theme_minimal()

print(p_sif_diff)

# 统计特征
cat("\n【3. 统计特征对比】\n\n")

stats_compare <- tibble(
  方法 = c("原始数据", "线性去趋势", "差分去趋势"),
  温度_均值 = c(
    mean(station_original$temp_anomaly, na.rm = TRUE),
    mean(station_linear$temp_anomaly_detrended, na.rm = TRUE),
    mean(station_diff$temp_anomaly_diff, na.rm = TRUE)
  ),
  温度_标准差 = c(
    sd(station_original$temp_anomaly, na.rm = TRUE),
    sd(station_linear$temp_anomaly_detrended, na.rm = TRUE),
    sd(station_diff$temp_anomaly_diff, na.rm = TRUE)
  ),
  SIF_均值 = c(
    mean(station_original$sif, na.rm = TRUE),
    mean(station_linear$sif_detrended, na.rm = TRUE),
    mean(station_diff$sif_diff, na.rm = TRUE)
  ),
  SIF_标准差 = c(
    sd(station_original$sif, na.rm = TRUE),
    sd(station_linear$sif_detrended, na.rm = TRUE),
    sd(station_diff$sif_diff, na.rm = TRUE)
  )
) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

print(stats_compare)

# 相关性对比
cat("\n【4. 温度-SIF相关性对比】\n\n")

station_linear_lag <- station_linear %>%
  mutate(sif_lag1 = lead(sif_detrended, 1))

cor_linear <- cor.test(station_linear_lag$temp_anomaly_detrended,
                       station_linear_lag$sif_lag1)

station_diff_lag <- station_diff %>%
  mutate(sif_lag1 = lead(sif_diff, 1))

cor_diff <- cor.test(station_diff_lag$temp_anomaly_diff,
                     station_diff_lag$sif_lag1)

cat("线性去趋势:\n")
cat("  r =", round(cor_linear$estimate, 3), "\n")
cat("  p =", format.pval(cor_linear$p.value), "\n\n")

cat("差分去趋势:\n")
cat("  r =", round(cor_diff$estimate, 3), "\n")
cat("  p =", format.pval(cor_diff$p.value), "\n\n")

# 散点图
p_cor_linear <- ggplot(station_linear_lag, 
                       aes(x = temp_anomaly_detrended, y = sif_lag1)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "线性去趋势：t月温度 vs t+1月SIF",
    subtitle = paste0("r = ", round(cor_linear$estimate, 3), 
                      ", p = ", format.pval(cor_linear$p.value, digits = 2)),
    x = "t月温度异常",
    y = "t+1月SIF"
  ) +
  theme_minimal()

print(p_cor_linear)

p_cor_diff <- ggplot(station_diff_lag,
                     aes(x = temp_anomaly_diff, y = sif_lag1)) +
  geom_point(alpha = 0.5, color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "差分去趋势：t月温度变化 vs t+1月SIF变化",
    subtitle = paste0("r = ", round(cor_diff$estimate, 3), 
                      ", p = ", format.pval(cor_diff$p.value, digits = 2)),
    x = "t月温度异常变化量",
    y = "t+1月SIF变化量"
  ) +
  theme_minimal()

print(p_cor_diff)

# ============================================================================
# 2. 完整CCM分析对比（单站点）
# ============================================================================

cat("\n【5. 完整CCM分析对比】\n\n")

library(rEDM)

# 准备线性去趋势的CCM数据
ccm_data_linear <- data.frame(
  time = 1:nrow(station_linear),
  temp_anomaly = station_linear$temp_anomaly_detrended,
  sif = station_linear$sif_detrended
)

# 准备差分去趋势的CCM数据
ccm_data_diff <- data.frame(
  time = 1:nrow(station_diff),
  temp_anomaly = station_diff$temp_anomaly_diff,
  sif = station_diff$sif_diff
)

cat("--- 线性去趋势 CCM 分析 ---\n")

# 确定最优E（线性）
n_linear <- nrow(ccm_data_linear)
lib_pred_linear <- min(100, n_linear)

embed_temp_linear <- EmbedDimension(
  dataFrame = ccm_data_linear,
  lib = paste("1", lib_pred_linear),
  pred = paste("1", lib_pred_linear),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  maxE = 10,
  showPlot = FALSE
)

embed_sif_linear <- EmbedDimension(
  dataFrame = ccm_data_linear,
  lib = paste("1", lib_pred_linear),
  pred = paste("1", lib_pred_linear),
  columns = "sif",
  target = "sif",
  maxE = 10,
  showPlot = FALSE
)

best_E_temp_linear <- embed_temp_linear$E[which.max(embed_temp_linear$rho)]
best_E_sif_linear <- embed_sif_linear$E[which.max(embed_sif_linear$rho)]
E_ccm_linear <- max(best_E_temp_linear, best_E_sif_linear)

cat("最优 E =", E_ccm_linear, "\n")

# CCM分析（线性）
tau <- 1
embedding_loss_linear <- (E_ccm_linear - 1) * tau
tp <- 1
max_lib_linear <- n_linear - embedding_loss_linear - tp

lib_start_linear <- max(E_ccm_linear + 2, 10)
lib_end_linear <- max_lib_linear
lib_step_linear <- max(2, floor((lib_end_linear - lib_start_linear) / 15))
libSizes_str_linear <- paste(lib_start_linear, lib_end_linear, lib_step_linear)

ccm_linear_result <- CCM(
  dataFrame = ccm_data_linear,
  E = E_ccm_linear,
  Tp = 1,
  columns = "temp_anomaly",
  target = "sif",
  libSizes = libSizes_str_linear,
  sample = 100,
  random = TRUE,
  showPlot = FALSE
)

ccm_linear_summary <- ccm_linear_result %>%
  rename(lib_size = LibSize) %>%
  group_by(lib_size) %>%
  summarise(
    rho_temp_to_sif = mean(`temp_anomaly:sif`, na.rm = TRUE),
    rho_sif_to_temp = mean(`sif:temp_anomaly`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method = "线性去趋势")

cat("最终 ρ (Temp→SIF) =", 
    round(ccm_linear_summary$rho_temp_to_sif[nrow(ccm_linear_summary)], 3), "\n")
cat("最终 ρ (SIF→Temp) =", 
    round(ccm_linear_summary$rho_sif_to_temp[nrow(ccm_linear_summary)], 3), "\n\n")

cat("--- 差分去趋势 CCM 分析 ---\n")

# 确定最优E（差分）
n_diff <- nrow(ccm_data_diff)
lib_pred_diff <- min(100, n_diff)

embed_temp_diff <- EmbedDimension(
  dataFrame = ccm_data_diff,
  lib = paste("1", lib_pred_diff),
  pred = paste("1", lib_pred_diff),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  maxE = 10,
  showPlot = FALSE
)

embed_sif_diff <- EmbedDimension(
  dataFrame = ccm_data_diff,
  lib = paste("1", lib_pred_diff),
  pred = paste("1", lib_pred_diff),
  columns = "sif",
  target = "sif",
  maxE = 10,
  showPlot = FALSE
)

best_E_temp_diff <- embed_temp_diff$E[which.max(embed_temp_diff$rho)]
best_E_sif_diff <- embed_sif_diff$E[which.max(embed_sif_diff$rho)]
E_ccm_diff <- max(best_E_temp_diff, best_E_sif_diff)

cat("最优 E =", E_ccm_diff, "\n")

# CCM分析（差分）
embedding_loss_diff <- (E_ccm_diff - 1) * tau
max_lib_diff <- n_diff - embedding_loss_diff - tp

lib_start_diff <- max(E_ccm_diff + 2, 10)
lib_end_diff <- max_lib_diff
lib_step_diff <- max(2, floor((lib_end_diff - lib_start_diff) / 15))
libSizes_str_diff <- paste(lib_start_diff, lib_end_diff, lib_step_diff)

ccm_diff_result <- CCM(
  dataFrame = ccm_data_diff,
  E = E_ccm_diff,
  Tp = 1,
  columns = "temp_anomaly",
  target = "sif",
  libSizes = libSizes_str_diff,
  sample = 100,
  random = TRUE,
  showPlot = FALSE
)

ccm_diff_summary <- ccm_diff_result %>%
  rename(lib_size = LibSize) %>%
  group_by(lib_size) %>%
  summarise(
    rho_temp_to_sif = mean(`temp_anomaly:sif`, na.rm = TRUE),
    rho_sif_to_temp = mean(`sif:temp_anomaly`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method = "差分去趋势")

cat("最终 ρ (Temp→SIF) =", 
    round(ccm_diff_summary$rho_temp_to_sif[nrow(ccm_diff_summary)], 3), "\n")
cat("最终 ρ (SIF→Temp) =", 
    round(ccm_diff_summary$rho_sif_to_temp[nrow(ccm_diff_summary)], 3), "\n\n")

# CCM对比图
ccm_combined <- bind_rows(ccm_linear_summary, ccm_diff_summary)

p_ccm_compare <- ggplot(ccm_combined, 
                        aes(x = lib_size, y = rho_temp_to_sif, 
                            color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", 
             color = "gray", alpha = 0.5) +
  scale_color_manual(values = c("线性去趋势" = "blue", 
                                "差分去趋势" = "red")) +
  labs(
    title = paste("CCM对比: 温度→SIF - 站点", example_station),
    subtitle = "不同去趋势方法的因果强度",
    x = "Library Size",
    y = "Cross Map Skill (ρ)",
    color = "去趋势方法",
    linetype = "去趋势方法"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_ccm_compare)

# ============================================================================
# 3. 批量分析：线性去趋势（主分析）
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("                开始批量CCM分析（线性去趋势）\n")
cat(rep("=", 70), "\n\n", sep = "")

# CCM分析函数
perform_ccm_with_correlation <- function(station_id, data, min_points = 30) {
  
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif) %>%
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    station_data_with_lag <- station_data %>%
      mutate(sif_lag1 = lead(sif, 1))
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif
    )
    
    lib_pred_size <- min(100, n_data)
    
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
    
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    tp <- 1
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) {
      return(NULL)
    }
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    ccm_result <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = 1,
      columns = "temp_anomaly",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
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
    
    final_rho_temp_to_sif <- ccm_summary$rho_temp_to_sif_mean[nrow(ccm_summary)]
    final_rho_sif_to_temp <- ccm_summary$rho_sif_to_temp_mean[nrow(ccm_summary)]
    
    trend_temp_to_sif <- cor(ccm_summary$lib_size, 
                             ccm_summary$rho_temp_to_sif_mean)
    trend_sif_to_temp <- cor(ccm_summary$lib_size, 
                             ccm_summary$rho_sif_to_temp_mean)
    
    cor_result <- cor.test(station_data_with_lag$temp_anomaly, 
                           station_data_with_lag$sif_lag1)
    
    correlation_r <- cor_result$estimate
    correlation_p <- cor_result$p.value
    
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
    
    effect_direction <- case_when(
      correlation_p >= 0.05 ~ "不显著",
      correlation_r < 0 ~ "负向(温度↑ SIF↓)",
      correlation_r > 0 ~ "正向(温度↑ SIF↑)",
      TRUE ~ "无关"
    )
    
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

# 批量分析
all_stations <- data_for_ccm_linear %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("\n准备分析", length(all_stations), "个站点（线性去趋势）...\n\n")

ccm_results_linear <- map_dfr(seq_along(all_stations), function(i) {
  if (i %% 50 == 0) {
    cat("已完成:", i, "/", length(all_stations), "\n")
  }
  perform_ccm_with_correlation(all_stations[i], data_for_ccm_linear, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_linear), "/", length(all_stations), "\n")

# ============================================================================
# 4. 批量分析：差分去趋势（对比分析）
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("                开始批量CCM分析（差分去趋势）\n")
cat(rep("=", 70), "\n\n", sep = "")

all_stations_diff <- data_for_ccm_diff %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("准备分析", length(all_stations_diff), "个站点（差分去趋势）...\n\n")

ccm_results_diff <- map_dfr(seq_along(all_stations_diff), function(i) {
  if (i %% 50 == 0) {
    cat("已完成:", i, "/", length(all_stations_diff), "\n")
  }
  perform_ccm_with_correlation(all_stations_diff[i], data_for_ccm_diff, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_diff), "/", length(all_stations_diff), "\n")

# ============================================================================
# 5. 两种方法的批量结果对比
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("              两种去趋势方法的批量结果对比\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("=== 线性去趋势结果 ===\n\n")

causality_stats_linear <- table(ccm_results_linear$causality_type)
print(causality_stats_linear)

cat("\n各因果类型占比:\n")
print(round(prop.table(causality_stats_linear) * 100, 1))

effect_stats_linear <- ccm_results_linear %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  count(effect_direction) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\n温度→SIF的效应方向统计:\n")
print(effect_stats_linear)

cat("\n=== 差分去趋势结果 ===\n\n")

causality_stats_diff <- table(ccm_results_diff$causality_type)
print(causality_stats_diff)

cat("\n各因果类型占比:\n")
print(round(prop.table(causality_stats_diff) * 100, 1))

effect_stats_diff <- ccm_results_diff %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  count(effect_direction) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\n温度→SIF的效应方向统计:\n")
print(effect_stats_diff)

# 对比汇总
cat("\n=== 综合对比 ===\n\n")

comparison_summary <- tibble(
  去趋势方法 = c("线性去趋势", "差分去趋势"),
  总站点数 = c(nrow(ccm_results_linear), nrow(ccm_results_diff)),
  有Temp到SIF因果 = c(
    sum(ccm_results_linear$causality_type %in% c("Temp → SIF", "双向因果")),
    sum(ccm_results_diff$causality_type %in% c("Temp → SIF", "双向因果"))
  ),
  其中负向影响 = c(
    sum(ccm_results_linear$causality_type %in% c("Temp → SIF", "双向因果") & 
          ccm_results_linear$correlation_r < 0 & 
          ccm_results_linear$correlation_p < 0.05),
    sum(ccm_results_diff$causality_type %in% c("Temp → SIF", "双向因果") & 
          ccm_results_diff$correlation_r < 0 & 
          ccm_results_diff$correlation_p < 0.05)
  ),
  其中正向影响 = c(
    sum(ccm_results_linear$causality_type %in% c("Temp → SIF", "双向因果") & 
          ccm_results_linear$correlation_r > 0 & 
          ccm_results_linear$correlation_p < 0.05),
    sum(ccm_results_diff$causality_type %in% c("Temp → SIF", "双向因果") & 
          ccm_results_diff$correlation_r > 0 & 
          ccm_results_diff$correlation_p < 0.05)
  ),
  平均CCM强度 = c(
    mean(ccm_results_linear$rho_temp_to_sif[
      ccm_results_linear$causality_type %in% c("Temp → SIF", "双向因果")
    ], na.rm = TRUE),
    mean(ccm_results_diff$rho_temp_to_sif[
      ccm_results_diff$causality_type %in% c("Temp → SIF", "双向因果")
    ], na.rm = TRUE)
  ),
  平均相关系数 = c(
    mean(ccm_results_linear$correlation_r[
      ccm_results_linear$causality_type %in% c("Temp → SIF", "双向因果")
    ], na.rm = TRUE),
    mean(ccm_results_diff$correlation_r[
      ccm_results_diff$causality_type %in% c("Temp → SIF", "双向因果")
    ], na.rm = TRUE)
  )
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(comparison_summary)

# 可视化对比
library(patchwork)

# 因果类型分布对比
causality_compare_data <- bind_rows(
  ccm_results_linear %>% 
    count(causality_type) %>% 
    mutate(method = "线性去趋势"),
  ccm_results_diff %>% 
    count(causality_type) %>% 
    mutate(method = "差分去趋势")
)

p_causality_compare <- ggplot(causality_compare_data, 
                              aes(x = causality_type, y = n, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("线性去趋势" = "#377EB8", 
                               "差分去趋势" = "#E41A1C")) +
  labs(
    title = "因果类型分布对比",
    x = "因果类型",
    y = "站点数",
    fill = "去趋势方法"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

print(p_causality_compare)

# CCM强度分布对比
ccm_strength_compare <- bind_rows(
  ccm_results_linear %>% 
    filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
    select(meteo_stat_id, rho_temp_to_sif) %>%
    mutate(method = "线性去趋势"),
  ccm_results_diff %>% 
    filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
    select(meteo_stat_id, rho_temp_to_sif) %>%
    mutate(method = "差分去趋势")
)

p_rho_compare <- ggplot(ccm_strength_compare, 
                        aes(x = rho_temp_to_sif, fill = method)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("线性去趋势" = "#377EB8", 
                               "差分去趋势" = "#E41A1C")) +
  labs(
    title = "CCM因果强度分布对比 (Temp→SIF)",
    x = "ρ (Cross Map Skill)",
    y = "站点数",
    fill = "去趋势方法"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_rho_compare)

# 相关系数分布对比
cor_compare <- bind_rows(
  ccm_results_linear %>% 
    filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
    select(meteo_stat_id, correlation_r) %>%
    mutate(method = "线性去趋势"),
  ccm_results_diff %>% 
    filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
    select(meteo_stat_id, correlation_r) %>%
    mutate(method = "差分去趋势")
)

p_cor_compare_hist <- ggplot(cor_compare, 
                             aes(x = correlation_r, fill = method)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("线性去趋势" = "#377EB8", 
                               "差分去趋势" = "#E41A1C")) +
  labs(
    title = "相关系数分布对比 (Temp→SIF)",
    x = "Pearson r",
    y = "站点数",
    fill = "去趋势方法"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_cor_compare_hist)

# 散点图：两种方法的结果对比（共同站点）
common_stations <- intersect(ccm_results_linear$meteo_stat_id,
                             ccm_results_diff$meteo_stat_id)

if (length(common_stations) > 0) {
  
  cat("\n共同分析的站点数:", length(common_stations), "\n")
  
  comparison_data <- ccm_results_linear %>%
    filter(meteo_stat_id %in% common_stations) %>%
    select(meteo_stat_id, rho_linear = rho_temp_to_sif, 
           r_linear = correlation_r) %>%
    left_join(
      ccm_results_diff %>%
        filter(meteo_stat_id %in% common_stations) %>%
        select(meteo_stat_id, rho_diff = rho_temp_to_sif, 
               r_diff = correlation_r),
      by = "meteo_stat_id"
    )
  
  # CCM强度对比
  p_rho_scatter <- ggplot(comparison_data, 
                          aes(x = rho_linear, y = rho_diff)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(
      title = "CCM强度对比：线性 vs 差分",
      x = "线性去趋势 ρ",
      y = "差分去趋势 ρ",
      caption = paste("相关系数:", 
                      round(cor(comparison_data$rho_linear, 
                                comparison_data$rho_diff, 
                                use = "complete.obs"), 3))
    ) +
    theme_minimal()
  
  print(p_rho_scatter)
  
  # 相关系数对比
  p_r_scatter <- ggplot(comparison_data, 
                        aes(x = r_linear, y = r_diff)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(
      title = "相关系数对比：线性 vs 差分",
      x = "线性去趋势 r",
      y = "差分去趋势 r",
      caption = paste("相关系数:", 
                      round(cor(comparison_data$r_linear, 
                                comparison_data$r_diff, 
                                use = "complete.obs"), 3))
    ) +
    theme_minimal()
  
  print(p_r_scatter)
}

# ============================================================================
# 6. 继续使用线性去趋势进行后续分析
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("            继续使用线性去趋势数据进行后续分析\n")
cat(rep("=", 70), "\n\n", sep = "")

# 使用线性去趋势的结果
ccm_results_all <- ccm_results_linear

# 因果方向 × 相关性质的交叉统计
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
# 7. 空间分布分析
# ============================================================================

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

cat("=== 准备空间分析 ===\n")

station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

cat("读取到", nrow(station_coords), "个站点的坐标\n")

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

china_map <- ne_countries(country = "china", scale = "medium", 
                          returnclass = "sf")

bbox <- st_bbox(c(
  xmin = min(temp_to_sif_stations$longitude) - 2,
  xmax = max(temp_to_sif_stations$longitude) + 2,
  ymin = min(temp_to_sif_stations$lat) - 2,
  ymax = max(temp_to_sif_stations$lat) + 2
))

p_map_main <- ggplot() +
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
                      ") | Tp=1 (滞后1月) | 线性去趋势"),
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

print(p_map_main)

temp_to_sif_sig <- temp_to_sif_stations %>%
  filter(effect_category != "无显著相关")

if (nrow(temp_to_sif_sig) > 0) {
  p_map_facet <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(data = temp_to_sif_sig,
               aes(x = longitude, y = latitude, 
                   color = effect_category,
                   size = effect_strength),
               alpha = 0.7) +
    scale_color_manual(
      values = c(
        "显著负相关" = "#d73027",
        "显著正相关" = "#4575b4"
      ),
      name = "效应类型"
    ) +
    scale_size_continuous(
      range = c(2, 6),
      name = "|r|"
    ) +
    facet_wrap(~ effect_category, ncol = 2) +
    labs(
      title = "正向与负向影响的空间分布对比",
      subtitle = paste0("显著相关的站点 (n=", nrow(temp_to_sif_sig), ")"),
      x = "经度",
      y = "纬度"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.background = element_rect(fill = "aliceblue", color = NA)
    )
  
  print(p_map_facet)
}

temp_to_sif_stations_region <- temp_to_sif_stations %>%
  mutate(
    region = case_when(
      longitude < 105 ~ "西部",
      longitude >= 105 & longitude < 115 ~ "中部",
      TRUE ~ "东部"
    ),
    lat_zone = case_when(
      lat < 30 ~ "南方",
      lat >= 30 & lat < 40 ~ "中纬度",
      TRUE ~ "北方"
    )
  )


# ============================================================================
# 空间分布对比：线性 vs 差分去趋势
# ============================================================================

cat("\n=== 空间分布对比：两种去趋势方法 ===\n\n")

# 线性去趋势的空间数据
temp_to_sif_linear <- ccm_results_linear %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    effect_category = case_when(
      correlation_r < 0 & correlation_p < 0.05 ~ "显著负相关",
      correlation_r > 0 & correlation_p < 0.05 ~ "显著正相关",
      TRUE ~ "无显著相关"
    ),
    method = "线性去趋势"
  )

# 差分去趋势的空间数据
temp_to_sif_diff <- ccm_results_diff %>%
  filter(causality_type %in% c("Temp → SIF", "双向因果")) %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    effect_category = case_when(
      correlation_r < 0 & correlation_p < 0.05 ~ "显著负相关",
      correlation_r > 0 & correlation_p < 0.05 ~ "显著正相关",
      TRUE ~ "无显著相关"
    ),
    method = "差分去趋势"
  )

cat("线性去趋势 - 空间站点数:", nrow(temp_to_sif_linear), "\n")
cat("  负向:", sum(temp_to_sif_linear$effect_category == "显著负相关"), "\n")
cat("  正向:", sum(temp_to_sif_linear$effect_category == "显著正相关"), "\n\n")

cat("差分去趋势 - 空间站点数:", nrow(temp_to_sif_diff), "\n")
cat("  负向:", sum(temp_to_sif_diff$effect_category == "显著负相关"), "\n")
cat("  正向:", sum(temp_to_sif_diff$effect_category == "显著正相关"), "\n\n")

# 合并两种方法的数据
temp_to_sif_both <- bind_rows(temp_to_sif_linear, temp_to_sif_diff)

# 并排对比地图
p_spatial_compare <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = temp_to_sif_both,
             aes(x = longitude, y = latitude, 
                 color = effect_category,
                 size = abs(correlation_r)),
             alpha = 0.6) +
  scale_color_manual(
    values = c(
      "显著负相关" = "#d73027",
      "显著正相关" = "#4575b4",
      "无显著相关" = "#999999"
    ),
    name = "效应方向"
  ) +
  scale_size_continuous(range = c(1, 3), name = "|r|") +
  facet_wrap(~ method, ncol = 2) +
  labs(
    title = "空间分布对比：线性去趋势 vs 差分去趋势",
    subtitle = "温度→SIF 因果关系的空间格局",
    x = "经度", y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )
print(p_spatial_compare)
