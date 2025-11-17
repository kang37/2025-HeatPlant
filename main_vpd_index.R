# ============================================================================
# 基于每日VPD的热事件分析
# ============================================================================

pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, 
               tidyr, showtext, rEDM, sf, rnaturalearth, rnaturalearthdata)
showtext_auto()

cat("\n", rep("=", 70), "\n", sep = "")
cat("              VPD热事件分析\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================================
# 1. 读取每日气象数据
# ============================================================================

cat("【1. 读取每日数据】\n")

read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1) %>%
    return()
}

meteo_file_list <- list.files(
  path = "data_raw/meteo_data_1961-2023",
  pattern = "\\.txt$",
  full.names = TRUE
) %>%
  .[!grepl("sta_lonlat_china.txt", .)]

# 读取SIF站点列表
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

target_stations <- unique(meteo_sif_data$meteo_stat_id)

cat("目标站点数:", length(target_stations), "\n")

# 读取每日数据
meteo_data_daily <- map(
  meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations],
  read_one_meteo_daily
) %>%
  list_rbind() %>%
  rename_with(~ tolower(.x)) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    day = day(date)
  ) %>%
  filter(month %in% 5:9)  # 扩展到5-9月

cat("读取完成，数据行数:", nrow(meteo_data_daily), "\n\n")

# ============================================================================
# 2. 计算每日VPD
# ============================================================================

cat("【2. 计算每日VPD】\n")

meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    # 饱和水汽压 (kPa)
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    
    # 实际水汽压 (kPa)
    avp = (rh / 100) * svp,
    
    # VPD (kPa)
    vpd = svp - avp
  ) %>%
  # 数据质量控制
  mutate(
    vpd = case_when(
      vpd < 0 ~ NA_real_,
      vpd > 10 ~ NA_real_,
      is.na(tavg) | is.na(rh) ~ NA_real_,
      TRUE ~ vpd
    )
  ) %>%
  filter(!is.na(vpd))

# VPD统计
vpd_stats <- meteo_data_daily_vpd %>%
  summarise(
    mean = mean(vpd, na.rm = TRUE),
    sd = sd(vpd, na.rm = TRUE),
    min = min(vpd, na.rm = TRUE),
    q25 = quantile(vpd, 0.25, na.rm = TRUE),
    median = quantile(vpd, 0.5, na.rm = TRUE),
    q75 = quantile(vpd, 0.75, na.rm = TRUE),
    q90 = quantile(vpd, 0.90, na.rm = TRUE),
    q95 = quantile(vpd, 0.95, na.rm = TRUE),
    max = max(vpd, na.rm = TRUE)
  )

cat("\nVPD统计信息 (kPa):\n")
print(vpd_stats)

# VPD分布图
p_vpd_dist <- ggplot(meteo_data_daily_vpd, aes(x = vpd)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = vpd_stats$median, 
             linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = vpd_stats$q90, 
             linetype = "dashed", color = "orange", linewidth = 1) +
  geom_vline(xintercept = vpd_stats$q95, 
             linetype = "dashed", color = "darkred", linewidth = 1) +
  annotate("text", x = vpd_stats$median, y = Inf, 
           label = paste0("中位数 = ", round(vpd_stats$median, 2), " kPa"),
           vjust = 2, hjust = -0.1, color = "red") +
  annotate("text", x = vpd_stats$q90, y = Inf, 
           label = paste0("P90 = ", round(vpd_stats$q90, 2), " kPa"),
           vjust = 4, hjust = -0.1, color = "orange") +
  annotate("text", x = vpd_stats$q95, y = Inf, 
           label = paste0("P95 = ", round(vpd_stats$q95, 2), " kPa"),
           vjust = 6, hjust = -0.1, color = "darkred") +
  labs(
    title = "每日VPD分布（5-9月）",
    x = "VPD (kPa)",
    y = "频数"
  ) +
  theme_minimal()

print(p_vpd_dist)

# ============================================================================
# 3. 定义热事件阈值
# ============================================================================

cat("\n【3. 定义热事件阈值】\n")

# 方法1：基于百分位数（推荐）
vpd_threshold_p90 <- vpd_stats$q90  # 90th百分位
vpd_threshold_p95 <- vpd_stats$q95  # 95th百分位

# 方法2：基于文献的绝对阈值
vpd_threshold_abs <- 2.0  # 2.0 kPa是常用的高VPD阈值

# 选择使用的阈值
VPD_THRESHOLD <- vpd_threshold_p90  # 使用90th百分位

cat("使用的VPD热事件阈值:", round(VPD_THRESHOLD, 3), "kPa\n")
cat("  (对应第", 90, "百分位)\n\n")

cat("其他参考阈值:\n")
cat("  P95 =", round(vpd_threshold_p95, 3), "kPa\n")
cat("  绝对阈值 =", vpd_threshold_abs, "kPa\n\n")

# ============================================================================
# 4. 识别每日热事件并构建月度指标
# ============================================================================

cat("【4. 构建月度VPD热事件指标】\n")

# 每日标记热事件
meteo_data_daily_events <- meteo_data_daily_vpd %>%
  mutate(
    # 是否为热事件
    is_heat_event = vpd > VPD_THRESHOLD,
    
    # 热事件强度（超过阈值的部分）
    heat_intensity = pmax(vpd - VPD_THRESHOLD, 0)
  )

# 月度汇总
monthly_heat_metrics <- meteo_data_daily_events %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(
    # 基础统计
    n_days = n(),
    vpd_mean = mean(vpd, na.rm = TRUE),
    vpd_max = max(vpd, na.rm = TRUE),
    vpd_sd = sd(vpd, na.rm = TRUE),
    
    # 热事件频次
    heat_event_days = sum(is_heat_event, na.rm = TRUE),
    heat_event_freq = heat_event_days / n_days,  # 热事件日数比例
    
    # 热事件强度
    heat_intensity_sum = sum(heat_intensity, na.rm = TRUE),  # 累积强度
    heat_intensity_mean = mean(heat_intensity[is_heat_event], na.rm = TRUE),  # 平均强度
    
    # 组合指标1：频次 × 强度
    heat_index_v1 = heat_event_freq * heat_intensity_sum,
    
    # 组合指标2：事件日数 × 平均强度
    heat_index_v2 = heat_event_days * ifelse(is.finite(heat_intensity_mean), 
                                             heat_intensity_mean, 0),
    
    # 组合指标3：加权综合指数
    # HI = α × (标准化频次) + β × (标准化强度)
    .groups = "drop"
  ) %>%
  # 标准化处理（按站点）
  group_by(meteo_stat_id) %>%
  mutate(
    # Z-score标准化
    heat_freq_z = scale(heat_event_freq)[,1],
    heat_intensity_z = scale(heat_intensity_sum)[,1],
    
    # 综合热事件指数（权重可调整）
    heat_index_composite = 0.5 * heat_freq_z + 0.5 * heat_intensity_z
  ) %>%
  ungroup()

cat("月度热事件指标构建完成\n")
cat("站点数:", length(unique(monthly_heat_metrics$meteo_stat_id)), "\n")
cat("数据行数:", nrow(monthly_heat_metrics), "\n\n")

# 热事件统计
heat_event_summary <- monthly_heat_metrics %>%
  summarise(
    总月数 = n(),
    有热事件的月份 = sum(heat_event_days > 0),
    平均热事件天数 = mean(heat_event_days),
    最大热事件天数 = max(heat_event_days),
    平均热事件频率 = mean(heat_event_freq) * 100
  )

cat("热事件月度统计:\n")
print(heat_event_summary)
cat("\n")

# ============================================================================
# 5. 合并SIF数据
# ============================================================================

cat("【5. 合并SIF数据】\n")

# 月度SIF数据
sif_monthly <- meteo_sif_data %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")

# 合并
data_heat_sif <- monthly_heat_metrics %>%
  inner_join(sif_monthly, by = c("meteo_stat_id", "year", "month")) %>%
  filter(!is.na(sif))

cat("合并后数据行数:", nrow(data_heat_sif), "\n")
cat("站点数:", length(unique(data_heat_sif$meteo_stat_id)), "\n\n")

# ============================================================================
# 6. 去趋势处理
# ============================================================================

cat("【6. 去趋势处理】\n")

# 安全的去趋势函数
safe_detrend <- function(x, time_idx) {
  # 检查是否有足够的非NA值
  valid_data <- !is.na(x)
  n_valid <- sum(valid_data)
  
  if (n_valid < 3) {
    # 数据点太少，返回原值或NA
    return(rep(NA_real_, length(x)))
  }
  
  tryCatch({
    # 拟合线性模型
    model <- lm(x ~ time_idx)
    # 返回残差
    residuals(model)
  }, error = function(e) {
    # 如果出错，返回NA
    rep(NA_real_, length(x))
  })
}

data_heat_sif_detrended <- data_heat_sif %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    
    # 安全去趋势
    sif_detrended = safe_detrend(sif, time_idx),
    heat_index_composite_detrended = safe_detrend(heat_index_composite, time_idx),
    heat_event_freq_detrended = safe_detrend(heat_event_freq, time_idx),
    heat_intensity_sum_detrended = safe_detrend(heat_intensity_sum, time_idx),
    vpd_mean_detrended = safe_detrend(vpd_mean, time_idx)
  ) %>%
  ungroup() %>%
  # 过滤掉去趋势失败的站点
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(sif_detrended)) >= 10) %>%  # 至少10个有效数据点
  ungroup()

cat("去趋势完成\n")
cat("保留的站点数:", length(unique(data_heat_sif_detrended$meteo_stat_id)), "\n")
cat("数据行数:", nrow(data_heat_sif_detrended), "\n\n")

# ============================================================================
# 7. CCM分析：热事件指标 → SIF
# ============================================================================

cat("【7. CCM分析：热事件指标 → SIF】\n\n")

perform_ccm_heat_sif <- function(station_id, data, min_points = 30) {
  
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(
      sif = sif_detrended,
      heat_index = heat_index_composite_detrended,
      heat_freq = heat_event_freq_detrended,
      heat_intensity = heat_intensity_sum_detrended,
      vpd_mean = vpd_mean_detrended
    ) %>%
    filter(!is.na(sif), !is.na(heat_index)) %>%  # 移除NA
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      sif = station_data$sif,
      heat_index = station_data$heat_index,
      heat_freq = station_data$heat_freq,
      heat_intensity = station_data$heat_intensity,
      vpd_mean = station_data$vpd_mean
    )
    
    lib_pred_size <- min(100, n_data)
    
    # 确定最优E
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    embed_heat <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "heat_index",
      target = "heat_index",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    best_E_heat <- embed_heat$E[which.max(embed_heat$rho)]
    E_ccm <- max(best_E_sif, best_E_heat)
    
    # CCM参数
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    tp <- 0  # 同步关系
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) {
      return(NULL)
    }
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # CCM: 热事件指数 → SIF
    ccm_heat_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "heat_index",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # CCM: 平均VPD → SIF（对比）
    ccm_vpd_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "vpd_mean",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 汇总
    ccm_summary_heat <- ccm_heat_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`heat_index:sif`, na.rm = TRUE),
        rho_sd = sd(`heat_index:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    ccm_summary_vpd <- ccm_vpd_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`vpd_mean:sif`, na.rm = TRUE),
        rho_sd = sd(`vpd_mean:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    final_rho_heat <- ccm_summary_heat$rho_mean[nrow(ccm_summary_heat)]
    final_rho_vpd <- ccm_summary_vpd$rho_mean[nrow(ccm_summary_vpd)]
    
    trend_heat <- cor(ccm_summary_heat$lib_size, ccm_summary_heat$rho_mean)
    trend_vpd <- cor(ccm_summary_vpd$lib_size, ccm_summary_vpd$rho_mean)
    
    # S-map分析
    smap_heat <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "heat_index",
      target = "sif",
      embedded = FALSE
    )
    
    smap_coeffs <- smap_heat$coefficients
    
    if (!is.null(smap_coeffs) && ncol(smap_coeffs) >= 2) {
      heat_coef_col <- which(grepl("heat", colnames(smap_coeffs), ignore.case = TRUE))
      
      if (length(heat_coef_col) > 0) {
        smap_coef_heat <- smap_coeffs[, heat_coef_col[1]]
      } else {
        smap_coef_heat <- smap_coeffs[, 2]
      }
      
      smap_coef_heat <- smap_coef_heat[!is.na(smap_coef_heat)]
      
      mean_smap_coef <- mean(smap_coef_heat, na.rm = TRUE)
      sd_smap_coef <- sd(smap_coef_heat, na.rm = TRUE)
      
    } else {
      mean_smap_coef <- NA
      sd_smap_coef <- NA
      smap_coef_heat <- NA
    }
    
    # 判断因果
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    heat_causes_sif <- (final_rho_heat > ccm_threshold_rho & 
                          trend_heat > ccm_threshold_trend)
    vpd_causes_sif <- (final_rho_vpd > ccm_threshold_rho & 
                         trend_vpd > ccm_threshold_trend)
    
    effect_type <- case_when(
      !heat_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0.001 ~ "促进效应(+)",
      mean_smap_coef < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    effect_stability <- case_when(
      is.na(mean_smap_coef) | is.na(sd_smap_coef) ~ "未知",
      abs(mean_smap_coef) > sd_smap_coef ~ "稳定",
      TRUE ~ "波动大"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,
      
      # 热事件指数 → SIF
      rho_heat_to_sif = final_rho_heat,
      trend_heat_to_sif = trend_heat,
      heat_causes_sif = heat_causes_sif,
      effect_index_heat = mean_smap_coef,
      effect_index_sd_heat = sd_smap_coef,
      effect_type_heat = effect_type,
      effect_stability_heat = effect_stability,
      
      # VPD均值 → SIF（对比）
      rho_vpd_to_sif = final_rho_vpd,
      trend_vpd_to_sif = trend_vpd,
      vpd_causes_sif = vpd_causes_sif,
      
      # 保存数据
      smap_coefficients_heat = list(smap_coef_heat),
      ccm_summary_heat = list(ccm_summary_heat),
      ccm_summary_vpd = list(ccm_summary_vpd)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# 批量分析
all_stations_heat <- data_heat_sif_detrended %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("准备分析", length(all_stations_heat), "个站点...\n\n")

ccm_results_heat <- map_dfr(seq_along(all_stations_heat), function(i) {
  if (i %% 10 == 0) {
    cat("已完成:", i, "/", length(all_stations_heat), "\n")
  }
  perform_ccm_heat_sif(all_stations_heat[i], data_heat_sif_detrended, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_heat), "/", length(all_stations_heat), "\n\n")

# ============================================================================
# 8. 结果统计
# ============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("                   统计结果\n")
cat(rep("=", 70), "\n\n", sep = "")

# 热事件指数 → SIF
heat_to_sif_results <- ccm_results_heat %>%
  filter(heat_causes_sif)

cat("【热事件指数 → SIF】\n")
cat("有因果关系的站点数:", nrow(heat_to_sif_results), "\n")
cat("效应类型分布:\n")
print(table(heat_to_sif_results$effect_type_heat))
cat("效应稳定性:\n")
print(table(heat_to_sif_results$effect_stability_heat))
cat("\n")

# 对比VPD均值
vpd_to_sif_results <- ccm_results_heat %>%
  filter(vpd_causes_sif)

cat("【VPD均值 → SIF (对比)】\n")
cat("有因果关系的站点数:", nrow(vpd_to_sif_results), "\n\n")

# 效应指数统计
if (nrow(heat_to_sif_results) > 0) {
  heat_index_stats <- heat_to_sif_results %>%
    filter(!is.na(effect_index_heat)) %>%
    summarise(
      n = n(),
      mean = mean(effect_index_heat),
      sd = sd(effect_index_heat),
      median = median(effect_index_heat),
      min = min(effect_index_heat),
      max = max(effect_index_heat),
      n_negative = sum(effect_index_heat < -0.001),
      pct_negative = round(sum(effect_index_heat < -0.001) / n() * 100, 1)
    )
  
  cat("热事件效应指数统计:\n")
  print(heat_index_stats)
  cat("\n")
}

# 对比表
comparison_heat_vpd <- tibble(
  指标 = c("有因果站点数", "抑制效应站点", "促进效应站点"),
  热事件指数 = c(
    nrow(heat_to_sif_results),
    sum(heat_to_sif_results$effect_type_heat == "抑制效应(-)", na.rm = TRUE),
    sum(heat_to_sif_results$effect_type_heat == "促进效应(+)", na.rm = TRUE)
  ),
  VPD均值 = c(
    nrow(vpd_to_sif_results),
    NA,
    NA
  )
)

cat("【热事件指数 vs VPD均值对比】\n")
print(comparison_heat_vpd)
cat("\n")

# ============================================================================
# 9. 空间可视化
# ============================================================================

cat("【9. 空间可视化】\n")

# 读取坐标
station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并
spatial_heat_sif <- ccm_results_heat %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    effect_strength = abs(effect_index_heat),
    combined_label = case_when(
      !heat_causes_sif ~ "无因果",
      effect_type_heat == "促进效应(+)" ~ "热事件促进SIF",
      effect_type_heat == "抑制效应(-)" ~ "热事件抑制SIF",
      TRUE ~ "效应不明"
    )
  )

cat("有坐标的站点数:", nrow(spatial_heat_sif), "\n\n")

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 主地图
p_heat_spatial <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_heat_sif,
    aes(x = longitude, y = latitude,
        color = combined_label,
        size = effect_strength),
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "热事件促进SIF" = "#4575b4",
      "热事件抑制SIF" = "#d73027",
      "无因果" = "#999999",
      "效应不明" = "#fee090"
    ),
    name = "因果效应"
  ) +
  scale_size_continuous(range = c(1, 8), name = "效应强度\n|∂|") +
  labs(
    title = "VPD热事件指数对SIF的因果影响",
    subtitle = paste0(
      "抑制效应: ", sum(spatial_heat_sif$combined_label == "热事件抑制SIF"), " | ",
      "促进效应: ", sum(spatial_heat_sif$combined_label == "热事件促进SIF"), " | ",
      "VPD阈值 = ", round(VPD_THRESHOLD, 2), " kPa"
    ),
    x = "经度",
    y = "纬度",
    caption = "热事件指数 = 频次 × 强度的综合指标"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p_heat_spatial)

# 分面地图
spatial_heat_sif_sig <- spatial_heat_sif %>%
  filter(combined_label %in% c("热事件促进SIF", "热事件抑制SIF"))

if (nrow(spatial_heat_sif_sig) > 0) {
  p_heat_facet <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(
      data = spatial_heat_sif_sig,
      aes(x = longitude, y = latitude,
          color = effect_stability_heat,
          size = effect_strength),
      alpha = 0.8
    ) +
    scale_color_manual(
      values = c("稳定" = "#1b9e77", "波动大" = "#d95f02"),
      name = "效应稳定性"
    ) +
    scale_size_continuous(range = c(1, 3), name = "|∂|") +
    facet_wrap(~ combined_label, ncol = 2) +
    labs(
      title = "热事件对SIF的效应：促进 vs 抑制",
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
  
  print(p_heat_facet)
}

# ============================================================================
# 10. 对比分析：热事件指数 vs VPD均值
# ============================================================================

cat("\n【10. 热事件指数 vs VPD均值对比】\n\n")

# 计算相关性
comparison_cor <- ccm_results_heat %>%
  filter(!is.na(rho_heat_to_sif), !is.na(rho_vpd_to_sif)) %>%
  summarise(
    cor_rho = cor(rho_heat_to_sif, rho_vpd_to_sif),
    n_both_causal = sum(heat_causes_sif & vpd_causes_sif),
    n_only_heat = sum(heat_causes_sif & !vpd_causes_sif),
    n_only_vpd = sum(!heat_causes_sif & vpd_causes_sif),
    mean_rho_heat = mean(rho_heat_to_sif[heat_causes_sif]),
    mean_rho_vpd = mean(rho_vpd_to_sif[vpd_causes_sif])
  )

cat("CCM强度相关性（热事件指数 vs VPD均值）:", 
    round(comparison_cor$cor_rho, 3), "\n\n")

cat("因果关系对比:\n")
cat("  两者都有因果:", comparison_cor$n_both_causal, "站点\n")
cat("  仅热事件指数有因果:", comparison_cor$n_only_heat, "站点\n")
cat("  仅VPD均值有因果:", comparison_cor$n_only_vpd, "站点\n\n")

cat("平均CCM强度:\n")
cat("  热事件指数:", round(comparison_cor$mean_rho_heat, 3), "\n")
cat("  VPD均值:", round(comparison_cor$mean_rho_vpd, 3), "\n\n")

# 散点图对比
p_comparison_scatter <- ggplot(ccm_results_heat, 
                               aes(x = rho_vpd_to_sif, y = rho_heat_to_sif)) +
  geom_point(aes(color = heat_causes_sif & vpd_causes_sif), 
             alpha = 0.6, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  scale_color_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "#999999"),
    labels = c("单方向因果", "双向因果"),
    name = ""
  ) +
  labs(
    title = "热事件指数 vs VPD均值的CCM强度对比",
    subtitle = paste0("相关系数 r = ", round(comparison_cor$cor_rho, 3)),
    x = "CCM强度：VPD均值 → SIF",
    y = "CCM强度：热事件指数 → SIF",
    caption = "红线 = 1:1参考线；蓝线 = 拟合线"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "bottom"
  )

print(p_comparison_scatter)


