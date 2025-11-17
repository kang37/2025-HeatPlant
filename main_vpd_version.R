pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, 
               tidyr, showtext, rEDM, sf, rnaturalearth, rnaturalearthdata)
showtext_auto()

# ============================================================================
# 数据读取（保持你原有的代码）
# ============================================================================

meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop") %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

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

# ============================================================================
# 计算VPD（参考Yuan et al. 2019）
# 根据论文方法计算VPD
meteo_data_with_vpd <- meteo_data_lin %>%
  mutate(
    vpd = 6.112 * exp((17.67 * tavg) / (tavg + 243.5)) * (1 - rh / 100)
  ) %>%
  # 检查VPD的合理性
  mutate(
    vpd = case_when(
      vpd < 0 ~ NA_real_,      # VPD不能为负
      vpd > 10 ~ NA_real_,     # 异常高值
      TRUE ~ vpd
    )
  )
hist(meteo_data_with_vpd$vpd)
cat("VPD平均值:", round(mean(meteo_data_with_vpd$vpd, na.rm = TRUE), 3), "kPa\n\n")

# 合并所有数据
meteo_data <- meteo_sif_data %>%
  inner_join(meteo_data_with_vpd, by = c("meteo_stat_id", "year", "month")) %>% 
  drop_na() %>% 
  arrange(meteo_stat_id, year, month)

# ============================================================================
# 夏季数据处理
# ============================================================================

meteo_data_summer <- meteo_data %>%
  filter(month %in% 5:10) %>%
  group_by(meteo_stat_id, month) %>%
  mutate(
    temp_q90 = quantile(tavg, 0.9, na.rm = TRUE),
    temp_anomaly = tavg - temp_q90,
    is_hot_month = tavg > temp_q90
  ) %>%
  ungroup() %>%
  arrange(meteo_stat_id, year, month)

# ============================================================================
# 去趋势处理
# ============================================================================

meteo_data_detrended <- meteo_data_summer %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    temp_anomaly_detrended = residuals(lm(temp_anomaly ~ time_idx)),
    sif_detrended = residuals(lm(sif ~ time_idx)),
    vpd_detrended = residuals(lm(vpd ~ time_idx))  # VPD也去趋势
  ) %>%
  ungroup()

data_for_ccm <- meteo_data_detrended %>%
  select(meteo_stat_id, year, month, 
         temp_anomaly = temp_anomaly_detrended,
         sif = sif_detrended,
         vpd = vpd_detrended,  # 添加VPD
         is_hot_month)

# ============================================================================
# CCM + S-map 分析函数（增强版：包含VPD分析）
# ============================================================================

perform_ccm_vpd_analysis <- function(station_id, data, min_points = 30) {
  
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif, vpd) %>%  # 包含VPD
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif,
      vpd = station_data$vpd  # VPD
    )
    
    lib_pred_size <- min(100, n_data)
    
    # ========================================================================
    # 1. 确定最优嵌入维度 E
    # ========================================================================
    
    # Temp
    embed_temp <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "temp_anomaly",
      target = "temp_anomaly",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    # SIF
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    # VPD
    embed_vpd <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "vpd",
      target = "vpd",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    best_E_vpd <- embed_vpd$E[which.max(embed_vpd$rho)]
    
    E_ccm <- max(best_E_temp, best_E_sif, best_E_vpd)
    
    # ========================================================================
    # 2. CCM 分析
    # ========================================================================
    
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
    
    # --- Temp → SIF ---
    ccm_temp_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "temp_anomaly",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # --- VPD → SIF ---
    ccm_vpd_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "vpd",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # CCM汇总
    ccm_summary_temp_sif <- ccm_temp_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    ccm_summary_vpd_sif <- ccm_vpd_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`vpd:sif`, na.rm = TRUE),
        rho_sd = sd(`vpd:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    final_rho_temp_to_sif <- ccm_summary_temp_sif$rho_mean[nrow(ccm_summary_temp_sif)]
    final_rho_vpd_to_sif <- ccm_summary_vpd_sif$rho_mean[nrow(ccm_summary_vpd_sif)]
    
    trend_temp_to_sif <- cor(ccm_summary_temp_sif$lib_size, 
                             ccm_summary_temp_sif$rho_mean)
    trend_vpd_to_sif <- cor(ccm_summary_vpd_sif$lib_size, 
                            ccm_summary_vpd_sif$rho_mean)
    
    # ========================================================================
    # 3. S-map 分析
    # ========================================================================
    
    # --- Temp → SIF ---
    smap_temp_to_sif <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "temp_anomaly",
      target = "sif",
      embedded = FALSE
    )
    
    smap_coeffs_temp <- smap_temp_to_sif$coefficients
    
    if (!is.null(smap_coeffs_temp) && ncol(smap_coeffs_temp) >= 2) {
      temp_coef_col <- which(grepl("temp", colnames(smap_coeffs_temp), ignore.case = TRUE))
      
      if (length(temp_coef_col) > 0) {
        smap_coef_temp <- smap_coeffs_temp[, temp_coef_col[1]]
      } else {
        smap_coef_temp <- smap_coeffs_temp[, 2]
      }
      
      smap_coef_temp <- smap_coef_temp[!is.na(smap_coef_temp)]
      
      mean_smap_coef_temp <- mean(smap_coef_temp, na.rm = TRUE)
      sd_smap_coef_temp <- sd(smap_coef_temp, na.rm = TRUE)
      
    } else {
      mean_smap_coef_temp <- NA
      sd_smap_coef_temp <- NA
      smap_coef_temp <- NA
    }
    
    # --- VPD → SIF ---
    smap_vpd_to_sif <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "vpd",
      target = "sif",
      embedded = FALSE
    )
    
    smap_coeffs_vpd <- smap_vpd_to_sif$coefficients
    
    if (!is.null(smap_coeffs_vpd) && ncol(smap_coeffs_vpd) >= 2) {
      vpd_coef_col <- which(grepl("vpd", colnames(smap_coeffs_vpd), ignore.case = TRUE))
      
      if (length(vpd_coef_col) > 0) {
        smap_coef_vpd <- smap_coeffs_vpd[, vpd_coef_col[1]]
      } else {
        smap_coef_vpd <- smap_coeffs_vpd[, 2]
      }
      
      smap_coef_vpd <- smap_coef_vpd[!is.na(smap_coef_vpd)]
      
      mean_smap_coef_vpd <- mean(smap_coef_vpd, na.rm = TRUE)
      sd_smap_coef_vpd <- sd(smap_coef_vpd, na.rm = TRUE)
      
    } else {
      mean_smap_coef_vpd <- NA
      sd_smap_coef_vpd <- NA
      smap_coef_vpd <- NA
    }
    
    # ========================================================================
    # 4. 判断因果和效应
    # ========================================================================
    
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    # Temp → SIF
    temp_causes_sif <- (final_rho_temp_to_sif > ccm_threshold_rho & 
                          trend_temp_to_sif > ccm_threshold_trend)
    
    effect_type_temp_sif <- case_when(
      !temp_causes_sif ~ "无因果",
      is.na(mean_smap_coef_temp) ~ "S-map失败",
      mean_smap_coef_temp > 0.001 ~ "促进效应(+)",
      mean_smap_coef_temp < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    # 检查效应稳定性
    effect_stability_temp <- case_when(
      is.na(mean_smap_coef_temp) | is.na(sd_smap_coef_temp) ~ "未知",
      abs(mean_smap_coef_temp) > sd_smap_coef_temp ~ "稳定",
      TRUE ~ "波动大"
    )
    
    # VPD → SIF
    vpd_causes_sif <- (final_rho_vpd_to_sif > ccm_threshold_rho & 
                         trend_vpd_to_sif > ccm_threshold_trend)
    
    effect_type_vpd_sif <- case_when(
      !vpd_causes_sif ~ "无因果",
      is.na(mean_smap_coef_vpd) ~ "S-map失败",
      mean_smap_coef_vpd > 0.001 ~ "促进效应(+)",
      mean_smap_coef_vpd < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    effect_stability_vpd <- case_when(
      is.na(mean_smap_coef_vpd) | is.na(sd_smap_coef_vpd) ~ "未知",
      abs(mean_smap_coef_vpd) > sd_smap_coef_vpd ~ "稳定",
      TRUE ~ "波动大"
    )
    
    # ========================================================================
    # 5. 返回结果
    # ========================================================================
    
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,
      
      # Temp → SIF
      rho_temp_to_sif = final_rho_temp_to_sif,
      trend_temp_to_sif = trend_temp_to_sif,
      temp_causes_sif = temp_causes_sif,
      effect_index_temp_sif = mean_smap_coef_temp,
      effect_index_sd_temp_sif = sd_smap_coef_temp,
      effect_type_temp_sif = effect_type_temp_sif,
      effect_stability_temp = effect_stability_temp,
      
      # VPD → SIF
      rho_vpd_to_sif = final_rho_vpd_to_sif,
      trend_vpd_to_sif = trend_vpd_to_sif,
      vpd_causes_sif = vpd_causes_sif,
      effect_index_vpd_sif = mean_smap_coef_vpd,
      effect_index_sd_vpd_sif = sd_smap_coef_vpd,
      effect_type_vpd_sif = effect_type_vpd_sif,
      effect_stability_vpd = effect_stability_vpd,
      
      # 保存数据
      smap_coefficients_temp = list(smap_coef_temp),
      smap_coefficients_vpd = list(smap_coef_vpd),
      ccm_summary_temp = list(ccm_summary_temp_sif),
      ccm_summary_vpd = list(ccm_summary_vpd_sif)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# ============================================================================
# 批量分析
# ============================================================================

cat("\n=== 开始批量CCM+S-map分析（Temp和VPD） ===\n\n")

all_stations <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("准备分析", length(all_stations), "个站点...\n\n")

ccm_results_all <- map_dfr(seq_along(all_stations), function(i) {
  if (i %% 10 == 0) {
    cat("已完成:", i, "/", length(all_stations), "\n")
  }
  perform_ccm_vpd_analysis(all_stations[i], data_for_ccm, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_all), "/", length(all_stations), "\n\n")

# ============================================================================
# 统计结果

# Temp → SIF
cat("【Temp → SIF】\n")
temp_to_sif_results <- ccm_results_all %>%
  filter(temp_causes_sif)

cat("有因果关系的站点数:", nrow(temp_to_sif_results), "\n")
print(table(temp_to_sif_results$effect_type_temp_sif))
cat("效应稳定性:\n")
print(table(temp_to_sif_results$effect_stability_temp))
cat("\n")

# VPD → SIF
cat("【VPD → SIF】\n")
vpd_to_sif_results <- ccm_results_all %>%
  filter(vpd_causes_sif)

cat("有因果关系的站点数:", nrow(vpd_to_sif_results), "\n")
print(table(vpd_to_sif_results$effect_type_vpd_sif))
cat("效应稳定性:\n")
print(table(vpd_to_sif_results$effect_stability_vpd))
cat("\n")

# 对比
cat("【Temp vs VPD 对比】\n")
comparison <- tibble(
  指标 = c("有因果站点数", "抑制效应站点", "促进效应站点", "稳定效应站点"),
  Temp到SIF = c(
    nrow(temp_to_sif_results),
    sum(temp_to_sif_results$effect_type_temp_sif == "抑制效应(-)", na.rm = TRUE),
    sum(temp_to_sif_results$effect_type_temp_sif == "促进效应(+)", na.rm = TRUE),
    sum(temp_to_sif_results$effect_stability_temp == "稳定", na.rm = TRUE)
  ),
  VPD到SIF = c(
    nrow(vpd_to_sif_results),
    sum(vpd_to_sif_results$effect_type_vpd_sif == "抑制效应(-)", na.rm = TRUE),
    sum(vpd_to_sif_results$effect_type_vpd_sif == "促进效应(+)", na.rm = TRUE),
    sum(vpd_to_sif_results$effect_stability_vpd == "稳定", na.rm = TRUE)
  )
)

print(comparison)

cat("\n分析完成！\n")

# ============================================================================
# 空间分布可视化：Temp和VPD对SIF的影响
# ============================================================================

cat("\n【准备空间可视化】\n")

# 读取站点坐标
station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并结果和坐标
spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude))

cat("有坐标的站点数:", nrow(spatial_data), "\n\n")

# 获取中国地图
china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# ============================================================================
# 1. Temp → SIF 空间分布
# ============================================================================

cat("【1. Temp → SIF 空间分布】\n")

spatial_temp_sif <- spatial_data %>%
  filter(temp_causes_sif) %>%
  mutate(
    effect_strength = abs(effect_index_temp_sif)
  )

# 地图边界
bbox <- st_bbox(c(
  xmin = min(spatial_temp_sif$longitude, na.rm = TRUE) - 2,
  xmax = max(spatial_temp_sif$longitude, na.rm = TRUE) + 2,
  ymin = min(spatial_temp_sif$latitude, na.rm = TRUE) - 2,
  ymax = max(spatial_temp_sif$latitude, na.rm = TRUE) + 2
))

# 分面地图：按效应类型分组
spatial_temp_sif_sig <- spatial_temp_sif %>%
  filter(effect_type_temp_sif %in% c("促进效应(+)", "抑制效应(-)"))

if (nrow(spatial_temp_sif_sig) > 0) {
  p_spatial_temp_facet <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(data = spatial_temp_sif_sig,
               aes(x = longitude, y = latitude,
                   color = effect_stability_temp,
                   size = effect_strength),
               alpha = 0.7) +
    scale_color_manual(
      values = c("稳定" = "#1b9e77", "波动大" = "#d95f02"),
      name = "效应稳定性"
    ) +
    scale_size_continuous(range = c(2, 6), name = "|∂|") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"])) +
    facet_wrap(~ effect_type_temp_sif, ncol = 2) +
    labs(
      title = "Temp → SIF: 促进效应 vs 抑制效应的空间分布",
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
  
  print(p_spatial_temp_facet)
}

# 效应指数的空间分布（连续色标）
p_spatial_temp_continuous <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = spatial_temp_sif %>% filter(!is.na(effect_index_temp_sif)),
             aes(x = longitude, y = latitude,
                 color = effect_index_temp_sif,
                 size = abs(effect_index_temp_sif)),
             alpha = 0.7) +
  scale_color_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    name = "效应指数\n∂"
  ) +
  scale_size_continuous(range = c(2, 8), name = "|∂|") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  labs(
    title = "Temp → SIF: 效应指数（S-map系数）的空间分布",
    subtitle = "蓝色=促进效应，红色=抑制效应",
    x = "经度",
    y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p_spatial_temp_continuous)

# ============================================================================
# 2. VPD → SIF 空间分布
# ============================================================================

cat("\n【2. VPD → SIF 空间分布】\n")

spatial_vpd_sif <- spatial_data %>%
  filter(vpd_causes_sif) %>%
  mutate(
    effect_strength = abs(effect_index_vpd_sif)
  )

cat("VPD → SIF 有因果关系的站点数:", nrow(spatial_vpd_sif), "\n")

# 更新地图边界（如果VPD站点分布不同）
bbox_vpd <- st_bbox(c(
  xmin = min(spatial_vpd_sif$longitude, na.rm = TRUE) - 2,
  xmax = max(spatial_vpd_sif$longitude, na.rm = TRUE) + 2,
  ymin = min(spatial_vpd_sif$latitude, na.rm = TRUE) - 2,
  ymax = max(spatial_vpd_sif$latitude, na.rm = TRUE) + 2
))

# 分面地图：按效应类型分组
spatial_vpd_sif_sig <- spatial_vpd_sif %>%
  filter(effect_type_vpd_sif %in% c("促进效应(+)", "抑制效应(-)"))

if (nrow(spatial_vpd_sif_sig) > 0) {
  p_spatial_vpd_facet <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(data = spatial_vpd_sif_sig,
               aes(x = longitude, y = latitude,
                   color = effect_stability_vpd,
                   size = effect_strength),
               alpha = 0.7) +
    scale_color_manual(
      values = c("稳定" = "#1b9e77", "波动大" = "#d95f02"),
      name = "效应稳定性"
    ) +
    scale_size_continuous(range = c(2, 6), name = "|∂|") +
    coord_sf(xlim = c(bbox_vpd["xmin"], bbox_vpd["xmax"]),
             ylim = c(bbox_vpd["ymin"], bbox_vpd["ymax"])) +
    facet_wrap(~ effect_type_vpd_sif, ncol = 2) +
    labs(
      title = "VPD → SIF: 促进效应 vs 抑制效应的空间分布",
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
  
  print(p_spatial_vpd_facet)
}

# 效应指数的空间分布（连续色标）
p_spatial_vpd_continuous <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = spatial_vpd_sif %>% filter(!is.na(effect_index_vpd_sif)),
             aes(x = longitude, y = latitude,
                 color = abs(effect_index_vpd_sif),
                 size = abs(effect_index_vpd_sif)),
             alpha = 0.7) +
  scale_color_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    name = "效应指数\n∂"
  ) +
  scale_size_continuous(range = c(1, 3), name = "|∂|") +
  coord_sf(xlim = c(bbox_vpd["xmin"], bbox_vpd["xmax"]),
           ylim = c(bbox_vpd["ymin"], bbox_vpd["ymax"])) +
  labs(
    title = "VPD → SIF: 效应指数（S-map系数）的空间分布",
    subtitle = "蓝色=促进效应，红色=抑制效应",
    x = "经度",
    y = "纬度"
  ) +
  facet_wrap(.~ c(effect_index_vpd_sif > 0)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p_spatial_vpd_continuous)

# ============================================================================
# 3. Temp vs VPD 对比地图
# ============================================================================

cat("\n【3. Temp vs VPD 对比】\n")

# 准备对比数据
spatial_comparison <- spatial_data %>%
  filter(temp_causes_sif | vpd_causes_sif) %>%
  mutate(
    # 创建组合变量
    driver = case_when(
      temp_causes_sif & !vpd_causes_sif ~ "仅Temp",
      !temp_causes_sif & vpd_causes_sif ~ "仅VPD",
      temp_causes_sif & vpd_causes_sif ~ "Temp+VPD",
      TRUE ~ "无因果"
    ),
    # 主导因子（基于CCM强度）
    dominant = case_when(
      !temp_causes_sif & !vpd_causes_sif ~ "无",
      !vpd_causes_sif ~ "Temp主导",
      !temp_causes_sif ~ "VPD主导",
      abs(rho_temp_to_sif) > abs(rho_vpd_to_sif) ~ "Temp主导",
      TRUE ~ "VPD主导"
    )
  )

# 并排对比图
p_comparison_sidebyside <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = spatial_comparison,
             aes(x = longitude, y = latitude,
                 color = driver,
                 size = pmax(abs(effect_index_temp_sif), 
                             abs(effect_index_vpd_sif), na.rm = TRUE)),
             alpha = 0.7) +
  scale_color_manual(
    values = c(
      "仅Temp" = "#e41a1c",
      "仅VPD" = "#377eb8",
      "Temp+VPD" = "#4daf4a",
      "无因果" = "#999999"
    ),
    name = "因果驱动"
  ) +
  scale_size_continuous(range = c(2, 8), name = "效应强度\n|∂|") +
  labs(
    title = "Temp vs VPD: 对SIF因果影响的空间对比",
    subtitle = paste0("仅Temp: ", sum(spatial_comparison$driver == "仅Temp"), 
                      " | 仅VPD: ", sum(spatial_comparison$driver == "仅VPD"),
                      " | 两者都有: ", sum(spatial_comparison$driver == "Temp+VPD")),
    x = "经度",
    y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p_comparison_sidebyside)

# ============================================================================
# 4. 经济数据关联分析
# ============================================================================

cat("\n【4. VPD-SIF因果关系与经济发展指标的关系】\n\n")

# ----------------------------------------------------------------------------
# 4.1 读取并准备经济数据
# ----------------------------------------------------------------------------

cat("读取经济和空间数据...\n")

# 读取城市shapefile
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326)

# 读取GDP数据
gdp_data <- readxl::read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  select(
    city_name = 城市,
    gdp_per_capita = "人均地区生产总值(元)",
    gdp_total = "地区生产总值(万元)"
  ) %>%
  mutate(
    gdp_per_capita = as.numeric(gdp_per_capita),
    gdp_total = as.numeric(gdp_total)
  )

# 空间匹配站点到城市
stations_sf <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_city <- stations_sf %>%
  st_join(china_cities_shp, join = st_within) %>%
  as.data.frame() %>%
  select(-geometry) %>%
  rename(
    city_code = ct_adcode,
    city_name = ct_name,
    province_name = pr_name
  )

# 合并所有数据
economic_data <- station_coords %>%
  left_join(stations_with_city, by = "meteo_stat_id") %>%
  left_join(gdp_data, by = "city_name") %>%
  select(meteo_stat_id, longitude, latitude, city_name, province_name,
         gdp_per_capita, gdp_total)

# 合并CCM结果和经济数据
spatial_econ <- spatial_data %>%
  left_join(economic_data, by = c("meteo_stat_id", "longitude", "latitude"))

cat("有GDP数据的站点:", sum(!is.na(spatial_econ$gdp_per_capita)), "\n\n")

# ----------------------------------------------------------------------------
# 4.2 创建综合因果标签
# ----------------------------------------------------------------------------

spatial_econ <- spatial_econ %>%
  mutate(
    # VPD因果标签
    vpd_causality = case_when(
      !vpd_causes_sif ~ "无因果",
      effect_type_vpd_sif == "促进效应(+)" ~ "VPD促进SIF",
      effect_type_vpd_sif == "抑制效应(-)" ~ "VPD抑制SIF",
      TRUE ~ "效应不明"
    ),
    
    # Temp因果标签
    temp_causality = case_when(
      !temp_causes_sif ~ "无因果",
      effect_type_temp_sif == "促进效应(+)" ~ "Temp促进SIF",
      effect_type_temp_sif == "抑制效应(-)" ~ "Temp抑制SIF",
      TRUE ~ "效应不明"
    ),
    
    # 综合驱动因子
    driver_type = case_when(
      temp_causes_sif & vpd_causes_sif ~ "双因子驱动",
      temp_causes_sif ~ "仅温度驱动",
      vpd_causes_sif ~ "仅VPD驱动",
      TRUE ~ "无驱动"
    )
  )

# ----------------------------------------------------------------------------
# 4.3 通用可视化函数
# ----------------------------------------------------------------------------

# GDP分级函数
categorize_gdp <- function(data, gdp_column, breaks, labels) {
  data %>%
    mutate(
      gdp_category = cut(
        .data[[gdp_column]],
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE
      )
    )
}

# 堆积条形图
plot_stacked_bar <- function(data, causality_col, title_text, gdp_type) {
  data %>%
    filter(.data[[causality_col]] != "效应不明", !is.na(gdp_category)) %>%
    count(gdp_category, .data[[causality_col]]) %>%
    ggplot(aes(x = gdp_category, y = n, fill = .data[[causality_col]])) +
    geom_col(position = "fill") +
    scale_fill_manual(
      values = c(
        "VPD促进SIF" = "#4575b4",
        "VPD抑制SIF" = "#d73027",
        "Temp促进SIF" = "#91cf60",
        "Temp抑制SIF" = "#fc8d59",
        "无因果" = "#999999"
      ),
      name = "因果类型"
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = title_text,
      subtitle = paste0("不同", gdp_type, "水平下的因果类型分布"),
      x = paste0(gdp_type, "类别"),
      y = "比例 (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )
}

# 箱线图
plot_boxplot_econ <- function(data, causality_col, gdp_column, title_text, y_label) {
  data %>%
    filter(.data[[causality_col]] %in% c("VPD促进SIF", "VPD抑制SIF", "Temp促进SIF", "Temp抑制SIF", "无因果"),
           !is.na(.data[[gdp_column]])) %>%
    ggplot(aes(x = reorder(.data[[causality_col]], -.data[[gdp_column]], FUN = median), 
               y = .data[[gdp_column]],
               fill = .data[[causality_col]])) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
    scale_fill_manual(
      values = c(
        "VPD促进SIF" = "#4575b4",
        "VPD抑制SIF" = "#d73027",
        "Temp促进SIF" = "#91cf60",
        "Temp抑制SIF" = "#fc8d59",
        "无因果" = "#999999"
      ),
      guide = "none"
    ) +
    labs(
      title = title_text,
      x = "因果类型",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 9)
    )
}

# ----------------------------------------------------------------------------
# 4.4 VPD-SIF关系与经济分析
# ----------------------------------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("【VPD-SIF因果关系与经济指标分析】\n")
cat(rep("=", 70), "\n\n", sep = "")

# 过滤有效数据
spatial_econ_valid <- spatial_econ %>%
  filter(!is.na(gdp_per_capita) | !is.na(gdp_total))

if (nrow(spatial_econ_valid) == 0) {
  cat("警告: 无有效经济数据\n\n")
} else {
  
  # === A. 人均GDP分析 ===
  cat("【A. 人均GDP与VPD-SIF因果关系】\n\n")
  
  spatial_econ_percapita <- spatial_econ_valid %>%
    filter(!is.na(gdp_per_capita))
  
  if (nrow(spatial_econ_percapita) > 0) {
    
    # 统计表
    cat("1. VPD因果类型的经济统计:\n")
    vpd_gdp_stats <- spatial_econ_percapita %>%
      group_by(vpd_causality) %>%
      summarise(
        站点数 = n(),
        平均GDP_万元 = round(mean(gdp_per_capita) / 10000, 2),
        中位GDP_万元 = round(median(gdp_per_capita) / 10000, 2),
        GDP标准差_万元 = round(sd(gdp_per_capita) / 10000, 2),
        .groups = "drop"
      ) %>%
      arrange(desc(站点数))
    
    print(vpd_gdp_stats)
    cat("\n")
    
    # 对比Temp因果类型
    cat("2. Temp因果类型的经济统计:\n")
    temp_gdp_stats <- spatial_econ_percapita %>%
      group_by(temp_causality) %>%
      summarise(
        站点数 = n(),
        平均GDP_万元 = round(mean(gdp_per_capita) / 10000, 2),
        中位GDP_万元 = round(median(gdp_per_capita) / 10000, 2),
        .groups = "drop"
      ) %>%
      arrange(desc(站点数))
    
    print(temp_gdp_stats)
    cat("\n")
    
    # GDP分级
    spatial_econ_percapita_cat <- categorize_gdp(
      spatial_econ_percapita,
      "gdp_per_capita",
      breaks = c(0, 40000, 80000, 120000, Inf),
      labels = c("低(<4万)", "中(4-8万)", "较高(8-12万)", "高(>12万)")
    )
    
    # 图1: VPD堆积条形图
    p1_vpd <- plot_stacked_bar(
      spatial_econ_percapita_cat,
      "vpd_causality",
      "VPD-SIF因果关系",
      "人均GDP"
    )
    print(p1_vpd)
    
    # 图2: Temp堆积条形图（对比）
    p1_temp <- plot_stacked_bar(
      spatial_econ_percapita_cat,
      "temp_causality",
      "Temp-SIF因果关系",
      "人均GDP"
    )
    print(p1_temp)
    
    # 图3: VPD箱线图
    spatial_econ_percapita_plot <- spatial_econ_percapita %>%
      mutate(gdp_per_capita_万元 = gdp_per_capita / 10000)
    
    p2_vpd <- plot_boxplot_econ(
      spatial_econ_percapita_plot,
      "vpd_causality",
      "gdp_per_capita_万元",
      "VPD-SIF因果类型与人均GDP",
      "人均GDP (万元)"
    )
    print(p2_vpd)
    
    # 图4: Temp箱线图（对比）
    p2_temp <- plot_boxplot_econ(
      spatial_econ_percapita_plot,
      "temp_causality",
      "gdp_per_capita_万元",
      "Temp-SIF因果类型与人均GDP",
      "人均GDP (万元)"
    )
    print(p2_temp)
    
    # 统计检验
    cat("\n3. 统计检验:\n")
    
    # VPD检验
    test_vpd <- spatial_econ_percapita %>%
      filter(vpd_causality %in% c("VPD促进SIF", "VPD抑制SIF", "无因果"))
    
    if (nrow(test_vpd) > 10 && length(unique(test_vpd$vpd_causality)) >= 2) {
      kw_vpd <- kruskal.test(gdp_per_capita ~ vpd_causality, data = test_vpd)
      cat("  VPD因果: χ² =", round(kw_vpd$statistic, 3), 
          ", p =", format.pval(kw_vpd$p.value, digits = 3), "\n")
    }
    
    # Temp检验
    test_temp <- spatial_econ_percapita %>%
      filter(temp_causality %in% c("Temp促进SIF", "Temp抑制SIF", "无因果"))
    
    if (nrow(test_temp) > 10 && length(unique(test_temp$temp_causality)) >= 2) {
      kw_temp <- kruskal.test(gdp_per_capita ~ temp_causality, data = test_temp)
      cat("  Temp因果: χ² =", round(kw_temp$statistic, 3), 
          ", p =", format.pval(kw_temp$p.value, digits = 3), "\n")
    }
    cat("\n")
  }
  
  # === B. 总GDP分析 ===
  cat("【B. 总GDP与VPD-SIF因果关系】\n\n")
  
  spatial_econ_total <- spatial_econ_valid %>%
    filter(!is.na(gdp_total))
  
  if (nrow(spatial_econ_total) > 0) {
    
    # 统计表
    cat("1. VPD因果类型的总GDP统计:\n")
    vpd_gdp_total_stats <- spatial_econ_total %>%
      group_by(vpd_causality) %>%
      summarise(
        站点数 = n(),
        平均GDP_亿元 = round(mean(gdp_total), 2),
        中位GDP_亿元 = round(median(gdp_total), 2),
        .groups = "drop"
      ) %>%
      arrange(desc(站点数))
    
    print(vpd_gdp_total_stats)
    cat("\n")
    
    # GDP分级（四分位数）
    gdp_total_quantiles <- quantile(spatial_econ_total$gdp_total, 
                                    probs = c(0, 0.25, 0.5, 0.75, 1), 
                                    na.rm = TRUE)
    
    spatial_econ_total_cat <- categorize_gdp(
      spatial_econ_total,
      "gdp_total",
      breaks = gdp_total_quantiles,
      labels = c("低(Q1)", "中(Q2)", "较高(Q3)", "高(Q4)")
    )
    
    cat("总GDP分级 (按四分位数):\n")
    cat("  Q1: <", round(gdp_total_quantiles[2], 1), "亿元\n", sep = "")
    cat("  Q2: ", round(gdp_total_quantiles[2], 1), "-", 
        round(gdp_total_quantiles[3], 1), "亿元\n", sep = "")
    cat("  Q3: ", round(gdp_total_quantiles[3], 1), "-", 
        round(gdp_total_quantiles[4], 1), "亿元\n", sep = "")
    cat("  Q4: >", round(gdp_total_quantiles[4], 1), "亿元\n\n", sep = "")
    
    # 图5: VPD堆积条形图
    p3_vpd <- plot_stacked_bar(
      spatial_econ_total_cat,
      "vpd_causality",
      "VPD-SIF因果关系",
      "总GDP"
    )
    print(p3_vpd)
    
    # 图6: Temp堆积条形图
    p3_temp <- plot_stacked_bar(
      spatial_econ_total_cat,
      "temp_causality",
      "Temp-SIF因果关系",
      "总GDP"
    )
    print(p3_temp)
    
    # 图7: VPD箱线图
    p4_vpd <- plot_boxplot_econ(
      spatial_econ_total,
      "vpd_causality",
      "gdp_total",
      "VPD-SIF因果类型与总GDP",
      "总GDP (亿元)"
    )
    print(p4_vpd)
    
    # 图8: Temp箱线图
    p4_temp <- plot_boxplot_econ(
      spatial_econ_total,
      "temp_causality",
      "gdp_total",
      "Temp-SIF因果类型与总GDP",
      "总GDP (亿元)"
    )
    print(p4_temp)
    
    # 统计检验
    cat("\n2. 统计检验:\n")
    
    test_vpd_total <- spatial_econ_total %>%
      filter(vpd_causality %in% c("VPD促进SIF", "VPD抑制SIF", "无因果"))
    
    if (nrow(test_vpd_total) > 10) {
      kw_vpd_total <- kruskal.test(gdp_total ~ vpd_causality, data = test_vpd_total)
      cat("  VPD因果: χ² =", round(kw_vpd_total$statistic, 3), 
          ", p =", format.pval(kw_vpd_total$p.value, digits = 3), "\n")
    }
    cat("\n")
  }
  
  # === C. 驱动因子对比分析 ===
  cat("【C. 双因子驱动分析】\n\n")
  
  driver_gdp_stats <- spatial_econ_percapita %>%
    group_by(driver_type) %>%
    summarise(
      站点数 = n(),
      平均GDP_万元 = round(mean(gdp_per_capita) / 10000, 2),
      中位GDP_万元 = round(median(gdp_per_capita) / 10000, 2),
      .groups = "drop"
    ) %>%
    arrange(desc(站点数))
  
  cat("驱动因子类型的经济统计:\n")
  print(driver_gdp_stats)
  cat("\n")
  
  # 图9: 驱动因子箱线图
  p5_driver <- spatial_econ_percapita %>%
    mutate(gdp_per_capita_万元 = gdp_per_capita / 10000) %>%
    filter(driver_type != "无驱动") %>%
    ggplot(aes(x = reorder(driver_type, -gdp_per_capita_万元, FUN = median),
               y = gdp_per_capita_万元,
               fill = driver_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    scale_fill_manual(
      values = c(
        "仅温度驱动" = "#fc8d59",
        "仅VPD驱动" = "#4575b4",
        "双因子驱动" = "#4daf4a"
      ),
      guide = "none"
    ) +
    labs(
      title = "驱动因子类型与人均GDP的关系",
      subtitle = "比较单一驱动 vs 双因子驱动",
      x = "驱动因子类型",
      y = "人均GDP (万元)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  print(p5_driver)
  
  # === D. 空间地图：叠加GDP ===
  cat("\n【D. GDP空间分布图】\n")
  
  p_gdp_spatial_vpd <- ggplot() +
    geom_sf(data = china_cities_shp, fill = "gray95", color = "gray70", 
            linewidth = 0.2, alpha = 0.3) +
    geom_point(
      data = spatial_econ_percapita,
      aes(x = longitude, y = latitude,
          size = gdp_per_capita / 10000,
          color = vpd_causality),
      alpha = 0.7
    ) +
    scale_color_manual(
      values = c(
        "VPD促进SIF" = "#4575b4",
        "VPD抑制SIF" = "#d73027",
        "无因果" = "#999999",
        "效应不明" = "#fee090"
      ),
      name = "VPD因果类型"
    ) +
    scale_size_continuous(
      range = c(1, 10), 
      name = "人均GDP\n(万元)"
    ) +
    labs(
      title = "VPD-SIF因果关系与GDP的空间分布",
      subtitle = "点的大小=GDP水平，颜色=因果类型",
      x = "经度",
      y = "纬度"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      panel.background = element_rect(fill = "aliceblue", color = NA)
    )
  
  print(p_gdp_spatial_vpd)
}

