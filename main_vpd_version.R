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
    # 1. 计算饱和水汽压 SVP (kPa)
    # Buck方程: SVP = 0.61121 * exp((17.67*T)/(T+243.5))
    p_mst = 1013.25 * ((tavg + 273.16) / (tavg + 273.16 + 0.0065 * 0))^5.625, 
    f_w = 1 + 7 * 10^(-4) + 3.46 * 10^(-6) * p_mst, 
    svp = 6.112 * f_w * exp((17.67 * tavg) / (tavg + 243.5)),
    
    # 2. 计算实际水汽压 AVP (kPa)
    # AVP = (RH/100) * SVP
    avp = (rh / 100) * svp,
    
    # 3. 计算VPD (kPa)
    # VPD = SVP - AVP
    vpd = svp - avp
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
  filter(month %in% 6:8) %>%
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

