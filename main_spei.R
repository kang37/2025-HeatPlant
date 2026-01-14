pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, 
               ncdf4, rEDM, sf, rnaturalearth)

# ============================================================================
# 1. 读取站点SIF数据
# ============================================================================

meteo_sif_data <- fread("data_raw/meteo_stat_SIF_data.csv") %>% 
  rename_with(tolower) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop") %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 读取站点坐标
station_coords <- fread("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(tolower) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# ============================================================================
# 2. 批量读取ERA5干旱指数（优化版）
# ============================================================================

cat("\n读取ERA5干旱指数...\n")

# 获取所有nc文件
nc_files <- list.files("data_raw/era5_drought", pattern = "SPEI1.*\\.nc$", 
                       full.names = TRUE)

cat("找到", length(nc_files), "个文件\n")

# 快速提取单个nc文件的数据
extract_nc_fast <- function(nc_file, stations) {
  
  # 从文件名提取年月
  ym <- str_extract(basename(nc_file), "\\d{6}(?=\\.area)")
  if (is.na(ym)) return(NULL)
  
  year <- as.integer(substr(ym, 1, 4))
  month <- as.integer(substr(ym, 5, 6))
  
  # 打开nc文件
  nc <- nc_open(nc_file)
  
  # 读取坐标和数据
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  spei <- ncvar_get(nc, "SPEI1")  # 变量名是大写SPEI1
  
  nc_close(nc)
  
  # 为每个站点找最近网格
  station_spei <- stations %>%
    rowwise() %>%
    mutate(
      lon_idx = which.min(abs(lon - longitude)),
      lat_idx = which.min(abs(lat - latitude)),
      spei = spei[lon_idx, lat_idx],
      year = year,
      month = month
    ) %>%
    ungroup() %>%
    select(meteo_stat_id, year, month, spei)
  
  return(station_spei)
}

# 批量处理（仅处理2000-2023年）
drought_data <- map_dfr(nc_files, ~extract_nc_fast(.x, station_coords), 
                        .progress = TRUE) %>%
  filter(!is.na(spei), !is.infinite(spei))

cat("提取完成! 数据行数:", nrow(drought_data), "\n")
cat("SPEI统计: 均值=", round(mean(drought_data$spei, na.rm=TRUE), 3),
    ", 范围=[", round(min(drought_data$spei, na.rm=TRUE), 2), ",",
    round(max(drought_data$spei, na.rm=TRUE), 2), "]\n\n")

# ============================================================================
# 3. 读取温度数据（只读必要变量）
# ============================================================================

cat("读取温度数据...\n")

read_one_meteo <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  fread(path, skip = 1, na.strings = c("", "NA")) %>%
    filter(across(where(is.numeric), ~.x < 999990)) %>%
    mutate(meteo_stat_id = station_id) %>%
    select(meteo_stat_id, date, tavg) %>%  # 小写列名
    return()
}

meteo_file_lin <- list.files("data_raw/meteo_data_1961-2023",
                             pattern = "\\.txt$", full.names = TRUE) %>%
  .[!grepl("sta_lonlat_china.txt", .)]

# 只读取有SIF数据的站点
target_stations <- unique(meteo_sif_data$meteo_stat_id)

meteo_data_lin <- map_dfr(
  meteo_file_lin[str_extract(basename(meteo_file_lin), "\\d+") %in% target_stations],
  read_one_meteo
) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date)
  ) %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(tavg = mean(tavg, na.rm = TRUE), .groups = "drop")

# ============================================================================
# 4. 合并所有数据
# ============================================================================

meteo_data <- meteo_sif_data %>%
  inner_join(meteo_data_lin, by = c("meteo_stat_id", "year", "month")) %>%
  inner_join(drought_data, by = c("meteo_stat_id", "year", "month")) %>%
  drop_na() %>%
  arrange(meteo_stat_id, year, month)

cat("合并后数据:", nrow(meteo_data), "行,", 
    length(unique(meteo_data$meteo_stat_id)), "个站点\n\n")

# ============================================================================
# 5. 夏季数据 + 去趋势
# ============================================================================

data_for_ccm <- meteo_data %>%
  filter(month %in% 5:10) %>%
  group_by(meteo_stat_id, month) %>%
  mutate(
    temp_q90 = quantile(tavg, 0.9, na.rm = TRUE),
    temp_anomaly = tavg - temp_q90
  ) %>%
  ungroup() %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    temp_anomaly = residuals(lm(temp_anomaly ~ time_idx)),
    sif = residuals(lm(sif ~ time_idx)),
    spei = residuals(lm(spei ~ time_idx))  # SPEI去趋势
  ) %>%
  ungroup() %>%
  select(meteo_stat_id, year, month, temp_anomaly, sif, spei)

# ============================================================================
# 6. CCM分析函数（简化版）
# ============================================================================

perform_ccm_drought <- function(station_id, data, min_points = 30) {
  
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif, spei) %>%
    as.data.frame()
  
  n_data <- nrow(station_data)
  if (n_data < min_points) return(NULL)
  
  tryCatch({
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif,
      spei = station_data$spei
    )
    
    lib_pred_size <- min(100, n_data)
    
    # 确定最优E
    embed_spei <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "spei",
      target = "spei",
      maxE = min(6, floor(n_data/10)),
      showPlot = FALSE
    )
    
    E_ccm <- embed_spei$E[which.max(embed_spei$rho)]
    
    # CCM设置
    tau <- 1
    tp <- 1
    embedding_loss <- (E_ccm - 1) * tau
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) return(NULL)
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # CCM: SPEI → SIF
    ccm_spei_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "spei",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    ccm_summary <- ccm_spei_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(rho_mean = mean(`spei:sif`, na.rm = TRUE), .groups = "drop")
    
    final_rho <- ccm_summary$rho_mean[nrow(ccm_summary)]
    trend_rho <- cor(ccm_summary$lib_size, ccm_summary$rho_mean)
    
    # S-map: SPEI → SIF
    smap_spei_to_sif <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "spei",
      target = "sif",
      embedded = FALSE
    )
    
    smap_coeffs <- smap_spei_to_sif$coefficients
    
    if (!is.null(smap_coeffs) && ncol(smap_coeffs) >= 2) {
      spei_coef <- smap_coeffs[, 2][!is.na(smap_coeffs[, 2])]
      mean_smap_coef <- mean(spei_coef, na.rm = TRUE)
      sd_smap_coef <- sd(spei_coef, na.rm = TRUE)
    } else {
      mean_smap_coef <- NA
      sd_smap_coef <- NA
      spei_coef <- NA
    }
    
    # 判断因果
    spei_causes_sif <- (final_rho > 0.1 & trend_rho > 0)
    
    effect_type <- case_when(
      !spei_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0.001 ~ "促进效应(+)",
      mean_smap_coef < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,
      rho_spei_to_sif = final_rho,
      trend_spei_to_sif = trend_rho,
      spei_causes_sif = spei_causes_sif,
      effect_index = mean_smap_coef,
      effect_index_sd = sd_smap_coef,
      effect_type = effect_type,
      smap_coefficients = list(spei_coef),
      ccm_summary = list(ccm_summary)
    )
    
  }, error = function(e) {
    return(NULL)
  })
}

# ============================================================================
# 7. 批量分析
# ============================================================================

cat("\n=== 开始CCM分析 ===\n")

all_stations <- data_for_ccm %>%
  count(meteo_stat_id) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("分析站点数:", length(all_stations), "\n\n")

ccm_results_all <- map_dfr(all_stations, 
                           ~perform_ccm_drought(.x, data_for_ccm, 30),
                           .progress = TRUE)

cat("\n成功分析:", nrow(ccm_results_all), "个站点\n\n")

# ============================================================================
# 8. 结果统计
# ============================================================================

cat("【SPEI → SIF 因果分析结果】\n\n")

spei_to_sif <- ccm_results_all %>% filter(spei_causes_sif)

cat("有因果关系的站点:", nrow(spei_to_sif), "/", nrow(ccm_results_all), "\n")
cat("因果比例:", round(100 * nrow(spei_to_sif) / nrow(ccm_results_all), 1), "%\n\n")

cat("效应类型分布:\n")
print(table(spei_to_sif$effect_type))
cat("\n")

cat("效应统计:\n")
summary_stats <- spei_to_sif %>%
  filter(effect_type %in% c("促进效应(+)", "抑制效应(-)")) %>%
  summarise(
    平均效应指数 = mean(effect_index, na.rm = TRUE),
    中位效应指数 = median(effect_index, na.rm = TRUE),
    效应标准差 = sd(effect_index, na.rm = TRUE)
  )
print(summary_stats)
cat("\n")

# ============================================================================
# 9. 空间可视化
# ============================================================================

cat("【空间可视化】\n")

spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude))

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 地图1: 因果关系空间分布
p1 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_data %>% filter(spei_causes_sif),
    aes(x = longitude, y = latitude, 
        color = effect_type,
        size = abs(effect_index)),
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "促进效应(+)" = "#4575b4",
      "抑制效应(-)" = "#d73027",
      "效应极弱" = "#fee090"
    ),
    name = "效应类型"
  ) +
  scale_size_continuous(range = c(1, 5), name = "|效应指数|") +
  labs(
    title = "SPEI → SIF 因果关系的空间分布",
    subtitle = paste0("有因果站点: ", nrow(spei_to_sif), " / 总站点: ", 
                      nrow(ccm_results_all)),
    x = "经度", y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p1)

# 地图2: 效应指数连续分布
p2 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_data %>% filter(spei_causes_sif, !is.na(effect_index)),
    aes(x = longitude, y = latitude,
        color = effect_index,
        size = abs(effect_index)),
    alpha = 0.7
  ) +
  scale_color_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    name = "效应指数"
  ) +
  scale_size_continuous(range = c(1, 5), name = "|效应|") +
  labs(
    title = "SPEI → SIF 效应指数的空间梯度",
    subtitle = "蓝色=干旱促进SIF，红色=干旱抑制SIF",
    x = "经度", y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p2)

cat("\n分析完成!\n")

# 保存结果
# fwrite(ccm_results_all, "ccm_results_spei_sif.csv")
# fwrite(spatial_data, "ccm_results_spatial.csv")