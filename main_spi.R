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
nc_files <- list.files("data_raw/era5_drought", pattern = "SPI1.*\\.nc$", 
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
  SPI <- ncvar_get(nc, "SPI1")  # 变量名是大写SPI1
  
  nc_close(nc)
  
  # 为每个站点找最近网格
  station_SPI <- stations %>%
    rowwise() %>%
    mutate(
      lon_idx = which.min(abs(lon - longitude)),
      lat_idx = which.min(abs(lat - latitude)),
      SPI = SPI[lon_idx, lat_idx],
      year = year,
      month = month
    ) %>%
    ungroup() %>%
    select(meteo_stat_id, year, month, SPI)
  
  return(station_SPI)
}

# 批量处理（仅处理2000-2023年）
drought_data <- map_dfr(nc_files, ~extract_nc_fast(.x, station_coords), 
                        .progress = TRUE) %>%
  filter(!is.na(SPI), !is.infinite(SPI))

cat("提取完成! 数据行数:", nrow(drought_data), "\n")
cat("SPI统计: 均值=", round(mean(drought_data$SPI, na.rm=TRUE), 3),
    ", 范围=[", round(min(drought_data$SPI, na.rm=TRUE), 2), ",",
    round(max(drought_data$SPI, na.rm=TRUE), 2), "]\n\n")

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
    SPI = residuals(lm(SPI ~ time_idx))  # SPI去趋势
  ) %>%
  ungroup() %>%
  select(meteo_stat_id, year, month, temp_anomaly, sif, SPI)

# ============================================================================
# 6. CCM分析函数（双向因果版）
# ============================================================================

perform_ccm_drought <- function(station_id, data, min_points = 30) {

  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif, SPI) %>%
    as.data.frame()

  n_data <- nrow(station_data)
  if (n_data < min_points) return(NULL)

  tryCatch({

    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif,
      SPI = station_data$SPI
    )

    lib_pred_size <- min(100, n_data)

    # 确定最优E（分别为SPI和SIF）
    embed_SPI <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "SPI",
      target = "SPI",
      maxE = min(6, floor(n_data/10)),
      showPlot = FALSE
    )

    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(6, floor(n_data/10)),
      showPlot = FALSE
    )

    best_E_SPI <- embed_SPI$E[which.max(embed_SPI$rho)]
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    E_ccm <- max(best_E_SPI, best_E_sif)

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

    # ========== CCM: SPI → SIF ==========
    ccm_SPI_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "SPI",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )

    ccm_summary_SPI_sif <- ccm_SPI_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(rho_mean = mean(`SPI:sif`, na.rm = TRUE), .groups = "drop")

    final_rho_SPI_to_sif <- ccm_summary_SPI_sif$rho_mean[nrow(ccm_summary_SPI_sif)]
    trend_SPI_to_sif <- cor(ccm_summary_SPI_sif$lib_size, ccm_summary_SPI_sif$rho_mean)

    # ========== CCM: SIF → SPI ==========
    ccm_sif_to_SPI <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "sif",
      target = "SPI",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )

    ccm_summary_sif_SPI <- ccm_sif_to_SPI %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(rho_mean = mean(`sif:SPI`, na.rm = TRUE), .groups = "drop")

    final_rho_sif_to_SPI <- ccm_summary_sif_SPI$rho_mean[nrow(ccm_summary_sif_SPI)]
    trend_sif_to_SPI <- cor(ccm_summary_sif_SPI$lib_size, ccm_summary_sif_SPI$rho_mean)

    # ========== S-map: SPI → SIF ==========
    smap_SPI_to_sif <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "SPI",
      target = "sif",
      embedded = FALSE
    )

    smap_coeffs_SPI <- smap_SPI_to_sif$coefficients

    if (!is.null(smap_coeffs_SPI) && ncol(smap_coeffs_SPI) >= 2) {
      SPI_coef <- smap_coeffs_SPI[, 2][!is.na(smap_coeffs_SPI[, 2])]
      mean_smap_coef_SPI_sif <- mean(SPI_coef, na.rm = TRUE)
      sd_smap_coef_SPI_sif <- sd(SPI_coef, na.rm = TRUE)
    } else {
      mean_smap_coef_SPI_sif <- NA
      sd_smap_coef_SPI_sif <- NA
      SPI_coef <- NA
    }

    # ========== S-map: SIF → SPI ==========
    smap_sif_to_SPI <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "sif",
      target = "SPI",
      embedded = FALSE
    )

    smap_coeffs_sif <- smap_sif_to_SPI$coefficients

    if (!is.null(smap_coeffs_sif) && ncol(smap_coeffs_sif) >= 2) {
      sif_coef <- smap_coeffs_sif[, 2][!is.na(smap_coeffs_sif[, 2])]
      mean_smap_coef_sif_SPI <- mean(sif_coef, na.rm = TRUE)
      sd_smap_coef_sif_SPI <- sd(sif_coef, na.rm = TRUE)
    } else {
      mean_smap_coef_sif_SPI <- NA
      sd_smap_coef_sif_SPI <- NA
      sif_coef <- NA
    }

    # ========== 判断因果方向 ==========
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0

    SPI_causes_sif <- (final_rho_SPI_to_sif > ccm_threshold_rho &
                         trend_SPI_to_sif > ccm_threshold_trend)
    sif_causes_SPI <- (final_rho_sif_to_SPI > ccm_threshold_rho &
                         trend_sif_to_SPI > ccm_threshold_trend)

    # 因果方向分类
    causality_direction <- case_when(
      SPI_causes_sif & !sif_causes_SPI ~ "SPI → SIF",
      !SPI_causes_sif & sif_causes_SPI ~ "SIF → SPI",
      SPI_causes_sif & sif_causes_SPI ~ "双向因果",
      TRUE ~ "无显著因果"
    )

    # 效应类型（SPI → SIF）
    effect_type_SPI_sif <- case_when(
      !SPI_causes_sif ~ "无因果",
      is.na(mean_smap_coef_SPI_sif) ~ "S-map失败",
      mean_smap_coef_SPI_sif > 0.001 ~ "促进效应(+)",
      mean_smap_coef_SPI_sif < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )

    # 效应类型（SIF → SPI）
    effect_type_sif_SPI <- case_when(
      !sif_causes_SPI ~ "无因果",
      is.na(mean_smap_coef_sif_SPI) ~ "S-map失败",
      mean_smap_coef_sif_SPI > 0.001 ~ "促进效应(+)",
      mean_smap_coef_sif_SPI < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )

    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,

      # CCM结果
      rho_SPI_to_sif = final_rho_SPI_to_sif,
      trend_SPI_to_sif = trend_SPI_to_sif,
      rho_sif_to_SPI = final_rho_sif_to_SPI,
      trend_sif_to_SPI = trend_sif_to_SPI,

      # 因果判断
      SPI_causes_sif = SPI_causes_sif,
      sif_causes_SPI = sif_causes_SPI,
      causality_direction = causality_direction,

      # S-map效应指数（SPI → SIF）
      effect_index_SPI_sif = mean_smap_coef_SPI_sif,
      effect_index_sd_SPI_sif = sd_smap_coef_SPI_sif,
      effect_type_SPI_sif = effect_type_SPI_sif,

      # S-map效应指数（SIF → SPI）
      effect_index_sif_SPI = mean_smap_coef_sif_SPI,
      effect_index_sd_sif_SPI = sd_smap_coef_sif_SPI,
      effect_type_sif_SPI = effect_type_sif_SPI,

      # 保存详细数据
      smap_coefficients_SPI_sif = list(SPI_coef),
      smap_coefficients_sif_SPI = list(sif_coef)
    )

  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
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

cat("\n", rep("=", 70), "\n", sep = "")
cat("                   统计结果\n")
cat(rep("=", 70), "\n\n", sep = "")

# 1. 因果方向统计
cat("【1. 因果方向统计】\n")
causality_stats <- table(ccm_results_all$causality_direction)
print(causality_stats)
cat("\n百分比:\n")
print(round(prop.table(causality_stats) * 100, 1))
cat("\n")

# 2. SPI → SIF 的效应类型统计
cat("【2. SPI → SIF 效应类型统计】\n")
SPI_to_sif_results <- ccm_results_all %>%
  filter(causality_direction %in% c("SPI → SIF", "双向因果"))

if (nrow(SPI_to_sif_results) > 0) {
  effect_stats_SPI_sif <- table(SPI_to_sif_results$effect_type_SPI_sif)
  print(effect_stats_SPI_sif)
  cat("\n百分比:\n")
  print(round(prop.table(effect_stats_SPI_sif) * 100, 1))
} else {
  cat("无 SPI → SIF 因果关系的站点\n")
}
cat("\n")

# 3. SIF → SPI 的效应类型统计
cat("【3. SIF → SPI 效应类型统计】\n")
sif_to_SPI_results <- ccm_results_all %>%
  filter(causality_direction %in% c("SIF → SPI", "双向因果"))

if (nrow(sif_to_SPI_results) > 0) {
  effect_stats_sif_SPI <- table(sif_to_SPI_results$effect_type_sif_SPI)
  print(effect_stats_sif_SPI)
  cat("\n百分比:\n")
  print(round(prop.table(effect_stats_sif_SPI) * 100, 1))
} else {
  cat("无 SIF → SPI 因果关系的站点\n")
}
cat("\n")

# 4. 效应指数统计
cat("【4. 效应指数（S-map系数）统计】\n\n")

cat("SPI → SIF 方向:\n")
SPI_sif_index_stats <- SPI_to_sif_results %>%
  filter(!is.na(effect_index_SPI_sif)) %>%
  summarise(
    n = n(),
    mean = mean(effect_index_SPI_sif),
    sd = sd(effect_index_SPI_sif),
    median = median(effect_index_SPI_sif),
    min = min(effect_index_SPI_sif),
    max = max(effect_index_SPI_sif),
    n_positive = sum(effect_index_SPI_sif > 0.001),
    n_negative = sum(effect_index_SPI_sif < -0.001),
    pct_negative = round(sum(effect_index_SPI_sif < -0.001) / n() * 100, 1)
  )
print(SPI_sif_index_stats)
cat("\n")

if (nrow(sif_to_SPI_results) > 0) {
  cat("SIF → SPI 方向:\n")
  sif_SPI_index_stats <- sif_to_SPI_results %>%
    filter(!is.na(effect_index_sif_SPI)) %>%
    summarise(
      n = n(),
      mean = mean(effect_index_sif_SPI),
      sd = sd(effect_index_sif_SPI),
      median = median(effect_index_sif_SPI),
      min = min(effect_index_sif_SPI),
      max = max(effect_index_sif_SPI)
    )
  print(sif_SPI_index_stats)
}
cat("\n")

# ============================================================================
# 9. 空间可视化
# ============================================================================

cat("【空间可视化】\n")

spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude))

cat("有坐标的站点数:", nrow(spatial_data), "\n\n")

# 准备SPI → SIF因果关系的数据
spatial_SPI_sif <- spatial_data %>%
  filter(causality_direction %in% c("SPI → SIF", "双向因果")) %>%
  mutate(
    effect_strength = abs(effect_index_SPI_sif)
  )

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 地图1: 因果方向 + 效应类型（SPI → SIF）
p1 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_SPI_sif,
    aes(x = longitude, y = latitude,
        color = effect_type_SPI_sif,
        shape = causality_direction,
        size = effect_strength),
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "促进效应(+)" = "#4575b4",
      "抑制效应(-)" = "#d73027",
      "效应极弱" = "#fdae61",
      "S-map失败" = "#999999"
    ),
    name = "效应类型"
  ) +
  scale_shape_manual(
    values = c("SPI → SIF" = 16, "双向因果" = 17),
    name = "因果方向"
  ) +
  scale_size_continuous(range = c(1, 5), name = "效应强度\n|∂|") +
  labs(
    title = "SPI对SIF的因果效应空间分布",
    subtitle = paste0("n=", nrow(spatial_SPI_sif),
                      " | 抑制效应: ", sum(spatial_SPI_sif$effect_type_SPI_sif == "抑制效应(-)", na.rm = TRUE),
                      " | 促进效应: ", sum(spatial_SPI_sif$effect_type_SPI_sif == "促进效应(+)", na.rm = TRUE),
                      " | 双向因果: ", sum(spatial_SPI_sif$causality_direction == "双向因果")),
    x = "经度", y = "纬度",
    caption = "圆形=单向因果（SPI→SIF），三角形=双向因果"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p1)

# 地图2: 效应指数连续分布（SPI → SIF）
p2 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_SPI_sif %>% filter(!is.na(effect_index_SPI_sif)),
    aes(x = longitude, y = latitude,
        color = effect_index_SPI_sif,
        shape = causality_direction,
        size = abs(effect_index_SPI_sif)),
    alpha = 0.7
  ) +
  scale_color_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = 0,
    name = "效应指数\n∂"
  ) +
  scale_shape_manual(
    values = c("SPI → SIF" = 16, "双向因果" = 17),
    name = "因果方向"
  ) +
  scale_size_continuous(range = c(1, 5), name = "|∂|") +
  labs(
    title = "SPI → SIF 效应指数的空间梯度",
    subtitle = "蓝色=SPI促进SIF（湿润有利），红色=SPI抑制SIF（干旱有害）",
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

# 地图3: 所有因果方向的分布
p3 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_data %>% filter(causality_direction != "无显著因果"),
    aes(x = longitude, y = latitude,
        color = causality_direction),
    alpha = 0.7,
    size = 2.5
  ) +
  scale_color_manual(
    values = c(
      "SPI → SIF" = "#1b9e77",
      "SIF → SPI" = "#d95f02",
      "双向因果" = "#7570b3"
    ),
    name = "因果方向"
  ) +
  labs(
    title = "SPI与SIF因果方向的空间分布",
    subtitle = paste0("SPI→SIF: ", sum(spatial_data$causality_direction == "SPI → SIF"),
                      " | SIF→SPI: ", sum(spatial_data$causality_direction == "SIF → SPI"),
                      " | 双向: ", sum(spatial_data$causality_direction == "双向因果"),
                      " | 无因果: ", sum(spatial_data$causality_direction == "无显著因果")),
    x = "经度", y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p3)

# ============================================================================
# 最终总结表
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("                   最终总结\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_table <- tibble(
  分类 = c(
    "总站点数",
    "有SPI→SIF因果",
    "  - 单向因果",
    "  - 双向因果",
    "有SIF→SPI因果",
    "  - 单向因果",
    "  - 双向因果",
    "无显著因果",
    "",
    "SPI→SIF效应类型:",
    "  - 抑制效应(-)",
    "  - 促进效应(+)",
    "  - 效应极弱"
  ),
  数量 = c(
    nrow(ccm_results_all),
    nrow(SPI_to_sif_results),
    sum(ccm_results_all$causality_direction == "SPI → SIF"),
    sum(ccm_results_all$causality_direction == "双向因果"),
    nrow(sif_to_SPI_results),
    sum(ccm_results_all$causality_direction == "SIF → SPI"),
    sum(ccm_results_all$causality_direction == "双向因果"),
    sum(ccm_results_all$causality_direction == "无显著因果"),
    "",
    "",
    sum(SPI_to_sif_results$effect_type_SPI_sif == "抑制效应(-)", na.rm = TRUE),
    sum(SPI_to_sif_results$effect_type_SPI_sif == "促进效应(+)", na.rm = TRUE),
    sum(SPI_to_sif_results$effect_type_SPI_sif == "效应极弱", na.rm = TRUE)
  )
)

print(summary_table)

cat("\n分析完成!\n")

# 保存结果
# fwrite(ccm_results_all, "data_proc/ccm_results_SPI_sif.csv")
# fwrite(spatial_data, "data_proc/ccm_results_SPI_spatial.csv")
