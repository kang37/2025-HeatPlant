# ============================================================================
# VPD热事件综合分析（改进版）
# 基于文献阈值（2.0 kPa），双向因果分析，经济关联
# ============================================================================

pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr,
               tidyr, showtext, rEDM, sf, rnaturalearth, rnaturalearthdata,
               patchwork, readxl)
showtext_auto()

cat("\n", rep("=", 70), "\n", sep = "")
cat("        VPD热事件综合分析（文献阈值 + 双向因果）\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================================
# 1. 读取站点SIF数据
# ============================================================================

cat("【1. 读取数据】\n")

meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  tibble() %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 读取站点坐标
station_coords <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

target_stations <- unique(meteo_sif_data$meteo_stat_id)
cat("目标站点数:", length(target_stations), "\n")

# ============================================================================
# 2. 读取每日气象数据
# ============================================================================

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
  filter(month %in% 5:9)

cat("读取完成，数据行数:", nrow(meteo_data_daily), "\n\n")

# ============================================================================
# 3. 计算每日VPD
# ============================================================================

cat("【2. 计算每日VPD】\n")

meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    avp = (rh / 100) * svp,
    vpd = svp - avp
  ) %>%
  mutate(
    vpd = case_when(
      vpd < 0 ~ NA_real_,
      vpd > 10 ~ NA_real_,
      is.na(tavg) | is.na(rh) ~ NA_real_,
      TRUE ~ vpd
    )
  ) %>%
  filter(!is.na(vpd))

cat("VPD数据行数:", nrow(meteo_data_daily_vpd), "\n")
cat("VPD统计: 均值=", round(mean(meteo_data_daily_vpd$vpd, na.rm = TRUE), 3),
    " kPa, 中位=", round(median(meteo_data_daily_vpd$vpd, na.rm = TRUE), 3),
    " kPa\n\n")

# ============================================================================
# 4. 基于文献阈值构建月度VPD热事件指标
# ============================================================================

cat("【3. 构建月度VPD热事件指标】\n")

# 基于文献的绝对阈值（Yuan et al. 2019等）
VPD_THRESHOLD <- 2.0  # kPa
cat("VPD热事件阈值:", VPD_THRESHOLD, "kPa（基于文献）\n")

monthly_heat <- meteo_data_daily_vpd %>%
  mutate(
    vpd_excess = pmax(vpd - VPD_THRESHOLD, 0),
    is_heat_day = vpd > VPD_THRESHOLD
  ) %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(
    n_days = n(),
    vpd_mean = mean(vpd, na.rm = TRUE),
    vpd_max = max(vpd, na.rm = TRUE),
    heat_days = sum(is_heat_day, na.rm = TRUE),
    cumulative_heat_impact = sum(vpd_excess, na.rm = TRUE),
    mean_heat_intensity = mean(vpd_excess[is_heat_day], na.rm = TRUE),
    heat_freq = heat_days / n_days,
    .groups = "drop"
  )

heat_summary <- monthly_heat %>%
  summarise(
    总月数 = n(),
    有热影响月份 = sum(cumulative_heat_impact > 0),
    比例 = round(sum(cumulative_heat_impact > 0) / n() * 100, 1),
    平均累积热影响 = round(mean(cumulative_heat_impact[cumulative_heat_impact > 0]), 3)
  )
print(heat_summary)
cat("\n")

# ============================================================================
# 5. 合并SIF数据并去趋势
# ============================================================================

cat("【4. 合并SIF数据并去趋势】\n")

sif_monthly <- meteo_sif_data %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")

safe_detrend <- function(x, time_idx) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch({
    residuals(lm(x ~ time_idx))
  }, error = function(e) rep(NA_real_, length(x)))
}

data_for_ccm <- monthly_heat %>%
  inner_join(sif_monthly, by = c("meteo_stat_id", "year", "month")) %>%
  filter(!is.na(sif)) %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    sif_detrended = safe_detrend(sif, time_idx),
    heat_impact_detrended = safe_detrend(cumulative_heat_impact, time_idx)
  ) %>%
  ungroup() %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(sif_detrended)) >= 10,
         sum(!is.na(heat_impact_detrended)) >= 10) %>%
  ungroup()

cat("可用站点:", length(unique(data_for_ccm$meteo_stat_id)), "\n")
cat("数据行数:", nrow(data_for_ccm), "\n\n")

# ============================================================================
# 6. CCM分析函数（双向因果版）
# ============================================================================

perform_ccm_heat_sif <- function(station_id, data, min_points = 30) {

  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(
      sif = sif_detrended,
      heat_impact = heat_impact_detrended
    ) %>%
    filter(!is.na(sif), !is.na(heat_impact)) %>%
    as.data.frame()

  n_data <- nrow(station_data)
  if (n_data < min_points) return(NULL)

  tryCatch({

    station_data_ccm <- data.frame(
      time = 1:n_data,
      sif = station_data$sif,
      heat_impact = station_data$heat_impact
    )

    lib_pred_size <- min(100, n_data)

    # 确定最优E
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data / 10)),
      showPlot = FALSE
    )

    embed_heat <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "heat_impact",
      target = "heat_impact",
      maxE = min(8, floor(n_data / 10)),
      showPlot = FALSE
    )

    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    best_E_heat <- embed_heat$E[which.max(embed_heat$rho)]
    E_ccm <- max(best_E_sif, best_E_heat)

    # CCM参数
    tau <- 1
    tp <- 0
    embedding_loss <- (E_ccm - 1) * tau
    max_available_lib <- n_data - embedding_loss - tp

    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib

    if (lib_end <= lib_start || lib_end < 15) return(NULL)

    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)

    # ========== CCM: 热影响 → SIF ==========
    ccm_heat_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm, Tp = tp,
      columns = "heat_impact", target = "sif",
      libSizes = libSizes_str,
      sample = 50, random = TRUE, showPlot = FALSE
    )

    ccm_summary_heat_sif <- ccm_heat_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(rho_mean = mean(`heat_impact:sif`, na.rm = TRUE), .groups = "drop")

    final_rho_heat_sif <- ccm_summary_heat_sif$rho_mean[nrow(ccm_summary_heat_sif)]
    trend_heat_sif <- cor(ccm_summary_heat_sif$lib_size, ccm_summary_heat_sif$rho_mean)

    # ========== CCM: SIF → 热影响 ==========
    ccm_sif_to_heat <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm, Tp = tp,
      columns = "sif", target = "heat_impact",
      libSizes = libSizes_str,
      sample = 50, random = TRUE, showPlot = FALSE
    )

    ccm_summary_sif_heat <- ccm_sif_to_heat %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(rho_mean = mean(`sif:heat_impact`, na.rm = TRUE), .groups = "drop")

    final_rho_sif_heat <- ccm_summary_sif_heat$rho_mean[nrow(ccm_summary_sif_heat)]
    trend_sif_heat <- cor(ccm_summary_sif_heat$lib_size, ccm_summary_sif_heat$rho_mean)

    # ========== S-map: 热影响 → SIF ==========
    smap_heat_sif <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm, theta = 2,
      columns = "heat_impact", target = "sif",
      embedded = FALSE
    )

    smap_coeffs_heat <- smap_heat_sif$coefficients

    if (!is.null(smap_coeffs_heat) && ncol(smap_coeffs_heat) >= 2) {
      heat_coef_col <- which(grepl("heat", colnames(smap_coeffs_heat), ignore.case = TRUE))
      if (length(heat_coef_col) > 0) {
        smap_coef_heat <- smap_coeffs_heat[, heat_coef_col[1]]
      } else {
        smap_coef_heat <- smap_coeffs_heat[, 2]
      }
      smap_coef_heat <- smap_coef_heat[!is.na(smap_coef_heat)]
      mean_smap_coef_heat_sif <- mean(smap_coef_heat, na.rm = TRUE)
      sd_smap_coef_heat_sif <- sd(smap_coef_heat, na.rm = TRUE)
    } else {
      mean_smap_coef_heat_sif <- NA
      sd_smap_coef_heat_sif <- NA
      smap_coef_heat <- NA
    }

    # ========== S-map: SIF → 热影响 ==========
    smap_sif_heat <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm, theta = 2,
      columns = "sif", target = "heat_impact",
      embedded = FALSE
    )

    smap_coeffs_sif <- smap_sif_heat$coefficients

    if (!is.null(smap_coeffs_sif) && ncol(smap_coeffs_sif) >= 2) {
      sif_coef_col <- which(grepl("sif", colnames(smap_coeffs_sif), ignore.case = TRUE))
      if (length(sif_coef_col) > 0) {
        smap_coef_sif <- smap_coeffs_sif[, sif_coef_col[1]]
      } else {
        smap_coef_sif <- smap_coeffs_sif[, 2]
      }
      smap_coef_sif <- smap_coef_sif[!is.na(smap_coef_sif)]
      mean_smap_coef_sif_heat <- mean(smap_coef_sif, na.rm = TRUE)
      sd_smap_coef_sif_heat <- sd(smap_coef_sif, na.rm = TRUE)
    } else {
      mean_smap_coef_sif_heat <- NA
      sd_smap_coef_sif_heat <- NA
      smap_coef_sif <- NA
    }

    # ========== 判断因果方向 ==========
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0

    heat_causes_sif <- (final_rho_heat_sif > ccm_threshold_rho &
                          trend_heat_sif > ccm_threshold_trend)
    sif_causes_heat <- (final_rho_sif_heat > ccm_threshold_rho &
                          trend_sif_heat > ccm_threshold_trend)

    causality_direction <- case_when(
      heat_causes_sif & !sif_causes_heat ~ "VPD热影响 → SIF",
      !heat_causes_sif & sif_causes_heat ~ "SIF → VPD热影响",
      heat_causes_sif & sif_causes_heat ~ "双向因果",
      TRUE ~ "无显著因果"
    )

    # 效应类型（VPD热影响 → SIF）
    effect_type_heat_sif <- case_when(
      !heat_causes_sif ~ "无因果",
      is.na(mean_smap_coef_heat_sif) ~ "S-map失败",
      mean_smap_coef_heat_sif > 0.001 ~ "促进效应(+)",
      mean_smap_coef_heat_sif < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )

    # 效应类型（SIF → VPD热影响）
    effect_type_sif_heat <- case_when(
      !sif_causes_heat ~ "无因果",
      is.na(mean_smap_coef_sif_heat) ~ "S-map失败",
      mean_smap_coef_sif_heat > 0.001 ~ "促进效应(+)",
      mean_smap_coef_sif_heat < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )

    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,

      # CCM结果
      rho_heat_to_sif = final_rho_heat_sif,
      trend_heat_to_sif = trend_heat_sif,
      rho_sif_to_heat = final_rho_sif_heat,
      trend_sif_to_heat = trend_sif_heat,

      # 因果判断
      heat_causes_sif = heat_causes_sif,
      sif_causes_heat = sif_causes_heat,
      causality_direction = causality_direction,

      # S-map效应指数（VPD热影响 → SIF）
      effect_index_heat_sif = mean_smap_coef_heat_sif,
      effect_index_sd_heat_sif = sd_smap_coef_heat_sif,
      effect_type_heat_sif = effect_type_heat_sif,

      # S-map效应指数（SIF → VPD热影响）
      effect_index_sif_heat = mean_smap_coef_sif_heat,
      effect_index_sd_sif_heat = sd_smap_coef_sif_heat,
      effect_type_sif_heat = effect_type_sif_heat,

      # 保存详细数据
      smap_coefficients_heat_sif = list(smap_coef_heat),
      smap_coefficients_sif_heat = list(smap_coef_sif)
    )

  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# ============================================================================
# 7. 批量分析
# ============================================================================

cat("\n【5. 开始CCM分析】\n")

all_stations <- data_for_ccm %>%
  count(meteo_stat_id) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("分析站点数:", length(all_stations), "\n\n")

ccm_results_all <- map_dfr(seq_along(all_stations), function(i) {
  if (i %% 20 == 0) cat("已完成:", i, "/", length(all_stations), "\n")
  perform_ccm_heat_sif(all_stations[i], data_for_ccm, min_points = 30)
})

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

# 2. VPD热影响 → SIF 的效应类型统计
cat("【2. VPD热影响 → SIF 效应类型统计】\n")
heat_to_sif_results <- ccm_results_all %>%
  filter(causality_direction %in% c("VPD热影响 → SIF", "双向因果"))

if (nrow(heat_to_sif_results) > 0) {
  effect_stats_heat_sif <- table(heat_to_sif_results$effect_type_heat_sif)
  print(effect_stats_heat_sif)
  cat("\n百分比:\n")
  print(round(prop.table(effect_stats_heat_sif) * 100, 1))
} else {
  cat("无 VPD热影响 → SIF 因果关系的站点\n")
}
cat("\n")

# 3. SIF → VPD热影响 的效应类型统计
cat("【3. SIF → VPD热影响 效应类型统计】\n")
sif_to_heat_results <- ccm_results_all %>%
  filter(causality_direction %in% c("SIF → VPD热影响", "双向因果"))

if (nrow(sif_to_heat_results) > 0) {
  effect_stats_sif_heat <- table(sif_to_heat_results$effect_type_sif_heat)
  print(effect_stats_sif_heat)
  cat("\n百分比:\n")
  print(round(prop.table(effect_stats_sif_heat) * 100, 1))
} else {
  cat("无 SIF → VPD热影响 因果关系的站点\n")
}
cat("\n")

# 4. 效应指数统计
cat("【4. 效应指数统计】\n\n")

cat("VPD热影响 → SIF 方向:\n")
heat_sif_index_stats <- heat_to_sif_results %>%
  filter(!is.na(effect_index_heat_sif)) %>%
  summarise(
    n = n(),
    mean = mean(effect_index_heat_sif),
    sd = sd(effect_index_heat_sif),
    median = median(effect_index_heat_sif),
    min = min(effect_index_heat_sif),
    max = max(effect_index_heat_sif),
    n_positive = sum(effect_index_heat_sif > 0.001),
    n_negative = sum(effect_index_heat_sif < -0.001),
    pct_negative = round(sum(effect_index_heat_sif < -0.001) / n() * 100, 1)
  )
print(heat_sif_index_stats)
cat("\n")

if (nrow(sif_to_heat_results) > 0) {
  cat("SIF → VPD热影响 方向:\n")
  sif_heat_index_stats <- sif_to_heat_results %>%
    filter(!is.na(effect_index_sif_heat)) %>%
    summarise(
      n = n(),
      mean = mean(effect_index_sif_heat),
      sd = sd(effect_index_sif_heat),
      median = median(effect_index_sif_heat),
      min = min(effect_index_sif_heat),
      max = max(effect_index_sif_heat)
    )
  print(sif_heat_index_stats)
}
cat("\n")

# ============================================================================
# 9. 空间可视化
# ============================================================================

cat("【6. 空间可视化】\n")

spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude))

cat("有坐标的站点数:", nrow(spatial_data), "\n\n")

# 准备VPD热影响 → SIF因果关系的数据
spatial_heat_sif <- spatial_data %>%
  filter(causality_direction %in% c("VPD热影响 → SIF", "双向因果")) %>%
  mutate(effect_strength = abs(effect_index_heat_sif))

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 地图1: VPD热影响对SIF的因果效应（圆形=单向，三角形=双向）
p1 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_heat_sif,
    aes(x = longitude, y = latitude,
        color = effect_type_heat_sif,
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
    values = c("VPD热影响 → SIF" = 16, "双向因果" = 17),
    name = "因果方向"
  ) +
  scale_size_continuous(range = c(1, 5), name = "效应强度\n|∂|") +
  labs(
    title = "VPD热事件对SIF的因果效应空间分布",
    subtitle = paste0("VPD阈值 = ", VPD_THRESHOLD, " kPa | n=", nrow(spatial_heat_sif),
                      " | 抑制: ", sum(spatial_heat_sif$effect_type_heat_sif == "抑制效应(-)", na.rm = TRUE),
                      " | 促进: ", sum(spatial_heat_sif$effect_type_heat_sif == "促进效应(+)", na.rm = TRUE),
                      " | 双向: ", sum(spatial_heat_sif$causality_direction == "双向因果")),
    x = "经度", y = "纬度",
    caption = "圆形=单向因果（VPD热影响→SIF），三角形=双向因果"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p1)

# 地图2: 效应指数连续色标
p2 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_heat_sif %>% filter(!is.na(effect_index_heat_sif)),
    aes(x = longitude, y = latitude,
        color = effect_index_heat_sif,
        shape = causality_direction,
        size = abs(effect_index_heat_sif)),
    alpha = 0.7
  ) +
  scale_color_gradient2(
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = 0, name = "效应指数\n∂"
  ) +
  scale_shape_manual(
    values = c("VPD热影响 → SIF" = 16, "双向因果" = 17),
    name = "因果方向"
  ) +
  scale_size_continuous(range = c(1, 5), name = "|∂|") +
  labs(
    title = "VPD热影响 → SIF 效应指数的空间梯度",
    subtitle = "蓝色=热事件促进SIF，红色=热事件抑制SIF",
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

# 地图3: 所有因果方向分布
p3 <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = spatial_data %>% filter(causality_direction != "无显著因果"),
    aes(x = longitude, y = latitude,
        color = causality_direction),
    alpha = 0.7, size = 2.5
  ) +
  scale_color_manual(
    values = c(
      "VPD热影响 → SIF" = "#1b9e77",
      "SIF → VPD热影响" = "#d95f02",
      "双向因果" = "#7570b3"
    ),
    name = "因果方向"
  ) +
  labs(
    title = "VPD热事件与SIF因果方向的空间分布",
    subtitle = paste0(
      "VPD热影响→SIF: ", sum(spatial_data$causality_direction == "VPD热影响 → SIF"),
      " | SIF→VPD热影响: ", sum(spatial_data$causality_direction == "SIF → VPD热影响"),
      " | 双向: ", sum(spatial_data$causality_direction == "双向因果"),
      " | 无因果: ", sum(spatial_data$causality_direction == "无显著因果")),
    x = "经度", y = "纬度"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p3)

# ============================================================================
# 10. 经济数据关联分析
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("    VPD热事件因果关系与经济发展指标的关系\n")
cat(rep("=", 70), "\n\n", sep = "")

# 10.1 读取经济数据
cat("读取经济和空间数据...\n")

china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326)

gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
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

economic_data <- station_coords %>%
  left_join(stations_with_city, by = "meteo_stat_id") %>%
  left_join(gdp_data, by = "city_name") %>%
  select(meteo_stat_id, longitude, latitude, city_name, province_name,
         gdp_per_capita, gdp_total)

spatial_econ <- spatial_data %>%
  left_join(economic_data, by = c("meteo_stat_id", "longitude", "latitude"))

cat("有GDP数据的站点:", sum(!is.na(spatial_econ$gdp_per_capita)), "\n\n")

# 10.2 创建因果标签（聚焦VPD热影响 → SIF方向）
spatial_econ <- spatial_econ %>%
  mutate(
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      causality_direction == "SIF → VPD热影响" ~ "SIF→VPD热影响",
      # 以下两类都包含 heat→SIF 方向
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "VPD热影响 → SIF" ~ "热事件促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "VPD热影响 → SIF" ~ "热事件抑制SIF",
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "双向因果" ~ "双向-促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "双向因果" ~ "双向-抑制SIF",
      TRUE ~ "其他"
    )
  )

# 10.3 人均GDP分析
cat("【人均GDP分析】\n\n")

spatial_econ_valid <- spatial_econ %>%
  filter(!is.na(gdp_per_capita))

if (nrow(spatial_econ_valid) > 0) {

  # 统计表
  cat("因果类型的人均GDP统计:\n")
  gdp_stats <- spatial_econ_valid %>%
    group_by(heat_causality) %>%
    summarise(
      站点数 = n(),
      平均GDP_万元 = round(mean(gdp_per_capita) / 10000, 2),
      中位GDP_万元 = round(median(gdp_per_capita) / 10000, 2),
      GDP标准差_万元 = round(sd(gdp_per_capita) / 10000, 2),
      .groups = "drop"
    ) %>%
    arrange(desc(站点数))

  print(gdp_stats)
  cat("\n")

  # GDP分级
  spatial_econ_cat <- spatial_econ_valid %>%
    mutate(
      gdp_category = cut(
        gdp_per_capita,
        breaks = c(0, 40000, 80000, 120000, Inf),
        labels = c("低(<4万)", "中(4-8万)", "较高(8-12万)", "高(>12万)"),
        include.lowest = TRUE
      )
    )

  # 堆积条形图
  p_bar_percapita <- spatial_econ_cat %>%
    filter(heat_causality != "其他", !is.na(gdp_category)) %>%
    count(gdp_category, heat_causality) %>%
    ggplot(aes(x = gdp_category, y = n, fill = heat_causality)) +
    geom_col(position = "fill") +
    scale_fill_manual(
      values = c(
        "热事件促进SIF" = "#4575b4",
        "热事件抑制SIF" = "#d73027",
        "双向-促进SIF" = "#91bfdb",
        "双向-抑制SIF" = "#fc8d59",
        "SIF→VPD热影响" = "#984ea3",
        "无因果" = "#999999"
      ),
      name = "因果类型"
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "不同人均GDP水平下的VPD热事件因果类型分布",
      x = "人均GDP类别", y = "比例 (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )

  print(p_bar_percapita)

  # 箱线图
  p_box_percapita <- spatial_econ_valid %>%
    filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                                  "双向-促进SIF", "双向-抑制SIF", "无因果")) %>%
    mutate(gdp_per_capita_万元 = gdp_per_capita / 10000) %>%
    ggplot(aes(x = reorder(heat_causality, -gdp_per_capita_万元, FUN = median),
               y = gdp_per_capita_万元,
               fill = heat_causality)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
    scale_fill_manual(
      values = c(
        "热事件促进SIF" = "#4575b4",
        "热事件抑制SIF" = "#d73027",
        "双向-促进SIF" = "#91bfdb",
        "双向-抑制SIF" = "#fc8d59",
        "无因果" = "#999999"
      ),
      guide = "none"
    ) +
    labs(
      title = "VPD热事件因果类型与人均GDP的关系",
      x = "因果类型", y = "人均GDP (万元)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 9)
    )

  print(p_box_percapita)

  # 统计检验
  cat("统计检验（人均GDP）:\n")
  test_data <- spatial_econ_valid %>%
    filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF", "无因果"))

  if (nrow(test_data) > 10 && length(unique(test_data$heat_causality)) >= 2) {
    kw_result <- kruskal.test(gdp_per_capita ~ heat_causality, data = test_data)
    cat("  Kruskal-Wallis检验: χ² =", round(kw_result$statistic, 3),
        ", p =", format.pval(kw_result$p.value, digits = 3), "\n")

    # 两两比较（如有3组以上）
    if (length(unique(test_data$heat_causality)) >= 2) {
      pw_result <- pairwise.wilcox.test(
        test_data$gdp_per_capita,
        test_data$heat_causality,
        p.adjust.method = "BH"
      )
      cat("  两两比较（Wilcoxon, BH校正）:\n")
      print(pw_result$p.value)
    }
  }
  cat("\n")
}

# 10.4 总GDP分析
cat("【总GDP分析】\n\n")

spatial_econ_total <- spatial_econ %>%
  filter(!is.na(gdp_total))

if (nrow(spatial_econ_total) > 0) {

  # 统计表
  cat("因果类型的总GDP统计:\n")
  gdp_total_stats <- spatial_econ_total %>%
    group_by(heat_causality) %>%
    summarise(
      站点数 = n(),
      平均GDP_亿元 = round(mean(gdp_total), 2),
      中位GDP_亿元 = round(median(gdp_total), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(站点数))

  print(gdp_total_stats)
  cat("\n")

  # GDP分级（四分位数）
  gdp_total_quantiles <- quantile(spatial_econ_total$gdp_total,
                                  probs = c(0, 0.25, 0.5, 0.75, 1),
                                  na.rm = TRUE)

  spatial_econ_total_cat <- spatial_econ_total %>%
    mutate(
      gdp_category = cut(
        gdp_total,
        breaks = gdp_total_quantiles,
        labels = c("低(Q1)", "中(Q2)", "较高(Q3)", "高(Q4)"),
        include.lowest = TRUE
      )
    )

  # 堆积条形图
  p_bar_total <- spatial_econ_total_cat %>%
    filter(heat_causality != "其他", !is.na(gdp_category)) %>%
    count(gdp_category, heat_causality) %>%
    ggplot(aes(x = gdp_category, y = n, fill = heat_causality)) +
    geom_col(position = "fill") +
    scale_fill_manual(
      values = c(
        "热事件促进SIF" = "#4575b4",
        "热事件抑制SIF" = "#d73027",
        "双向-促进SIF" = "#91bfdb",
        "双向-抑制SIF" = "#fc8d59",
        "SIF→VPD热影响" = "#984ea3",
        "无因果" = "#999999"
      ),
      name = "因果类型"
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "不同总GDP水平下的VPD热事件因果类型分布",
      x = "总GDP类别（四分位数）", y = "比例 (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )

  print(p_bar_total)

  # 箱线图
  p_box_total <- spatial_econ_total %>%
    filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                                  "双向-促进SIF", "双向-抑制SIF", "无因果")) %>%
    ggplot(aes(x = reorder(heat_causality, -gdp_total, FUN = median),
               y = gdp_total,
               fill = heat_causality)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
    scale_fill_manual(
      values = c(
        "热事件促进SIF" = "#4575b4",
        "热事件抑制SIF" = "#d73027",
        "双向-促进SIF" = "#91bfdb",
        "双向-抑制SIF" = "#fc8d59",
        "无因果" = "#999999"
      ),
      guide = "none"
    ) +
    labs(
      title = "VPD热事件因果类型与总GDP的关系",
      x = "因果类型", y = "总GDP (亿元)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 9)
    )

  print(p_box_total)

  # 统计检验
  cat("统计检验（总GDP）:\n")
  test_data_total <- spatial_econ_total %>%
    filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF", "无因果"))

  if (nrow(test_data_total) > 10 && length(unique(test_data_total$heat_causality)) >= 2) {
    kw_result_total <- kruskal.test(gdp_total ~ heat_causality, data = test_data_total)
    cat("  Kruskal-Wallis检验: χ² =", round(kw_result_total$statistic, 3),
        ", p =", format.pval(kw_result_total$p.value, digits = 3), "\n")

    if (length(unique(test_data_total$heat_causality)) >= 2) {
      pw_result_total <- pairwise.wilcox.test(
        test_data_total$gdp_total,
        test_data_total$heat_causality,
        p.adjust.method = "BH"
      )
      cat("  两两比较（Wilcoxon, BH校正）:\n")
      print(pw_result_total$p.value)
    }
  }
  cat("\n")
}

# 10.5 GDP空间地图
cat("【GDP空间分布图】\n")

if (nrow(spatial_econ_valid) > 0) {
  p_gdp_spatial <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70",
            linewidth = 0.3) +
    geom_point(
      data = spatial_econ_valid %>%
        filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                                      "双向-促进SIF", "双向-抑制SIF")),
      aes(x = longitude, y = latitude,
          size = gdp_per_capita / 10000,
          color = heat_causality),
      alpha = 0.7
    ) +
    scale_color_manual(
      values = c(
        "热事件促进SIF" = "#4575b4",
        "热事件抑制SIF" = "#d73027",
        "双向-促进SIF" = "#91bfdb",
        "双向-抑制SIF" = "#fc8d59"
      ),
      name = "因果类型"
    ) +
    scale_size_continuous(range = c(1, 10), name = "人均GDP\n(万元)") +
    labs(
      title = "VPD热事件因果关系与GDP的空间分布",
      subtitle = "点的大小=GDP水平，颜色=因果类型",
      x = "经度", y = "纬度"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      panel.background = element_rect(fill = "aliceblue", color = NA)
    )

  print(p_gdp_spatial)
}

# ============================================================================
# 最终总结表
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("                   最终总结\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_table <- tibble(
  分类 = c(
    "总站点数",
    "有VPD热影响→SIF因果",
    "  - 单向因果",
    "  - 双向因果",
    "有SIF→VPD热影响因果",
    "  - 单向因果",
    "  - 双向因果",
    "无显著因果",
    "",
    "VPD热影响→SIF效应类型:",
    "  - 抑制效应(-)",
    "  - 促进效应(+)",
    "  - 效应极弱"
  ),
  数量 = c(
    nrow(ccm_results_all),
    nrow(heat_to_sif_results),
    sum(ccm_results_all$causality_direction == "VPD热影响 → SIF"),
    sum(ccm_results_all$causality_direction == "双向因果"),
    nrow(sif_to_heat_results),
    sum(ccm_results_all$causality_direction == "SIF → VPD热影响"),
    sum(ccm_results_all$causality_direction == "双向因果"),
    sum(ccm_results_all$causality_direction == "无显著因果"),
    "",
    "",
    sum(heat_to_sif_results$effect_type_heat_sif == "抑制效应(-)", na.rm = TRUE),
    sum(heat_to_sif_results$effect_type_heat_sif == "促进效应(+)", na.rm = TRUE),
    sum(heat_to_sif_results$effect_type_heat_sif == "效应极弱", na.rm = TRUE)
  )
)

print(summary_table)

cat("\n分析完成!\n")

# 保存结果
# fwrite(ccm_results_all %>% select(-starts_with("smap_coefficients")),
#        "data_proc/ccm_results_vpd_heat_sif.csv")
