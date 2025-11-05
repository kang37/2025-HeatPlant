pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, tidyr, showtext, rEDM)
showtext_auto()

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

# CCM ----
pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, 
               tidyr, rEDM, sf, rnaturalearth, rnaturalearthdata)

# [之前的数据读取代码保持不变...]

# ============================================================================
# 完整的 CCM + S-map 分析
# ============================================================================

perform_ccm_with_smap <- function(station_id, data, min_points = 30) {
  
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
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif
    )
    
    lib_pred_size <- min(100, n_data)
    
    # 1. 确定最优 E
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
    
    # 2. CCM 分析
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
    
    # Temp → SIF 方向
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
    
    # SIF → Temp 方向
    ccm_sif_to_temp <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "sif",
      target = "temp_anomaly",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # CCM 结果汇总
    ccm_summary_temp_sif <- ccm_temp_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(direction = "Temp → SIF")
    
    ccm_summary_sif_temp <- ccm_sif_to_temp %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`sif:temp_anomaly`, na.rm = TRUE),
        rho_sd = sd(`sif:temp_anomaly`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(direction = "SIF → Temp")
    
    final_rho_temp_to_sif <- ccm_summary_temp_sif$rho_mean[nrow(ccm_summary_temp_sif)]
    final_rho_sif_to_temp <- ccm_summary_sif_temp$rho_mean[nrow(ccm_summary_sif_temp)]
    
    trend_temp_to_sif <- cor(ccm_summary_temp_sif$lib_size, 
                             ccm_summary_temp_sif$rho_mean)
    trend_sif_to_temp <- cor(ccm_summary_sif_temp$lib_size, 
                             ccm_summary_sif_temp$rho_mean)
    
    # 3. S-map 分析（Temp → SIF）
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
    
    # 提取 S-map 系数
    smap_coeffs <- smap_temp_to_sif$coefficients
    
    if (!is.null(smap_coeffs) && ncol(smap_coeffs) >= 2) {
      # 找到温度对应的列
      temp_coef_col <- which(grepl("temp", colnames(smap_coeffs), ignore.case = TRUE))
      
      if (length(temp_coef_col) > 0) {
        smap_coef_temp <- smap_coeffs[, temp_coef_col[1]]
      } else {
        # 假设第2列是自变量系数（第1列通常是截距）
        smap_coef_temp <- smap_coeffs[, 2]
      }
      
      smap_coef_temp <- smap_coef_temp[!is.na(smap_coef_temp)]
      
      mean_smap_coef <- mean(smap_coef_temp, na.rm = TRUE)
      sd_smap_coef <- sd(smap_coef_temp, na.rm = TRUE)
      median_smap_coef <- median(smap_coef_temp, na.rm = TRUE)
      
    } else {
      mean_smap_coef <- NA
      sd_smap_coef <- NA
      median_smap_coef <- NA
      smap_coef_temp <- NA
    }
    
    # S-map 分析（SIF → Temp）
    smap_sif_to_temp <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "sif",
      target = "temp_anomaly",
      embedded = FALSE
    )
    
    smap_coeffs_sif <- smap_sif_to_temp$coefficients
    
    if (!is.null(smap_coeffs_sif) && ncol(smap_coeffs_sif) >= 2) {
      sif_coef_col <- which(grepl("sif", colnames(smap_coeffs_sif), ignore.case = TRUE))
      
      if (length(sif_coef_col) > 0) {
        smap_coef_sif <- smap_coeffs_sif[, sif_coef_col[1]]
      } else {
        smap_coef_sif <- smap_coeffs_sif[, 2]
      }
      
      smap_coef_sif <- smap_coef_sif[!is.na(smap_coef_sif)]
      
      mean_smap_coef_sif <- mean(smap_coef_sif, na.rm = TRUE)
      sd_smap_coef_sif <- sd(smap_coef_sif, na.rm = TRUE)
      
    } else {
      mean_smap_coef_sif <- NA
      sd_smap_coef_sif <- NA
      smap_coef_sif <- NA
    }
    
    # 4. 判断因果关系
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    temp_causes_sif <- (final_rho_temp_to_sif > ccm_threshold_rho & 
                          trend_temp_to_sif > ccm_threshold_trend)
    sif_causes_temp <- (final_rho_sif_to_temp > ccm_threshold_rho & 
                          trend_sif_to_temp > ccm_threshold_trend)
    
    # 因果方向
    causality_direction <- case_when(
      temp_causes_sif & !sif_causes_temp ~ "Temp → SIF",
      !temp_causes_sif & sif_causes_temp ~ "SIF → Temp",
      temp_causes_sif & sif_causes_temp ~ "双向因果",
      TRUE ~ "无显著因果"
    )
    
    # 效应类型（仅针对 Temp → SIF）
    effect_type_temp_sif <- case_when(
      !temp_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0.001 ~ "促进效应(+)",
      mean_smap_coef < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    # 效应类型（仅针对 SIF → Temp）
    effect_type_sif_temp <- case_when(
      !sif_causes_temp ~ "无因果",
      is.na(mean_smap_coef_sif) ~ "S-map失败",
      mean_smap_coef_sif > 0.001 ~ "促进效应(+)",
      mean_smap_coef_sif < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      E = E_ccm,
      
      # CCM 结果
      rho_temp_to_sif = final_rho_temp_to_sif,
      trend_temp_to_sif = trend_temp_to_sif,
      rho_sif_to_temp = final_rho_sif_to_temp,
      trend_sif_to_temp = trend_sif_to_temp,
      
      # 因果判断
      temp_causes_sif = temp_causes_sif,
      sif_causes_temp = sif_causes_temp,
      causality_direction = causality_direction,
      
      # S-map 效应指数（Temp → SIF）
      effect_index_temp_sif = mean_smap_coef,      # 效应指数
      effect_index_sd_temp_sif = sd_smap_coef,
      effect_index_median_temp_sif = median_smap_coef,
      effect_type_temp_sif = effect_type_temp_sif,
      
      # S-map 效应指数（SIF → Temp）
      effect_index_sif_temp = mean_smap_coef_sif,
      effect_index_sd_sif_temp = sd_smap_coef_sif,
      effect_type_sif_temp = effect_type_sif_temp,
      
      # 保存详细数据
      smap_coefficients_temp_sif = list(smap_coef_temp),
      smap_coefficients_sif_temp = list(smap_coef_sif),
      ccm_summary_data = list(bind_rows(ccm_summary_temp_sif, ccm_summary_sif_temp))
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# ============================================================================
# 批量分析
# ============================================================================

cat("\n开始批量CCM+S-map分析...\n\n")

all_stations <- data_for_ccm_linear %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("准备分析", length(all_stations), "个站点...\n\n")

ccm_results_all <- map_dfr(seq_along(all_stations), function(i) {
  if (i %% 10 == 0) {
    cat("已完成:", i, "/", length(all_stations), "\n")
  }
  perform_ccm_with_smap(all_stations[i], data_for_ccm_linear, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_all), "/", length(all_stations), "\n\n")

# ============================================================================
# 详细统计
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

# 2. Temp → SIF 的效应类型统计
cat("【2. Temp → SIF 效应类型统计】\n")
temp_to_sif_results <- ccm_results_all %>%
  filter(causality_direction %in% c("Temp → SIF", "双向因果"))

effect_stats_temp_sif <- table(temp_to_sif_results$effect_type_temp_sif)
print(effect_stats_temp_sif)
cat("\n百分比:\n")
print(round(prop.table(effect_stats_temp_sif) * 100, 1))
cat("\n")

# 3. SIF → Temp 的效应类型统计
cat("【3. SIF → Temp 效应类型统计】\n")
sif_to_temp_results <- ccm_results_all %>%
  filter(causality_direction %in% c("SIF → Temp", "双向因果"))

if (nrow(sif_to_temp_results) > 0) {
  effect_stats_sif_temp <- table(sif_to_temp_results$effect_type_sif_temp)
  print(effect_stats_sif_temp)
  cat("\n百分比:\n")
  print(round(prop.table(effect_stats_sif_temp) * 100, 1))
} else {
  cat("无 SIF → Temp 因果关系的站点\n")
}
cat("\n")

# 4. 因果方向 × 效应类型交叉表
cat("【4. 因果方向 × 效应类型交叉表（Temp → SIF）】\n")
cross_table <- temp_to_sif_results %>%
  count(causality_direction, effect_type_temp_sif) %>%
  pivot_wider(names_from = effect_type_temp_sif, 
              values_from = n, 
              values_fill = 0)

print(cross_table)
cat("\n")

# 5. 效应指数统计
cat("【5. 效应指数（S-map系数）统计】\n\n")

cat("Temp → SIF 方向:\n")
temp_sif_index_stats <- temp_to_sif_results %>%
  filter(!is.na(effect_index_temp_sif)) %>%
  summarise(
    n = n(),
    mean = mean(effect_index_temp_sif),
    sd = sd(effect_index_temp_sif),
    median = median(effect_index_temp_sif),
    min = min(effect_index_temp_sif),
    max = max(effect_index_temp_sif),
    n_positive = sum(effect_index_temp_sif > 0.001),
    n_negative = sum(effect_index_temp_sif < -0.001),
    pct_negative = round(sum(effect_index_temp_sif < -0.001) / n() * 100, 1)
  )

print(temp_sif_index_stats)
cat("\n")

if (nrow(sif_to_temp_results) > 0) {
  cat("SIF → Temp 方向:\n")
  sif_temp_index_stats <- sif_to_temp_results %>%
    filter(!is.na(effect_index_sif_temp)) %>%
    summarise(
      n = n(),
      mean = mean(effect_index_sif_temp),
      sd = sd(effect_index_sif_temp),
      median = median(effect_index_sif_temp),
      min = min(effect_index_sif_temp),
      max = max(effect_index_sif_temp)
    )
  
  print(sif_temp_index_stats)
}
cat("\n")

# ============================================================================
# 可视化
# ============================================================================

# 效应指数分布图
p_effect_index <- ggplot(temp_to_sif_results %>% filter(!is.na(effect_index_temp_sif)),
                         aes(x = effect_index_temp_sif, fill = effect_type_temp_sif)) +
  geom_histogram(bins = 40, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = mean(temp_to_sif_results$effect_index_temp_sif, na.rm = TRUE),
             linetype = "solid", color = "red", linewidth = 1) +
  scale_fill_manual(
    values = c("促进效应(+)" = "#4575b4",
               "抑制效应(-)" = "#d73027",
               "效应极弱" = "#999999",
               "S-map失败" = "#fee090")
  ) +
  labs(
    title = "效应指数分布（Temp → SIF）",
    subtitle = paste0("平均效应指数 = ", 
                      round(temp_sif_index_stats$mean, 4),
                      " | 抑制效应占比 = ",
                      temp_sif_index_stats$pct_negative, "%"),
    x = "效应指数（S-map系数 ∂）",
    y = "站点数",
    fill = "效应类型",
    caption = "负值=抑制效应（温度↑ SIF↓），正值=促进效应（温度↑ SIF↑）"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_effect_index)

# ============================================================================
# 空间分布地图
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

# Temp → SIF 因果关系的站点
spatial_temp_sif <- spatial_data %>%
  filter(causality_direction %in% c("Temp → SIF", "双向因果")) %>%
  mutate(
    # 组合标签
    combined_label = paste0(
      case_when(
        causality_direction == "Temp → SIF" ~ "单向",
        causality_direction == "双向因果" ~ "双向",
        TRUE ~ ""
      ),
      " - ",
      effect_type_temp_sif
    ),
    # 效应强度（用于点的大小）
    effect_strength = abs(effect_index_temp_sif)
  )

# 获取中国地图
china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 地图边界
bbox <- st_bbox(c(
  xmin = min(spatial_temp_sif$lon, na.rm = TRUE) - 2,
  xmax = max(spatial_temp_sif$lon, na.rm = TRUE) + 2,
  ymin = min(spatial_temp_sif$lat, na.rm = TRUE) - 2,
  ymax = max(spatial_temp_sif$lat, na.rm = TRUE) + 2
))

# 主地图：因果方向 + 效应类型
p_spatial_main <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = spatial_temp_sif,
             aes(x = longitude, y = latitude,
                 color = effect_type_temp_sif,
                 shape = causality_direction,
                 size = effect_strength),
             alpha = 0.7) +
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
    values = c("Temp → SIF" = 16, "双向因果" = 17),
    name = "因果方向"
  ) +
  scale_size_continuous(
    range = c(1, 5),
    name = "效应强度\n|∂|"
  ) +
  labs(
    title = "热月温度对SIF的因果效应空间分布",
    subtitle = paste0("n=", nrow(spatial_temp_sif), 
                      " | 抑制效应: ", sum(spatial_temp_sif$effect_type_temp_sif == "抑制效应(-)", na.rm = TRUE),
                      " | 促进效应: ", sum(spatial_temp_sif$effect_type_temp_sif == "促进效应(+)", na.rm = TRUE)),
    x = "经度",
    y = "纬度",
    caption = "圆形=单向因果，三角形=双向因果"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )

print(p_spatial_main)

# 分面地图：按效应类型分组
spatial_temp_sif_sig <- spatial_temp_sif %>%
  filter(effect_type_temp_sif %in% c("促进效应(+)", "抑制效应(-)"))

if (nrow(spatial_temp_sif_sig) > 0) {
  p_spatial_facet <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(data = spatial_temp_sif_sig,
               aes(x = longitude, y = latitude,
                   color = causality_direction,
                   size = effect_strength),
               alpha = 0.7) +
    scale_color_manual(
      values = c("Temp → SIF" = "#1b9e77", "双向因果" = "#d95f02"),
      name = "因果方向"
    ) +
    scale_size_continuous(range = c(2, 6), name = "|∂|") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"])) +
    facet_wrap(~ effect_type_temp_sif, ncol = 2) +
    labs(
      title = "促进效应 vs 抑制效应的空间分布",
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
  
  print(p_spatial_facet)
}

# 效应指数的空间分布（连续色标）
p_spatial_continuous <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = spatial_temp_sif %>% filter(!is.na(effect_index_temp_sif)),
             aes(x = lon, y = lat,
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
    title = "效应指数（S-map系数）的空间分布",
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

print(p_spatial_continuous)

# ============================================================================
# 最终总结表
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("                   最终总结\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_table <- tibble(
  分类 = c(
    "总站点数",
    "有Temp→SIF因果",
    "  - 单向因果",
    "  - 双向因果",
    "有SIF→Temp因果",
    "无显著因果",
    "",
    "Temp→SIF效应类型:",
    "  - 抑制效应(-)",
    "  - 促进效应(+)",
    "  - 效应极弱",
    "  - S-map失败",
    "",
    "平均效应指数",
    "抑制效应比例"),
  数量 = c(
    nrow(ccm_results_all),
    nrow(temp_to_sif_results),
    sum(temp_to_sif_results$causality_direction == "Temp → SIF"),
    sum(temp_to_sif_results$causality_direction == "双向因果"),
    nrow(sif_to_temp_results),
    sum(ccm_results_all$causality_direction == "无显著因果"),
    "",
    "",
    sum(temp_to_sif_results$effect_type_temp_sif == "抑制效应(-)", na.rm = TRUE),
    sum(temp_to_sif_results$effect_type_temp_sif == "促进效应(+)", na.rm = TRUE),
    sum(temp_to_sif_results$effect_type_temp_sif == "效应极弱", na.rm = TRUE),
    sum(temp_to_sif_results$effect_type_temp_sif == "S-map失败", na.rm = TRUE),
    "",
    round(temp_sif_index_stats$mean, 4),
    paste0(temp_sif_index_stats$pct_negative, "%")
  )
)

print(summary_table)

