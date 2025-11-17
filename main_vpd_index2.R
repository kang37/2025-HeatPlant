# ============================================================================
# VPD热事件综合分析（改进版）
# ============================================================================

pacman::p_load(dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, 
               tidyr, showtext, rEDM, sf, rnaturalearth, rnaturalearthdata,
               patchwork)
showtext_auto()

cat("\n", rep("=", 70), "\n", sep = "")
cat("              VPD热事件综合分析\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================================
# 1. 读取每日气象数据（与之前相同，这里省略）
# ============================================================================

# [前面读取数据的代码保持不变...]

# ============================================================================
# 2. 定义VPD阈值（两个版本）
# ============================================================================

cat("【2. 定义VPD阈值】\n")

# 全局VPD统计
vpd_global_stats <- meteo_data_daily_vpd %>%
  summarise(
    mean = mean(vpd, na.rm = TRUE),
    sd = sd(vpd, na.rm = TRUE),
    q80 = quantile(vpd, 0.80, na.rm = TRUE),
    q85 = quantile(vpd, 0.85, na.rm = TRUE),
    q90 = quantile(vpd, 0.90, na.rm = TRUE),
    q95 = quantile(vpd, 0.95, na.rm = TRUE)
  )

# 版本1：基于百分位（80%）
VPD_THRESHOLD_P80 <- vpd_global_stats$q80

# 版本2：基于文献的绝对阈值
# 参考 Yuan et al. (2019) 和其他文献
VPD_THRESHOLD_ABS <- 2.0  # kPa，常用的高VPD阈值

cat("\n阈值设定:\n")
cat("  版本1 (P80): ", round(VPD_THRESHOLD_P80, 3), " kPa\n")
cat("  版本2 (文献): ", VPD_THRESHOLD_ABS, " kPa\n\n")

# ============================================================================
# 3. 构建月度VPD热影响指标
# ============================================================================

cat("【3. 构建月度VPD热影响指标】\n")

build_monthly_heat_index <- function(daily_data, vpd_threshold, version_name) {
  
  daily_data %>%
    mutate(
      # 超过阈值的部分
      vpd_excess = pmax(vpd - vpd_threshold, 0),
      
      # 是否为热事件日
      is_heat_day = vpd > vpd_threshold
    ) %>%
    group_by(meteo_stat_id, year, month) %>%
    summarise(
      # 基础统计
      n_days = n(),
      vpd_mean = mean(vpd, na.rm = TRUE),
      vpd_max = max(vpd, na.rm = TRUE),
      
      # 热事件统计
      heat_days = sum(is_heat_day, na.rm = TRUE),
      
      # 【关键指标】累积热影响 = Σ(VPD超过阈值的部分)
      cumulative_heat_impact = sum(vpd_excess, na.rm = TRUE),
      
      # 平均热强度
      mean_heat_intensity = mean(vpd_excess[is_heat_day], na.rm = TRUE),
      
      # 其他辅助指标
      heat_freq = heat_days / n_days,
      max_heat_intensity = max(vpd_excess, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      version = version_name,
      threshold = vpd_threshold
    )
}

# 版本1：P80阈值
monthly_heat_p80 <- build_monthly_heat_index(
  meteo_data_daily_vpd, 
  VPD_THRESHOLD_P80, 
  "P80"
)

# 版本2：文献阈值
monthly_heat_abs <- build_monthly_heat_index(
  meteo_data_daily_vpd, 
  VPD_THRESHOLD_ABS, 
  "Literature"
)

cat("版本1 (P80) - 热影响统计:\n")
monthly_heat_p80 %>%
  summarise(
    总月数 = n(),
    有热影响月份 = sum(cumulative_heat_impact > 0),
    比例 = round(sum(cumulative_heat_impact > 0) / n() * 100, 1),
    平均累积热影响 = mean(cumulative_heat_impact[cumulative_heat_impact > 0])
  ) %>%
  print()

cat("\n版本2 (文献) - 热影响统计:\n")
monthly_heat_abs %>%
  summarise(
    总月数 = n(),
    有热影响月份 = sum(cumulative_heat_impact > 0),
    比例 = round(sum(cumulative_heat_impact > 0) / n() * 100, 1),
    平均累积热影响 = mean(cumulative_heat_impact[cumulative_heat_impact > 0])
  ) %>%
  print()

cat("\n")

# ============================================================================
# 4. 合并SIF数据并去趋势
# ============================================================================

cat("【4. 合并SIF数据】\n")

# SIF数据
sif_monthly <- meteo_sif_data %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")

# 安全去趋势函数
safe_detrend <- function(x, time_idx) {
  valid_data <- !is.na(x)
  n_valid <- sum(valid_data)
  
  if (n_valid < 3) {
    return(rep(NA_real_, length(x)))
  }
  
  tryCatch({
    model <- lm(x ~ time_idx)
    residuals(model)
  }, error = function(e) {
    rep(NA_real_, length(x))
  })
}

# 处理两个版本
prepare_data_for_ccm <- function(monthly_heat, sif_data) {
  monthly_heat %>%
    inner_join(sif_data, by = c("meteo_stat_id", "year", "month")) %>%
    filter(!is.na(sif)) %>%
    group_by(meteo_stat_id) %>%
    arrange(year, month) %>%
    mutate(
      time_idx = row_number(),
      sif_detrended = safe_detrend(sif, time_idx),
      heat_impact_detrended = safe_detrend(cumulative_heat_impact, time_idx),
      heat_days_detrended = safe_detrend(heat_days, time_idx)
    ) %>%
    ungroup() %>%
    group_by(meteo_stat_id) %>%
    filter(sum(!is.na(sif_detrended)) >= 10,
           sum(!is.na(heat_impact_detrended)) >= 10) %>%
    ungroup()
}

data_ccm_p80 <- prepare_data_for_ccm(monthly_heat_p80, sif_monthly)
data_ccm_abs <- prepare_data_for_ccm(monthly_heat_abs, sif_monthly)

cat("版本1 (P80) - 可用站点:", length(unique(data_ccm_p80$meteo_stat_id)), "\n")
cat("版本2 (文献) - 可用站点:", length(unique(data_ccm_abs$meteo_stat_id)), "\n\n")

# ============================================================================
# 5. CCM分析函数（增强版）
# ============================================================================

cat("【5. CCM分析】\n\n")

perform_ccm_heat_sif_v2 <- function(station_id, data, version_name, min_points = 30) {
  
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(
      sif = sif_detrended,
      heat_impact = heat_impact_detrended,
      heat_days = heat_days_detrended
    ) %>%
    filter(!is.na(sif), !is.na(heat_impact)) %>%
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    
    station_data_ccm <- data.frame(
      time = 1:n_data,
      sif = station_data$sif,
      heat_impact = station_data$heat_impact,
      heat_days = station_data$heat_days
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
      columns = "heat_impact",
      target = "heat_impact",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    best_E_heat <- embed_heat$E[which.max(embed_heat$rho)]
    E_ccm <- max(best_E_sif, best_E_heat)
    
    # CCM参数
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    tp <- 0
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) {
      return(NULL)
    }
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # CCM: 热影响 → SIF
    ccm_heat_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "heat_impact",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # CCM: SIF → 热影响（反向）
    ccm_sif_to_heat <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "sif",
      target = "heat_impact",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 汇总
    ccm_summary_heat_sif <- ccm_heat_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`heat_impact:sif`, na.rm = TRUE),
        rho_sd = sd(`heat_impact:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    ccm_summary_sif_heat <- ccm_sif_to_heat %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`sif:heat_impact`, na.rm = TRUE),
        rho_sd = sd(`sif:heat_impact`, na.rm = TRUE),
        .groups = "drop"
      )
    
    final_rho_heat_sif <- ccm_summary_heat_sif$rho_mean[nrow(ccm_summary_heat_sif)]
    final_rho_sif_heat <- ccm_summary_sif_heat$rho_mean[nrow(ccm_summary_sif_heat)]
    
    trend_heat_sif <- cor(ccm_summary_heat_sif$lib_size, ccm_summary_heat_sif$rho_mean)
    trend_sif_heat <- cor(ccm_summary_sif_heat$lib_size, ccm_summary_sif_heat$rho_mean)
    
    # S-map分析
    smap_heat <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "heat_impact",
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
    
    # 判断因果关系
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    heat_causes_sif <- (final_rho_heat_sif > ccm_threshold_rho & 
                          trend_heat_sif > ccm_threshold_trend)
    sif_causes_heat <- (final_rho_sif_heat > ccm_threshold_rho & 
                          trend_sif_heat > ccm_threshold_trend)
    
    # 因果方向
    causality_direction <- case_when(
      heat_causes_sif & !sif_causes_heat ~ "热影响→SIF",
      !heat_causes_sif & sif_causes_heat ~ "SIF→热影响",
      heat_causes_sif & sif_causes_heat ~ "双向因果",
      TRUE ~ "无因果"
    )
    
    # 效应类型
    effect_type <- case_when(
      !heat_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0.001 ~ "促进效应(+)",
      mean_smap_coef < -0.001 ~ "抑制效应(-)",
      TRUE ~ "效应极弱"
    )
    
    # 组合标签
    combined_label <- case_when(
      causality_direction == "无因果" ~ "无因果",
      causality_direction == "热影响→SIF" & effect_type == "促进效应(+)" ~ "热影响促进SIF",
      causality_direction == "热影响→SIF" & effect_type == "抑制效应(-)" ~ "热影响抑制SIF",
      causality_direction == "SIF→热影响" ~ "SIF→热影响",
      causality_direction == "双向因果" & effect_type == "促进效应(+)" ~ "双向-促进",
      causality_direction == "双向因果" & effect_type == "抑制效应(-)" ~ "双向-抑制",
      TRUE ~ "其他"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      version = version_name,
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
      
      # 效应
      effect_index = mean_smap_coef,
      effect_index_sd = sd_smap_coef,
      effect_type = effect_type,
      combined_label = combined_label,
      
      # 保存数据
      smap_coefficients = list(smap_coef_heat)
    )
    
  }, error = function(e) {
    return(NULL)
  })
}

# 分析两个版本
cat("开始CCM分析...\n")

# 版本1
all_stations_p80 <- unique(data_ccm_p80$meteo_stat_id)
cat("版本1 (P80) - 分析", length(all_stations_p80), "个站点\n")

ccm_results_p80 <- map_dfr(seq_along(all_stations_p80), function(i) {
  if (i %% 20 == 0) cat("  P80已完成:", i, "/", length(all_stations_p80), "\n")
  perform_ccm_heat_sif_v2(all_stations_p80[i], data_ccm_p80, "P80", min_points = 30)
})

cat("版本1完成，成功分析:", nrow(ccm_results_p80), "个站点\n\n")

# 版本2
all_stations_abs <- unique(data_ccm_abs$meteo_stat_id)
cat("版本2 (文献) - 分析", length(all_stations_abs), "个站点\n")

ccm_results_abs <- map_dfr(seq_along(all_stations_abs), function(i) {
  if (i %% 20 == 0) cat("  文献已完成:", i, "/", length(all_stations_abs), "\n")
  perform_ccm_heat_sif_v2(all_stations_abs[i], data_ccm_abs, "Literature", min_points = 30)
})

cat("版本2完成，成功分析:", nrow(ccm_results_abs), "个站点\n\n")

# 合并结果
ccm_results_all <- bind_rows(ccm_results_p80, ccm_results_abs)

# ============================================================================
# 6. 统计汇总
# ============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("                   统计结果\n")
cat(rep("=", 70), "\n\n", sep = "")

# 按版本统计
for (ver in c("P80", "Literature")) {
  
  cat("【", ver, "版本】\n", sep = "")
  
  ver_results <- ccm_results_all %>% filter(version == ver)
  
  # 因果方向统计
  causality_stats <- table(ver_results$causality_direction)
  cat("\n因果方向分布:\n")
  print(causality_stats)
  cat("百分比:\n")
  print(round(prop.table(causality_stats) * 100, 1))
  
  # 效应类型统计（仅热影响→SIF）
  heat_to_sif <- ver_results %>% filter(heat_causes_sif)
  
  if (nrow(heat_to_sif) > 0) {
    cat("\n热影响→SIF的效应类型:\n")
    effect_stats <- table(heat_to_sif$effect_type)
    print(effect_stats)
    cat("百分比:\n")
    print(round(prop.table(effect_stats) * 100, 1))
  }
  
  # 组合标签统计
  cat("\n综合分类统计:\n")
  combined_stats <- ver_results %>%
    count(combined_label) %>%
    mutate(percentage = round(n / sum(n) * 100, 1)) %>%
    arrange(desc(n))
  
  print(combined_stats)
  
  cat("\n", rep("-", 70), "\n\n", sep = "")
}

# 对比表
comparison_summary <- ccm_results_all %>%
  group_by(version) %>%
  summarise(
    总站点数 = n(),
    有因果站点 = sum(causality_direction != "无因果"),
    热影响促进SIF = sum(combined_label == "热影响促进SIF"),
    热影响抑制SIF = sum(combined_label == "热影响抑制SIF"),
    SIF到热影响 = sum(combined_label == "SIF→热影响"),
    双向因果 = sum(causality_direction == "双向因果"),
    无因果 = sum(causality_direction == "无因果"),
    .groups = "drop"
  ) %>%
  pivot_longer(-version, names_to = "指标", values_to = "数量")

cat("【版本对比】\n")
comparison_summary %>%
  pivot_wider(names_from = version, values_from = 数量) %>%
  print()

cat("\n")

# ============================================================================
# 7. 空间可视化
# ============================================================================

cat("【7. 空间可视化】\n\n")

# 读取坐标
station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并坐标
spatial_results <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    effect_strength = abs(effect_index)
  )

cat("有坐标的站点数:", nrow(spatial_results), "\n\n")

# 中国地图
china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 为两个版本分别制图
for (ver in c("P80", "Literature")) {
  
  ver_spatial <- spatial_results %>% filter(version == ver)
  
  threshold_val <- ifelse(ver == "P80", VPD_THRESHOLD_P80, VPD_THRESHOLD_ABS)
  
  # 主地图：所有类型
  p_main <- ggplot() +
    geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
    geom_point(
      data = ver_spatial,
      aes(x = longitude, y = latitude,
          color = combined_label,
          size = effect_strength),
      alpha = 0.7
    ) +
    scale_color_manual(
      values = c(
        "热影响促进SIF" = "#4575b4",
        "热影响抑制SIF" = "#d73027",
        "SIF→热影响" = "#984ea3",
        "双向-促进" = "#377eb8",
        "双向-抑制" = "#e41a1c",
        "无因果" = "#999999",
        "其他" = "#fee090"
      ),
      name = "因果类型"
    ) +
    scale_size_continuous(range = c(1, 8), name = "效应强度\n|∂|") +
    labs(
      title = paste0("VPD热影响对SIF的因果关系空间分布 (", ver, ")"),
      subtitle = paste0(
        "VPD阈值 = ", round(threshold_val, 2), " kPa | ",
        "有因果站点: ", sum(ver_spatial$causality_direction != "无因果"), " / ",
        nrow(ver_spatial)
      ),
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
  
  print(p_main)
  
  # 分面图：按因果方向分组
  ver_spatial_sig <- ver_spatial %>%
    filter(causality_direction %in% c("热影响→SIF", "SIF→热影响", "双向因果"))
  
  if (nrow(ver_spatial_sig) > 0) {
    p_facet <- ggplot() +
      geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
      geom_point(
        data = ver_spatial_sig,
        aes(x = longitude, y = latitude,
            color = effect_type,
            size = effect_strength),
        alpha = 0.8
      ) +
      scale_color_manual(
        values = c(
          "促进效应(+)" = "#4575b4",
          "抑制效应(-)" = "#d73027",
          "效应极弱" = "#fdae61"
        ),
        name = "效应类型"
      ) +
      scale_size_continuous(range = c(1, 3), name = "|∂|") +
      facet_wrap(~ causality_direction, ncol = 3) +
      labs(
        title = paste0("因果关系类型空间分布 (", ver, ")"),
        x = "经度",
        y = "纬度"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        panel.background = element_rect(fill = "aliceblue", color = NA)
      )
    
    print(p_facet)
  }
}

# ============================================================================
# 8. 与经济发展指标的关系
# ============================================================================

cat("\n【8. 因果关系与经济发展指标的关系】\n\n")

# 读取GDP/经济数据（需要你提供数据文件）
# 这里我提供一个示例框架

# 示例：创建模拟的经济数据
# 实际使用时，替换为真实数据
economic_data <- station_coords %>%
  mutate(
    # 这里应该从真实数据源读取
    # 示例：根据经纬度估算（仅用于演示）
    gdp_per_capita = case_when(
      longitude > 110 & latitude > 30 ~ rnorm(n(), 80000, 20000),  # 东部发达
      longitude > 105 ~ rnorm(n(), 50000, 15000),  # 中部
      TRUE ~ rnorm(n(), 30000, 10000)  # 西部
    ),
    urbanization_rate = case_when(
      longitude > 110 & latitude > 30 ~ runif(n(), 0.6, 0.9),
      longitude > 105 ~ runif(n(), 0.4, 0.7),
      TRUE ~ runif(n(), 0.2, 0.5)
    )
  )

# 如果有真实数据，使用以下代码读取：
# economic_data <- read_csv("data_raw/economic_indicators.csv") %>%
#   rename_with(~tolower(.x)) %>%
#   mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并数据
spatial_econ <- spatial_results %>%
  left_join(economic_data, by = c("meteo_stat_id", "longitude", "latitude"))

# 为两个版本分别分析
for (ver in c("P80", "Literature")) {
  
  cat("【", ver, "版本 - 经济指标分析】\n\n", sep = "")
  
  ver_econ <- spatial_econ %>% filter(version == ver)
  
  # 按因果类型统计经济指标
  econ_by_causality <- ver_econ %>%
    filter(!is.na(gdp_per_capita)) %>%
    group_by(combined_label) %>%
    summarise(
      站点数 = n(),
      平均GDP = round(mean(gdp_per_capita, na.rm = TRUE), 0),
      GDP标准差 = round(sd(gdp_per_capita, na.rm = TRUE), 0),
      平均城镇化率 = round(mean(urbanization_rate, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(站点数))
  
  cat("因果类型与经济指标:\n")
  print(econ_by_causality)
  cat("\n")
  
  # 箱线图：GDP vs 因果类型
  p_gdp_box <- ver_econ %>%
    filter(!is.na(gdp_per_capita),
           combined_label %in% c("热影响促进SIF", "热影响抑制SIF", "无因果")) %>%
    ggplot(aes(x = combined_label, y = gdp_per_capita, fill = combined_label)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(
      values = c(
        "热影响促进SIF" = "#4575b4",
        "热影响抑制SIF" = "#d73027",
        "无因果" = "#999999"
      ),
      guide = "none"
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste0("因果类型与人均GDP的关系 (", ver, ")"),
      x = "因果类型",
      y = "人均GDP (元)",
      caption = "方框 = 四分位数；点 = 单个站点"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  print(p_gdp_box)
  散点图：GDP vs 效应强度
  p_gdp_scatter <- ver_econ %>%
    filter(!is.na(gdp_per_capita), !is.na(effect_index),
           combined_label %in% c("热影响促进SIF", "热影响抑制SIF")) %>%
    ggplot(aes(x = gdp_per_capita, y = effect_index, color = combined_label)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    scale_color_manual(
      values = c(
        "热影响促进SIF" = "#4575b4",
        "热影响抑制SIF" = "#d73027"
      ),
      name = "因果类型"
    ) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = paste0("GDP与热影响效应强度的关系 (", ver, ")"),
      x = "人均GDP (元)",
      y = "效应指数 (S-map系数)",
      caption = paste0("相关系数: ",
                       "促进 r=", round(cor(
                         ver_econgdppercapita[verecongdp_per_capita[ver_econ
                                                                    gdpp​erc​apita[vere​concombined_label == "热影响促进SIF"],
                                                                    ver_econeffectindex[vereconeffect_index[ver_econ
                                                                                                            effecti​ndex[vere​concombined_label == "热影响促进SIF"],
                                                                                                            use = "complete.obs"
                       ), 3),
                       " | 抑制 r=", round(cor(
                         ver_econgdppercapita[verecongdp_per_capita[ver_econ
                                                                    gdpp​erc​apita[vere​concombined_label == "热影响抑制SIF"],
                                                                    ver_econeffectindex[vereconeffect_index[ver_econ
                                                                                                            effecti​ndex[vere​concombined_label == "热影响抑制SIF"],
                                                                                                            use = "complete.obs"
                       ), 3))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom"
    )
  
  print(p_gdp_scatter)
  城镇化率分析
  p_urban_box <- ver_econ %>%
    filter(!is.na(urbanization_rate),
           combined_label %in% c("热影响促进SIF", "热影响抑制SIF", "无因果")) %>%
    ggplot(aes(x = combined_label, y = urbanization_rate * 100, fill = combined_label)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(
      values = c(
        "热影响促进SIF" = "#4575b4",
        "热影响抑制SIF" = "#d73027",
        "无因果" = "#999999"
      ),
      guide = "none"
    ) +
    labs(
      title = paste0("因果类型与城镇化率的关系 (", ver, ")"),
      x = "因果类型",
      y = "城镇化率 (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  print(p_urban_box)
  统计检验
  if (sum(!is.na(ver_econ$gdp_per_capita)) > 0) {
    cat("ANOVA检验 - GDP差异:\n")
    anova_result <- aov(gdp_per_capita ~ combined_label,
                        data = ver_econ %>%
                          filter(combined_label %in% c("热影响促进SIF", "热影响抑制SIF", "无因果")))
    print(summary(anova_result))
    cat("\n")