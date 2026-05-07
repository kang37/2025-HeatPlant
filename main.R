# Statement ----
# 基于每日VPD的热事件分析

# Preparation ----
pacman::p_load(
  dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr, tidyr, 
  showtext, rEDM, sf, rnaturalearth, rnaturalearthdata, knitr, ggalluvial,
  patchwork, readxl
)
showtext_auto()

# Data ----
# 函数：读取每日气象数据。
read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1) %>%
    return()
}

# 气象数据文件列表。
meteo_file_list <- list.files(
  path = "data_raw/meteo_data_1961-2023",
  pattern = "\\.txt$",
  full.names = TRUE
) %>%
  .[!grepl("sta_lonlat_china.txt", .)]

# 读取SIF站点列表。
# Bug: 为何SIF数据只有2000-2023年？LHSIF是否覆盖更多年份？
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>% 
  filter(!is.na(sif)) %>% 
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 目标站点。
target_stations <- unique(meteo_sif_data$meteo_stat_id)
cat("目标站点数:", length(target_stations), "\n")

# 地图。
china_map <- 
  ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 站点坐标。
# Bug：可以合并它和target_stations？
station_coords <- read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 读取每日数据。
meteo_data_daily <- map(
  meteo_file_list[
    str_extract(basename(meteo_file_list), "\\d+") %in% target_stations
  ],
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
  # Bug：仅分析每年的5-9月是否合理？
  filter(year %in% c(unique(meteo_sif_data$year)), month %in% 5:9)

# 计算每日VPD。
meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    # 饱和水汽压 (kPa)
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    # 实际水汽压 (kPa)
    avp = (rh / 100) * svp,
    # VPD (kPa)
    vpd = svp - avp
  ) %>%
  # 数据质量控制。
  mutate(
    vpd = case_when(
      vpd < 0 ~ NA_real_,
      vpd > 10 ~ NA_real_,
      is.na(tavg) | is.na(rh) ~ NA_real_,
      TRUE ~ vpd
    )
  )

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

# 定义热事件阈值：基于文献的绝对阈值。
# Bug：是否需要基于不同阈值做敏感性分析？
# 2.0 kPa是常用的高VPD阈值。
vpd_threshold <- 2.0  

# 识别每日热事件并构建月度指标。
meteo_data_daily_events <- meteo_data_daily_vpd %>%
  mutate(
    # 是否为热事件
    is_heat_event = vpd > vpd_threshold,
    # 热事件强度（超过阈值的部分）。
    heat_intensity = pmax(vpd - vpd_threshold, 0)
  )

# 月度汇总：每个月的胁迫指标和SIF。
monthly_heat_metrics <- meteo_data_daily_events %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(
    # 基础统计
    n_days = n(),
    vpd_mean = mean(vpd, na.rm = TRUE),
    vpd_max = max(vpd, na.rm = TRUE),
    vpd_sd = sd(vpd, na.rm = TRUE),
    
    # 热事件频次和频率。
    heat_event_days = sum(is_heat_event, na.rm = TRUE),
    heat_event_freq = heat_event_days / n_days, 
    
    # 热事件强度：计算超过阈值的累计强度。
    heat_over_sum = sum(heat_intensity, na.rm = TRUE),  
    .groups = "drop"
  ) %>%
  # 标准化处理（按站点）
  group_by(meteo_stat_id) %>%
  mutate(
    # Z-score标准化
    heat_freq_z = scale(heat_event_freq)[,1],
    heat_intensity_z = scale(heat_over_sum)[,1],
    # 综合热事件指数（Bug：权重可调整）。
    heat_index_composite = 0.5 * heat_freq_z + 0.5 * heat_intensity_z
  ) %>%
  filter(n() > 30) %>% 
  ungroup()

# 月度SIF数据
sif_monthly <- meteo_sif_data %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")

# 函数：安全去趋势
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

# 合并热胁迫和SIF数据并去趋势。
data_heat_sif <- monthly_heat_metrics %>%
  inner_join(sif_monthly, by = c("meteo_stat_id", "year", "month")) %>%
  filter(!is.na(sif)) %>% 
  # 去趋势。
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    
    # 安全去趋势
    sif_detrended = safe_detrend(sif, time_idx),
    heat_index_composite_detrended = safe_detrend(heat_index_composite, time_idx),
    heat_event_freq_detrended = safe_detrend(heat_event_freq, time_idx),
    heat_intensity_sum_detrended = safe_detrend(heat_over_sum, time_idx)
  ) %>%
  ungroup() %>%
  # 过滤掉去趋势失败的站点
  group_by(meteo_stat_id) %>%
  # Bug: 保留多少个有效数据点？
  filter(sum(!is.na(sif_detrended)) >= 30) %>%  
  ungroup()

cat("合并后数据行数:", nrow(data_heat_sif), "\n")
cat("站点数:", length(unique(data_heat_sif$meteo_stat_id)), "\n\n")

tar_make()
tar_load(data_heat_sif)

# 各站点热事件与SIF相关图 ----
# Bug：没有时间滞后情况下的相关性。
# 准备输出目录。
if (!dir.exists("data_proc/station_plots")) {
  dir.create("data_proc/station_plots", recursive = TRUE)
}

# 获取所有站点并分组 (每页25个站点)
all_stations <- unique(data_heat_sif$meteo_stat_id)
stations_per_page <- 25
num_pages <- ceiling(length(all_stations) / stations_per_page)

cat("共有站点:", length(all_stations), "，将分", num_pages, "页输出...\n")

for (i in 1:num_pages) {
  start_idx <- (i - 1) * stations_per_page + 1
  end_idx <- min(i * stations_per_page, length(all_stations))
  current_stations <- all_stations[start_idx:end_idx]
  
  # 提取当前页数据
  plot_df <- data_heat_sif %>%
    filter(meteo_stat_id %in% current_stations)
  
  # 绘图
  p <- ggplot(plot_df, aes(x = heat_index_composite_detrended, y = sif_detrended)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1) +
    geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 0.8) +
    facet_wrap(~meteo_stat_id, scales = "free", ncol = 5) +
    labs(
      title = paste0("Heat Index vs SIF (Page ", i, "/", num_pages, ")"),
      subtitle = "X: heat_index_composite_detrended | Y: sif_detrended",
      x = "Heat Impact (Detrended)",
      y = "SIF (Detrended)"
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 8))
  
  # 保存
  file_name <- sprintf("data_proc/station_plots/station_scatters_page_%02d.png", i)
  ggsave(file_name, p, width = 15, height = 12, dpi = 150)
  
  if (i %% 5 == 0) cat("已完成", i, "页...\n")
}

# CCM ----
# 函数：进行CCM分析。
perform_ccm_heat_sif <- function(
    station_id, data, min_points = 30, 
    sif_col = "sif_detrended", heat_col = "heat_index_composite_detrended", 
    tp_x = 0) {
  # 提取特定站点的数据。
  station_data_ccm <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    group_by(year) %>% 
    mutate(heat_index = dplyr::lag(!!sym(heat_col), n = tp_x)) %>% 
    ungroup() %>% 
    select(sif = all_of(sif_col), heat_index) %>%
    filter(!is.na(sif), !is.na(heat_index)) %>%
    mutate(time = row_number(), .before = 1) %>%
    as.data.frame()
  
  # 如果样本不足，返回空值。
  n_data <- nrow(station_data_ccm)
  if (n_data < min_points) return(NULL)
  
  tryCatch({
    # 确定最优E
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    embed_heat <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
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
    # Tp设置为0为同步，否则为滞后关系。
    tp <- tp_x
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) return(NULL)
    
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
    
    # 汇总
    ccm_summary_heat <- ccm_heat_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`heat_index:sif`, na.rm = TRUE),
        rho_sd = sd(`heat_index:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    final_rho_heat <- ccm_summary_heat$rho_mean[nrow(ccm_summary_heat)]
    trend_heat <- cor(ccm_summary_heat$lib_size, ccm_summary_heat$rho_mean)
    
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
      heat_coef_col <- 
        which(grepl("heat", colnames(smap_coeffs), ignore.case = TRUE))
      
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
    
    heat_causes_sif <- (
      final_rho_heat > ccm_threshold_rho & 
        trend_heat > ccm_threshold_trend
    )
    
    effect_type <- case_when(
      !heat_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0 ~ "促进",
      mean_smap_coef < 0 ~ "抑制"
    )
    
    effect_stability <- case_when(
      is.na(mean_smap_coef) | is.na(sd_smap_coef) ~ "未知",
      abs(mean_smap_coef) > sd_smap_coef ~ "稳定",
      TRUE ~ "波动大"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      tp = tp_x, 
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
      
      # 保存数据
      smap_coefficients_heat = list(smap_coef_heat),
      ccm_summary_heat = list(ccm_summary_heat)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# 批量分析。
all_stations_heat <- data_heat_sif %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("准备分析", length(all_stations_heat), "个站点...\n\n")

# 不同Tp下各站点CCM结果。
ccm_results_heat <- map_dfr(c(0:4), function(current_tp) {
  # 内部循环遍历所有站点
  map_dfr(seq_along(all_stations_heat), function(i) {
    perform_ccm_heat_sif(
      all_stations_heat[i], data_heat_sif, tp_x = current_tp
    )
  })
}) %>% 
  # 加入站点的经纬度数据。
  left_join(station_coords, by = "meteo_stat_id") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(effect_strength = abs(effect_index_heat))

# Station causal effect change Sankey ----
# 提取所有Tp对。
available_tps <- sort(unique(ccm_results_heat$tp))
tp_pairs <- map(
  1:(length(available_tps) - 1), ~ c(available_tps[.x], available_tps[.x+1])
)

# 函数：作桑基图，展示各站点因果关系变化。
plot_sankey_pair <- function(data, pair) {
  # 定义作图样式。
  lvls <- c("促进", "抑制", "无因果", "S-map失败")
  cols <- c("促进" = "#E41A1C", "抑制" = "#377EB8", "无因果" = "#999999")
  
  # 提取Tp对。
  tp_start <- pair[1]
  tp_end   <- pair[2]
  
  # 提取这对 tp 的数据并转宽
  pair_data <- data %>%
    filter(tp %in% pair) %>%
    select(meteo_stat_id, tp, effect_type_heat) %>%
    pivot_wider(
      id_cols = meteo_stat_id, names_from = tp, values_from = effect_type_heat
    ) %>%
    # 仅保留在该对中均有数据的站点
    drop_na() %>%
    # 转换为因子
    mutate(
      start_val = 
        factor(as.character(get(as.character(tp_start))), levels = lvls),
      end_val = 
        factor(as.character(get(as.character(tp_end))),   levels = lvls)
    )
  
  # 统计频数
  pair_summary <- pair_data %>%
    group_by(start_val, end_val) %>%
    summarise(n = n(), .groups = "drop")
  
  # 绘图
  ggplot(pair_summary,
         aes(axis1 = start_val, axis2 = end_val, y = n)) +
    geom_alluvium(aes(fill = start_val), width = 1/12, alpha = 0.6) +
    geom_stratum(width = 1/6, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c(paste0("tp = ", tp_start), paste0("tp = ", tp_end)), expand = c(.1, .1)) +
    scale_fill_manual(values = cols) +
    labs(subtitle = paste0("流转: tp=", tp_start, " -> tp=", tp_end),
         y = "站点数", fill = "起始性质") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.y = element_blank())
}

# 批量生成并组合。
all_plots <- map(tp_pairs, ~plot_sankey_pair(ccm_results_heat, .x))

# 使用patchwork组合，每行放2个图。
combined_plot <- wrap_plots(all_plots, ncol = 2) + 
  plot_annotation(
    title = "各站点因果性质演变",
    theme = theme(plot.title = element_text(size = 18))
  )

# 结果统计 ----
# 因果类型数量。
print(table(ccm_results_heat$effect_type_heat))
print(table(ccm_results_heat$effect_stability_heat))

# 空间可视化
# 主地图
p_heat_spatial <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(
    data = ccm_results_heat,
    aes(x = longitude, y = latitude, color = effect_type_heat, size = effect_strength),
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("促进" = "#4575b4", "抑制" = "#d73027", "无因果" = "#999999"),
    name = "因果效应"
  ) +
  scale_size_continuous(range = c(1, 8), name = "效应强度\n|∂|") +
  labs(
    title = "VPD热事件指数对SIF的因果影响",
    subtitle = paste0(
      "抑制效应: ", sum(ccm_results_heat$effect_type_heat == "抑制"), " | ",
      "促进效应: ", sum(ccm_results_heat$effect_type_heat == "促进"), " | ",
      "VPD阈值 = ", round(vpd_threshold, 2), " kPa"
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
  ) + 
  facet_wrap(.~ tp)

# ============================================================================
# 全维度整合分析 (气候带 + 投资) ----
# ============================================================================

cat("【正在构建全维度 ccm_results_heat 对象 (自包含版)...】\n")

# 1. 变量富集：气候带匹配 ----
if (!require("terra")) install.packages("terra")
koppen_raster <- terra::rast("data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")
koppen_lookup <- c(
  "1" = "Af", "2" = "Am", "3" = "As", "4" = "Aw", "5" = "BWh", "6" = "BWk", "7" = "BSh", "8" = "BSk",
  "9" = "Csa", "10" = "Csb", "11" = "Csc", "12" = "Cwa", "13" = "Cwb", "14" = "Cwc", "15" = "Cfa",
  "16" = "Cfb", "17" = "Cfc", "18" = "Dsa", "19" = "Dsb", "20" = "Dsc", "21" = "Dsd", "22" = "Dwa",
  "23" = "Dwb", "24" = "Dwc", "25" = "Dwd", "26" = "Dfa", "27" = "Dfb", "28" = "Dfc", "29" = "Dfd", "30" = "ET", "31" = "EF"
)
pts <- terra::vect(ccm_results_heat, geom = c("longitude", "latitude"), crs = "EPSG:4326")
ccm_results_heat$koppen_code <- as.character(terra::extract(koppen_raster, pts)[,2])
ccm_results_heat$koppen_class <- koppen_lookup[ccm_results_heat$koppen_code]

# 2. 变量富集：绿化投资指标匹配 (独立集成) ----
cat("正在读取经济与绿化数据...\n")

# A. 2020年绿化投资 (万元)
green_invest_2020 <- read_csv("data_raw/green_invest/中国城市数据.csv", show_col_types = FALSE) %>%
  dplyr::select(city_name = 城市名称, invest = 园林绿化_2020) %>%
  mutate(invest = as.numeric(invest),
         city_name = ifelse(str_detect(city_name, "市$"), city_name, paste0(city_name, "市")))

# B. 2020年GDP (万元)
gdp_data_2020 <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, gdp = "地区生产总值(万元)") %>%
  mutate(gdp = as.numeric(gdp))

# C. 2020年绿地面积 (公顷)
# 请注意：此处使用了您之前修正过的列名，如果仍有空格请微调
green_area_2020 <- read_excel("data_raw/城市绿地面积数据_2003-2023.xlsx", sheet = "绿地面积_明细数据") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, area_green = `绿地面积(公顷)`) %>%
  mutate(area_green = as.numeric(area_green),
         city_name = ifelse(str_detect(city_name, "市$"), city_name, paste0(city_name, "市")))

# D. 空间匹配站点到城市
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326)
pts_sf <- st_as_sf(ccm_results_heat %>% distinct(meteo_stat_id, longitude, latitude), 
                   coords = c("longitude", "latitude"), crs = 4326)
station_city_map <- pts_sf %>% 
  st_join(china_cities_shp, join = st_within) %>% 
  as.data.frame() %>%
  dplyr::select(meteo_stat_id, city_name = ct_name)

# E. 指标整合
invest_metrics <- station_city_map %>%
  left_join(green_invest_2020, by = "city_name") %>%
  left_join(gdp_data_2020, by = "city_name") %>%
  left_join(green_area_2020, by = "city_name") %>%
  mutate(
    invest_ratio = invest / gdp * 100,
    intensity_green = invest / area_green
  ) %>%
  dplyr::select(meteo_stat_id, invest_ratio, intensity_green)

ccm_results_heat <- ccm_results_heat %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  mutate(tp_label = paste0("Lag ", tp))

# 3. 简洁可视化 ----
cat("生成整合分析图表...\n")

# A. 气候带分析
p_climate <- ggplot(ccm_results_heat %>% filter(!is.na(koppen_class), effect_type_heat %in% c("促进", "抑制")),
                    aes(x = koppen_class, y = rho_heat_to_sif, fill = effect_type_heat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~tp_label) +
  scale_fill_manual(values = c("促进" = "#E41A1C", "抑制" = "#377EB8")) +
  labs(title = "不同气候带下的因果强度分布", x = "柯本气候分类", y = "因果强度 (Rho)", fill = "性质") +
  theme_minimal(base_size = 20) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# B. 投资强度分析
p_invest <- ggplot(ccm_results_heat %>% filter(effect_type_heat %in% c("促进", "抑制")),
                   aes(x = effect_type_heat, y = intensity_green, fill = effect_type_heat)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~tp_label) +
  scale_fill_manual(values = c("促进" = "#E41A1C", "抑制" = "#377EB8")) +
  labs(title = "投资强度与因果性质的关系", x = "因果性质", y = "单位绿地投资 (万元/公顷, log)", fill = "性质") +
  theme_minimal(base_size = 20)

# 保存 (使用 png 确保字大)
png("data_proc/integrated_climate_invest_analysis.png", width = 2400, height = 3000, res = 300)
print(p_climate / p_invest)
dev.off()

# 因果效应强度分布 ----
## 根据时间延迟 ----
# 预处理绘图数据
plot_data_rho <- ccm_results_heat %>%
  filter(!is.na(effect_type_heat)) %>%
  mutate(
    tp_label = paste0("Time Lag = ", tp, "月"),
    effect_type = factor(effect_type_heat, 
                         levels = c("促进", "抑制", "无因果", "S-map失败"))
  )

# 我们关注有因果关系的站点 (促进和抑制)
p_rho_dist <- ggplot(plot_data_rho %>% filter(effect_type %in% c("促进", "抑制")), 
                     aes(x = rho_heat_to_sif, fill = effect_type)) +
  geom_density(alpha = 0.6, color = "white") +
  facet_wrap(~tp_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("促进" = "#E41A1C", "抑制" = "#377EB8")) +
  labs(
    title = "不同延迟下的因果强度 (CCM Rho) 密度分布",
    subtitle = "对比分析促进效应与抑制效应在不同滞后时间下的强度特征",
    x = "因果强度 (Rho)",
    y = "密度",
    fill = "因果性质"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

# 保存图表
png("data_proc/ccm_rho_density_distribution.png", width = 2000, height = 2000, res = 300)
p_rho_dist
dev.off()

## 加入纬度 ----
p_lat_dist <- ggplot(plot_data_rho %>% filter(effect_type %in% c("促进", "抑制")), 
                     aes(x = latitude, y = rho_heat_to_sif, color = effect_type)) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, span = 0.6, linewidth = 1.5) +
  facet_wrap(~tp_label, ncol = 1) +
  scale_color_manual(values = c("促进" = "#377EB8", "抑制" = "#E41A1C")) +
  labs(
    title = "因果强度随纬度的地理分布趋势",
    subtitle = "展示不同延迟下，因果信号强度的南北分异 (LOESS平滑)",
    x = "纬度 (Latitude)",
    y = "因果强度 (CCM Rho)",
    color = "因果性质"
  ) +
  theme_minimal(base_size = 100) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        legend.position = "bottom")
ggsave("data_proc/ccm_rho_latitude_gradient.png", p_lat_dist, width = 20, height = 32, dpi = 300)

