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
tar_make()
tar_load(data_heat_sif)
tar_load(ccm_results_heat)

# 地图。
china_map <- 
  ne_countries(country = "china", scale = "medium", returnclass = "sf")

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
combined_plot <- wrap_plots(all_plots, nrow = 1) + 
  plot_annotation(
    title = "各站点因果性质演变",
    theme = theme(plot.title = element_text(size = 18))
  )
combined_plot

# 结果统计 ----
# 因果类型数量。
print(table(ccm_results_heat$effect_type_heat))
print(table(ccm_results_heat$effect_stability_heat))

# 空间可视化
# Bug: 应该放在哪里？
vpd_threshold <- 2
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
p_heat_spatial

# 全维度整合分析 (气候带 + 投资) ----
# 1. 变量富集：气候带匹配 ----
koppen_raster <- terra::rast("data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")

koppen_lookup <- c(
  "1" = "Af", "2" = "Am", "3" = "As", "4" = "Aw", "5" = "BWh", "6" = "BWk", "7" = "BSh", "8" = "BSk",
  "9" = "Csa", "10" = "Csb", "11" = "Csc", "12" = "Cwa", "13" = "Cwb", "14" = "Cwc", "15" = "Cfa",
  "16" = "Cfb", "17" = "Cfc", "18" = "Dsa", "19" = "Dsb", "20" = "Dsc", "21" = "Dsd", "22" = "Dwa",
  "23" = "Dwb", "24" = "Dwc", "25" = "Dwd", "26" = "Dfa", "27" = "Dfb", "28" = "Dfc", "29" = "Dfd", "30" = "ET", "31" = "EF"
)
pts <- 
  terra::vect(
    ccm_results_heat, geom = c("longitude", "latitude"), crs = "EPSG:4326"
  )
ccm_results_heat$koppen_code <- 
  as.character(terra::extract(koppen_raster, pts)[,2])
ccm_results_heat$koppen_class <- 
  koppen_lookup[ccm_results_heat$koppen_code]

# 变量富集：绿化投资指标匹配 (独立集成) ----
ccm_results_heat_var <- ccm_results_heat %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  mutate(tp_label = paste0("Lag ", tp))

# 3. 简洁可视化 ----
# A. 气候带分析
p_climate <- ggplot(
  ccm_results_heat_var %>% 
    filter(!is.na(koppen_class), effect_type_heat %in% c("促进", "抑制")),
  aes(x = koppen_class, y = rho_heat_to_sif, fill = effect_type_heat)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~tp_label, nrow = 1) +
  scale_fill_manual(values = c("促进" = "#E41A1C", "抑制" = "#377EB8")) +
  labs(title = "不同气候带下的因果强度分布", x = "柯本气候分类", y = "因果强度 (Rho)", fill = "性质") +
  theme_minimal(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# B. 投资强度分析
p_invest <- ggplot(
  ccm_results_heat_var %>% filter(effect_type_heat %in% c("促进", "抑制")),
  aes(x = effect_type_heat, y = invest_pa_park, fill = effect_type_heat)
) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~tp_label, nrow = 1) +
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
plot_data_rho <- ccm_results_heat_var %>%
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

