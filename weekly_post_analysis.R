# weekly_post_analysis.R
# 仿照 main.R，对周度 CCM 结果进行多维度后续分析 (增强大字号版)

# Preparation ----
pacman::p_load(
  dplyr, ggplot2, purrr, tidyr, showtext, sf, rnaturalearth, 
  patchwork, readxl, terra, targets, ggalluvial
)
showtext_auto()

# 设定超大字号基准 (5倍放大)
BASE_FONT_SIZE <- 60 

# 1. 数据加载与匹配 ----
cat("【1. 加载周度数据并匹配投资与气候带】\n")

# 加载 CCM 结果
if (!file.exists("data_proc/results_weekly_test_20.rds")) {
  stop("找不到周度测试结果文件，请先运行 main_weekly.R")
}
results_weekly <- readRDS("data_proc/results_weekly_test_20.rds") %>%
  mutate(tp_label = paste0("Lag ", tp, " Week"))

# 匹配气候带 (Koppen)
cat("- 匹配气候带...\n")
koppen_raster <- terra::rast("data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")
koppen_lookup <- c(
  "1"="Af","2"="Am","3"="As","4"="Aw","5"="BWh","6"="BWk","7"="BSh","8"="BSk",
  "9"="Csa","10"="Csb","11"="Csc","12"="Cwa","13"="Cwb","14"="Cwc","15"="Cfa",
  "16"="Cfb","17"="Cfc","18"="Dsa","19"="Dsb","20"="Dsc","21"="Dsd","22"="Dwa",
  "23"="Dwb","24"="Dwc","25"="Dwd","26"="Dfa","27"="Dfb","28"="Dfc","29"="Dfd","30"="ET","31"="EF"
)

# 提取经纬度
tar_load(station_coords)
results_weekly <- results_weekly %>% left_join(station_coords, by = "meteo_stat_id")

pts <- terra::vect(results_weekly, geom = c("longitude", "latitude"), crs = "EPSG:4326")
results_weekly$koppen_class <- koppen_lookup[as.character(terra::extract(koppen_raster, pts)[,2])]
results_weekly$koppen_group <- substr(results_weekly$koppen_class, 1, 1)

# 匹配投资指标
cat("- 匹配投资与省份城市信息...\n")
tar_load(invest_metrics)
results_weekly_var <- results_weekly %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  mutate(abs_effect = abs(mean_coef))

# ============================================================================
# 2. 桑基图分析 (因果性质演变) ----
# ============================================================================
cat("【2. 因果性质演变桑基图】\n")

available_tps <- sort(unique(results_weekly_var$tp))
tp_pairs <- map(1:(length(available_tps) - 1), ~ c(available_tps[.x], available_tps[.x+1]))

plot_sankey_pair_weekly <- function(data, pair) {
  lvls <- c("促进", "抑制", "无因果", "未知")
  cols <- c("促进" = "#E41A1C", "抑制" = "#377EB8", "无因果" = "#999999", "未知" = "#FF7F00")
  tp_start <- pair[1]; tp_end   <- pair[2]
  
  pair_data <- data %>% filter(tp %in% pair) %>% select(meteo_stat_id, tp, effect_type) %>%
    pivot_wider(id_cols = meteo_stat_id, names_from = tp, values_from = effect_type) %>% drop_na() %>%
    mutate(start_val = factor(as.character(get(as.character(tp_start))), levels = lvls),
           end_val = factor(as.character(get(as.character(tp_end))), levels = lvls))
  
  pair_summary <- pair_data %>% group_by(start_val, end_val) %>% summarise(n = n(), .groups = "drop")
  
  ggplot(pair_summary, aes(axis1 = start_val, axis2 = end_val, y = n)) +
    geom_alluvium(aes(fill = start_val), width = 1/12, alpha = 0.6) +
    geom_stratum(width = 1/6, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = BASE_FONT_SIZE/4) +
    scale_x_discrete(limits = c(paste0("tp=", tp_start), paste0("tp=", tp_end)), expand = c(.1, .1)) +
    scale_fill_manual(values = cols) + 
    theme_minimal(base_size = BASE_FONT_SIZE) + 
    theme(legend.position = "none", panel.grid = element_blank(),
          axis.text.y = element_blank(), axis.title = element_blank())
}

if (length(tp_pairs) > 0) {
  all_sankey_plots <- map(tp_pairs, ~plot_sankey_pair_weekly(results_weekly_var, .x))
  combined_sankey <- wrap_plots(all_sankey_plots, nrow = 1) + 
    plot_annotation(title = "周度因果性质演变", theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5)))
  ggsave("data_proc/weekly_post_causal_sankey.png", combined_sankey, width = 30, height = 20, dpi = 150)
}

# ============================================================================
# 3. 投资区间分布分析 (Combined Boxplot) ----
# ============================================================================
cat("【3. 投资区间分布综合图 (Raw & Log)】\n")

plot_causal_dist_weekly <- function(df_input, tp_val, use_log = FALSE) {
  df_long <- df_input %>%
    filter(effect_type %in% c("促进", "抑制")) %>%
    select(effect_type, abs_effect, invest_pa_tot, invest_pa_built, invest_pa_park) %>%
    pivot_longer(cols = starts_with("invest_pa_"), names_to = "invest_var", values_to = "invest_val") %>%
    mutate(invest_label = case_when(
      invest_var == "invest_pa_tot" ~ "总绿地",
      invest_var == "invest_pa_built" ~ "建成区",
      invest_var == "invest_pa_park" ~ "公园"
    )) %>%
    drop_na(invest_val) %>%
    group_by(invest_label) %>%
    mutate(invest_bin = factor(paste0("Q", cut(invest_val, 
                                               breaks = unique(quantile(invest_val, probs = seq(0, 1, 0.2), na.rm = TRUE)),
                                               include.lowest = TRUE, labels = FALSE)))) %>%
    ungroup()

  p <- ggplot(df_long, aes(x = invest_bin, y = abs_effect, fill = effect_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 2) +
    geom_jitter(aes(color = effect_type), width = 0.2, alpha = 0.3, size = 3) +
    facet_grid(effect_type ~ invest_label, scales = "free_y") +
    scale_fill_manual(values = c("促进"="#E41A1C", "抑制"="#377EB8")) +
    scale_color_manual(values = c("促进"="#E41A1C", "抑制"="#377EB8")) +
    labs(subtitle = paste0("Lag ", tp_val, " Week"),
         x = "投资强度区间", y = ifelse(use_log, "log10(因果强度)", "因果强度")) +
    theme_minimal(base_size = BASE_FONT_SIZE) + 
    theme(legend.position = "none", strip.text = element_text(face="bold"),
          axis.title = element_text(face="bold"), plot.subtitle = element_text(hjust = 0.5, face="bold"))
  
  if (use_log) p <- p + scale_y_log10()
  return(p)
}

# 组合 Tp 0-2 的图 (原始尺度)
p0 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 0), 0, FALSE)
p1 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 1), 1, FALSE)
p2 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 2), 2, FALSE)

combined_dist_raw <- (p0 / p1 / p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "周度因果强度随投资区间分布 (原始尺度)", 
                  theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5), legend.position = "bottom"))
ggsave("data_proc/weekly_post_dist_combined_raw.png", combined_dist_raw, width = 30, height = 45, dpi = 100)

# 组合 Tp 0-2 的图 (对数尺度)
pl0 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 0), 0, TRUE)
pl1 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 1), 1, TRUE)
pl2 <- plot_causal_dist_weekly(results_weekly_var %>% filter(tp == 2), 2, TRUE)

combined_dist_log <- (pl0 / pl1 / pl2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "周度因果强度随投资区间分布 (对数尺度)", 
                  theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5), legend.position = "bottom"))
ggsave("data_proc/weekly_post_dist_combined_log.png", combined_dist_log, width = 30, height = 45, dpi = 100)

# ============================================================================
# 4. 空间分布分析 (China Map) ----
# ============================================================================
cat("【4. 空间分布分析 (优化版)】\n")

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

p_spatial <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 1) +
  geom_point(data = results_weekly_var, 
             aes(x = longitude, y = latitude, color = effect_type), 
             size = 4, alpha = 0.8) + # 减小点大小 (相对大字号而言)
  facet_wrap(~ tp_label, ncol = 3) +
  scale_color_manual(values = c("促进"="#E41A1C", "抑制"="#377EB8", "无因果"="#999999")) +
  labs(title = "中国城市周度热事件因果性质空间分布",
       color = "因果性质") +
  theme_minimal(base_size = BASE_FONT_SIZE) + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "bottom")

ggsave("data_proc/weekly_post_spatial_causality.png", p_spatial, width = 40, height = 15, dpi = 100)

cat("\n所有周度后续分析已完成！请查看 data_proc/ 下的图片。\n")
