# weekly_post_analysis.R
# 仿照 main.R，对周度 CCM 结果进行多维度后续分析 (增强大字号版)

# Preparation ----
pacman::p_load(
  dplyr, ggplot2, purrr, tidyr, showtext, sf, rnaturalearth, 
  patchwork, readxl, terra, targets, ggalluvial, segmented
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
  
  pair_data <- data %>% filter(tp %in% pair) %>% dplyr::select(meteo_stat_id, tp, effect_type) %>%
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
    dplyr::select(effect_type, abs_effect, invest_pa_tot, invest_pa_built, invest_pa_park) %>%
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

# 组合所有可用 Tp 的图 (动态生成)
unique_lags <- sort(unique(results_weekly_var$tp))

cat("- 正在生成原始尺度综合图...\n")
plots_raw_list <- map(unique_lags, ~plot_causal_dist_weekly(results_weekly_var %>% filter(tp == .x), .x, FALSE))
combined_dist_raw <- wrap_plots(plots_raw_list, ncol = 1) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "周度因果强度随投资区间分布 (原始尺度)", 
                  theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5), legend.position = "bottom"))

# 动态计算高度：每个 lag 约 15 英寸
target_height <- length(unique_lags) * 15
ggsave("data_proc/weekly_post_dist_combined_raw.png", combined_dist_raw, width = 30, height = target_height, dpi = 100, limitsize = FALSE)

cat("- 正在生成对数尺度综合图...\n")
plots_log_list <- map(unique_lags, ~plot_causal_dist_weekly(results_weekly_var %>% filter(tp == .x), .x, TRUE))
combined_dist_log <- wrap_plots(plots_log_list, ncol = 1) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "周度因果强度随投资区间分布 (对数尺度)", 
                  theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5), legend.position = "bottom"))

ggsave("data_proc/weekly_post_dist_combined_log.png", combined_dist_log, width = 30, height = target_height, dpi = 100, limitsize = FALSE)

# ============================================================================
# 4. 空间分布分析 (China Map) ----
# ============================================================================
cat("【4. 空间分布分析 (优化版)】\n")

china_map <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

p_spatial <- ggplot() +
  geom_sf(data = china_map, fill = "gray95", color = "gray70", linewidth = 1) +
  geom_point(data = results_weekly_var, 
             aes(x = longitude, y = latitude, color = effect_type), 
             size = 0.3, alpha = 0.3) + # 减小点大小 (相对大字号而言)
  facet_wrap(~ tp_label, ncol = 3) +
  scale_color_manual(values = c("促进"="#E41A1C", "抑制"="#377EB8", "无因果"="#999999")) +
  labs(title = "中国城市周度热事件因果性质空间分布",
       color = "因果性质") +
  theme_minimal(base_size = BASE_FONT_SIZE) + 
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "bottom")

ggsave("data_proc/weekly_post_spatial_causality.png", p_spatial, width = 40, height = 15, dpi = 100)

# ============================================================================
# 5. 断点回归：识别投资阈值 (仅限 Lag 0 促进效应) ----
# ============================================================================
cat("【5. 断点回归：Lag 0 投资阈值识别 (Log-Log)】\n")

# 筛选 Lag 0 的促进作用数据
df_lag0_promote <- results_weekly_var %>%
  filter(tp == 0, effect_type == "促进", abs_effect > 0)

plot_breakpoint_weekly <- function(df_input, v_name, label) {
  # 准备数据并取 Log
  df_plot <- df_input %>%
    dplyr::select(x = !!sym(v_name), y = abs_effect) %>%
    filter(!is.na(x), x > 0, y > 0) %>%
    mutate(x_log = log10(x), y_log = log10(y)) %>%
    arrange(x_log)
  
  n_samples <- nrow(df_plot)
  if (n_samples < 15) {
    return(ggplot() + theme_void() + labs(subtitle = paste0(label, "\n(样本数不足: ", n_samples, ")")))
  }
  
  # 1. 基础线性模型
  fit_lm <- lm(y_log ~ x_log, data = df_plot)
  
  # 2. 尝试断点回归
  fit_seg <- try(segmented(fit_lm, seg.Z = ~x_log, npsi = 1), silent = TRUE)
  
  p <- ggplot(df_plot, aes(x = x_log, y = y_log)) +
    geom_point(alpha = 0.4, color = "gray50", size = 4) +
    labs(title = label, x = "log10(投资强度)", y = "log10(促进强度)") +
    theme_minimal(base_size = BASE_FONT_SIZE) +
    theme(plot.title = element_text(face="bold", hjust = 0.5),
          axis.title = element_text(face="bold"))

  if (!inherits(fit_seg, "try-error")) {
    # 提取断点
    psi_log <- fit_seg$psi[1, "Est."]
    threshold_raw <- 10^psi_log
    df_plot$fitted <- predict(fit_seg)
    
    p <- p + 
      geom_line(data = df_plot, aes(y = fitted), color = "blue", linewidth = 3) +
      geom_vline(xintercept = psi_log, color = "#d73027", linetype = "dashed", linewidth = 2) +
      annotate("label", x = psi_log, y = max(df_plot$y_log, na.rm=T), 
               label = paste0("T=", round(threshold_raw, 1)), 
               color="#d73027", fontface="bold", size = BASE_FONT_SIZE/4, fill = "white")
  } else {
    # 拟合失败：提取并输出失败信息
    err_msg <- attr(fit_seg, "condition")$message
    cat(paste0("\n[断点回归失败诊断 - ", label, "]\n"))
    cat(paste0("  - 错误信息: ", err_msg, "\n"))
    cat(paste0("  - 有效样本量: ", n_samples, "\n"))
    cat("  - 线性模型摘要:\n")
    print(summary(fit_lm))
    
    # 在图中显示失败原因
    p <- p + geom_smooth(method = "lm", color = "black", linetype = "dotted", linewidth = 2, se = FALSE) +
      labs(caption = paste0("断点回归未收敛: ", err_msg, " (n=", n_samples, ")")) +
      theme(plot.caption = element_text(size = BASE_FONT_SIZE/2, color = "red", hjust = 0))
  }
  
  return(p)
}

cat("- 正在计算各投资维度的断点...\n")
p_thr_tot   <- plot_breakpoint_weekly(df_lag0_promote, "invest_pa_tot",   "单位总绿地投资 (Lag 0)")
p_thr_built <- plot_breakpoint_weekly(df_lag0_promote, "invest_pa_built", "单位建成区投资 (Lag 0)")
p_thr_park  <- plot_breakpoint_weekly(df_lag0_promote, "invest_pa_park",  "单位公园投资 (Lag 0)")

# 组合三张图
combined_thresholds <- (p_thr_tot / p_thr_built / p_thr_park) +
  plot_annotation(title = "Lag 0 促进效应随投资强度的突变阈值分析",
                  theme = theme(plot.title = element_text(size = BASE_FONT_SIZE*1.2, face="bold", hjust=0.5)))

ggsave("data_proc/weekly_post_threshold_lag0.png", combined_thresholds, width = 25, height = 40, dpi = 100)

cat("\n所有周度后续分析（包括断点回归）已完成！请查看 data_proc/ 下的图片。\n")
