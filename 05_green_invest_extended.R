# ============================================================================
# Script: 05_green_invest_extended.R
# Description: 绿化投资扩展分析：时间变化、与SIF关系、阈值分析
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, readr, tidyr, stringr, sf, readxl, showtext,
  segmented, mgcv, patchwork
)
showtext_auto()

# ============================================================================
# 1. 加载数据
# ============================================================================

cat("【1. 加载数据】\n")

# CCM结果
ccm_results_path <- "data_proc/ccm_results_vpd_heat_sif.rds"
if (!file.exists(ccm_results_path)) {
  ccm_csv_path <- "data_proc/ccm_results_vpd_heat_sif.csv"
  if (file.exists(ccm_csv_path)) {
    ccm_results_all <- read_csv(ccm_csv_path, show_col_types = FALSE)
  } else {
    stop("请先运行 main_vpd_index2.R 生成CCM分析结果")
  }
} else {
  ccm_results_all <- readRDS(ccm_results_path)
}

# 站点坐标
station_coords <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  dplyr::select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# SIF数据（用于计算站点平均SIF）
sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# ============================================================================
# 2. 读取完整的绿化投资时间序列数据
# ============================================================================

cat("\n【2. 读取绿化投资时间序列】\n")

green_invest_raw <- read_csv(
  "data_raw/green_invest/中国城市数据.csv",
  show_col_types = FALSE
)

# 转换为长格式
green_invest_long <- green_invest_raw %>%
  pivot_longer(
    cols = starts_with("园林绿化_"),
    names_to = "year",
    values_to = "green_invest"
  ) %>%
  mutate(
    year = as.integer(str_extract(year, "\\d+")),
    green_invest = as.numeric(green_invest),
    city_name = paste0(城市名称, "市")
  ) %>%
  dplyr::select(city_name, year, green_invest) %>%
  filter(!is.na(green_invest), green_invest > 0)

cat("绿化投资数据年份范围:", min(green_invest_long$year), "-", max(green_invest_long$year), "\n")
cat("城市数:", n_distinct(green_invest_long$city_name), "\n")

# ============================================================================
# 3. 绘制绿化投资时间变化图 (分体量等级)
# ============================================================================

cat("\n【3. 绘制绿化投资时间变化图】\n")

# 1. 预处理：剔除2017年，筛选数据较全的城市
# 计算每个城市拥有的有效年份数
city_year_counts <- green_invest_long %>%
  filter(year != 2017) %>%
  group_by(city_name) %>%
  summarise(n_years = n(), .groups = "drop")

# 筛选年份数超过10年的城市（约占总年份的一半以上）
keep_cities <- city_year_counts %>%
  filter(n_years >= 10) %>%
  pull(city_name)

green_invest_filtered <- green_invest_long %>%
  filter(year != 2017, city_name %in% keep_cities)

# 2. 按投资体量分级 (基于城市平均投资额)
city_tiers <- green_invest_filtered %>%
  group_by(city_name) %>%
  summarise(avg_invest = mean(green_invest, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    tier = cut(
      avg_invest,
      breaks = quantile(avg_invest, probs = seq(0, 1, 0.05), na.rm = TRUE), 
      include.lowest = TRUE
    )
  )

# 合并等级信息
plot_data <- green_invest_filtered %>%
  inner_join(city_tiers, by = "city_name") %>%
  mutate(green_invest_亿元 = green_invest / 10000)

# 3. 绘图：分等级展示趋势
p_invest_tiers <- ggplot(plot_data, aes(x = year, y = green_invest_亿元)) +
  geom_line(aes(group = city_name), alpha = 0.2, color = "#4575b4") +
  # 添加等级平均线
  stat_summary(fun = mean, geom = "line", color = "#d73027", linewidth = 1.2) +
  facet_wrap(~tier, scales = "free_y", ncol = 4) +
  labs(
    title = "不同体量等级城市的园林绿化投资趋势",
    subtitle = "已剔除2017年数据，仅保留有效年份 >= 10年的城市 (红线为等级平均值)",
    x = "年份",
    y = "园林绿化投资 (亿元)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, size = 8), 
    panel.spacing = unit(1, "lines")
  )

print(p_invest_tiers)
ggsave("data_proc/green_invest_trend_tiers.png", p_invest_tiers,
       width = 14, height = 10, dpi = 300, scale = 0.4)
cat("等级趋势图已保存: data_proc/green_invest_trend_tiers.png\n")

# 计算全国年度汇总（同样剔除2017）
national_invest <- green_invest_filtered %>%
  group_by(year) %>%
  summarise(
    total_invest = sum(green_invest, na.rm = TRUE) / 10000,
    mean_invest = mean(green_invest, na.rm = TRUE) / 10000,
    n_cities = n(),
    .groups = "drop"
  )

p_invest_trend <- ggplot(national_invest, aes(x = year, y = total_invest)) +
  geom_line(linewidth = 1.2, color = "#2166ac") +
  geom_point(size = 3, color = "#2166ac") +
  labs(
    title = "中国城市园林绿化投资时间变化 (汇总)",
    subtitle = "2002-2024年全国汇总 (已剔除2017年)",
    x = "年份",
    y = "园林绿化投资总额 (亿元)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(p_invest_trend)
ggsave("data_proc/green_invest_trend.png", p_invest_trend,
       width = 12, height = 8, dpi = 300, scale = 0.4)
cat("图1已保存: data_proc/green_invest_trend.png\n")

# ============================================================================
# 4. 空间匹配站点到城市
# ============================================================================

cat("\n【4. 空间匹配站点到城市】\n")

china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326)

stations_sf <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_city <- stations_sf %>%
  st_join(china_cities_shp, join = st_within) %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  rename(city_name = ct_name) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# ============================================================================
# 5. 计算站点平均SIF并与投资合并
# ============================================================================

cat("\n【5. 计算站点SIF与投资关系】\n")

# 计算每个站点的平均SIF
station_sif_mean <- sif_data %>%
  group_by(meteo_stat_id) %>%
  summarise(
    sif_mean = mean(sif, na.rm = TRUE),
    sif_sd = sd(sif, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

# 读取2020年GDP用于计算比例
gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, gdp_total = "地区生产总值(万元)") %>%
  mutate(
    gdp_total = as.numeric(gdp_total)
  )

# 2020年绿化投资
green_invest_2020 <- green_invest_long %>%
  filter(year == 2020)

# 合并数据
sif_invest_data <- station_sif_mean %>%
  left_join(stations_with_city %>% dplyr::select(meteo_stat_id, city_name), by = "meteo_stat_id") %>%
  left_join(green_invest_2020, by = "city_name") %>%
  left_join(gdp_data, by = "city_name") %>%
  filter(!is.na(green_invest), !is.na(sif_mean)) %>%
  mutate(
    green_ratio = green_invest / gdp_total * 100,
    green_invest_亿元 = green_invest / 10000
  ) %>%
  filter(green_ratio > 0, green_ratio < 10)  # 过滤异常值

cat("有效站点数:", nrow(sif_invest_data), "\n")

# ============================================================================
# 6. 绘制绿化投资与SIF关系图
# ============================================================================

cat("\n【6. 绘制投资与SIF关系图】\n")

p_invest_sif <- ggplot(sif_invest_data, aes(x = green_invest_亿元, y = sif_mean)) +
  geom_point(alpha = 0.5, size = 2.5, color = "#4575b4") +
  geom_smooth(method = "loess", color = "#d73027", se = TRUE, linewidth = 1.2) +
  labs(
    title = "城市绿化投资与站点平均SIF的关系",
    subtitle = "2020年绿化投资 vs 站点多年平均SIF",
    x = "园林绿化投资 (亿元)",
    y = "平均SIF"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(p_invest_sif)
ggsave("data_proc/green_invest_sif_scatter.png", p_invest_sif,
       width = 12, height = 9, dpi = 300, scale = 0.4)
cat("图2已保存: data_proc/green_invest_sif_scatter.png\n")

# 绿化占比与SIF
# ============================================================================
# 7. 投资与因果效应强度的阈值分析
# ============================================================================

cat("\n【7. 投资-因果效应强度阈值分析】\n")

# 合并CCM结果并计算基础指标
threshold_base <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id)) %>%
  left_join(stations_with_city %>% dplyr::select(meteo_stat_id, city_name), by = "meteo_stat_id") %>%
  left_join(green_invest_2020, by = "city_name") %>%
  left_join(gdp_data, by = "city_name") %>%
  filter(!is.na(green_invest), !is.na(effect_index_heat_sif)) %>%
  mutate(
    green_ratio = green_invest / gdp_total * 100,
    green_invest_亿元 = green_invest / 10000,
    log_green_ratio = log10(green_ratio),
    log_green_invest = log10(green_invest_亿元),
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      effect_type_heat_sif == "促进效应(+)" ~ "促进",
      effect_type_heat_sif == "抑制效应(-)" ~ "抑制",
      TRUE ~ "其他"
    )
  ) %>%
  filter(green_ratio > 0, green_ratio < 10)

# --- 7.1 促进效应 (Promoting Effect) 深度分析 ---
cat("\n7.1 促进效应分析\n")
promote_data <- threshold_base %>% filter(heat_causality == "促进")

# 7.1.1 分组比较 (Quantile) - 占比
promote_grouped_ratio <- promote_data %>%
  mutate(invest_group = cut(green_ratio, breaks = quantile(green_ratio, probs = seq(0, 1, 0.2), na.rm = TRUE),
                            labels = c("Q1", "Q2", "Q3", "Q4", "Q5"), include.lowest = TRUE)) %>%
  group_by(invest_group) %>%
  summarise(n = n(), mean_eff = mean(effect_index_heat_sif, na.rm = TRUE),
            se_eff = sd(effect_index_heat_sif, na.rm = TRUE) / sqrt(n()), .groups = "drop")

p_quant_p_ratio <- ggplot(promote_grouped_ratio, aes(x = invest_group, y = mean_eff)) +
  geom_col(fill = "#4575b4", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_eff - se_eff, ymax = mean_eff + se_eff), width = 0.3, linewidth = 2) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 20, fontface = "bold") +
  labs(title = "不同投资占比下的平均促进效应强度", x = "占比分位数", y = "平均促进效应强度") +
  theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
ggsave("data_proc/threshold_quantile_promote.png", p_quant_p_ratio, width = 20, height = 15, dpi = 300)

# 7.1.1b 分组比较 (Quantile) - 绝对值
promote_grouped_abs <- promote_data %>%
  mutate(invest_group = cut(green_invest_亿元, breaks = quantile(green_invest_亿元, probs = seq(0, 1, 0.2), na.rm = TRUE),
                            labels = c("Q1", "Q2", "Q3", "Q4", "Q5"), include.lowest = TRUE)) %>%
  group_by(invest_group) %>%
  summarise(n = n(), mean_eff = mean(effect_index_heat_sif, na.rm = TRUE),
            se_eff = sd(effect_index_heat_sif, na.rm = TRUE) / sqrt(n()), .groups = "drop")

p_quant_p_abs <- ggplot(promote_grouped_abs, aes(x = invest_group, y = mean_eff)) +
  geom_col(fill = "#4575b4", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_eff - se_eff, ymax = mean_eff + se_eff), width = 0.3, linewidth = 2) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 20, fontface = "bold") +
  labs(title = "不同投资额(亿元)下的平均促进效应强度", x = "金额分位数", y = "平均促进效应强度") +
  theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
ggsave("data_proc/threshold_quantile_promote_abs.png", p_quant_p_abs, width = 20, height = 15, dpi = 300)

# 7.1.2 对数化分段回归 (占比)
try({
  lm_p <- lm(effect_index_heat_sif ~ log_green_ratio, data = promote_data)
  seg_p <- segmented(lm_p, seg.Z = ~log_green_ratio, npsi = 1)
  bp_p <- 10^seg_p$psi[1, "Est."]
  promote_data$fitted_log <- predict(seg_p)
  p_seg_log_p <- ggplot(promote_data, aes(x = log_green_ratio, y = effect_index_heat_sif)) +
    geom_point(alpha = 0.4, color = "#4575b4", size = 8) + geom_line(aes(y = fitted_log), color = "black", linewidth = 4) +
    geom_vline(xintercept = log10(bp_p), linetype = "dashed", color = "#d73027", linewidth = 3) +
    annotate("text", x = log10(bp_p) + 0.1, y = max(promote_data$effect_index_heat_sif)*0.8, label = paste0("阈值 ≈ ", round(bp_p, 3), "%"), color = "#d73027", size = 25, fontface = "bold") +
    labs(title = "促进效应阈值 (对数占比尺度)", x = "log10(占比 %)", y = "促进效应强度") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_segmented_log_promote.png", p_seg_log_p, width = 20, height = 15, dpi = 300)
})

# 7.1.2c 线性分段回归 (占比) - 补充分析
try({
  lm_p_lin <- lm(effect_index_heat_sif ~ green_ratio, data = promote_data)
  seg_p_lin <- segmented(lm_p_lin, seg.Z = ~green_ratio, npsi = 1)
  bp_p_lin <- seg_p_lin$psi[1, "Est."]
  promote_data$fitted_lin <- predict(seg_p_lin)
  p_seg_lin_p <- ggplot(promote_data, aes(x = green_ratio, y = effect_index_heat_sif)) +
    geom_point(alpha = 0.4, color = "#4575b4", size = 8) + geom_line(aes(y = fitted_lin), color = "black", linewidth = 4) +
    geom_vline(xintercept = bp_p_lin, linetype = "dashed", color = "#d73027", linewidth = 3) +
    annotate("text", x = bp_p_lin + 0.5, y = max(promote_data$effect_index_heat_sif)*0.8, label = paste0("阈值 ≈ ", round(bp_p_lin, 2), "%"), color = "#d73027", size = 25, fontface = "bold") +
    labs(title = "促进效应阈值 (原始占比尺度)", x = "绿化投资占比 (%)", y = "促进效应强度") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_segmented_linear_promote.png", p_seg_lin_p, width = 20, height = 15, dpi = 300)
})

# 7.1.2b 对数化分段回归 (金额)
try({
  lm_p_abs <- lm(effect_index_heat_sif ~ log_green_invest, data = promote_data)
  seg_p_abs <- segmented(lm_p_abs, seg.Z = ~log_green_invest, npsi = 1)
  bp_p_abs <- 10^seg_p_abs$psi[1, "Est."]
  promote_data$fitted_log_abs <- predict(seg_p_abs)
  p_seg_log_p_abs <- ggplot(promote_data, aes(x = log_green_invest, y = effect_index_heat_sif)) +
    geom_point(alpha = 0.4, color = "#4575b4", size = 8) + geom_line(aes(y = fitted_log_abs), color = "black", linewidth = 4) +
    geom_vline(xintercept = log10(bp_p_abs), linetype = "dashed", color = "#d73027", linewidth = 3) +
    annotate("text", x = log10(bp_p_abs) + 0.1, y = max(promote_data$effect_index_heat_sif)*0.8, label = paste0("阈值 ≈ ", round(bp_p_abs, 2), "亿"), color = "#d73027", size = 25, fontface = "bold") +
    labs(title = "促进效应阈值 (对数金额尺度)", x = "log10(金额 亿元)", y = "促进效应强度") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_segmented_log_promote_abs.png", p_seg_log_p_abs, width = 20, height = 15, dpi = 300)
})

# 7.1.2d 线性分段回归 (金额) - 补充分析
try({
  lm_p_abs_lin <- lm(effect_index_heat_sif ~ green_invest_亿元, data = promote_data)
  seg_p_abs_lin <- segmented(lm_p_abs_lin, seg.Z = ~green_invest_亿元, npsi = 1)
  bp_p_abs_lin <- seg_p_abs_lin$psi[1, "Est."]
  promote_data$fitted_lin_abs <- predict(seg_p_abs_lin)
  p_seg_lin_p_abs <- ggplot(promote_data, aes(x = green_invest_亿元, y = effect_index_heat_sif)) +
    geom_point(alpha = 0.4, color = "#4575b4", size = 8) + geom_line(aes(y = fitted_lin_abs), color = "black", linewidth = 4) +
    geom_vline(xintercept = bp_p_abs_lin, linetype = "dashed", color = "#d73027", linewidth = 3) +
    annotate("text", x = bp_p_abs_lin + 2, y = max(promote_data$effect_index_heat_sif)*0.8, label = paste0("阈值 ≈ ", round(bp_p_abs_lin, 1), "亿"), color = "#d73027", size = 25, fontface = "bold") +
    labs(title = "促进效应阈值 (原始金额尺度)", x = "绿化投资额 (亿元)", y = "促进效应强度") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_segmented_linear_promote_abs.png", p_seg_lin_p_abs, width = 20, height = 15, dpi = 300)
})

# --- 7.2 抑制效应 (Inhibiting Effect) 深度分析 ---
cat("\n7.2 抑制效应分析\n")
inhibit_data <- threshold_base %>% filter(heat_causality == "抑制") %>% mutate(abs_effect = abs(effect_index_heat_sif))

# 7.2.1 分组比较 (Quantile) - 占比
# ... (unchanged)
# 7.2.2 对数化GAM (抑制效应 - 占比)
if (nrow(inhibit_data) > 30) {
  gam_i <- gam(abs_effect ~ s(log_green_ratio, k = 5), data = inhibit_data)
  pred_i <- data.frame(log_green_ratio = seq(min(inhibit_data$log_green_ratio), max(inhibit_data$log_green_ratio), length.out = 100))
  pred_i$fitted <- predict(gam_i, newdata = pred_i)
  p_gam_log_i <- ggplot(inhibit_data, aes(x = log_green_ratio, y = abs_effect)) +
    geom_point(alpha = 0.4, color = "#d73027", size = 8) + geom_line(data = pred_i, aes(y = fitted), color = "black", linewidth = 4) +
    labs(title = "抑制强度随对数占比的变化 (GAM)", x = "log10(占比 %)", y = "抑制强度 |效应指数|") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_gam_log_inhibit.png", p_gam_log_i, width = 20, height = 15, dpi = 300)
  
  # 7.2.2c 线性GAM (占比) - 补充分析
  gam_i_lin <- gam(abs_effect ~ s(green_ratio, k = 5), data = inhibit_data)
  pred_i_lin <- data.frame(green_ratio = seq(min(inhibit_data$green_ratio), max(inhibit_data$green_ratio), length.out = 100))
  pred_i_lin$fitted <- predict(gam_i_lin, newdata = pred_i_lin)
  p_gam_lin_i <- ggplot(inhibit_data, aes(x = green_ratio, y = abs_effect)) +
    geom_point(alpha = 0.4, color = "#d73027", size = 8) + geom_line(data = pred_i_lin, aes(y = fitted), color = "black", linewidth = 4) +
    labs(title = "抑制强度随原始占比的变化 (GAM)", x = "绿化投资占比 (%)", y = "抑制强度 |效应指数|") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_gam_linear_inhibit.png", p_gam_lin_i, width = 20, height = 15, dpi = 300)
}

# 7.2.2b 对数化GAM (抑制效应 - 金额)
if (nrow(inhibit_data) > 30) {
  # ... (existing log GAM for absolute value)
  # 7.2.2d 线性GAM (金额) - 补充分析
  gam_i_abs_lin <- gam(abs_effect ~ s(green_invest_亿元, k = 5), data = inhibit_data)
  pred_i_abs_lin <- data.frame(green_invest_亿元 = seq(min(inhibit_data$green_invest_亿元), max(inhibit_data$green_invest_亿元), length.out = 100))
  pred_i_abs_lin$fitted <- predict(gam_i_abs_lin, newdata = pred_i_abs_lin)
  p_gam_lin_i_abs <- ggplot(inhibit_data, aes(x = green_invest_亿元, y = abs_effect)) +
    geom_point(alpha = 0.4, color = "#d73027", size = 8) + geom_line(data = pred_i_abs_lin, aes(y = fitted), color = "black", linewidth = 4) +
    labs(title = "抑制强度随原始金额的变化 (GAM)", x = "绿化投资额 (亿元)", y = "抑制强度 |效应指数|") +
    theme_minimal(base_size = 60) + theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5), axis.title = element_text(face = "bold", size = 60), axis.text = element_text(size = 50))
  ggsave("data_proc/threshold_gam_linear_inhibit_abs.png", p_gam_lin_i_abs, width = 20, height = 15, dpi = 300)
}

cat("\n分析完成!\n")
