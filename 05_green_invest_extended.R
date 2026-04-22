# ============================================================================
# Script: 05_green_invest_extended.R
# Description: 绿化投资扩展分析：时间变化、与SIF关系、阈值分析
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, readr, tidyr, stringr, sf, readxl, showtext,
  segmented, mgcv, patchwork, purrr
)
showtext_auto()

# ============================================================================
# 1. 加载数据
# ============================================================================

cat("【1. 加载数据】\n")

# CCM结果 (恢复原始逻辑)
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

# SIF数据
sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# ============================================================================
# 2. 读取绿化投资时间序列数据
# ============================================================================

cat("\n【2. 读取绿化投资时间序列】\n")

green_invest_raw <- read_csv(
  "data_raw/green_invest/中国城市数据.csv",
  show_col_types = FALSE
)

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

# ============================================================================
# 3. 绘制绿化投资时间变化图 (分体量等级)
# ============================================================================

cat("\n【3. 绘制绿化投资时间变化图】\n")

city_year_counts <- green_invest_long %>%
  filter(year != 2017) %>%
  group_by(city_name) %>%
  summarise(n_years = n(), .groups = "drop")

keep_cities <- city_year_counts %>% filter(n_years >= 10) %>% pull(city_name)
green_invest_filtered <- green_invest_long %>% filter(year != 2017, city_name %in% keep_cities)

city_tiers <- green_invest_filtered %>%
  group_by(city_name) %>%
  summarise(avg_invest = mean(green_invest, na.rm = TRUE), .groups = "drop") %>%
  mutate(tier = cut(avg_invest, breaks = quantile(avg_invest, probs = seq(0, 1, 0.05), na.rm = TRUE), include.lowest = TRUE))

plot_data <- green_invest_filtered %>%
  inner_join(city_tiers, by = "city_name") %>%
  mutate(green_invest_亿元 = green_invest / 10000)

p_invest_tiers <- ggplot(plot_data, aes(x = year, y = green_invest_亿元)) +
  geom_line(aes(group = city_name), alpha = 0.2, color = "#4575b4") +
  stat_summary(fun = mean, geom = "line", color = "#d73027", linewidth = 1.2) +
  facet_wrap(~tier, scales = "free_y", ncol = 4) +
  labs(title = "不同体量等级城市的园林绿化投资趋势", x = "年份", y = "投资 (亿元)") +
  theme_minimal(base_size = 16)

ggsave("data_proc/green_invest_trend_tiers.png", p_invest_tiers, width = 14, height = 10, dpi = 300, scale = 0.4)

# 全国年度汇总
national_invest <- green_invest_filtered %>%
  group_by(year) %>%
  summarise(total_invest = sum(green_invest, na.rm = TRUE) / 10000, .groups = "drop")

p_invest_trend <- ggplot(national_invest, aes(x = year, y = total_invest)) +
  geom_line(linewidth = 1.2, color = "#2166ac") +
  geom_point(size = 3, color = "#2166ac") +
  labs(title = "中国城市园林绿化投资时间变化 (汇总)", x = "年份", y = "总额 (亿元)") +
  theme_minimal(base_size = 18)

ggsave("data_proc/green_invest_trend.png", p_invest_trend, width = 12, height = 8, dpi = 300, scale = 0.4)

# ============================================================================
# 4. 空间匹配站点到城市
# ============================================================================

cat("\n【4. 空间匹配站点到城市】\n")
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326)
stations_sf <- station_coords %>% filter(!is.na(longitude), !is.na(latitude)) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
stations_with_city <- stations_sf %>% st_join(china_cities_shp, join = st_within) %>% as.data.frame() %>% 
  dplyr::select(meteo_stat_id, city_name = ct_name) %>% mutate(meteo_stat_id = as.character(meteo_stat_id))

# ============================================================================
# 5. 计算站点平均SIF并整合多维度投入 (2020年)
# ============================================================================

cat("\n【5. 整合多维度投入指标】\n")

station_sif_mean <- sif_data %>%
  group_by(meteo_stat_id) %>%
  summarise(sif_mean = mean(sif, na.rm = TRUE), .groups = "drop")

gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, gdp_total = "地区生产总值(万元)") %>%
  mutate(gdp_total = as.numeric(gdp_total))

# 读取2020年三类绿地面积 (新维度)
green_area_2020 <- read_excel("data_raw/城市绿地面积数据_2003-2023.xlsx", sheet = "绿地面积_明细数据") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, area_green = `绿地面积(公顷)`, 
                area_built = `建成区绿地面积(公顷)`, area_park = `公园绿地面积(公顷)`) %>%
  mutate(across(starts_with("area_"), as.numeric),
         city_name = ifelse(str_detect(city_name, "市$"), city_name, paste0(city_name, "市")))

green_invest_2020 <- green_invest_long %>% filter(year == 2020)

# 合并底表
analysis_base <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id)) %>%
  left_join(stations_with_city, by = "meteo_stat_id") %>%
  left_join(green_invest_2020, by = "city_name") %>%
  left_join(gdp_data, by = "city_name") %>%
  left_join(green_area_2020, by = "city_name") %>%
  filter(!is.na(green_invest), !is.na(effect_index_heat_sif)) %>%
  mutate(
    green_ratio = green_invest / gdp_total * 100,
    green_invest_亿元 = green_invest / 10000,
    intensity_green = green_invest / area_green,
    intensity_built = green_invest / area_built,
    intensity_park = green_invest / area_park,
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      effect_type_heat_sif == "促进效应(+)" ~ "促进",
      effect_type_heat_sif == "抑制效应(-)" ~ "抑制",
      TRUE ~ "其他"
    )
  ) %>%
  filter(green_ratio > 0, green_ratio < 10)

# ============================================================================
# 6. 绘制多维度投资与SIF关系图
# ============================================================================

cat("\n【6. 绘制关系散点图】\n")
sif_invest_data <- station_sif_mean %>% inner_join(analysis_base, by = "meteo_stat_id")

p_invest_sif <- ggplot(sif_invest_data, aes(x = green_invest_亿元, y = sif_mean)) +
  geom_point(alpha = 0.5, size = 2.5, color = "#4575b4") +
  geom_smooth(method = "loess", color = "#d73027", se = TRUE, linewidth = 1.2) +
  labs(title = "总投资 vs SIF", x = "总投资 (亿元)", y = "平均SIF") + theme_minimal()

p_intensity_sif <- ggplot(sif_invest_data %>% filter(!is.na(intensity_green)), aes(x = intensity_green, y = sif_mean)) +
  geom_point(alpha = 0.5, size = 2.5, color = "#4daf4a") +
  geom_smooth(method = "loess", color = "#d73027", se = TRUE, linewidth = 1.2) + scale_x_log10() +
  labs(title = "投资强度 vs SIF", x = "万元/公顷 (log)", y = "平均SIF") + theme_minimal()

ggsave("data_proc/green_invest_sif_extended.png", p_invest_sif + p_intensity_sif, width = 12, height = 6, dpi = 300)

# ============================================================================
# 7. 全维度阈值分析 (超大字号版)
# ============================================================================

cat("\n【7. 阈值分析：多维度对比 (超大字号)】\n")
promote_data <- analysis_base %>% filter(heat_causality == "促进", effect_index_heat_sif > 0)

plot_threshold_massive <- function(data, x_var, x_label, is_log = FALSE) {
  d <- data %>% filter(!is.na(.data[[x_var]]), .data[[x_var]] > 0)
  if (is_log) {
    d$x_p <- log10(d[[x_var]]); d$y_p <- log10(d$effect_index_heat_sif)
    lab_x <- paste0("log10(", x_label, ")"); lab_y <- "log10(促进强度)"
  } else {
    d$x_p <- d[[x_var]]; d$y_p <- d$effect_index_heat_sif
    lab_x <- x_label; lab_y <- "促进强度"
  }
  
  p <- ggplot(d, aes(x = x_p, y = y_p)) +
    geom_point(alpha = 0.3, color = "gray50", size = 6) +
    labs(title = x_label, x = lab_x, y = lab_y) +
    theme_minimal(base_size = 60) + 
    theme(plot.title = element_text(face = "bold", size = 80, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 65), axis.text = element_text(size = 50))
  
  try({
    fit_seg <- segmented(lm(y_p ~ x_p, data = d), seg.Z = ~x_p, npsi = 1)
    d$fitted <- predict(fit_seg); psi_val <- fit_seg$psi[1, "Est."]; thresh <- if(is_log) 10^psi_val else psi_val
    p <- p + geom_line(data = d, aes(y = fitted), color = "black", linewidth = 4) +
      geom_vline(xintercept = psi_val, linetype = "dashed", color = "#d73027", linewidth = 3) +
      annotate("text", x = psi_val, y = max(d$y_p)*0.9, label = paste0("T=", round(thresh, 1)), color = "#d73027", fontface = "bold", size = 25)
  }, silent = TRUE)
  return(p)
}

vars <- c("green_invest_亿元", "green_ratio", "intensity_green", "intensity_built", "intensity_park")
labs <- c("绿化投资", "投资占比", "单位绿地投资", "单位建成区投资", "单位公园投资")

cat("生成组合图...\n")
p_log_all <- wrap_plots(map2(vars, labs, ~plot_threshold_massive(promote_data, .x, .y, is_log = TRUE)), ncol = 3) + 
  plot_annotation(title = "多维度投入阈值分析 (Log-Log)", theme = theme(plot.title = element_text(size = 100, face = "bold", hjust = 0.5)))
ggsave("data_proc/threshold_segmented_loglog_all.png", p_log_all, width = 30, height = 20, dpi = 150)

p_lin_all <- wrap_plots(map2(vars, labs, ~plot_threshold_massive(promote_data, .x, .y, is_log = FALSE)), ncol = 3) + 
  plot_annotation(title = "多维度投入阈值分析 (线性)", theme = theme(plot.title = element_text(size = 100, face = "bold", hjust = 0.5)))
ggsave("data_proc/threshold_segmented_linear_all.png", p_lin_all, width = 30, height = 20, dpi = 150)

# ============================================================================
# 8. 抑制效应分析 (Section 7.2)
# ============================================================================

cat("\n【8. 抑制效应分析】\n")
inhibit_data <- analysis_base %>% filter(heat_causality == "抑制") %>% mutate(abs_effect = abs(effect_index_heat_sif))
if (nrow(inhibit_data) > 30) {
  gam_i <- gam(abs_effect ~ s(green_ratio, k = 5), data = inhibit_data)
  png("data_proc/threshold_gam_linear_inhibit.png", width = 800, height = 600)
  plot(gam_i, main = "抑制强度随投资占比的变化 (GAM)")
  dev.off()
}

cat("\n分析完成！逻辑已还原，阈值图字体已放大。\n")
