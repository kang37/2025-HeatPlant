# ============================================================================
# Script: 05_green_invest_extended.R
# Description: 绿化投资全维度分析：多维度投入强度、SIF关系及阈值分析
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

# 站点与SIF数据
station_data_raw <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat))

station_coords <- station_data_raw %>%
  dplyr::select(meteo_stat_id, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE)

station_sif_mean <- station_data_raw %>%
  filter(!is.na(sif)) %>%
  group_by(meteo_stat_id) %>%
  summarise(sif_mean = mean(sif, na.rm = TRUE), .groups = "drop")

# ============================================================================
# 2. 读取绿化投资与多维度绿地面积
# ============================================================================

cat("\n【2. 加载投入数据 (2020年)】\n")

# 2.1 绿化投资时间序列
green_invest_raw <- read_csv("data_raw/green_invest/中国城市数据.csv", show_col_types = FALSE)
green_invest_long <- green_invest_raw %>%
  pivot_longer(cols = starts_with("园林绿化_"), names_to = "year", values_to = "invest") %>%
  mutate(year = as.integer(str_extract(year, "\\d+")), invest = as.numeric(invest),
         city_name = paste0(城市名称, "市")) %>%
  filter(!is.na(invest), invest > 0)

green_invest_2020 <- green_invest_long %>% filter(year == 2020)

# 2.2 多维度绿地面积
green_area_2020 <- read_excel("data_raw/城市绿地面积数据_2003-2023.xlsx", sheet = "绿地面积_明细数据") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, area_green = `绿地面积(公顷)`, 
                area_built = `建成区绿地面积(公顷)`, area_park = `公园绿地面积(公顷)`) %>%
  mutate(across(starts_with("area_"), as.numeric),
         city_name = ifelse(str_detect(city_name, "市$"), city_name, paste0(city_name, "市")))

# 2.3 GDP (用于占比计算)
gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  dplyr::select(city_name = 城市, gdp_total = "地区生产总值(万元)") %>%
  mutate(gdp_total = as.numeric(gdp_total))

# ============================================================================
# 3. 空间匹配与指标整合
# ============================================================================

cat("\n【3. 整合指标】\n")

china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326)
stations_sf <- station_coords %>% filter(!is.na(longitude), !is.na(latitude)) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
stations_with_city <- stations_sf %>% st_join(china_cities_shp, join = st_within) %>% as.data.frame() %>% 
  dplyr::select(meteo_stat_id, city_name = ct_name)

# 最终分析底表
analysis_base <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id)) %>%
  inner_join(stations_with_city, by = "meteo_stat_id") %>%
  left_join(green_invest_2020, by = "city_name") %>%
  left_join(gdp_data, by = "city_name") %>%
  left_join(green_area_2020, by = "city_name") %>%
  mutate(
    invest_亿元 = invest / 10000,
    invest_ratio = invest / gdp_total * 100,
    intensity_green = invest / area_green,
    intensity_built = invest / area_built,
    intensity_park = invest / area_park,
    is_promote = ifelse(effect_type_heat_sif == "促进效应(+)", TRUE, FALSE)
  )

# ============================================================================
# 4. 阈值分析绘图逻辑
# ============================================================================

cat("\n【4. 执行全维度阈值分析】\n")

promote_data <- analysis_base %>% filter(is_promote, effect_index_heat_sif > 0)

plot_threshold_all <- function(data, x_var, x_label, is_log = FALSE) {
  d <- data %>% filter(!is.na(.data[[x_var]]), .data[[x_var]] > 0)
  
  if (is_log) {
    d$x_p <- log10(d[[x_var]])
    d$y_p <- log10(d$effect_index_heat_sif)
    lab_x <- paste0("log10(", x_label, ")")
    lab_y <- "log10(促进强度)"
  } else {
    d$x_p <- d[[x_var]]
    d$y_p <- d$effect_index_heat_sif
    lab_x <- x_label
    lab_y <- "促进效应强度"
  }
  
  p <- ggplot(d, aes(x = x_p, y = y_p)) +
    geom_point(alpha = 0.2, color = "gray40", size = 0.8) +
    labs(title = x_label, x = lab_x, y = lab_y) +
    theme_minimal(base_size = 10)
  
  try({
    fit_seg <- segmented(lm(y_p ~ x_p, data = d), seg.Z = ~x_p, npsi = 1)
    d$fitted <- predict(fit_seg)
    psi_val <- fit_seg$psi[1, "Est."]
    thresh <- if(is_log) 10^psi_val else psi_val
    
    p <- p + geom_line(data = d, aes(y = fitted), color = "black", linewidth = 1) +
      geom_vline(xintercept = psi_val, linetype = "dashed", color = "#d73027") +
      annotate("text", x = psi_val, y = max(d$y_p)*0.9, label = paste0("T=", round(thresh, 2)), 
               color = "#d73027", fontface = "bold", size = 3)
  }, silent = TRUE)
  return(p)
}

vars <- c("invest_亿元", "invest_ratio", "intensity_green", "intensity_built", "intensity_park")
labs <- c("绿化投资(亿元)", "投资占比(%)", "单位绿地投资", "单位建成区投资", "单位公园投资")

# 生成并保存 Log-Log 组合图
cat("生成 Log-Log 组合图...\n")
plots_log <- map2(vars, labs, ~plot_threshold_all(promote_data, .x, .y, is_log = TRUE))
p_log_final <- wrap_plots(plots_log, ncol = 3) + plot_annotation(title = "多维度投入指标阈值分析 (Log-Log 尺度)")
ggsave("data_proc/threshold_segmented_loglog_all.png", p_log_final, width = 15, height = 10, dpi = 300)

# 生成并保存 线性 组合图
cat("生成线性组合图...\n")
plots_lin <- map2(vars, labs, ~plot_threshold_all(promote_data, .x, .y, is_log = FALSE))
p_lin_final <- wrap_plots(plots_lin, ncol = 3) + plot_annotation(title = "多维度投入指标阈值分析 (线性尺度)")
ggsave("data_proc/threshold_segmented_linear_all.png", p_lin_final, width = 15, height = 10, dpi = 300)

cat("\n分析完成！组合图已保存至 data_proc 目录。\n")
