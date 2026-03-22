# ============================================================================
# Script: 03_green_ratio_spatial_scatter.R
# Description: 绿化投资占GDP比例与经纬度的散点图，展示因果类型和强度
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, readr, tidyr, stringr, sf, readxl, showtext, patchwork
)
showtext_auto()

# ============================================================================
# 1. 加载数据（复用02脚本的数据准备逻辑）
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
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 绿化投资
green_invest <- read_csv(
  "data_raw/green_invest/中国城市数据.csv",
  show_col_types = FALSE
) %>%
  select(city_name = 城市名称, green_invest_2020 = 园林绿化_2020) %>%
  mutate(
    green_invest_2020 = as.numeric(green_invest_2020),
    city_name = paste0(city_name, "市")
  )

# GDP数据
gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  select(city_name = 城市, gdp_total = "地区生产总值(万元)") %>%
  mutate(
    gdp_total = as.numeric(gdp_total)
  )

# 合并绿化和GDP
green_gdp <- green_invest %>%
  inner_join(gdp_data, by = "city_name") %>%
  mutate(green_ratio = green_invest_2020 / gdp_total * 100) %>%
  filter(!is.na(green_invest_2020), !is.na(gdp_total), gdp_total > 0)

# 空间匹配
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326)

stations_sf <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_city <- stations_sf %>%
  st_join(china_cities_shp, join = st_within) %>%
  as.data.frame() %>%
  select(-geometry) %>%
  rename(city_name = ct_name) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并所有数据
ccm_results_all <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  left_join(stations_with_city %>% select(meteo_stat_id, city_name),
            by = "meteo_stat_id") %>%
  left_join(green_gdp, by = "city_name") %>%
  filter(!is.na(longitude), !is.na(latitude))

# 创建因果类型标签和效应强度
plot_data <- spatial_data %>%
  mutate(
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      causality_direction == "SIF → VPD热影响" ~ "SIF->VPD热影响",
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "VPD热影响 → SIF" ~ "热事件促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "VPD热影响 → SIF" ~ "热事件抑制SIF",
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "双向因果" ~ "双向-促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "双向因果" ~ "双向-抑制SIF",
      TRUE ~ "其他"
    ),
    # 因果强度：使用效应指数的绝对值
    effect_strength = abs(effect_index_heat_sif)
  ) %>%
  filter(
    !is.na(green_ratio),
    green_ratio > 0, green_ratio < 10,  # 过滤异常值
    heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                          "双向-促进SIF", "双向-抑制SIF", "无因果"),
    !is.na(effect_strength)
  )

cat("用于绑图的站点数:", nrow(plot_data), "\n")

# ============================================================================
# 2. 定义颜色
# ============================================================================

causality_colors <- c(
  "热事件促进SIF" = "#4575b4",
  "热事件抑制SIF" = "#d73027",
  "双向-促进SIF" = "#91bfdb",
  "双向-抑制SIF" = "#fc8d59",
  "无因果" = "#999999"
)

# ============================================================================
# 3. 绘图1: 绿化占比 vs 纬度
# ============================================================================

cat("\n【2. 绑制散点图】\n")

p_lat <- plot_data %>%
  ggplot(aes(
    x = green_ratio,
    y = latitude,
    color = heat_causality,
    size = effect_strength
  )) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = causality_colors, name = "因果类型") +
  scale_size_continuous(range = c(1, 6), name = "因果强度\n|效应指数|") +
  labs(
    title = "绿化投资占GDP比例与纬度的关系",
    subtitle = "点大小表示因果效应强度",
    x = "绿化投资占GDP比例 (%)",
    y = "纬度 (°N)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_lat)
ggsave("data_proc/scatter_green_ratio_latitude.png", p_lat,
       width = 10, height = 7, dpi = 300)
cat("图1已保存: data_proc/scatter_green_ratio_latitude.png\n")

# ============================================================================
# 4. 绘图2: 绿化占比 vs 经度
# ============================================================================

p_lon <- plot_data %>%
  ggplot(aes(
    x = green_ratio,
    y = longitude,
    color = heat_causality,
    size = effect_strength
  )) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = causality_colors, name = "因果类型") +
  scale_size_continuous(range = c(1, 6), name = "因果强度\n|效应指数|") +
  labs(
    title = "绿化投资占GDP比例与经度的关系",
    subtitle = "点大小表示因果效应强度",
    x = "绿化投资占GDP比例 (%)",
    y = "经度 (°E)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_lon)
ggsave("data_proc/scatter_green_ratio_longitude.png", p_lon,
       width = 10, height = 7, dpi = 300)
cat("图2已保存: data_proc/scatter_green_ratio_longitude.png\n")

# ============================================================================
# 5. 组合图
# ============================================================================

p_combined <- p_lat + p_lon +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("data_proc/scatter_green_ratio_spatial_combined.png", p_combined,
       width = 14, height = 7, dpi = 300)
cat("组合图已保存: data_proc/scatter_green_ratio_spatial_combined.png\n")

cat("\n分析完成!\n")
