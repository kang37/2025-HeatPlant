# ============================================================================
# Script: 02_green_invest_analysis.R
# Description: 分析不同绿化投资水平下VPD热事件因果关系的差异
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, readr, tidyr, stringr, sf, readxl, showtext
)
showtext_auto()

# ============================================================================
# 1. 加载CCM分析结果
# ============================================================================

cat("【1. 加载数据】\n")

# 检查是否有已保存的CCM结果
ccm_results_path <- "data_proc/ccm_results_vpd_heat_sif.rds"
if (!file.exists(ccm_results_path)) {
  # 如果没有rds，尝试运行main_vpd_index2.R或使用已有csv
  ccm_csv_path <- "data_proc/ccm_results_vpd_heat_sif.csv"
  if (file.exists(ccm_csv_path)) {
    ccm_results_all <- read_csv(ccm_csv_path, show_col_types = FALSE)
  } else {
    stop("请先运行 main_vpd_index2.R 生成CCM分析结果")
  }
} else {
  ccm_results_all <- readRDS(ccm_results_path)
}

# 读取站点坐标
station_coords <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

cat("CCM结果站点数:", nrow(ccm_results_all), "\n")

# ============================================================================
# 2. 读取绿化投资数据
# ============================================================================

cat("\n【2. 读取绿化投资数据】\n")

green_invest <- read_csv(
  "data_raw/green_invest/中国城市数据.csv",
  show_col_types = FALSE
) %>%
  select(city_name = 城市名称, green_invest_2020 = 园林绿化_2020) %>%
  mutate(
    green_invest_2020 = as.numeric(green_invest_2020),
    # 添加"市"字以便与shapefile城市名匹配
    city_name = paste0(city_name, "市")
  )

cat("绿化投资数据城市数:", nrow(green_invest), "\n")
cat("有2020年数据的城市:", sum(!is.na(green_invest$green_invest_2020)), "\n")

# ============================================================================
# 3. 读取GDP数据（用于计算绿化投资占比）
# ============================================================================

cat("\n【3. 读取GDP数据】\n")

gdp_data <- read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
  filter(年份 == 2020) %>%
  select(
    city_name = 城市,
    gdp_total = "地区生产总值(万元)"
  ) %>%
  mutate(
    gdp_total = as.numeric(gdp_total)
  )

cat("GDP数据城市数:", nrow(gdp_data), "\n")

# 合并绿化投资和GDP
green_gdp <- green_invest %>%
  inner_join(gdp_data, by = "city_name") %>%
  mutate(
    # 绿化投资占GDP比例 (%)
    green_ratio = green_invest_2020 / gdp_total * 100
  ) %>%
  filter(!is.na(green_invest_2020), !is.na(gdp_total), gdp_total > 0)

cat("同时有绿化和GDP数据的城市:", nrow(green_gdp), "\n")

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
  select(-geometry) %>%
  rename(city_name = ct_name)

cat("匹配到城市的站点数:", sum(!is.na(stations_with_city$city_name)), "\n")

# ============================================================================
# 5. 合并所有数据
# ============================================================================

cat("\n【5. 合并数据】\n")

# 统一meteo_stat_id类型为character
ccm_results_all <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

stations_with_city <- stations_with_city %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 合并CCM结果与站点-城市映射
spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  left_join(stations_with_city %>% select(meteo_stat_id, city_name),
            by = "meteo_stat_id") %>%
  left_join(green_gdp, by = "city_name") %>%
  filter(!is.na(longitude), !is.na(latitude))

# 创建因果类型标签（与main_vpd_index2.R一致）
spatial_green <- spatial_data %>%
  mutate(
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      causality_direction == "SIF → VPD热影响" ~ "SIF→VPD热影响",
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "VPD热影响 → SIF" ~ "热事件促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "VPD热影响 → SIF" ~ "热事件抑制SIF",
      effect_type_heat_sif == "促进效应(+)" & causality_direction == "双向因果" ~ "双向-促进SIF",
      effect_type_heat_sif == "抑制效应(-)" & causality_direction == "双向因果" ~ "双向-抑制SIF",
      TRUE ~ "其他"
    )
  )

cat("有绿化投资数据的站点:", sum(!is.na(spatial_green$green_invest_2020)), "\n")
cat("有绿化/GDP比例数据的站点:", sum(!is.na(spatial_green$green_ratio)), "\n")

# ============================================================================
# 6. 绿化投资分析
# ============================================================================

cat("\n【6. 绿化投资分析】\n")

# 筛选有效数据
analysis_data <- spatial_green %>%
  filter(
    !is.na(green_invest_2020),
    heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                          "双向-促进SIF", "双向-抑制SIF", "无因果")
  )

cat("用于分析的站点数:", nrow(analysis_data), "\n\n")

# 统计表
cat("因果类型的绿化投资统计:\n")
green_stats <- analysis_data %>%
  group_by(heat_causality) %>%
  summarise(
    站点数 = n(),
    平均绿化投资_万元 = round(mean(green_invest_2020, na.rm = TRUE), 0),
    中位绿化投资_万元 = round(median(green_invest_2020, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  arrange(desc(站点数))
print(green_stats)

# 定义颜色
causality_colors <- c(
  "热事件促进SIF" = "#4575b4",
  "热事件抑制SIF" = "#d73027",
  "双向-促进SIF" = "#91bfdb",
  "双向-抑制SIF" = "#fc8d59",
  "无因果" = "#999999"
)

# ============================================================================
# 7. 绘图1: 绿化投资 vs 因果类型
# ============================================================================

cat("\n【7. 绘制箱线图】\n")

p1 <- analysis_data %>%
  mutate(green_invest_万元 = green_invest_2020 / 10000) %>%  # 转换为亿元便于显示
  ggplot(aes(
    x = reorder(heat_causality, -green_invest_万元, FUN = median),
    y = green_invest_万元,
    fill = heat_causality
  )) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, color = "gray30") +
  scale_fill_manual(values = causality_colors, guide = "none") +
  labs(
    title = "VPD热事件因果类型与城市绿化投资的关系",
    subtitle = "数据年份: 2020年",
    x = "因果类型",
    y = "园林绿化投资 (亿元)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    panel.grid.minor = element_blank()
  )

print(p1)
ggsave("data_proc/boxplot_green_invest_causality.png", p1,
       width = 10, height = 7, dpi = 300)
cat("图1已保存: data_proc/boxplot_green_invest_causality.png\n")

# ============================================================================
# 8. 绘图2: 绿化投资占GDP比例 vs 因果类型
# ============================================================================

# 筛选有比例数据的站点
analysis_data_ratio <- spatial_green %>%
  filter(
    !is.na(green_ratio),
    green_ratio > 0, green_ratio < 10,  # 过滤异常值
    heat_causality %in% c("热事件促进SIF", "热事件抑制SIF",
                          "双向-促进SIF", "双向-抑制SIF", "无因果")
  )

cat("\n绿化投资占GDP比例分析 (站点数:", nrow(analysis_data_ratio), ")\n")

# 统计表
cat("因果类型的绿化投资占比统计:\n")
ratio_stats <- analysis_data_ratio %>%
  group_by(heat_causality) %>%
  summarise(
    站点数 = n(),
    平均占比_百分比 = round(mean(green_ratio, na.rm = TRUE), 4),
    中位占比_百分比 = round(median(green_ratio, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(站点数))
print(ratio_stats)

p2 <- analysis_data_ratio %>%
  ggplot(aes(
    x = reorder(heat_causality, -green_ratio, FUN = median),
    y = green_ratio,
    fill = heat_causality
  )) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, color = "gray30") +
  scale_fill_manual(values = causality_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "VPD热事件因果类型与绿化投资占GDP比例的关系",
    subtitle = "数据年份: 2020年",
    x = "因果类型",
    y = "绿化投资占GDP比例 (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    panel.grid.minor = element_blank()
  )

print(p2)
ggsave("data_proc/boxplot_green_ratio_causality.png", p2,
       width = 10, height = 7, dpi = 300)
cat("图2已保存: data_proc/boxplot_green_ratio_causality.png\n")

# ============================================================================
# 9. 统计检验
# ============================================================================

cat("\n【8. 统计检验】\n")

# 绿化投资的Kruskal-Wallis检验
if (nrow(analysis_data) > 10 && length(unique(analysis_data$heat_causality)) >= 2) {
  kw_green <- kruskal.test(green_invest_2020 ~ heat_causality, data = analysis_data)
  cat("\n绿化投资 Kruskal-Wallis检验:\n")
  cat("  χ² =", round(kw_green$statistic, 3), ", p =", format.pval(kw_green$p.value, digits = 3), "\n")

  # 两两比较
  if (length(unique(analysis_data$heat_causality)) >= 2) {
    pw_green <- pairwise.wilcox.test(
      analysis_data$green_invest_2020,
      analysis_data$heat_causality,
      p.adjust.method = "BH"
    )
    cat("  两两比较 (Wilcoxon, BH校正):\n")
    print(pw_green$p.value)
  }
}

# 绿化占比的Kruskal-Wallis检验
if (nrow(analysis_data_ratio) > 10 && length(unique(analysis_data_ratio$heat_causality)) >= 2) {
  kw_ratio <- kruskal.test(green_ratio ~ heat_causality, data = analysis_data_ratio)
  cat("\n绿化占GDP比例 Kruskal-Wallis检验:\n")
  cat("  χ² =", round(kw_ratio$statistic, 3), ", p =", format.pval(kw_ratio$p.value, digits = 3), "\n")

  if (length(unique(analysis_data_ratio$heat_causality)) >= 2) {
    pw_ratio <- pairwise.wilcox.test(
      analysis_data_ratio$green_ratio,
      analysis_data_ratio$heat_causality,
      p.adjust.method = "BH"
    )
    cat("  两两比较 (Wilcoxon, BH校正):\n")
    print(pw_ratio$p.value)
  }
}

cat("\n分析完成!\n")
