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
  mutate(gdp_total = as.numeric(gdp_total))

cat("GDP数据城市数:", nrow(gdp_data), "\n")

# 合并绿化投资和GDP
green_gdp <- green_invest %>%
  inner_join(gdp_data, by = "city_name") %>%
  mutate(green_ratio = green_invest_2020 / gdp_total * 100) %>%
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

ccm_results_all <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

stations_with_city <- stations_with_city %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

spatial_data <- ccm_results_all %>%
  left_join(station_coords, by = "meteo_stat_id") %>%
  left_join(stations_with_city %>% select(meteo_stat_id, city_name),
            by = "meteo_stat_id") %>%
  left_join(green_gdp, by = "city_name") %>%
  filter(!is.na(longitude), !is.na(latitude))

# 简化因果类型：只分三类（促进、抑制、无因果）
spatial_green <- spatial_data %>%
  mutate(
    heat_causality = case_when(
      causality_direction == "无显著因果" ~ "无因果",
      # 促进效应（不区分单向/双向）
      effect_type_heat_sif == "促进效应(+)" ~ "促进",
      # 抑制效应（不区分单向/双向）
      effect_type_heat_sif == "抑制效应(-)" ~ "抑制",
      TRUE ~ "其他"
    )
  )

cat("有绿化投资数据的站点:", sum(!is.na(spatial_green$green_invest_2020)), "\n")
cat("有绿化/GDP比例数据的站点:", sum(!is.na(spatial_green$green_ratio)), "\n")

# ============================================================================
# 6. 绿化投资分析
# ============================================================================

cat("\n【6. 绿化投资分析】\n")

# 筛选有效数据（三类）
analysis_data <- spatial_green %>%
  filter(
    !is.na(green_invest_2020),
    heat_causality %in% c("促进", "抑制", "无因果")
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

# 定义颜色（三类）
causality_colors <- c(
  "促进" = "#4575b4",
  "抑制" = "#d73027",
  "无因果" = "#999999"
)

# 定义分段线性转换函数 (用于坐标轴压缩)
# 参数: x 数据, breakpoint 转折点, factor 压缩倍率
pw_trans <- function(x, bp, factor = 5) {
  ifelse(x <= bp, x, bp + (x - bp) / factor)
}

# 逆转换 (用于恢复标签)
pw_inv <- function(x, bp, factor = 5) {
  ifelse(x <= bp, x, bp + (x - bp) * factor)
}

# ============================================================================
# 7. 绘图1: 绿化投资 vs 因果类型 (坐标轴压缩)
# ============================================================================

cat("\n【7. 绘制箱线图 (投资额 - 20亿以上压缩)】\n")

bp_inv <- 20
p1 <- analysis_data %>%
  mutate(green_invest_亿元 = green_invest_2020 / 10000) %>%
  # 手动转换数据用于绘图
  mutate(y_trans = pw_trans(green_invest_亿元, bp_inv, 5)) %>%
  ggplot(aes(
    x = reorder(heat_causality, -green_invest_亿元, FUN = median),
    y = y_trans,
    fill = heat_causality
  )) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2, color = "gray30") +
  # 设置坐标轴标签
  scale_y_continuous(
    breaks = pw_trans(c(0, 5, 10, 15, 20, 40, 60, 80, 100), bp_inv, 5),
    labels = c(0, 5, 10, 15, 20, 40, 60, 80, 100)
  ) +
  geom_hline(yintercept = pw_trans(bp_inv, bp_inv), linetype = "dashed", color = "blue", alpha = 0.5) +
  scale_fill_manual(values = causality_colors, guide = "none") +
  labs(
    title = "因果类型与城市绿化投资的关系 (压缩尺度)",
    subtitle = paste0("虚线(20亿)以上部分按 1/5 比例压缩显示"),
    x = "因果类型",
    y = "园林绿化投资 (亿元)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(p1)
ggsave("data_proc/boxplot_green_invest_causality.png", p1, width = 11, height = 9, dpi = 300)

# ============================================================================
# 8. 绘图2: 绿化投资占GDP比例 vs 因果类型 (坐标轴压缩)
# ============================================================================

analysis_data_ratio <- spatial_green %>%
  filter(
    !is.na(green_ratio),
    green_ratio > 0, green_ratio < 10,
    heat_causality %in% c("促进", "抑制", "无因果")
  )

cat("\n【8. 绘制箱线图 (占比 - 0.5%以上压缩)】\n")

bp_ratio <- 0.5
p2 <- analysis_data_ratio %>%
  mutate(y_trans = pw_trans(green_ratio, bp_ratio, 5)) %>%
  ggplot(aes(
    x = reorder(heat_causality, -green_ratio, FUN = median),
    y = y_trans,
    fill = heat_causality
  )) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2, color = "gray30") +
  scale_y_continuous(
    breaks = pw_trans(c(0, 0.25, 0.5, 1, 1.5, 2, 3), bp_ratio, 5),
    labels = scales::percent_format(scale = 1)(c(0, 0.25, 0.5, 1, 1.5, 2, 3))
  ) +
  geom_hline(yintercept = pw_trans(bp_ratio, bp_ratio), linetype = "dashed", color = "blue", alpha = 0.5) +
  scale_fill_manual(values = causality_colors, guide = "none") +
  labs(
    title = "因果类型与绿化投资占比的关系 (压缩尺度)",
    subtitle = paste0("虚线(0.5%)以上部分按 1/5 比例压缩显示"),
    x = "因果类型",
    y = "绿化投资占GDP比例"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

print(p2)
ggsave("data_proc/boxplot_green_ratio_causality.png", p2, width = 11, height = 9, dpi = 300)

# ============================================================================
# 9. 统计检验
# ============================================================================

cat("\n【8. 统计检验】\n")

if (nrow(analysis_data) > 10 && length(unique(analysis_data$heat_causality)) >= 2) {
  kw_green <- kruskal.test(green_invest_2020 ~ heat_causality, data = analysis_data)
  cat("\n绿化投资 Kruskal-Wallis检验:\n")
  cat("  χ² =", round(kw_green$statistic, 3), ", p =", format.pval(kw_green$p.value, digits = 3), "\n")

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
