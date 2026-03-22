# ============================================================================
# Script: 04_koppen_climate_analysis.R
# Description: 分析不同柯本气候区与VPD热事件因果关系的联系
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, readr, tidyr, stringr, sf, terra, showtext
)
showtext_auto()

# ============================================================================
# 1. 定义柯本气候区编码表
# ============================================================================

koppen_legend <- tibble(
  code = 1:30,
  climate_code = c(
    "Af", "Am", "Aw",           # 1-3: 热带
    "BWh", "BWk", "BSh", "BSk", # 4-7: 干旱
    "Csa", "Csb", "Csc",        # 8-10: 温带干夏
    "Cwa", "Cwb", "Cwc",        # 11-13: 温带干冬
    "Cfa", "Cfb", "Cfc",        # 14-16: 温带无干季
    "Dsa", "Dsb", "Dsc", "Dsd", # 17-20: 冷干夏
    "Dwa", "Dwb", "Dwc", "Dwd", # 21-24: 冷干冬
    "Dfa", "Dfb", "Dfc", "Dfd", # 25-28: 冷无干季
    "ET", "EF"                   # 29-30: 极地
  ),
  climate_name = c(
    "热带雨林", "热带季风", "热带草原",
    "热沙漠", "冷沙漠", "热草原", "冷草原",
    "地中海热夏", "地中海暖夏", "地中海冷夏",
    "亚热带季风", "温带季风暖夏", "温带季风冷夏",
    "亚热带湿润", "温带海洋", "亚极地海洋",
    "冷干夏热夏", "冷干夏暖夏", "冷干夏冷夏", "冷干夏极冷冬",
    "冷干冬热夏", "冷干冬暖夏", "冷干冬冷夏", "冷干冬极冷冬",
    "湿润大陆热夏", "湿润大陆暖夏", "亚寒带", "极端大陆",
    "苔原", "冰盖"
  ),
  climate_group = c(
    rep("A 热带", 3),
    rep("B 干旱", 4),
    rep("C 温带", 9),
    rep("D 大陆", 12),
    rep("E 极地", 2)
  )
)

# ============================================================================
# 2. 加载CCM分析结果和站点坐标
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

cat("CCM结果站点数:", nrow(ccm_results_all), "\n")

# ============================================================================
# 3. 读取柯本气候区栅格数据
# ============================================================================

cat("\n【2. 读取柯本气候区数据】\n")

# 使用1991-2020时段的0.1度分辨率数据
koppen_raster <- rast("data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")
cat("栅格分辨率:", res(koppen_raster), "\n")
cat("栅格范围:", ext(koppen_raster)[1:4], "\n")

# ============================================================================
# 4. 提取每个站点的气候区
# ============================================================================

cat("\n【3. 提取站点气候区】\n")

# 创建站点的空间点
stations_vect <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  vect(geom = c("longitude", "latitude"), crs = "EPSG:4326")

# 提取气候区编码
climate_values <- extract(koppen_raster, stations_vect)

stations_climate <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    climate_code_num = climate_values[, 2]
  ) %>%
  left_join(koppen_legend, by = c("climate_code_num" = "code"))

cat("成功提取气候区的站点数:", sum(!is.na(stations_climate$climate_code)), "\n")

# 统计各气候区站点数
cat("\n各气候区站点分布:\n")
climate_dist <- stations_climate %>%
  filter(!is.na(climate_code)) %>%
  count(climate_group, climate_code, climate_name) %>%
  arrange(climate_group, climate_code)
print(as.data.frame(climate_dist))

# ============================================================================
# 5. 合并CCM结果与气候区数据
# ============================================================================

cat("\n【4. 合并数据】\n")

ccm_results_all <- ccm_results_all %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 创建因果类型标签
analysis_data <- ccm_results_all %>%
  left_join(stations_climate, by = "meteo_stat_id") %>%
  filter(!is.na(climate_code)) %>%
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
    effect_strength = abs(effect_index_heat_sif)
  ) %>%
  filter(heat_causality != "其他")

cat("用于分析的站点数:", nrow(analysis_data), "\n")

# ============================================================================
# 6. 统计分析
# ============================================================================

cat("\n【5. 统计分析】\n")

# 按气候大类统计因果类型分布
cat("\n各气候大类的因果类型分布:\n")
climate_causality_stats <- analysis_data %>%
  count(climate_group, heat_causality) %>%
  group_by(climate_group) %>%
  mutate(
    total = sum(n),
    pct = round(n / total * 100, 1)
  ) %>%
  ungroup()

print(as.data.frame(climate_causality_stats))

# ============================================================================
# 7. 可视化
# ============================================================================

cat("\n【6. 绑制图表】\n")

# 定义颜色
causality_colors <- c(
  "热事件促进SIF" = "#4575b4",
  "热事件抑制SIF" = "#d73027",
  "双向-促进SIF" = "#91bfdb",
  "双向-抑制SIF" = "#fc8d59",
  "SIF->VPD热影响" = "#984ea3",
  "无因果" = "#999999"
)

# --- 图1: 气候大类的因果类型堆积柱状图 ---
p1 <- analysis_data %>%
  count(climate_group, heat_causality) %>%
  ggplot(aes(x = climate_group, y = n, fill = heat_causality)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = causality_colors, name = "因果类型") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "不同柯本气候大类的VPD热事件因果类型分布",
    subtitle = "柱高表示各因果类型所占比例",
    x = "气候大类",
    y = "比例 (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 15, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p1)
ggsave("data_proc/koppen_causality_stacked_bar.png", p1,
       width = 10, height = 7, dpi = 300)
cat("图1已保存: data_proc/koppen_causality_stacked_bar.png\n")

# --- 图2: 详细气候类型的因果分布（热图） ---
# 只保留样本量>=5的气候类型
climate_counts <- analysis_data %>%
  count(climate_code, climate_name) %>%
  filter(n >= 5)

heatmap_data <- analysis_data %>%
  filter(climate_code %in% climate_counts$climate_code) %>%
  count(climate_code, climate_name, climate_group, heat_causality) %>%
  group_by(climate_code) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p2 <- heatmap_data %>%
  mutate(
    climate_label = paste0(climate_code, " ", climate_name),
    climate_label = fct_reorder(climate_label, as.numeric(factor(climate_group)))
  ) %>%
  ggplot(aes(x = heat_causality, y = climate_label, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct, 0)), size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "比例 (%)") +
  labs(
    title = "各柯本气候类型的因果关系分布热图",
    subtitle = "数字表示该因果类型在该气候区的占比 (%)",
    x = "因果类型",
    y = "气候类型"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p2)
ggsave("data_proc/koppen_causality_heatmap.png", p2,
       width = 12, height = 10, dpi = 300)
cat("图2已保存: data_proc/koppen_causality_heatmap.png\n")

# --- 图3: 气候大类的因果类型计数柱状图（分面） ---
p3 <- analysis_data %>%
  count(climate_group, heat_causality) %>%
  ggplot(aes(x = heat_causality, y = n, fill = heat_causality)) +
  geom_col(width = 0.7) +
  facet_wrap(~climate_group, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = causality_colors, guide = "none") +
  labs(
    title = "各柯本气候大类的因果类型站点数",
    x = "因果类型",
    y = "站点数"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

print(p3)
ggsave("data_proc/koppen_causality_facet_bar.png", p3,
       width = 12, height = 8, dpi = 300)
cat("图3已保存: data_proc/koppen_causality_facet_bar.png\n")

# ============================================================================
# 8. 卡方检验
# ============================================================================

cat("\n【7. 统计检验】\n")

# 气候大类与因果类型的独立性检验
contingency_table <- analysis_data %>%
  filter(heat_causality %in% c("热事件促进SIF", "热事件抑制SIF", "无因果")) %>%
  count(climate_group, heat_causality) %>%
  pivot_wider(names_from = heat_causality, values_from = n, values_fill = 0)

cat("\n列联表:\n")
print(contingency_table)

# 卡方检验
chi_matrix <- contingency_table %>%
  select(-climate_group) %>%
  as.matrix()
rownames(chi_matrix) <- contingency_table$climate_group

chi_test <- chisq.test(chi_matrix)
cat("\n卡方检验结果:\n")
cat("  χ² =", round(chi_test$statistic, 2), "\n")
cat("  df =", chi_test$parameter, "\n")
cat("  p =", format.pval(chi_test$p.value, digits = 4), "\n")

# ============================================================================
# 9. 中国柯本气候区地图
# ============================================================================

cat("\n【8. 绑制中国柯本气候区地图】\n")

# 加载中国边界
library(rnaturalearth)
library(rnaturalearthdata)

china_border <- ne_countries(country = "china", scale = "medium", returnclass = "sf")

# 裁剪栅格到中国范围（稍微扩大范围确保完整）
china_extent <- ext(73, 136, 18, 54)
koppen_china <- crop(koppen_raster, china_extent)

# 转换为数据框用于ggplot
koppen_df <- as.data.frame(koppen_china, xy = TRUE) %>%
  rename(climate_code_num = 3) %>%
  filter(!is.na(climate_code_num)) %>%
  left_join(koppen_legend, by = c("climate_code_num" = "code"))

cat("中国范围内的栅格点数:", nrow(koppen_df), "\n")

# 统计中国各气候区面积占比
cat("\n中国各气候区分布:\n")
china_climate_stats <- koppen_df %>%
  count(climate_group, climate_code, climate_name) %>%
  mutate(pct = round(n / sum(n) * 100, 2)) %>%
  arrange(desc(n))
print(as.data.frame(china_climate_stats))

# 定义柯本气候区颜色（基于官方配色）
koppen_colors <- c(
  "Af" = "#0000FF", "Am" = "#0078FF", "Aw" = "#46AAFA",

"BWh" = "#FF0000", "BWk" = "#FF9696", "BSh" = "#F5A500", "BSk" = "#FFDC64",
  "Csa" = "#FFFF00", "Csb" = "#C8C800", "Csc" = "#969600",
  "Cwa" = "#96FF96", "Cwb" = "#64C864", "Cwc" = "#329632",
  "Cfa" = "#C8FF50", "Cfb" = "#64FF50", "Cfc" = "#32C800",
  "Dsa" = "#FF00FF", "Dsb" = "#C800C8", "Dsc" = "#963296", "Dsd" = "#966496",
  "Dwa" = "#AAAFFF", "Dwb" = "#5A78DC", "Dwc" = "#4B50B4", "Dwd" = "#320087",
  "Dfa" = "#00FFFF", "Dfb" = "#37C8FF", "Dfc" = "#007D7D", "Dfd" = "#00465F",
  "ET" = "#B2B2B2", "EF" = "#666666"
)

# 绑制地图
p_map <- ggplot() +
  geom_raster(data = koppen_df, aes(x = x, y = y, fill = climate_code)) +
  geom_sf(data = china_border, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = koppen_colors,
    name = "气候类型",
    guide = guide_legend(ncol = 2)
  ) +
  coord_sf(xlim = c(73, 136), ylim = c(18, 54), expand = FALSE) +
  labs(
    title = "中国柯本气候区分布图",
    subtitle = "数据来源: Beck et al. (2023), 1991-2020时段",
    x = "经度", y = "纬度"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "gray90", linewidth = 0.3)
  )

print(p_map)
ggsave("data_proc/china_koppen_climate_map.png", p_map,
       width = 12, height = 10, dpi = 300)
cat("地图已保存: data_proc/china_koppen_climate_map.png\n")

# 绑制简化版（按气候大类）
p_map_group <- ggplot() +
  geom_raster(data = koppen_df, aes(x = x, y = y, fill = climate_group)) +
  geom_sf(data = china_border, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "A 热带" = "#0078FF",
      "B 干旱" = "#FF6400",
      "C 温带" = "#96FF50",
      "D 大陆" = "#00C8C8",
      "E 极地" = "#B2B2B2"
    ),
    name = "气候大类"
  ) +
  coord_sf(xlim = c(73, 136), ylim = c(18, 54), expand = FALSE) +
  labs(
    title = "中国柯本气候大类分布图",
    subtitle = "A=热带, B=干旱, C=温带, D=大陆, E=极地",
    x = "经度", y = "纬度"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    legend.position = "right",
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "gray90", linewidth = 0.3)
  )

print(p_map_group)
ggsave("data_proc/china_koppen_climate_group_map.png", p_map_group,
       width = 12, height = 10, dpi = 300)
cat("简化地图已保存: data_proc/china_koppen_climate_group_map.png\n")

cat("\n分析完成!\n")
