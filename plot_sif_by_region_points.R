# ============================================================================
# Script: plot_sif_by_region_points.R
# Description: 按地理区域以点图形式可视化SIF的季节性变化 (2000-2023)
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, ggplot2, lubridate, stringr, readr,
  tidyr, sf, showtext, patchwork
)
showtext_auto()

# --- 1. 定义省份到地理区域的映射 ---
province_to_region <- tibble(
  pr_name_raw = c(
    "辽宁省", "吉林省", "黑龙江省",
    "北京市", "天津市", "河北省", "山西省", "内蒙古自治区",
    "上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省",
    "河南省", "湖北省", "湖南省",
    "广东省", "广西壮族自治区", "海南省",
    "重庆市", "四川省", "贵州省", "云南省", "西藏自治区",
    "陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"
  ),
  region = c(
    rep("东北", 3), rep("华北", 5), rep("华东", 7),
    rep("华中", 3), rep("华南", 3), rep("西南", 5), rep("西北", 5)
  )
) %>%
  mutate(pr_name = str_remove(pr_name_raw, "省|市|自治区|回族|维吾尔|壮族"))

# --- 2. 读取源数据 ---
cat("--> 正在读取站点SIF数据和地理信息...\n")
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  tibble() %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

station_coords <- meteo_sif_data %>%
  select(meteo_stat_id, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE)

china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326) %>%
  mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))

# --- 3. 数据处理与空间关联 ---
cat("--> 正在进行空间匹配和数据聚合...\n")
stations_sf <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_province <- stations_sf %>%
  st_join(china_cities_shp, join = st_nearest_feature) %>%
  as.data.frame() %>%
  select(meteo_stat_id, pr_name)

# 筛选2000-2023年的数据，并关联区域信息
# 这是将要用于绘图的原始数据点
sif_regional_raw_points <- meteo_sif_data %>%
  filter(year >= 2000, year <= 2023) %>%
  left_join(stations_with_province, by = "meteo_stat_id") %>%
  left_join(province_to_region, by = "pr_name") %>%
  filter(!is.na(region))

# --- 4. 可视化 ---
cat("--> 正在生成图表...\n")

plot_sif_points <- ggplot(sif_regional_raw_points, aes(x = factor(month), y = sif)) +
  # 使用 geom_jitter 绘制点图，并设置低透明度和较小的点
  geom_jitter(aes(color = region), size = 0.5, alpha = 0.2, width = 0.3) +
  # 按区域分面
  facet_wrap(~region, scales = "free_y", ncol = 3) +
  # 添加标题和标签
  labs(
    title = "各区域原始SIF数据点季节性分布 (2000-2023)",
    subtitle = "每个点代表一个站点在一个月内的原始SIF值",
    x = "月份",
    y = "站点月 SIF 值"
  ) +
  scale_x_discrete(breaks = seq(1, 12, by = 1)) +
  # 设置一个简洁的主题
  theme_bw(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin=margin(b=15)),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

# --- 5. 保存图表 ---
ggsave("sif_seasonal_by_region_points.png", plot_sif_points, width = 15, height = 9, dpi = 300, bg = "white")

cat("--> 分析完成！图表已保存为 sif_seasonal_by_region_points.png\n")
