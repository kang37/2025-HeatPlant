# ============================================================================
# Script: plot_sif_by_region_timeseries.R
# Description: 按地理区域以时间序列折线图形式可视化SIF的变化 (2000-2023)
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, ggplot2, lubridate, stringr, readr,
  tidyr, sf, showtext, patchwork, scales
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

sif_regional_raw <- meteo_sif_data %>%
  filter(year >= 2000, year <= 2023) %>%
  left_join(stations_with_province, by = "meteo_stat_id") %>%
  left_join(province_to_region, by = "pr_name") %>%
  filter(!is.na(region))

# 按区域、年、月聚合，计算均值和标准差
sif_timeseries_agg <- sif_regional_raw %>%
  group_by(region, year, month) %>%
  summarise(
    sif_mean = mean(sif, na.rm = TRUE),
    sif_sd = sd(sif, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 创建一个日期列用于绘图
  mutate(
    date = as.Date(paste(year, month, "1", sep="-")),
    # 将NA值的sd替换为0，以避免ribbon出现问题
    sif_sd = ifelse(is.na(sif_sd), 0, sif_sd)
  )


# --- 4. 可视化 ---
cat("--> 正在生成图表...\n")

plot_sif_timeseries <- ggplot(sif_timeseries_agg, aes(x = date, y = sif_mean, group = 1)) +
  # 绘制阴影区域代表标准差范围
  geom_ribbon(aes(ymin = sif_mean - sif_sd, ymax = sif_mean + sif_sd, fill = region), alpha = 0.3) +
  # 绘制均值折线
  geom_line(aes(color = region), linewidth = 0.6) +
  # 按区域分面
  facet_wrap(~region, scales = "free_y", ncol = 1) + # 使用单列以更好地展示时间序列
  # 格式化X轴日期显示
  scale_x_date(
    date_breaks = "4 years",
    date_labels = "%Y"
  ) +
  # 添加标题和标签
  labs(
    title = "各区域SIF时间序列 (2000-2023)",
    subtitle = "实线为区域月平均SIF，阴影为区域内站点的标准差范围",
    x = "年份",
    y = "区域月平均 SIF"
  ) +
  theme_bw(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 100, hjust = 0.5),
    plot.subtitle = element_text(size = 100, hjust = 0.5, margin=margin(b=15)),
    strip.text = element_text(face = "bold", size = 100),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100),
    axis.title = element_text(size = 100),
    legend.position = "none"
  )

# --- 5. 保存图表 ---
ggsave("sif_timeseries_by_region.png", plot_sif_timeseries, width = 12, height = 20, dpi = 300, bg = "white")

cat("--> 分析完成！图表已保存为 sif_timeseries_by_region.png\n")
