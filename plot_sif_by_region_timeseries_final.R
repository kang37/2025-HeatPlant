# ============================================================================
# Script: plot_sif_by_region_timeseries_final.R
# Description: 重新绘制SIF时间序列图，修正中文显示并增大字体。
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, ggplot2, lubridate, stringr, readr,
  tidyr, sf, purrr, showtext, scales
)

# --- 1. 配置中文字体 ---
cat("--> 正在配置中文字体...\n")
font_add_google("Noto Sans SC", "notosans-sc")
showtext_auto()

# --- 2. 数据准备 ---
cat("--> 正在准备SIF数据...\n")

# 定义区域
province_to_region <- tibble(
  pr_name_raw = c(
    "辽宁省", "吉林省", "黑龙江省", "北京市", "天津市", "河北省", "山西省", "内蒙古自治区",
    "上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "河南省", "湖北省", "湖南省",
    "广东省", "广西壮族自治区", "海南省", "重庆市", "四川省", "贵州省", "云南省", "西藏自治区",
    "陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"
  ),
  region = c(
    rep("东北", 3), rep("华北", 5), rep("华东", 7), rep("华中", 3), rep("华南", 3), rep("西南", 5), rep("西北", 5)
  )
) %>%
  mutate(pr_name = str_remove(pr_name_raw, "省|市|自治区|回族|维吾尔|壮族"))

# 读取SIF数据
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  tibble() %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif), year >= 2000, year <= 2023) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 关联区域信息
station_coords <- meteo_sif_data %>% select(meteo_stat_id, longitude, latitude) %>% distinct(meteo_stat_id, .keep_all = TRUE)
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326) %>% mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))
stations_sf <- station_coords %>% filter(!is.na(longitude), !is.na(latitude)) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
stations_with_province <- stations_sf %>% st_join(china_cities_shp, join = st_nearest_feature) %>% as.data.frame() %>% select(meteo_stat_id, pr_name)

sif_regional_raw <- meteo_sif_data %>%
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
  mutate(
    date = as.Date(paste(year, month, "1", sep="-")),
    sif_sd = ifelse(is.na(sif_sd), 0, sif_sd)
  )

# --- 3. 可视化 (大字体、中文修正版) ---
cat("--> 正在生成最终的SIF图表...\n")

plot_sif_final <- ggplot(sif_timeseries_agg, aes(x = date, y = sif_mean, group = 1)) +
  geom_ribbon(aes(ymin = sif_mean - sif_sd, ymax = sif_mean + sif_sd, fill = region), alpha = 0.3) +
  geom_line(aes(color = region), linewidth = 0.8) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  labs(
    title = "各区域SIF时间序列 (2000-2023)",
    subtitle = "实线为区域月平均SIF，阴影为区域内站点的标准差范围",
    x = "年份",
    y = "区域月平均 SIF"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "notosans-sc"),
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

# --- 4. 保存图表 ---
ggsave("sif_timeseries_by_region_final.png", plot_sif_final, width = 14, height = 22, dpi = 300, bg = "white")

cat("--> 分析完成！最终图表已保存为 sif_timeseries_by_region_final.png\n")
