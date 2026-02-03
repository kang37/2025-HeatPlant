# ============================================================================
# Script: plot_sif_vs_vpd_scatter.R
# Description: 绘制SIF vs. 月平均VPD的散点图，按区域分面，按年份着色。
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, ggplot2, lubridate, stringr, readr,
  tidyr, sf, purrr, showtext, scales, viridis
)

# --- 1. 配置字体 ---
cat("--> 正在配置中文字体...\n")
font_add_google("Noto Sans SC", "notosans-sc")
showtext_auto()

# --- 2. 数据准备 ---
cat("--> 正在准备VPD指标数据... (此过程较慢)\n")

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

# 读取目标站点
target_stations <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% pull(meteo_stat) %>% unique() %>% as.character()

# 读取每日气象数据
read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1)
}
meteo_file_list <- list.files(path = "data_raw/meteo_data_1961-2023", pattern = "\\.txt$", full.names = TRUE) %>% .[!grepl("sta_lonlat_china.txt", .)]
meteo_data_daily <- map_dfr(meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations], read_one_meteo_daily) %>%
  rename_with(~ tolower(.x)) %>%
  mutate(date_d = as.Date(date), year = year(date_d), month = month(date_d)) %>%
  filter(year >= 2000, year <= 2023)

# 计算每日VPD
meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    avp = (rh / 100) * svp,
    vpd = svp - avp
  ) %>%
  mutate(vpd = if_else(vpd < 0 | vpd > 10, NA_real_, vpd)) %>%
  filter(!is.na(vpd))

# 计算月度VPD指标
VPD_THRESHOLD <- 2.0
monthly_vpd_indicators <- meteo_data_daily_vpd %>%
  mutate(
    vpd_excess = pmax(vpd - VPD_THRESHOLD, 0),
    is_heat_day = vpd > VPD_THRESHOLD
  ) %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(
    n_days = n(),
    vpd_mean = mean(vpd, na.rm = TRUE),
    vpd_max = max(vpd, na.rm = TRUE),
    heat_days = sum(is_heat_day, na.rm = TRUE),
    cumulative_heat_impact = sum(vpd_excess, na.rm = TRUE),
    .groups = "drop"
  )

# --- 3. 整合SIF、VPD和区域数据 ---
cat("--> 正在整合SIF、VPD和区域数据...\n")

# 读取SIF数据
meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  tibble() %>%
  rename_with(~tolower(.x)) %>%
  rename(meteo_stat_id = meteo_stat) %>%
  filter(!is.na(sif), year >= 2000, year <= 2023) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

# 获取区域信息
station_coords <- meteo_sif_data %>% select(meteo_stat_id, longitude, latitude) %>% distinct()
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326) %>% mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))
stations_sf <- station_coords %>% filter(!is.na(longitude), !is.na(latitude)) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
stations_with_province <- stations_sf %>% st_join(china_cities_shp, join = st_nearest_feature) %>% as.data.frame() %>% select(meteo_stat_id, pr_name)
stations_with_region <- stations_with_province %>% left_join(province_to_region, by = "pr_name") %>% select(meteo_stat_id, region)

# 合并成最终的数据框
final_data <- inner_join(
  meteo_sif_data,
  monthly_vpd_indicators,
  by = c("meteo_stat_id", "year", "month")
) %>%
  left_join(stations_with_region, by = "meteo_stat_id") %>%
  filter(!is.na(region))

# --- 4. 可视化 ---
cat("--> 正在生成图表...\n")

plot_scatter <- ggplot(final_data, aes(x = vpd_mean, y = sif)) +
  geom_point(alpha = 0.5, size = 1.5) +
  facet_grid(year ~ region, scales = "free") +
  # 使用viridis色板，对年份进行着色
  scale_color_viridis_c(option = "plasma", name = "年份") +
  labs(
    title = "各区域 SIF vs. 月平均VPD 关系图",
    subtitle = "每个点代表一个站点在一个月的数据，颜色代表年份",
    x = "月平均 VPD (kPa)",
    y = "SIF"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "notosans-sc"),
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# --- 5. 保存图表 ---
ggsave("sif_vs_vpd_scatter.png", plot_scatter, width = 18, height = 12, dpi = 300, bg = "white")

cat("--> 分析完成！图表已保存为 sif_vs_vpd_scatter.png\n")
