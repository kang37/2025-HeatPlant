# ============================================================================
# Script: process_vpd_data.R
# Description: 准备VPD时间序列数据并保存为RDS文件，以供快速绘图。
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr, lubridate, stringr, readr,
  tidyr, sf, purrr
)

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
cat("--> 正在读取站点、地理及每日气象数据... (此过程较慢，请稍候)\n")

target_stations <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  pull(meteo_stat) %>%
  unique() %>%
  as.character()

station_coords <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326) %>%
  mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))

read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1)
}

meteo_file_list <- list.files(
  path = "data_raw/meteo_data_1961-2023",
  pattern = "\\.txt$",
  full.names = TRUE
) %>%
  .[!grepl("sta_lonlat_china.txt", .)]

meteo_data_daily <- map_dfr(
  meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations],
  read_one_meteo_daily
) %>%
  rename_with(~ tolower(.x)) %>%
  mutate(
    date_d = as.Date(date),
    year = year(date_d),
    month = month(date_d)
  ) %>%
  filter(year >= 2000, year <= 2023)

# --- 3. 数据处理与计算VPD ---
cat("--> 正在计算VPD并进行空间聚合...\n")

meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    avp = (rh / 100) * svp,
    vpd = svp - avp
  ) %>%
  mutate(vpd = if_else(vpd < 0 | vpd > 10, NA_real_, vpd)) %>%
  filter(!is.na(vpd))

stations_sf <- station_coords %>%
  filter(meteo_stat_id %in% target_stations, !is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_province <- stations_sf %>%
  st_join(china_cities_shp, join = st_nearest_feature) %>%
  as.data.frame() %>%
  select(meteo_stat_id, pr_name)

vpd_regional_raw <- meteo_data_daily_vpd %>%
  left_join(stations_with_province, by = "meteo_stat_id") %>%
  left_join(province_to_region, by = "pr_name") %>%
  filter(!is.na(region))

vpd_timeseries_agg <- vpd_regional_raw %>%
  group_by(region, year, month) %>%
  summarise(
    vpd_mean = mean(vpd, na.rm = TRUE),
    vpd_sd = sd(vpd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    date = as.Date(paste(year, month, "1", sep="-")),
    vpd_sd = ifelse(is.na(vpd_sd), 0, vpd_sd)
  )

# --- 4. 保存处理好的数据 ---
cat("--> 正在保存聚合数据到 vpd_timeseries_agg.rds...\n")
saveRDS(vpd_timeseries_agg, file = "vpd_timeseries_agg.rds")

cat("--> 数据处理和保存完成！\n")
