# ============================================================================
# 各区域累积热影响指数时间序列可视化
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr,
  tidyr, showtext, sf
)
showtext_auto()

cat("\n", rep("=", 70), "\n", sep = "")
cat("        各区域累积热影响指数时间序列\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================================
# 1. 定义省份到地理区域的映射
# ============================================================================

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

# ============================================================================
# 2. 读取数据
# ============================================================================

cat("【1. 读取数据】\n")

# 读取站点坐标
station_coords <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>%
  rename_with(~tolower(.x)) %>%
  select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
  distinct(meteo_stat_id, .keep_all = TRUE) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

target_stations <- unique(station_coords$meteo_stat_id)
cat("目标站点数:", length(target_stations), "\n")

# 读取城市shapefile用于省份匹配
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>%
  st_transform(crs = 4326) %>%
  mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))

# 站点空间匹配到省份
stations_sf <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

stations_with_province <- stations_sf %>%
  st_join(china_cities_shp, join = st_nearest_feature) %>%
  as.data.frame() %>%
  select(meteo_stat_id, pr_name) %>%
  left_join(province_to_region, by = "pr_name") %>%
  filter(!is.na(region))

cat("有区域信息的站点:", nrow(stations_with_province), "\n")

# 读取每日气象数据
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

cat("读取每日气象数据...\n")
meteo_data_daily <- map(
  meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations],
  read_one_meteo_daily
) %>%
  list_rbind() %>%
  rename_with(~ tolower(.x)) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date)
  ) %>%
  filter(year >= 2000, year <= 2023, month %in% 5:9)

cat("数据行数:", nrow(meteo_data_daily), "\n\n")

# ============================================================================
# 3. 计算VPD和累积热影响指数
# ============================================================================

cat("【2. 计算VPD和累积热影响指数】\n")

VPD_THRESHOLD <- 2.0  # kPa

meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    avp = (rh / 100) * svp,
    vpd = svp - avp
  ) %>%
  mutate(
    vpd = case_when(
      vpd < 0 ~ NA_real_,
      vpd > 10 ~ NA_real_,
      is.na(tavg) | is.na(rh) ~ NA_real_,
      TRUE ~ vpd
    ),
    vpd_excess = pmax(vpd - VPD_THRESHOLD, 0, na.rm = TRUE)
  ) %>%
  filter(!is.na(vpd))

# 月度聚合
monthly_heat <- meteo_data_daily_vpd %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(
    n_days = n(),
    vpd_mean = mean(vpd, na.rm = TRUE),
    cumulative_heat_impact = sum(vpd_excess, na.rm = TRUE),
    heat_days = sum(vpd > VPD_THRESHOLD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(date = as.Date(paste(year, month, "15", sep = "-")))

cat("月度数据行数:", nrow(monthly_heat), "\n\n")

# ============================================================================
# 4. 按区域聚合
# ============================================================================

cat("【3. 按区域聚合】\n")

monthly_heat_regional <- monthly_heat %>%
  inner_join(stations_with_province, by = "meteo_stat_id") %>%
  group_by(region, year, month, date) %>%
  summarise(
    n_stations = n(),
    heat_impact_mean = mean(cumulative_heat_impact, na.rm = TRUE),
    heat_impact_sd = sd(cumulative_heat_impact, na.rm = TRUE),
    heat_impact_median = median(cumulative_heat_impact, na.rm = TRUE),
    vpd_mean = mean(vpd_mean, na.rm = TRUE),
    heat_days_mean = mean(heat_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    heat_impact_se = heat_impact_sd / sqrt(n_stations),
    heat_impact_lower = heat_impact_mean - heat_impact_sd,
    heat_impact_upper = heat_impact_mean + heat_impact_sd
  )

# 区域统计
region_stats <- monthly_heat_regional %>%
  group_by(region) %>%
  summarise(
    月数 = n(),
    平均站点数 = round(mean(n_stations), 1),
    平均累积热影响 = round(mean(heat_impact_mean), 3),
    最大累积热影响 = round(max(heat_impact_mean), 3),
    .groups = "drop"
  )

cat("区域统计:\n")
print(region_stats)
cat("\n")

# ============================================================================
# 5. 可视化
# ============================================================================

cat("【4. 生成图表】\n")

# 设置区域顺序（从北到南）
region_order <- c("东北", "华北", "西北", "华东", "华中", "西南", "华南")
monthly_heat_regional <- monthly_heat_regional %>%
  mutate(region = factor(region, levels = region_order))

# 图1: 累积热影响指数时间序列（分面图）
p1 <- ggplot(monthly_heat_regional, aes(x = date)) +
  geom_ribbon(
    aes(ymin = pmax(heat_impact_lower, 0), ymax = heat_impact_upper),
    fill = "#d73027", alpha = 0.2
  ) +
  geom_line(aes(y = heat_impact_mean), color = "#d73027", linewidth = 0.6) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = c(0.02, 0)
  ) +
  labs(
    title = "各区域累积热影响指数时间序列 (2000-2023)",
    subtitle = paste0("VPD阈值 = ", VPD_THRESHOLD, " kPa | 实线为区域月平均，阴影为±1标准差"),
    x = "年份",
    y = "累积热影响指数\n(kPa·days/月)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(p1)

# 图2: 各区域对比（同一坐标轴）
p2 <- ggplot(monthly_heat_regional, aes(x = date, y = heat_impact_mean, color = region)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  scale_color_brewer(palette = "Set1", name = "区域") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = c(0.02, 0)
  ) +
  labs(
    title = "各区域累积热影响指数对比",
    subtitle = paste0("VPD阈值 = ", VPD_THRESHOLD, " kPa"),
    x = "年份",
    y = "累积热影响指数 (kPa·days/月)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p2)

# 图3: 年度平均趋势
annual_heat_regional <- monthly_heat_regional %>%
  group_by(region, year) %>%
  summarise(
    heat_impact_annual = mean(heat_impact_mean, na.rm = TRUE),
    heat_impact_sd = sd(heat_impact_mean, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot(annual_heat_regional, aes(x = year, y = heat_impact_annual, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Set1", name = "区域") +
  scale_x_continuous(breaks = seq(2000, 2023, 4)) +
  labs(
    title = "各区域年均累积热影响指数趋势",
    subtitle = paste0("VPD阈值 = ", VPD_THRESHOLD, " kPa | 5-9月平均"),
    x = "年份",
    y = "年均累积热影响指数\n(kPa·days/月)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right"
  )

print(p3)

# 图4: 分面年度趋势
p4 <- ggplot(annual_heat_regional, aes(x = year, y = heat_impact_annual)) +
  geom_line(color = "#d73027", linewidth = 0.8) +
  geom_point(color = "#d73027", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#4575b4",
              linewidth = 0.5, linetype = "dashed", alpha = 0.2) +
  facet_wrap(~ region, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(2000, 2023, 5)) +
  labs(
    title = "各区域年均累积热影响指数趋势（含线性拟合）",
    subtitle = paste0("VPD阈值 = ", VPD_THRESHOLD, " kPa | 虚线为线性趋势"),
    x = "年份",
    y = "年均累积热影响指数\n(kPa·days/月)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(0.8, "lines")
  )

print(p4)

# ============================================================================
# 6. 趋势分析
# ============================================================================

cat("\n【5. 趋势分析】\n\n")

trend_results <- annual_heat_regional %>%
  group_by(region) %>%
  summarise(
    slope = coef(lm(heat_impact_annual ~ year))[2],
    p_value = summary(lm(heat_impact_annual ~ year))$coefficients[2, 4],
    r_squared = summary(lm(heat_impact_annual ~ year))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    trend = case_when(
      p_value < 0.05 & slope > 0 ~ "显著上升 ↑",
      p_value < 0.05 & slope < 0 ~ "显著下降 ↓",
      TRUE ~ "无显著趋势"
    ),
    slope = round(slope * 10, 4),  # 每10年变化
    p_value = round(p_value, 4),
    r_squared = round(r_squared, 3)
  ) %>%
  rename(
    区域 = region,
    "斜率(每10年)" = slope,
    "p值" = p_value,
    "R²" = r_squared,
    趋势 = trend
  )

print(trend_results)

cat("\n图表生成完成!\n")

# 保存数据（可选）
# saveRDS(monthly_heat_regional, "data_proc/monthly_heat_regional.rds")
# saveRDS(annual_heat_regional, "data_proc/annual_heat_regional.rds")
