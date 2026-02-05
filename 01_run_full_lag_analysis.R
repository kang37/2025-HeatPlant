# ============================================================================
# Script: run_full_lag_analysis.R (v2 - hotfix)
# Description: 统一执行VPD对SIF的滞后(lag=0,1,2)因果分析，
#              包含使用S-map确定效应性质，并将最终结果保存为CSV。
# ============================================================================

pacman::p_load(
  dplyr, ggplot2, lubridate, purrr, data.table, stringr, readr,
  tidyr, showtext, rEDM, sf
)
showtext_auto()

# --- 1. 数据加载与准备 ---
cat("--> 1. 加载和准备数据 (此过程较慢)...\n")

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

meteo_sif_data <- read.csv("data_raw/meteo_stat_SIF_data.csv") %>% tibble() %>% rename_with(~tolower(.x)) %>% rename(meteo_stat_id = meteo_stat) %>% filter(!is.na(sif)) %>% mutate(meteo_stat_id = as.character(meteo_stat_id))
station_coords <- meteo_sif_data %>% select(meteo_stat_id, longitude, latitude) %>% distinct(meteo_stat_id, .keep_all = TRUE)
target_stations <- unique(meteo_sif_data$meteo_stat_id)
china_cities_shp <- st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% st_transform(crs = 4326) %>% mutate(pr_name = str_remove(pr_name, "省|市|自治区|回族|维吾尔|壮族"))

read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1)
}
meteo_file_list <- list.files(path = "data_raw/meteo_data_1961-2023", pattern = "\\.txt$", full.names = TRUE) %>% .[!grepl("sta_lonlat_china.txt", .)]
meteo_data_daily <- map_dfr(meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations], read_one_meteo_daily) %>%
  rename_with(~ tolower(.x)) %>%
  mutate(date = as.Date(date), year = year(date), month = month(date)) %>%
  filter(month %in% 5:9)

meteo_data_daily_vpd <- meteo_data_daily %>%
  mutate(
    svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
    avp = (rh / 100) * svp,
    vpd = svp - avp
  ) %>%
  mutate(vpd = if_else(vpd < 0 | vpd > 10, NA_real_, vpd)) %>%
  filter(!is.na(vpd))

VPD_THRESHOLD <- 2.0
monthly_heat <- meteo_data_daily_vpd %>%
  mutate(vpd_excess = pmax(vpd - VPD_THRESHOLD, 0)) %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(cumulative_heat_impact = sum(vpd_excess, na.rm = TRUE), .groups = "drop")

sif_monthly <- meteo_sif_data %>%
  group_by(meteo_stat_id, year, month) %>%
  summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")

safe_detrend <- function(x, time_idx) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch({ residuals(lm(x ~ time_idx)) }, error = function(e) rep(NA_real_, length(x)))
}

data_for_ccm <- monthly_heat %>%
  inner_join(sif_monthly, by = c("meteo_stat_id", "year", "month")) %>%
  filter(!is.na(sif)) %>%
  group_by(meteo_stat_id) %>%
  arrange(year, month) %>%
  mutate(
    time_idx = row_number(),
    sif_detrended = safe_detrend(sif, time_idx),
    heat_impact_detrended = safe_detrend(cumulative_heat_impact, time_idx)
  ) %>%
  ungroup() %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(sif_detrended)) >= 10, sum(!is.na(heat_impact_detrended)) >= 10) %>%
  ungroup()
cat("--> 数据准备完成。\n")

# --- 2. 高级CCM/S-map分析函数 ---
perform_full_lag_analysis <- function(station_id, data, tp_lag = 0, min_points = 30) {
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(sif = sif_detrended, heat_impact = heat_impact_detrended) %>%
    filter(!is.na(sif), !is.na(heat_impact)) %>%
    as.data.frame()
  n_data <- nrow(station_data)
  if (n_data < min_points) return(NULL)
  tryCatch({
    station_data_ccm <- data.frame(time = 1:n_data, sif = station_data$sif, heat_impact = station_data$heat_impact)
    lib_pred_size <- min(100, n_data)
    embed_sif <- EmbedDimension(dataFrame = station_data_ccm, lib = paste("1", lib_pred_size), pred = paste("1", lib_pred_size), columns = "sif", target = "sif", maxE = min(8, floor(n_data / 10)), showPlot = FALSE)
    embed_heat <- EmbedDimension(dataFrame = station_data_ccm, lib = paste("1", lib_pred_size), pred = paste("1", lib_pred_size), columns = "heat_impact", target = "heat_impact", maxE = min(8, floor(n_data / 10)), showPlot = FALSE)
    E_ccm <- max(embed_sif$E[which.max(embed_sif$rho)], embed_heat$E[which.max(embed_heat$rho)])
    tau <- 1
    tp <- tp_lag
    embedding_loss <- (E_ccm - 1) * tau
    max_available_lib <- n_data - embedding_loss - tp
    if (max_available_lib < 15) return(NULL)
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    if (lib_end <= lib_start) return(NULL)
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    ccm_heat_to_sif <- CCM(dataFrame = station_data_ccm, E = E_ccm, Tp = tp, columns = "heat_impact", target = "sif", libSizes = libSizes_str, sample = 50, random = TRUE, showPlot = FALSE)
    ccm_summary <- ccm_heat_to_sif %>% group_by(LibSize) %>% summarise(rho_mean = mean(`heat_impact:sif`, na.rm = TRUE))
    final_rho <- ccm_summary$rho_mean[nrow(ccm_summary)]
    trend <- cor(ccm_summary$LibSize, ccm_summary$rho_mean)
    is_causal <- (final_rho > 0.1 & trend > 0)
    smap_coef <- NA
    if(is_causal) {
      smap_data <- station_data_ccm
      smap_data$heat_impact_lagged <- lag(smap_data$heat_impact, n = tp)
      smap_data <- smap_data %>% filter(!is.na(heat_impact_lagged))
      if(nrow(smap_data) > 15) {
        smap_res <- SMap(dataFrame = smap_data, E = E_ccm, theta = 2, columns = "heat_impact_lagged", target = "sif")
        smap_coeffs <- smap_res$coefficients
        if (!is.null(smap_coeffs) && ncol(smap_coeffs) >= 2) {
          smap_coef <- mean(smap_coeffs[, 2], na.rm = TRUE)
        }
      }
    }
    effect_type <- case_when(!is_causal ~ "无显著因果", is.na(smap_coef) ~ "S-map失败", smap_coef > 0.001 ~ "促进效应(+)", smap_coef < -0.001 ~ "抑制效应(-)", TRUE ~ "效应极弱")
    tibble(meteo_stat_id = station_id, tp_lag = tp_lag, is_causal = is_causal, effect_type = effect_type, rho = final_rho, smap_coef = smap_coef)
  }, error = function(e) { return(NULL) })
}

# --- 3. 批量执行 ---
cat("--> 2. 开始批量分析 (lags = 0, 1, 2)，这将非常耗时...\n")
all_stations <- data_for_ccm %>% count(meteo_stat_id) %>% filter(n >= 30) %>% pull(meteo_stat_id)
cat("分析站点数:", length(all_stations), "\n")

full_lag_results <- map_dfr(c(0, 1, 2), function(lag) {
  cat(sprintf("\n--- 正在分析 Lag = %d ---\n", lag))
  station_results <- map_dfr(seq_along(all_stations), function(i) {
    if (i %% 50 == 0) cat(sprintf("  Lag %d, 已完成: %d / %d\n", lag, i, length(all_stations)))
    perform_full_lag_analysis(all_stations[i], data_for_ccm, tp_lag = lag)
  })
  return(station_results)
})
cat("\n--> 分析完成！\n")

# --- 4. 保存最终结果 (修正部分) ---
cat("--> 3. 正在整合最终结果...\n")

# 1. 创建省份查找表
stations_with_province <- station_coords %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_join(china_cities_shp, join = st_nearest_feature) %>%
  as.data.frame() %>%
  select(meteo_stat_id, pr_name)

# 2. 创建包含所有站点信息的完整查找表
station_info <- station_coords %>%
    left_join(stations_with_province, by = "meteo_stat_id") %>%
    left_join(province_to_region, by = "pr_name")

# 3. 将分析结果与完整的站点信息表连接
final_results_to_save <- full_lag_results %>%
  left_join(station_info, by = "meteo_stat_id")

# 创建输出目录
if (!dir.exists("data_proc")) dir.create("data_proc")
fwrite(final_results_to_save, "data_proc/full_lag_analysis_results.csv")

cat("\n所有分析已完成，最终结果已保存到: data_proc/full_lag_analysis_results.csv\n")
cat("下一步，我们将使用这个文件来制作您需要的所有图表。\n")
