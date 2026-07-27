# =============================================================================
# 05_1_analysis.R
#
# CCM  ：VPD均值（vpd_mean_dt），读取已有结果，不重新跑
# 投资 ：pa_built_10y（单位建成区绿地面积投资，2011-2020均值）
# 回归 ：有序Logit（TRRI为因变量）
#         + 多项Logit（stype四类为因变量）
#         全国 + 分气候区；仅站点尺度
# 方差 ：varpart（投资 vs 地理+控制）
#
# 输出：
#   G1. 中国地图：Köppen气候区底色 + 站点CCM响应类型
#   G2. 热图：站点 × time lag（S-Map系数），气候区分面板
#   G3. 热图：城市 × stype站点占比
#   I.  有序Logit：全部变量系数表（全国+分气候区）
#   J.  多项Logit：全部变量系数表（全国+分气候区）
#   K.  方差分解：各分量独立+共享解释力（全国+分气候区）
# =============================================================================

pacman::p_load(
  dplyr, tidyr, purrr, ggplot2, stringr, readr,
  vegan, tibble, scales, showtext, sysfonts,
  targets, MASS, nnet, ggnewscale, patchwork, ggforce,
  terra, sf, rnaturalearth, rnaturalearthdata,
  geodata, osmextract
)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT     <- "data_proc/output_10y_built_up_05_01"
CCM_RDS <- "data_proc/output_10y_built_up_05_01/ccm_05_results.rds"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 14
theme_cn <- function(bs = BS) {
  theme_minimal(base_size = bs) +
    theme(text          = element_text(family = "heiti"),
          plot.title    = element_text(face = "bold", hjust = 0.5, size = bs + 4),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = bs),
          strip.text    = element_text(face = "bold", size = bs),
          axis.text     = element_text(size = bs - 2),
          axis.title    = element_text(size = bs))
}
theme_set(theme_cn())

stype_levels <- c("always_inhibit", "inhibit_promote",
                  "promote_inhibit", "always_promote")
stype_labels <- c("全程抑制", "抑制→促进", "促进→抑制", "全程促进")
stype_colors <- c(
  always_inhibit  = "#4575B4",
  inhibit_promote = "#74ADD1",
  promote_inhibit = "#FDAE61",
  always_promote  = "#D73027"
)
koppen_colors <- c(A = "#2CA25F", B = "#D95F0E", C = "#3182BD", D = "#756BB1")
koppen_labels <- c(A = "A 热带", B = "B 干旱", C = "C 温带", D = "D 大陆")

winsorize <- function(x, p = 0.99) {
  q <- quantile(x, p, na.rm = TRUE); ifelse(x > q, q, x)
}
sig_star <- function(p) case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "ns"
)


# =============================================================================
# B. 读取 CCM 结果
# =============================================================================

cat("\n=== B. 读取CCM结果 ===\n")
ccm_results <- readRDS(CCM_RDS)
cat(sprintf("站点数: %d | tp: %d..%d | 行数: %d\n",
            length(unique(ccm_results$meteo_stat_id)),
            min(ccm_results$tp), max(ccm_results$tp), nrow(ccm_results)))


# =============================================================================
# C. TRRI / stype 分类
# =============================================================================

cat("\n=== C. TRRI / stype 分类 ===\n")

N_TP     <- length(unique(ccm_results$tp))
MAX_TP   <- max(ccm_results$tp)
TRRI_MAX <- 2L * N_TP
cat(sprintf("N_TP=%d, TRRI范围: 1..%d\n", N_TP, TRRI_MAX))

find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

trri_df <- ccm_results %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq     = list(mean_coef),
    longitude    = first(longitude),
    latitude     = first(latitude),
    koppen_class = first(koppen_class),
    koppen_group = first(koppen_group),
    .groups      = "drop"
  ) %>%
  filter(map_lgl(coef_seq, ~length(.x) == N_TP)) %>%
  mutate(
    tp0_pos = map_lgl(coef_seq, ~!is.na(.x[[1]]) && .x[[1]] > 0),
    HTW     = map_int(coef_seq, ~find_transition(.x, to_neg = TRUE)),
    ITW     = map_int(coef_seq, ~find_transition(.x, to_neg = FALSE)),
    stype   = case_when(
      tp0_pos  & is.na(HTW)  ~ "always_promote",
      !tp0_pos & is.na(ITW)  ~ "always_inhibit",
      tp0_pos  & !is.na(HTW) ~ "promote_inhibit",
      !tp0_pos & !is.na(ITW) ~ "inhibit_promote",
      TRUE                   ~ "always_inhibit"
    ),
    TRRI = case_when(
      stype == "always_inhibit"  ~ 1L,
      stype == "inhibit_promote" ~ as.integer(1L + (N_TP - ITW)),
      stype == "promote_inhibit" ~ as.integer(N_TP + HTW),
      stype == "always_promote"  ~ as.integer(2L * N_TP)
    )
  )

cat("站点类型分布:\n"); print(count(trri_df, stype))
write_csv(trri_df %>% dplyr::select(-coef_seq),
          file.path(OUT, "trri_station.csv"))


# =============================================================================
# D. 投资变量（pa_built_10y）
# =============================================================================

cat("\n=== D. 投资变量构建 ===\n")

tar_load(green_invest_2020)
tar_load(green_area_2020)
tar_load(station_city_map)

invest_raw0 <- read.csv("data_raw/green_invest/city_invest_data.csv",
                        check.names = FALSE)
nc_inv <- ncol(invest_raw0)
names(invest_raw0) <- c("city_raw", paste0("inv_", 2002:(2001 + nc_inv - 1)))

invest_raw <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  dplyr::select(invest_raw0, -city_raw)
)
year_cols <- paste0("inv_", 2002:2030)
year_cols <- year_cols[year_cols %in% names(invest_raw)]
years_num <- as.integer(str_extract(year_cols, "\\d{4}$"))
cols_10y  <- year_cols[years_num >= 2011 & years_num <= 2020]
inv_10y   <- rowMeans(
  dplyr::select(invest_raw, all_of(cols_10y)) %>%
    mutate(across(everything(), as.numeric)),
  na.rm = TRUE
)

green_2020 <- green_area_2020 %>%
  dplyr::select(city_name, area_built = area_green_built) %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))

invest_tbl <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  tibble(inv_10y = inv_10y)
) %>%
  left_join(green_2020, by = "city_name") %>%
  mutate(pa_built_10y = inv_10y / area_built)

cat(sprintf("pa_built_10y 非NA城市: %d\n", sum(!is.na(invest_tbl$pa_built_10y))))


# =============================================================================
# E. 控制变量（人均GDP；CGI政府治理指数2021）
# =============================================================================

cat("\n=== E. 控制变量 ===\n")

city_db_raw <- readxl::read_excel(
  "data_raw/china_city_db.xlsx",
  sheet = 1, col_names = TRUE
) %>%
  rename(year = 1, city_name = 3, pgdp = 15, urban_rate = 26) %>%
  dplyr::select(year, city_name, pgdp, urban_rate) %>%
  filter(year >= 2011, year <= 2020) %>%
  mutate(across(c(pgdp, urban_rate), as.numeric))

city_ctrl <- city_db_raw %>%
  group_by(city_name) %>%
  summarise(pgdp_10y       = mean(pgdp,       na.rm = TRUE),
            urban_rate_10y = mean(urban_rate, na.rm = TRUE),
            .groups = "drop")

# CGI 政府治理与公共服务指数（2021）
cgi_raw <- readxl::read_excel(
  "data_raw/green_invest/cgi_2021.xlsx",
  col_names = TRUE
) %>%
  rename(cgi_rank = 1, city_name = 2, province = 3, cgi_score = 4) %>%
  dplyr::select(city_name, cgi_score) %>%
  mutate(
    # 统一城市名格式：末尾加"市"（部分可能已有）
    city_name = ifelse(str_detect(city_name, "市$"), city_name,
                       paste0(city_name, "市"))
  )

cat(sprintf("CGI数据：%d 城市，得分范围 %.1f ~ %.1f\n",
            nrow(cgi_raw), min(cgi_raw$cgi_score), max(cgi_raw$cgi_score)))

city_vars <- invest_tbl %>%
  left_join(city_ctrl, by = "city_name") %>%
  left_join(cgi_raw,   by = "city_name")

invest_station <- station_city_map %>% left_join(city_vars, by = "city_name")


# =============================================================================
# E2. 站点降水量（meteo_data_1961-2023，生长季4-10月，2011-2020均值）
# =============================================================================

cat("\n=== E2. 站点降水量 ===\n")

METEO_DIR   <- "data_raw/meteo_data_1961-2023"
MISSING_VAL <- 999990   # 气象数据缺失值阈值（≥此值视为缺失）

read_precip_station <- function(fpath) {
  # 第1行：station_id,lon,lat；第2行：列名；其余：数据
  header <- readLines(fpath, n = 1)
  stat_id <- str_split(header, ",")[[1]][1]
  d <- tryCatch(
    read_csv(fpath, skip = 1, col_types = cols(.default = "d", date = "c"),
             show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || !"precip" %in% names(d)) return(NULL)
  d %>%
    mutate(
      meteo_stat_id = stat_id,
      date          = as.Date(date),
      year          = as.integer(format(date, "%Y")),
      month         = as.integer(format(date, "%m")),
      precip_clean  = ifelse(precip >= MISSING_VAL, NA_real_, precip)
    ) %>%
    filter(year >= 2011, year <= 2020, month >= 4, month <= 10) %>%
    dplyr::select(meteo_stat_id, date, year, month, precip_clean)
}

PRECIP_CACHE <- "data_raw/precip_station_cache.rds"

if (file.exists(PRECIP_CACHE)) {
  cat("  [已缓存] 读取降水量站点数据...\n")
  precip_station <- readRDS(PRECIP_CACHE)
} else {
  meteo_files <- list.files(METEO_DIR, pattern = "\\.txt$", full.names = TRUE)
  cat(sprintf("  读取 %d 个气象站文件...\n", length(meteo_files)))

  precip_raw <- map_dfr(meteo_files, read_precip_station)
  cat(sprintf("  有效行数: %d\n", nrow(precip_raw)))

  # 站点级：2011-2020生长季日均降水（mm/d）
  precip_station <- precip_raw %>%
    group_by(meteo_stat_id) %>%
    summarise(
      precip_mean = mean(precip_clean, na.rm = TRUE),
      precip_n    = sum(!is.na(precip_clean)),
      .groups = "drop"
    ) %>%
    filter(precip_n >= 100)   # 至少100天有效观测

  saveRDS(precip_station, PRECIP_CACHE)
  cat(sprintf("  降水量已缓存 -> %s\n", PRECIP_CACHE))
}

cat(sprintf("  降水量站点数（有效观测≥100天）: %d\n", nrow(precip_station)))
cat(sprintf("  降水量范围: %.2f ~ %.2f mm/d\n",
            min(precip_station$precip_mean, na.rm = TRUE),
            max(precip_station$precip_mean, na.rm = TRUE)))


# =============================================================================
# E3. 土壤数据（SoilGrids via geodata，clay+sand，0-5cm，站点坐标提取）
# =============================================================================

cat("\n=== E3. 土壤数据（SoilGrids）===\n")

SOIL_DIR <- "data_raw/soilgrids"
dir.create(SOIL_DIR, recursive = TRUE, showWarnings = FALSE)

# 提取坐标（仅需要有CCM结果的站点）
station_coords <- trri_df %>%
  dplyr::select(meteo_stat_id, longitude, latitude) %>%
  filter(!is.na(longitude), !is.na(latitude))

# 下载中国范围的SoilGrids（clay + sand，0-5cm层，约20-50MB）
# geodata::soil_world 下载全球瓦片，按范围裁剪
china_ext <- terra::ext(72, 136, 17, 54)

get_soil_rast <- function(var, depth = 5, path = SOIL_DIR) {
  fname <- file.path(path, sprintf("soil_%s_%dcm.tif", var, depth))
  if (file.exists(fname)) {
    cat(sprintf("  [已缓存] %s\n", basename(fname)))
    return(terra::rast(fname))
  }
  cat(sprintf("  [下载] %s %dcm...\n", var, depth))
  r <- tryCatch(
    geodata::soil_world(var = var, depth = depth, stat = "mean", path = path),
    error = function(e) { cat(sprintf("  [soil下载失败: %s]\n", e$message)); NULL }
  )
  if (is.null(r)) return(NULL)
  r_crop <- terra::crop(r, china_ext)
  terra::writeRaster(r_crop, fname, overwrite = TRUE)
  r_crop
}

soil_clay <- get_soil_rast("clay", 5)
soil_sand <- get_soil_rast("sand", 5)

# 按站点坐标提取土壤值
extract_soil <- function(r, coords_df, var_name) {
  if (is.null(r)) {
    cat(sprintf("  [跳过 %s：栅格为NULL]\n", var_name))
    return(tibble(meteo_stat_id = coords_df$meteo_stat_id, !!var_name := NA_real_))
  }
  pts <- terra::vect(coords_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  vals <- terra::extract(r, pts)[, 2]
  tibble(meteo_stat_id = coords_df$meteo_stat_id, !!var_name := as.numeric(vals))
}

soil_df <- extract_soil(soil_clay, station_coords, "soil_clay") %>%
  left_join(extract_soil(soil_sand, station_coords, "soil_sand"),
            by = "meteo_stat_id")

cat(sprintf("  clay 非NA: %d / %d\n",
            sum(!is.na(soil_df$soil_clay)), nrow(soil_df)))
cat(sprintf("  sand 非NA: %d / %d\n",
            sum(!is.na(soil_df$soil_sand)), nrow(soil_df)))


# =============================================================================
# E4. 城市规模：常住人口（中国城市数据库1990-2023，2011-2020均值）
# =============================================================================

cat("\n=== E4. 城市常住人口 ===\n")

city_pop_raw <- readxl::read_excel(
  "data_raw/china_city_db2.xlsx",
  sheet = 1, col_names = TRUE
) %>%
  rename(year = 1, city_name = 3, pop_resident = 24) %>%   # col24=常住人口(万人)
  dplyr::select(year, city_name, pop_resident) %>%
  filter(year >= 2011, year <= 2020) %>%
  mutate(pop_resident = as.numeric(pop_resident))

city_pop <- city_pop_raw %>%
  group_by(city_name) %>%
  summarise(pop_10y = mean(pop_resident, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    city_name = ifelse(str_detect(city_name, "市$"), city_name,
                       paste0(city_name, "市"))
  )

cat(sprintf("  常住人口城市数: %d，范围: %.1f ~ %.1f 万人\n",
            sum(!is.na(city_pop$pop_10y)),
            min(city_pop$pop_10y, na.rm = TRUE),
            max(city_pop$pop_10y, na.rm = TRUE)))


# =============================================================================
# E5. 道路密度（OSM Overpass API，站点10km缓冲区内道路长度/面积）
# 逐站点查询 Overpass API，只下载每个站点周边的道路数据，无需整省下载
# =============================================================================

cat("\n=== E5. 道路密度（OSM Overpass API）===\n")

ROAD_CACHE <- "data_raw/road_density_station.rds"

if (file.exists(ROAD_CACHE)) {
  cat("  [已缓存] 读取道路密度...\n")
  road_df <- readRDS(ROAD_CACHE) %>%
    mutate(meteo_stat_id = as.character(meteo_stat_id))
} else {
  if (!requireNamespace("osmdata", quietly = TRUE)) {
    stop("请先安装 osmdata 包: install.packages('osmdata')")
  }

  main_road_types <- c("motorway","trunk","primary","secondary","tertiary",
                       "unclassified","residential",
                       "motorway_link","trunk_link","primary_link",
                       "secondary_link","tertiary_link")
  buf_radius_m  <- 10000          # 10 km 缓冲区半径
  buf_area_km2  <- pi * 10^2      # ≈ 314.16 km²
  # 10km半径对应的经纬度偏移量（粗略，用于Overpass bbox）
  deg_offset    <- 0.09           # ~10km ≈ 0.09°

  cat(sprintf("  共 %d 个站点，逐站查询 Overpass API...\n", nrow(station_coords)))

  query_road_density <- function(sid, lon, lat) {
    bbox <- c(lat - deg_offset, lon - deg_offset,
              lat + deg_offset, lon + deg_offset)
    q <- tryCatch(
      osmdata::opq(bbox = bbox, timeout = 60) %>%
        osmdata::add_osm_feature(key = "highway",
                                 value = main_road_types) %>%
        osmdata::osmdata_sf(quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(q) || is.null(q$osm_lines) || nrow(q$osm_lines) == 0) return(0)

    # 投影到 UTM zone 50N（中国中部，通用）
    roads_proj <- tryCatch(sf::st_transform(q$osm_lines, crs = 32650),
                           error = function(e) NULL)
    if (is.null(roads_proj)) return(0)

    center_sf <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326) %>%
      sf::st_transform(crs = 32650)
    buf <- sf::st_buffer(center_sf, dist = buf_radius_m)

    roads_clip <- tryCatch(sf::st_intersection(roads_proj, buf),
                           error = function(e) NULL)
    if (is.null(roads_clip) || nrow(roads_clip) == 0) return(0)

    as.numeric(sum(sf::st_length(roads_clip), na.rm = TRUE) / 1000) / buf_area_km2
  }

  # 逐站查询，失败则返回 NA，并每10站打印进度
  density_vals <- map_dbl(seq_len(nrow(station_coords)), function(i) {
    if (i %% 10 == 0)
      cat(sprintf("    进度: %d / %d\n", i, nrow(station_coords)))
    row <- station_coords[i, ]
    result <- tryCatch(
      query_road_density(row$meteo_stat_id, row$longitude, row$latitude),
      error = function(e) NA_real_
    )
    if (is.null(result)) NA_real_ else result
  })

  road_df <- station_coords %>%
    mutate(road_density = density_vals)

  n_ok <- sum(!is.na(road_df$road_density))
  cat(sprintf("  道路密度完成：%d / %d 站点有值\n", n_ok, nrow(road_df)))
  saveRDS(road_df, ROAD_CACHE)
  cat(sprintf("  道路密度已缓存 -> %s\n", ROAD_CACHE))
}

if (any(!is.na(road_df$road_density))) {
  cat(sprintf("  road_density 非NA: %d / %d，范围: %.3f ~ %.3f km/km²\n",
              sum(!is.na(road_df$road_density)), nrow(road_df),
              min(road_df$road_density, na.rm = TRUE),
              max(road_df$road_density, na.rm = TRUE)))
} else {
  cat(sprintf("  road_density 非NA: 0 / %d（全部为NA，等待OSM数据）\n", nrow(road_df)))
}


# =============================================================================
# E6. 建筑占地面积（王琳提供，站点30m缓冲区内建筑footprint，单位m²）
# =============================================================================

cat("\n=== E6. 建筑占地面积（王琳数据）===\n")

building_df <- read_csv(
  "data_raw/China_stations_buildings_with_coords30_local_rerun_filled.csv",
  show_col_types = FALSE
) %>%
  transmute(
    meteo_stat_id    = as.character(AirQualityStation),
    building_footprint = as.numeric(buildingFootprint)
  )

cat(sprintf("  建筑数据：%d 个站点，非NA: %d，正值: %d，零值: %d\n",
            nrow(building_df),
            sum(!is.na(building_df$building_footprint)),
            sum(building_df$building_footprint > 0, na.rm = TRUE),
            sum(building_df$building_footprint == 0, na.rm = TRUE)))
cat(sprintf("  范围: %.2f ~ %.2f m²，均值: %.2f m²\n",
            min(building_df$building_footprint, na.rm = TRUE),
            max(building_df$building_footprint, na.rm = TRUE),
            mean(building_df$building_footprint, na.rm = TRUE)))


# =============================================================================
# E7. 建筑平均高度（CMAB数据集，站点30m缓冲区内面积加权平均，单位m）
# 与E6的building_footprint结合计算建筑体积密度（方案A）
# =============================================================================

cat("\n=== E7. 建筑平均高度（CMAB）===\n")

HEIGHT_CACHE <- "data_raw/building_height_station.rds"
BUFFER_AREA_M2 <- pi * 30^2  # 30m缓冲区面积（m²）≈ 2827 m²

if (file.exists(HEIGHT_CACHE)) {
  height_df <- readRDS(HEIGHT_CACHE) %>%
    mutate(meteo_stat_id = as.character(meteo_stat_id))
  cat(sprintf("  建筑高度：%d 个站点，height > 0: %d\n",
              nrow(height_df), sum(height_df$mean_height > 0, na.rm = TRUE)))
  cat(sprintf("  高度范围: %.1f ~ %.1f m\n",
              min(height_df$mean_height[height_df$mean_height > 0], na.rm = TRUE),
              max(height_df$mean_height, na.rm = TRUE)))
} else {
  cat("  [警告] building_height_station.rds 不存在，请先运行 12_building_height.R\n")
  height_df <- tibble(meteo_stat_id = character(), mean_height = numeric(), n_buildings = integer())
}

# =============================================================================
# F. 合并分析数据框（仅站点级）
# =============================================================================

cat("\n=== F. 合并分析数据框 ===\n")

anal_df <- trri_df %>%
  left_join(invest_station,  by = "meteo_stat_id") %>%
  left_join(precip_station,  by = "meteo_stat_id") %>%
  left_join(soil_df,         by = "meteo_stat_id") %>%
  left_join(dplyr::select(road_df, meteo_stat_id, road_density),
                           by = "meteo_stat_id") %>%
  left_join(building_df,     by = "meteo_stat_id") %>%
  left_join(height_df,       by = "meteo_stat_id") %>%
  # 城市级背景变量（通过city_name连接）
  left_join(city_pop,        by = "city_name") %>%
  filter(!is.na(TRRI), !is.na(longitude), !is.na(latitude),
         !is.na(koppen_group)) %>%
  mutate(
    koppen_B         = as.integer(koppen_group == "B"),
    koppen_C         = as.integer(koppen_group == "C"),
    koppen_D         = as.integer(koppen_group == "D"),
    # 管理变量（截尾，不标准化）
    pa_built_10y_w   = winsorize(pa_built_10y),
    cgi_score_w      = winsorize(cgi_score),
    # Local变量
    precip_mean_w    = winsorize(precip_mean),
    soil_clay_w      = winsorize(soil_clay),
    soil_sand_w      = winsorize(soil_sand),
    # 背景变量
    pop_10y_w            = winsorize(pop_10y),
    road_density_w       = winsorize(road_density),
    building_footprint_w = winsorize(building_footprint),
    # 建筑体积密度（方案A）：footprint × 平均高度 / 缓冲区面积（m³/m²）
    building_volume      = replace_na(building_footprint, 0) *
                           replace_na(mean_height, 0),
    building_vol_density = building_volume / BUFFER_AREA_M2,
    building_vol_density_w = winsorize(building_vol_density),
    # 旧变量保留（兼容）
    pgdp_10y_w       = winsorize(pgdp_10y),
    TRRI_ord         = factor(TRRI, levels = 1:TRRI_MAX, ordered = TRUE),
    stype_fct        = factor(stype, levels = stype_levels, labels = stype_labels)
  )

cat(sprintf("站点总数: %d\n", nrow(anal_df)))
cat(sprintf("pa_built_10y 非NA: %d\n",    sum(!is.na(anal_df$pa_built_10y))))
cat(sprintf("cgi_score 非NA: %d\n",       sum(!is.na(anal_df$cgi_score))))
cat(sprintf("precip_mean 非NA: %d\n",     sum(!is.na(anal_df$precip_mean))))
cat(sprintf("soil_clay 非NA: %d\n",       sum(!is.na(anal_df$soil_clay))))
cat(sprintf("pop_10y 非NA: %d\n",         sum(!is.na(anal_df$pop_10y))))
cat(sprintf("road_density 非NA: %d\n",          sum(!is.na(anal_df$road_density))))
cat(sprintf("building_footprint 非NA: %d\n",   sum(!is.na(anal_df$building_footprint))))
cat(sprintf("mean_height > 0: %d\n",           sum(anal_df$mean_height > 0, na.rm = TRUE)))
cat(sprintf("building_vol_density > 0: %d\n",  sum(anal_df$building_vol_density > 0, na.rm = TRUE)))
cat("气候组分布:\n"); print(count(anal_df, koppen_group))

# 城市级（仅用于G3热图，不做回归）
city_agg <- anal_df %>%
  group_by(city_name) %>%
  summarise(
    n_stations        = n(),
    koppen_group      = first(koppen_group),
    n_always_inhibit  = sum(stype == "always_inhibit"),
    n_inhibit_promote = sum(stype == "inhibit_promote"),
    n_promote_inhibit = sum(stype == "promote_inhibit"),
    n_always_promote  = sum(stype == "always_promote"),
    .groups = "drop"
  )


# =============================================================================
# G. 可视化
# =============================================================================

cat("\n=== G. 可视化 ===\n")

# ---------- G1. 中国地图：Köppen底色（仅中国境内）+ 站点CCM类型 + 边际密度 ----------

koppen_tif  <- "data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif"
china_bbox  <- c(xmin = 72, xmax = 136, ymin = 17, ymax = 54)

koppen_to_group <- function(x) {
  case_when(
    x >= 1  & x <= 4  ~ "A",
    x >= 5  & x <= 9  ~ "B",
    x >= 10 & x <= 17 ~ "C",
    x >= 18 & x <= 28 ~ "D",
    TRUE ~ NA_character_
  )
}

# --- 底图矢量 ---
world_land  <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"),
                        error = function(e) NULL)
china_sf    <- tryCatch(
  ne_countries(country = c("China","Hong Kong S.A.R.","Macao S.A.R.",
                            "Taiwan"), scale = "medium", returnclass = "sf") %>%
    st_union(),
  error = function(e) NULL)
china_prov  <- tryCatch(ne_states(country = "China", returnclass = "sf"),
                        error = function(e) NULL)

# --- Köppen 栅格：仅保留中国境内像元 ---
koppen_rast <- tryCatch({
  r <- terra::rast(koppen_tif)
  terra::crop(r, terra::ext(china_bbox))
}, error = function(e) { cat("  [Köppen raster读取失败]\n"); NULL })

koppen_df_china <- NULL
if (!is.null(koppen_rast) && !is.null(china_sf)) {
  # 用中国边界 mask 栅格
  china_vect <- tryCatch(terra::vect(china_sf), error = function(e) NULL)
  if (!is.null(china_vect)) {
    kr_masked <- tryCatch(terra::mask(koppen_rast, china_vect),
                          error = function(e) koppen_rast)  # 退回无mask版
  } else {
    kr_masked <- koppen_rast
  }
  koppen_df_china <- terra::as.data.frame(kr_masked, xy = TRUE) %>%
    rename(koppen_code = 3) %>%
    mutate(koppen_grp = koppen_to_group(koppen_code)) %>%
    filter(!is.na(koppen_grp))
}

map_df <- trri_df %>%
  mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))

stype_colors_cn <- setNames(stype_colors, stype_labels)

# Köppen色系：绿/橙/蓝/紫（暖色系，区别于CCM的红蓝系）
koppen_colors_map <- c(A = "#5A8F76", B = "#EAD5A0", C = "#A3B86C", D = "#C2DFCD")
koppen_labels_map <- c(A = "A 热带/亚热带", B = "B 干旱", C = "C 温带", D = "D 大陆")

# --- 主地图 ---
theme_map <- theme_minimal(base_size = 13) +
  theme(
    text             = element_text(family = "heiti"),
    panel.background = element_rect(fill = "#D6EAF8", color = NA),  # 海洋浅蓝
    panel.grid.major = element_line(color = "white", linewidth = 0.3),
    axis.text        = element_text(size = 10),
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11, face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                    family = "heiti"),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11,
                                    family = "heiti")
  )

p_map_main <- ggplot()

# 世界陆地（浅灰）
if (!is.null(world_land))
  p_map_main <- p_map_main +
    geom_sf(data = world_land, fill = "grey88", color = "grey75",
            linewidth = 0.15, inherit.aes = FALSE)

# 中国境内 Köppen 底色
if (!is.null(koppen_df_china))
  p_map_main <- p_map_main +
    geom_raster(data = koppen_df_china,
                aes(x = x, y = y, fill = koppen_grp), alpha = 0.65) +
    scale_fill_manual(values = koppen_colors_map, labels = koppen_labels_map,
                      name = "Köppen气候区", na.value = "grey88")

# 省界
if (!is.null(china_prov))
  p_map_main <- p_map_main +
    geom_sf(data = china_prov, fill = NA, color = "grey55",
            linewidth = 0.2, inherit.aes = FALSE)

# 站点（CCM类型：红橙黄绿蓝系与Köppen色系区分）
stype_colors_point <- c(
  "全程抑制"  = "#1A6FBF",   # 深蓝
  "抑制→促进" = "#74C6E8",   # 浅蓝
  "促进→抑制" = "#F4A261",   # 橙
  "全程促进"  = "#C1121F"    # 深红
)

p_map_main <- p_map_main +
  ggnewscale::new_scale_color() +
  geom_point(data = map_df,
             aes(x = longitude, y = latitude, color = stype_label),
             size = 1.8, alpha = 0.88, shape = 16) +
  scale_color_manual(values = stype_colors_point, name = "CCM响应类型") +
  coord_sf(xlim = c(72, 136), ylim = c(17, 54), expand = FALSE) +
  labs(title    = "城市植被热响应类型（CCM）与Köppen气候区",
       subtitle = sprintf("X=VPD均值（去趋势）；tp=0..%d；共%d站", MAX_TP, nrow(map_df)),
       x = "经度", y = "纬度") +
  theme_map

# --- 保存主地图 ---
ggsave(file.path(OUT, "map_stype_koppen.png"),
       p_map_main, width = 14, height = 9, dpi = 300)
cat("-> map_stype_koppen.png\n")

# --- 经纬度 × CCM类型占比图（独立输出）---
make_prop_df <- function(data, coord_col, bin_size, min_n = 3) {
  data %>%
    filter(!is.na(.data[[coord_col]]), !is.na(stype_label)) %>%
    mutate(bin = floor(.data[[coord_col]] / bin_size) * bin_size + bin_size / 2) %>%
    group_by(bin, stype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(bin) %>%
    mutate(pct = n / sum(n), total = sum(n)) %>%
    ungroup() %>%
    filter(total >= min_n)
}

lon_prop <- make_prop_df(map_df, "longitude", 3)
lat_prop <- make_prop_df(map_df, "latitude",  2)

theme_prop <- theme_minimal(base_size = 13) +
  theme(text            = element_text(family = "heiti"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 11),
        axis.title       = element_text(size = 12),
        legend.text      = element_text(size = 11, family = "heiti"),
        legend.title     = element_text(size = 12, face = "bold", family = "heiti"),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                        family = "heiti"),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11,
                                        family = "heiti"))

# 图2：经度区间 × CCM占比（横轴=经度，纵轴=占比，堆叠条形）
n_lon_total <- n_distinct(map_df$meteo_stat_id[!is.na(map_df$longitude)])
p_lon_prop <- ggplot(lon_prop,
                     aes(x = factor(bin), y = pct, fill = stype_label)) +
  geom_col(position = "stack", alpha = 0.88, width = 0.85) +
  geom_text(data = lon_prop %>% group_by(bin) %>% slice(1),
            aes(x = factor(bin), y = 1.03, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3, color = "grey40", family = "heiti") +
  scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  scale_x_discrete(labels = function(x)
    paste0(as.numeric(x) - 1.5, "°–", as.numeric(x) + 1.5, "°E")) +
  labs(title    = "不同经度区间各类CCM响应类型的站点占比",
       subtitle = sprintf("每3°经度一组；共%d站；n=各组站点数", n_lon_total),
       x = "经度区间", y = "占比") +
  theme_prop +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(file.path(OUT, "ccm_prop_by_lon.png"),
       p_lon_prop, width = 13, height = 6, dpi = 300)
cat("-> ccm_prop_by_lon.png\n")

# 图3：纬度区间 × CCM占比（纵轴=纬度，横轴=占比，水平堆叠条形）
n_lat_total <- n_distinct(map_df$meteo_stat_id[!is.na(map_df$latitude)])
p_lat_prop <- ggplot(lat_prop,
                     aes(y = factor(bin), x = pct, fill = stype_label)) +
  geom_col(position = "stack", alpha = 0.88, width = 0.85) +
  geom_text(data = lat_prop %>% group_by(bin) %>% slice(1),
            aes(y = factor(bin), x = 1.03, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3, color = "grey40", family = "heiti",
            hjust = 0) +
  scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.12), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x)
    paste0(as.numeric(x) - 1, "°–", as.numeric(x) + 1, "°N")) +
  labs(title    = "不同纬度区间各类CCM响应类型的站点占比",
       subtitle = sprintf("每2°纬度一组；共%d站；n=各组站点数", n_lat_total),
       y = "纬度区间", x = "占比（累加至100%）") +
  theme_prop +
  theme(legend.position = "right")

ggsave(file.path(OUT, "ccm_prop_by_lat.png"),
       p_lat_prop, width = 9, height = 10, dpi = 300)
cat("-> ccm_prop_by_lat.png\n")

# 图3b：纬度 × 总数（条形）+ 各CCM比例（条形右侧饼图）
{
  # 每个纬度bin的总数和行位置（按纬度升序编号）
  bins_sorted <- sort(unique(lat_prop$bin))
  lat_total <- lat_prop %>%
    group_by(bin) %>%
    summarise(total = first(total), .groups = "drop") %>%
    arrange(bin) %>%
    mutate(y_pos = row_number())

  lat_prop_p <- lat_prop %>%
    left_join(dplyr::select(lat_total, bin, y_pos, total), by = c("bin","total"))

  n_bins  <- nrow(lat_total)
  max_n   <- max(lat_total$total)

  # 条形压短：BAR_MAX 控制最大宽度，PIE_X 控制饼图圆心位置
  BAR_MAX <- 3.5      # 条形最大宽度（缩放坐标，压短）
  PIE_X   <- 5.2      # 饼图圆心x坐标
  PIE_R   <- 0.40     # 饼图半径（y单位=1对应1个bin间距）
  BAR_H   <- 0.28     # 条形半高

  lat_total_s <- lat_total %>%
    mutate(x_scaled = total / max_n * BAR_MAX)
  lat_prop_s  <- lat_prop_p %>%
    mutate(x_scaled = total / max_n * BAR_MAX)

  y_labels <- paste0(bins_sorted - 1, "°–", bins_sorted + 1, "°N")

  p_lat_pie <- ggplot() +
    geom_rect(data = lat_total_s,
              aes(xmin = 0, xmax = x_scaled,
                  ymin = y_pos - BAR_H, ymax = y_pos + BAR_H),
              fill = "grey75", color = "grey55", linewidth = 0.3, alpha = 0.85) +
    geom_text(data = lat_total_s,
              aes(x = x_scaled + 0.05, y = y_pos, label = total),
              hjust = 0, size = 3, color = "grey30", family = "heiti") +
    ggforce::geom_arc_bar(
      data = lat_prop_s,
      aes(x0 = PIE_X, y0 = y_pos, r0 = 0, r = PIE_R,
          amount = n, fill = stype_label),
      stat = "pie", color = "white", linewidth = 0.35
    ) +
    scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
    scale_y_continuous(breaks = lat_total_s$y_pos, labels = y_labels,
                       expand = c(0.02, 0.02)) +
    scale_x_continuous(
      breaks = c(0, BAR_MAX/2, BAR_MAX),
      labels = c(0, round(max_n/2), max_n),
      limits = c(0, PIE_X + PIE_R + 0.4),
      expand = c(0, 0)
    ) +
    annotate("text", x = PIE_X, y = n_bins + 0.75,
             label = "CCM占比", size = 3.5, family = "heiti",
             color = "grey30", hjust = 0.5) +
    annotate("segment",
             x = PIE_X - PIE_R - 0.05, xend = PIE_X + PIE_R + 0.05,
             y = n_bins + 0.48, yend = n_bins + 0.48,
             color = "grey60", linewidth = 0.4) +
    coord_fixed() +
    labs(title    = "不同纬度区间各类CCM响应类型（条形=站点总数，饼图=各类占比）",
         subtitle = sprintf("每2°纬度一组；共%d站", n_lat_total),
         x = "站点数", y = "纬度区间") +
    theme_prop +
    theme(legend.position = "right", panel.grid.major.y = element_blank())

  ggsave(file.path(OUT, "ccm_prop_by_lat_pie.png"),
         p_lat_pie, width = 9, height = 10, dpi = 300)
  cat("-> ccm_prop_by_lat_pie.png\n")
}

# 图2b：经度区间 × 总数（竖向条形）+ 各CCM比例（条形上方饼图）
{
  bins_sorted_lon <- sort(unique(lon_prop$bin))
  lon_total <- lon_prop %>%
    group_by(bin) %>%
    summarise(total = first(total), .groups = "drop") %>%
    arrange(bin) %>%
    mutate(x_pos = row_number())

  lon_prop_p <- lon_prop %>%
    left_join(dplyr::select(lon_total, bin, x_pos, total), by = c("bin","total"))

  n_bins_lon <- nrow(lon_total)
  max_n_lon  <- max(lon_total$total)

  # 竖向：条形高度 = 站点数缩放；饼图在条形上方
  BAR_MAX_Y <- 3.5     # 条形最大高度（缩放坐标）
  PIE_Y     <- 5.2     # 饼图圆心y坐标
  PIE_R_L   <- 0.40    # 饼图半径（x单位=1对应1个bin间距）
  BAR_W     <- 0.28    # 条形半宽

  lon_total_s <- lon_total %>%
    mutate(y_scaled = total / max_n_lon * BAR_MAX_Y)
  lon_prop_s  <- lon_prop_p %>%
    mutate(y_scaled = total / max_n_lon * BAR_MAX_Y)

  x_labels_lon <- paste0(bins_sorted_lon - 1.5, "°–",
                          bins_sorted_lon + 1.5, "°E")

  p_lon_pie <- ggplot() +
    # 竖向条形
    geom_rect(data = lon_total_s,
              aes(xmin = x_pos - BAR_W, xmax = x_pos + BAR_W,
                  ymin = 0, ymax = y_scaled),
              fill = "grey75", color = "grey55", linewidth = 0.3, alpha = 0.85) +
    # 条形顶端标注总数
    geom_text(data = lon_total_s,
              aes(x = x_pos, y = y_scaled + 0.08, label = total),
              vjust = 0, size = 2.8, color = "grey30", family = "heiti") +
    # 饼图（条形上方）
    ggforce::geom_arc_bar(
      data = lon_prop_s,
      aes(x0 = x_pos, y0 = PIE_Y, r0 = 0, r = PIE_R_L,
          amount = n, fill = stype_label),
      stat = "pie", color = "white", linewidth = 0.35
    ) +
    scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
    scale_x_continuous(breaks = lon_total_s$x_pos, labels = x_labels_lon,
                       expand = c(0.02, 0.02)) +
    scale_y_continuous(
      breaks = c(0, BAR_MAX_Y/2, BAR_MAX_Y),
      labels = c(0, round(max_n_lon/2), max_n_lon),
      limits = c(-0.1, PIE_Y + PIE_R_L + 0.5),
      expand = c(0, 0)
    ) +
    annotate("text", x = n_bins_lon + 0.7, y = PIE_Y,
             label = "CCM\n占比", size = 3.2, family = "heiti",
             color = "grey30", hjust = 0, vjust = 0.5) +
    annotate("segment",
             x = n_bins_lon + 0.4, xend = n_bins_lon + 0.4,
             y = PIE_Y - PIE_R_L - 0.05, yend = PIE_Y + PIE_R_L + 0.05,
             color = "grey60", linewidth = 0.4) +
    coord_fixed() +
    labs(title    = "不同经度区间各类CCM响应类型（条形=站点总数，饼图=各类占比）",
         subtitle = sprintf("每3°经度一组；共%d站", n_lon_total),
         x = "经度区间", y = "站点数") +
    theme_prop +
    theme(legend.position  = "right",
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1))

  ggsave(file.path(OUT, "ccm_prop_by_lon_pie.png"),
         p_lon_pie, width = 22, height = 7, dpi = 300)
  cat("-> ccm_prop_by_lon_pie.png\n")
}

# ---------- G2. 热图：站点 × time lag，按stype+气候区 ----------

coef_long <- ccm_results %>%
  dplyr::select(meteo_stat_id, tp, mean_coef) %>%
  left_join(trri_df %>% dplyr::select(meteo_stat_id, stype, koppen_group, TRRI),
            by = "meteo_stat_id") %>%
  filter(!is.na(stype)) %>%
  mutate(coef_clip = pmax(pmin(mean_coef, 0.5), -0.5))

MAX_PER_ZONE <- 120
koppen_groups_plot <- sort(unique(coef_long$koppen_group))

p_heatmap_list <- map(koppen_groups_plot, function(grp) {
  d <- coef_long %>% filter(koppen_group == grp)
  stations <- d %>% dplyr::select(meteo_stat_id, stype, TRRI) %>%
    distinct() %>% arrange(factor(stype, levels = stype_levels), TRRI)
  if (nrow(stations) > MAX_PER_ZONE) {
    stations <- stations %>%
      group_by(stype) %>%
      slice_sample(prop = MAX_PER_ZONE / nrow(stations)) %>%
      ungroup() %>% arrange(factor(stype, levels = stype_levels), TRRI)
  }
  stations <- stations %>% mutate(y_rank = row_number())
  d <- d %>%
    filter(meteo_stat_id %in% stations$meteo_stat_id) %>%
    left_join(dplyr::select(stations, meteo_stat_id, y_rank, stype),
              by = c("meteo_stat_id", "stype"))

  # stype分隔线
  sep_lines <- stations %>% group_by(stype) %>%
    summarise(y_max = max(y_rank), .groups = "drop") %>%
    filter(y_max < max(stations$y_rank))

  # stype标注位置
  stype_anno <- stations %>% group_by(stype) %>%
    summarise(y_mid = mean(y_rank), .groups = "drop") %>%
    mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))

  koppen_name <- c(A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")[grp]

  ggplot(d, aes(x = factor(tp), y = factor(y_rank), fill = coef_clip)) +
    geom_tile() +
    geom_hline(data = sep_lines, aes(yintercept = y_max + 0.5),
               color = "white", linewidth = 1.5, inherit.aes = FALSE) +
    annotate("text",
             x    = rep(0.3, nrow(stype_anno)),
             y    = stype_anno$y_mid,
             label = as.character(stype_anno$stype_label),
             hjust = 0, size = 3.2, family = "heiti", color = "grey20") +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish,
      name = "S-Map系数"
    ) +
    scale_x_discrete(labels = paste0("tp", 0:MAX_TP)) +
    labs(
      title    = sprintf("%s气候区：站点热响应热图（n=%d）",
                         koppen_name, length(unique(d$meteo_stat_id))),
      subtitle = "行=站点（全程抑制→…→全程促进）；列=时间滞后；蓝=抑制 红=促进",
      x = "Time lag", y = NULL
    ) +
    theme_cn(bs = 12) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank())
})

iwalk(p_heatmap_list, function(p, i) {
  grp   <- koppen_groups_plot[i]
  fname <- file.path(OUT, sprintf("heatmap_coef_%s.png", grp))
  ggsave(fname, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("-> heatmap_coef_%s.png\n", grp))
})

# ---------- G3. 热图：城市 × stype站点占比 ----------

city_stype_long <- city_agg %>%
  pivot_longer(
    cols      = c(n_always_inhibit, n_inhibit_promote, n_promote_inhibit, n_always_promote),
    names_to  = "stype",
    values_to = "count"
  ) %>%
  mutate(
    stype       = str_remove(stype, "^n_"),
    stype_label = factor(stype, levels = stype_levels, labels = stype_labels),
    pct         = count / n_stations * 100
  )

city_order <- city_agg %>%
  arrange(koppen_group, desc(n_stations)) %>%
  mutate(city_rank = row_number())

city_stype_long <- city_stype_long %>%
  left_join(dplyr::select(city_order, city_name, city_rank), by = "city_name")

koppen_sep <- city_order %>%
  group_by(koppen_group) %>% summarise(y_max = max(city_rank), .groups = "drop") %>%
  filter(y_max < max(city_order$city_rank))

p_city_stype <- ggplot(city_stype_long,
                       aes(x = stype_label, y = factor(city_rank), fill = pct)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_hline(data = koppen_sep, aes(yintercept = y_max + 0.5),
             color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_gradient(low = "white", high = "#D73027",
                      name = "占比(%)", limits = c(0, 100)) +
  labs(
    title    = "城市各热响应类型站点占比",
    subtitle = sprintf("共%d城市，按气候区（A/B/C/D）和站点数降序排列",
                       nrow(city_agg)),
    x = "响应类型", y = "城市（↑占更多站点）"
  ) +
  theme_cn() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid  = element_blank())

ggsave(file.path(OUT, "heatmap_city_stype.png"),
       p_city_stype, width = 9, height = 14, dpi = 300)
cat("-> heatmap_city_stype.png\n")

write_csv(
  city_stype_long %>%
    dplyr::select(city_name, koppen_group, stype_label, count, pct) %>%
    arrange(koppen_group, city_name),
  file.path(OUT, "city_stype_counts.csv")
)
cat("-> city_stype_counts.csv\n")


# =============================================================================
# I. 有序Logit（TRRI ~ 所有变量）
#    仅站点级；全国 + 分气候区；单一控制变量集
#    全国：inv_w + cgi_score_w + lon + lat + Köppen哑变量
#    分区：inv_w + cgi_score_w + lon + lat
#    输出：所有回归变量的系数表
# =============================================================================

cat("\n=== I. 有序Logit（TRRI）===\n")

# 全国：地理+Köppen哑变量；分区内：仅地理（已按区分层，不需要Köppen哑变量）
ctrl_all   <- c("cgi_score_w", "longitude", "latitude",
                "koppen_B", "koppen_C", "koppen_D")
ctrl_inner <- c("cgi_score_w", "longitude", "latitude")

run_olr_full <- function(df, inv_var, ctrl_vars, grp_label) {
  d <- df %>%
    filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0,
           !is.na(TRRI),
           if_all(all_of(ctrl_vars[ctrl_vars %in% names(df)]), ~!is.na(.x)))
  if (nrow(d) < 15) {
    cat(sprintf("  [跳过 %s] n=%d\n", grp_label, nrow(d)))
    return(NULL)
  }
  d$inv_w  <- d[[inv_var]]   # 仅截尾，不标准化
  d$Y_ord  <- factor(d$TRRI, levels = sort(unique(d$TRRI)), ordered = TRUE)
  ctrl_use <- ctrl_vars[ctrl_vars %in% names(d)]
  ctrl_use <- ctrl_use[sapply(ctrl_use, function(v) var(d[[v]], na.rm = TRUE) > 0)]
  all_pred <- c("inv_w", ctrl_use)
  fml <- as.formula(paste("Y_ord ~", paste(all_pred, collapse = "+")))
  m <- tryCatch(
    MASS::polr(fml, data = d, Hess = TRUE, method = "logistic"),
    error = function(e) { cat(sprintf("  [polr error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s <- tryCatch(summary(m),
                error = function(e) { cat(sprintf("  [summary error: %s]\n", e$message)); NULL })
  if (is.null(s)) return(NULL)

  # 用模型实际保留的系数（rank-deficient 时 polr 会自动删掉部分变量）
  kept_pred <- intersect(all_pred, rownames(s$coefficients))
  coef_mat  <- s$coefficients[kept_pred, , drop = FALSE]
  z_vec     <- setNames(coef_mat[, "t value"], kept_pred)
  pval_vec  <- setNames(2 * pnorm(abs(z_vec), lower.tail = FALSE), kept_pred)

  m0  <- tryCatch(MASS::polr(Y_ord ~ 1, data = d, Hess = FALSE),
                  error = function(e) NULL)
  mcf <- if (!is.null(m0))
    round(1 - as.numeric(logLik(m)) / as.numeric(logLik(m0)), 4)
  else NA_real_

  get_v <- function(v, mat, pv) list(
    b = if (v %in% rownames(mat)) mat[v, "Value"] else NA_real_,
    p = if (v %in% names(pv))    pv[v]           else NA_real_,
    s = if (v %in% names(pv))    sig_star(pv[v]) else ""
  )
  iv  <- get_v("inv_w",      coef_mat, pval_vec)
  cgi <- get_v("cgi_score_w", coef_mat, pval_vec)
  cat(sprintf("  [%s] n=%d  inv_w: β=%.4f p=%.4f%s  cgi: β=%.4f p=%.4f%s  McF=%.4f\n",
              grp_label, nrow(d),
              iv$b, iv$p, iv$s, cgi$b, cgi$p, cgi$s, mcf))

  var_labels <- c(
    inv_w       = "投资强度（pa_built_10y）",
    cgi_score_w = "CGI政府治理指数（2021，截尾）",
    longitude   = "经度",
    latitude    = "纬度",
    koppen_B    = "Köppen B（干旱区）",
    koppen_C    = "Köppen C（温带）",
    koppen_D    = "Köppen D（大陆）"
  )

  tibble(
    koppen    = grp_label,
    n         = nrow(d),
    mcfadden  = mcf,
    variable  = rownames(coef_mat),
    var_label = var_labels[rownames(coef_mat)],
    coef      = round(coef_mat[, "Value"],      6),
    se        = round(coef_mat[, "Std. Error"], 6),
    z_val     = round(z_vec,    3),
    p_val     = round(pval_vec, 4),
    sig       = sig_star(pval_vec),
    direction = ifelse(coef_mat[, "Value"] > 0, "+", "-")
  )
}

koppen_all <- c("ALL", sort(unique(anal_df$koppen_group)))

olr_tbl <- map_dfr(koppen_all, function(grp) {
  df_g   <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
  run_olr_full(df_g, "pa_built_10y_w", ctrl_v, grp)
})

cat("\n=== I. 有序Logit系数汇总 ===\n")
print(olr_tbl %>%
        filter(variable %in% c("inv_w", "cgi_score_w")) %>%
        dplyr::select(koppen, variable, n, mcfadden, coef, se, p_val, sig, direction))

write_csv(olr_tbl, file.path(OUT, "olr_trri_full_coef.csv"))
cat("-> olr_trri_full_coef.csv（含所有变量系数）\n")

# 投资+CGI p值热图（气候区为y轴，变量为x轴）
if (nrow(olr_tbl) > 0) {
  koppen_nm <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  var_nm    <- c(inv_w = "投资（pa_built_10y）", cgi_score_w = "CGI治理指数")

  p_olr_pval <- olr_tbl %>%
    filter(variable %in% c("inv_w", "cgi_score_w")) %>%
    mutate(
      grp_label  = factor(koppen_nm[koppen],
                          levels = rev(koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)])),
      var_nm_lbl = factor(var_nm[variable], levels = var_nm),
      cell_label = sprintf("β=%.3f\np=%.3f%s\nn=%d", coef, p_val, sig, n),
      p_fill     = pmin(p_val, 0.5)
    ) %>%
    filter(!is.na(grp_label), !is.na(var_nm_lbl)) %>%
    ggplot(aes(x = var_nm_lbl, y = grp_label, fill = p_fill)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.8, family = "heiti") +
    scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#F0F0F0",
                         midpoint = 0.05, limits = c(0, 0.5),
                         name = "p值（红=显著）") +
    labs(title    = "有序Logit：管理变量对TRRI的效应",
         subtitle = "因变量=TRRI（有序因子，1-18）；管理=投资强度+CGI治理指数",
         x = NULL, y = NULL) +
    theme_cn() + theme(panel.grid = element_blank())

  ggsave(file.path(OUT, "olr_trri_invest_summary.png"),
         p_olr_pval, width = 8, height = 6, dpi = 300)
  cat("-> olr_trri_invest_summary.png\n")
}


# =============================================================================
# I2. 有序Logit 可视化
#     (a) 森林图：所有预测变量系数 ± 95%CI，分气候区
#     (b) 散点+箱线：投资/CGI 分位 × TRRI（原始数据关系）
#     (c) 预测概率曲线：固定其他变量，投资/CGI 沿范围变化时 P(TRRI=k) 变化
# =============================================================================

cat("\n=== I2. 有序Logit 可视化 ===\n")

if (nrow(olr_tbl) > 0) {
  koppen_nm  <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  koppen_col <- c("全部"="#555555","A(热带)"="#1B9E77","B(干旱)"="#D95F02",
                  "C(温带)"="#7570B3","D(大陆)"="#E7298A")

  # (a) 森林图：系数 ± 1.96*SE -------------------------------------------
  forest_df <- olr_tbl %>%
    mutate(
      grp_label = factor(koppen_nm[koppen],
                         levels = koppen_nm[c("ALL","B","C","D")]),
      ci_lo     = coef - 1.96 * se,
      ci_hi     = coef + 1.96 * se,
      var_short = case_when(
        variable == "inv_w"       ~ "单位面积投资",
        variable == "cgi_score_w" ~ "CGI治理指数",
        variable == "longitude"   ~ "经度",
        variable == "latitude"    ~ "纬度",
        variable == "koppen_B"    ~ "Köppen B",
        variable == "koppen_C"    ~ "Köppen C",
        variable == "koppen_D"    ~ "Köppen D",
        TRUE ~ variable
      ),
      var_short = factor(var_short,
                         levels = rev(c("单位面积投资","CGI治理指数","经度","纬度",
                                        "Köppen B","Köppen C","Köppen D")))
    ) %>%
    filter(!is.na(grp_label))

  p_forest <- ggplot(forest_df,
                     aes(x = coef, y = var_short, color = grp_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.6) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(0.6), linewidth = 0.8) +
    geom_point(aes(shape = sig %in% c("*","**","***",".")),
               position = position_dodge(0.6), size = 2.8) +
    scale_color_manual(values = koppen_col, name = "气候区") +
    scale_shape_manual(values = c("TRUE"=16, "FALSE"=1),
                       labels = c("TRUE"="显著(p<0.1)", "FALSE"="不显著"),
                       name = "") +
    labs(title    = "有序Logit各预测变量系数（森林图）",
         subtitle = "水平线=95%CI；实心圆=p<0.1显著；纵轴=预测变量；横轴=log-odds系数",
         x = "系数（Log-odds）", y = NULL) +
    theme_cn() +
    theme(legend.position = "right",
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  ggsave(file.path(OUT, "olr_forest.png"),
         p_forest, width = 10, height = 6, dpi = 300)
  cat("-> olr_forest.png\n")

  # (b) 散点+箱线：将连续投资/CGI 分成5分位组，展示各组 TRRI 分布 --------
  make_quantile_box <- function(df, xvar, xlabel) {
    df %>%
      filter(!is.na(.data[[xvar]]), .data[[xvar]] > 0,
             !is.na(TRRI), !is.na(koppen_group)) %>%
      mutate(
        q5     = cut(.data[[xvar]],
                     breaks = quantile(.data[[xvar]], probs = 0:5/5, na.rm = TRUE),
                     include.lowest = TRUE, labels = paste0("Q", 1:5)),
        kgroup = factor(koppen_nm[koppen_group],
                        levels = koppen_nm[c("B","C","D")])
      ) %>%
      filter(!is.na(q5), !is.na(kgroup)) %>%
      ggplot(aes(x = q5, y = TRRI, fill = kgroup)) +
      geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5,
                   linewidth = 0.5, alpha = 0.75) +
      geom_smooth(aes(x = as.integer(q5), group = kgroup, color = kgroup),
                  method = "loess", se = FALSE, linewidth = 1.2,
                  show.legend = FALSE) +
      scale_fill_manual(values  = c("B(干旱)"="#D95F02","C(温带)"="#7570B3",
                                    "D(大陆)"="#E7298A"), name = "气候区") +
      scale_color_manual(values = c("B(干旱)"="#D95F02","C(温带)"="#7570B3",
                                    "D(大陆)"="#E7298A")) +
      scale_y_continuous(breaks = seq(1, 18, 2)) +
      labs(title    = sprintf("%s 五分位组 × TRRI 分布", xlabel),
           subtitle = "箱线=各分位组TRRI分布；折线=loess趋势",
           x = sprintf("%s（Q1=最低，Q5=最高）", xlabel),
           y = "TRRI（1=全程抑制，18=全程促进）") +
      theme_cn() +
      theme(legend.position = "right")
  }

  p_box_inv <- make_quantile_box(anal_df, "pa_built_10y_w", "单位面积投资")
  p_box_cgi <- make_quantile_box(anal_df, "cgi_score_w",    "CGI治理指数")

  ggsave(file.path(OUT, "olr_box_invest.png"),
         p_box_inv, width = 10, height = 5, dpi = 300)
  cat("-> olr_box_invest.png\n")

  ggsave(file.path(OUT, "olr_box_cgi.png"),
         p_box_cgi, width = 10, height = 5, dpi = 300)
  cat("-> olr_box_cgi.png\n")

  # (c) 预测概率曲线：固定其他变量于中位数，沿自变量范围预测 TRRI --------
  # 用全国模型（ALL）预测；展示 predicted mean TRRI（加权期望值）
  make_pred_curve <- function(xvar, xlabel, n_pts = 50) {
    d_all <- anal_df %>%
      filter(!is.na(.data[[xvar]]), .data[[xvar]] > 0,
             !is.na(TRRI), !is.na(cgi_score_w),
             !is.na(longitude), !is.na(latitude),
             !is.na(koppen_B))
    if (nrow(d_all) < 30) return(NULL)

    # 重建全国模型（与 run_olr_full 一致，但用于预测）
    d_all$inv_w <- d_all[[xvar]]
    d_all$Y_ord <- factor(d_all$TRRI, levels = sort(unique(d_all$TRRI)), ordered = TRUE)
    ctrl_use <- ctrl_all[ctrl_all %in% names(d_all)]
    ctrl_use <- ctrl_use[sapply(ctrl_use, function(v) var(d_all[[v]], na.rm=TRUE) > 0)]
    # 投资变量名替换
    pred_vars <- c("inv_w", ctrl_use[ctrl_use != xvar])
    fml <- as.formula(paste("Y_ord ~", paste(pred_vars, collapse = "+")))
    m <- tryCatch(
      MASS::polr(fml, data = d_all, Hess = FALSE, method = "logistic"),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)

    # 预测数据框：xvar 沿范围变化，其余固定于中位数
    xseq   <- seq(min(d_all[[xvar]], na.rm=TRUE),
                  max(d_all[[xvar]], na.rm=TRUE), length.out = n_pts)
    newdat <- d_all[rep(1, n_pts), pred_vars, drop = FALSE]
    for (v in pred_vars) {
      if (v != "inv_w") newdat[[v]] <- median(d_all[[v]], na.rm = TRUE)
    }
    newdat$inv_w <- xseq
    # 强制转换 Y_ord levels 与训练集一致
    levels_ord <- levels(d_all$Y_ord)
    probs <- tryCatch(predict(m, newdat, type = "probs"), error = function(e) NULL)
    if (is.null(probs)) return(NULL)
    if (!is.matrix(probs)) probs <- matrix(probs, nrow = 1)

    # 计算预测 E[TRRI]
    trri_vals <- as.integer(levels_ord)
    pred_mean <- as.numeric(probs %*% trri_vals)

    # 95% CI via delta method（简化：用 ±1 TRRI 单位作为示意）
    tibble(x = xseq, pred_trri = pred_mean,
           xlabel = xlabel, xvar = xvar)
  }

  pred_inv <- make_pred_curve("pa_built_10y_w", "单位面积投资（万元/公顷）")
  pred_cgi <- make_pred_curve("cgi_score_w",    "CGI治理指数")

  pred_all <- bind_rows(pred_inv, pred_cgi) %>% filter(!is.na(pred_trri))

  if (nrow(pred_all) > 0) {
    # 同时展示两条预测曲线（分面）
    p_pred <- ggplot(pred_all, aes(x = x, y = pred_trri)) +
      geom_line(color = "#D73027", linewidth = 1.2) +
      geom_rug(data = anal_df %>%
                 pivot_longer(c(pa_built_10y_w, cgi_score_w),
                              names_to = "xvar", values_to = "x") %>%
                 filter(!is.na(x), x > 0) %>%
                 mutate(xlabel = ifelse(xvar == "pa_built_10y_w",
                                        "单位面积投资（万元/公顷）", "CGI治理指数")),
               aes(x = x), sides = "b", alpha = 0.15, color = "grey40",
               inherit.aes = FALSE) +
      facet_wrap(~xlabel, scales = "free_x") +
      labs(title    = "有序Logit：预测 TRRI 随管理变量的变化（全国模型）",
           subtitle = "其余变量固定于中位数；竖线密度=实际数据分布",
           x = "自变量取值", y = "预测 TRRI（期望值）") +
      theme_cn()

    ggsave(file.path(OUT, "olr_pred_curve.png"),
           p_pred, width = 11, height = 5, dpi = 300)
    cat("-> olr_pred_curve.png\n")
  }
}


# =============================================================================
# J. 多项Logit：stype（四类名义变量）~ pa_built_10y + cgi + 控制
#    仅站点级；全国 + 分气候区；输出所有变量系数
# =============================================================================

cat("\n=== J. 多项Logit（stype四类）===\n")

run_mnl_full <- function(df, inv_var, ctrl_vars, grp_label,
                         ref = "always_inhibit") {
  d <- df %>%
    filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0, !is.na(stype),
           if_all(all_of(ctrl_vars[ctrl_vars %in% names(df)]), ~!is.na(.x)))
  if (nrow(d) < 20) {
    cat(sprintf("  [跳过 %s] n=%d\n", grp_label, nrow(d)))
    return(NULL)
  }
  d$inv_w   <- d[[inv_var]]  # 仅截尾，不标准化
  d$stype_f <- relevel(factor(d$stype, levels = stype_levels), ref = ref)
  ctrl_use  <- ctrl_vars[ctrl_vars %in% names(d)]
  ctrl_use  <- ctrl_use[sapply(ctrl_use, function(v) var(d[[v]], na.rm = TRUE) > 0)]
  all_pred  <- c("inv_w", ctrl_use)
  fml <- as.formula(paste("stype_f ~", paste(all_pred, collapse = "+")))
  m <- tryCatch(
    suppressWarnings(nnet::multinom(fml, data = d, trace = FALSE, MaxNWts = 10000)),
    error = function(e) { cat(sprintf("  [multinom error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s      <- summary(m)
  coef_m <- s$coefficients
  se_m   <- s$standard.errors
  z_m    <- coef_m / se_m
  pval_m <- 2 * pnorm(abs(z_m), lower.tail = FALSE)

  cat(sprintf("  MNL [%s] n=%d  inv_w p: %s\n",
              grp_label, nrow(d),
              paste(sprintf("%s:%.3f%s", rownames(coef_m),
                            pval_m[, "inv_w"], sig_star(pval_m[, "inv_w"])),
                    collapse = " ")))

  var_labels <- c(
    inv_w       = "投资强度（pa_built_10y）",
    cgi_score_w = "CGI政府治理指数（2021，截尾）",
    longitude   = "经度",
    latitude    = "纬度",
    koppen_B    = "Köppen B（干旱区）",
    koppen_C    = "Köppen C（温带）",
    koppen_D    = "Köppen D（大陆）"
  )

  map_dfr(rownames(coef_m), function(cat) {
    vars <- colnames(coef_m)
    tibble(
      koppen      = grp_label,
      reference   = ref,
      outcome_cat = cat,
      n           = nrow(d),
      variable    = vars,
      var_label   = var_labels[vars],
      coef        = round(coef_m[cat, vars], 6),
      se          = round(se_m[cat,   vars], 6),
      z_val       = round(z_m[cat,    vars], 3),
      p_val       = round(pval_m[cat,  vars], 4),
      sig         = sig_star(pval_m[cat, vars]),
      direction   = ifelse(coef_m[cat, vars] > 0,
                           "↑更可能落入此类", "↓更不可能落入此类")
    )
  })
}

mnl_tbl <- map_dfr(koppen_all, function(grp) {
  df_g   <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
  run_mnl_full(df_g, "pa_built_10y_w", ctrl_v, grp)
})

cat("\n=== J. 多项Logit：管理变量系数汇总 ===\n")
print(mnl_tbl %>%
        filter(variable %in% c("inv_w", "cgi_score_w")) %>%
        dplyr::select(koppen, outcome_cat, variable, n, coef, se, p_val, sig))

write_csv(mnl_tbl, file.path(OUT, "mnl_stype_full_coef.csv"))
cat("-> mnl_stype_full_coef.csv（含所有变量系数）\n")

# 投资+CGI系数热图（气候区 × 响应类型，两个变量并列）
if (nrow(mnl_tbl) > 0) {
  koppen_nm  <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  stype_lmap <- setNames(stype_labels, stype_levels)
  var_nm     <- c(inv_w = "投资", cgi_score_w = "CGI治理")

  p_mnl <- mnl_tbl %>%
    filter(variable %in% c("inv_w", "cgi_score_w")) %>%
    mutate(
      grp_label  = factor(koppen_nm[koppen],
                          levels = koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)]),
      cat_label  = factor(stype_lmap[outcome_cat], levels = stype_labels),
      var_lbl    = factor(var_nm[variable], levels = var_nm),
      cell_label = sprintf("β=%.2f\np=%.3f%s", coef, p_val, sig),
      coef_clip  = pmax(pmin(coef, 1), -1)
    ) %>%
    filter(!is.na(grp_label), !is.na(cat_label)) %>%
    ggplot(aes(x = cat_label, y = grp_label, fill = coef_clip)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.2, family = "heiti") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                         midpoint = 0, name = "β（截断±1）") +
    facet_wrap(~var_lbl, nrow = 1) +
    labs(title    = "多项Logit：管理变量对响应类型归属的效应",
         subtitle = "参照类别=全程抑制；β>0表示该变量增加→更可能落入该类",
         x = "对比类别（vs 全程抑制）", y = NULL) +
    theme_cn() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          panel.grid  = element_blank())

  ggsave(file.path(OUT, "mnl_stype_invest_beta.png"),
         p_mnl, width = 14, height = 6, dpi = 300)
  cat("-> mnl_stype_invest_beta.png\n")
}


# =============================================================================
# K. 方差分解（varpart）
#    分组：[管理] = pa_built_10y_w + cgi_score_w
#          [地理] = lon/lat/Köppen/pgdp（基础版不含pgdp）
#    输出四分量：管理独立[a]、地理独立[b]、共享[c]、未解释[d]
#    仅站点级；全国 + 分气候区；基础控制 + 扩展控制（含pgdp）
# =============================================================================

cat("\n=== K. 方差分解（三组）===\n")
# 组1 管理：投资强度 + CGI
# 组2 Local：降水量 + 土壤（clay/sand）
# 组3 背景：常住人口 + 道路密度

mgmt_vars  <- c("pa_built_10y_w", "cgi_score_w")
local_vars <- c("precip_mean_w")
bg_vars    <- c("pop_10y_w", "road_density_w", "building_footprint_w")

run_vp3 <- function(df, outcome_var, mv, lv, bv, grp_label) {
  req_vars <- c(mv, lv, bv)
  d <- df %>%
    filter(!is.na(pa_built_10y), pa_built_10y > 0,
           !is.na(.data[[outcome_var]])) %>%
    filter(if_all(all_of(req_vars[req_vars %in% names(.)]), ~!is.na(.x)))

  if (nrow(d) < 20) {
    cat(sprintf("  [跳过 %s] n=%d（满足三组完整数据不足20）\n", grp_label, nrow(d)))
    return(NULL)
  }

  # 过滤零方差变量
  filter_nonzero <- function(vars) {
    vars[vars %in% names(d) & sapply(vars[vars %in% names(d)],
                                      function(v) var(d[[v]], na.rm=TRUE) > 0)]
  }
  mv_use <- filter_nonzero(mv)
  lv_use <- filter_nonzero(lv)
  bv_use <- filter_nonzero(bv)

  if (length(mv_use) == 0 || length(lv_use) == 0 || length(bv_use) == 0) {
    cat(sprintf("  [跳过 %s] 某组变量全为零方差\n", grp_label))
    return(NULL)
  }

  vp <- tryCatch(
    vegan::varpart(d[[outcome_var]],
                   dplyr::select(d, all_of(mv_use)),
                   dplyr::select(d, all_of(lv_use)),
                   dplyr::select(d, all_of(bv_use))),
    error = function(e) { cat(sprintf("  [varpart error: %s]\n", e$message)); NULL })
  if (is.null(vp)) return(NULL)

  # indfract 三分组时有8行（2^3），行序：
  # [a]=X1独立, [b]=X2独立, [c]=X3独立,
  # [d]=X1∩X2, [e]=X1∩X3, [f]=X2∩X3, [g]=X1∩X2∩X3, [h]=未解释
  fr <- vp$part$indfract$Adj.R.square

  cat(sprintf(
    "  [%s] n=%d  管理=%.2f%% Local=%.2f%% 背景=%.2f%% 未解释=%.2f%%\n",
    grp_label, nrow(d),
    fr[1]*100, fr[2]*100, fr[3]*100, fr[8]*100))

  tibble(
    koppen        = grp_label,
    n             = nrow(d),
    mgmt_vars_used  = paste(mv_use, collapse = "+"),
    local_vars_used = paste(lv_use, collapse = "+"),
    bg_vars_used    = paste(bv_use, collapse = "+"),
    mgmt_only     = round(fr[1]*100, 2),
    local_only    = round(fr[2]*100, 2),
    bg_only       = round(fr[3]*100, 2),
    mgmt_local    = round(fr[4]*100, 2),
    mgmt_bg       = round(fr[5]*100, 2),
    local_bg      = round(fr[6]*100, 2),
    triple        = round(fr[7]*100, 2),
    unexplained   = round(fr[8]*100, 2)
  )
}

vp_tbl <- map_dfr(koppen_all, function(grp) {
  df_g <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  run_vp3(df_g, "TRRI", mgmt_vars, local_vars, bg_vars, grp)
})

cat("\n=== K. 方差分解汇总 ===\n")
print(vp_tbl %>%
        dplyr::select(koppen, n, mgmt_only, local_only, bg_only, unexplained))

write_csv(vp_tbl, file.path(OUT, "varpart_results.csv"))
cat("-> varpart_results.csv\n")

koppen_nm <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")

if (nrow(vp_tbl) > 0) {
  vp_plot_df <- vp_tbl %>%
    mutate(grp_label = factor(koppen_nm[koppen],
                              levels = koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)])) %>%
    filter(!is.na(grp_label))

  # 图1：三组独立解释力堆叠条形图（仅独立分量，忽略共享）
  p_vp_stack <- vp_plot_df %>%
    pivot_longer(c(mgmt_only, local_only, bg_only, unexplained),
                 names_to = "component", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(component,
                          levels = c("mgmt_only","local_only","bg_only","unexplained"),
                          labels = c("[a] 管理（独立）","[b] Local（独立）",
                                     "[c] 背景（独立）","[d] 未解释/共享"))
    ) %>%
    ggplot(aes(x = grp_label, y = r2_show, fill = comp_label)) +
    geom_col(position = "stack", alpha = 0.88, width = 0.65) +
    geom_text(aes(label = ifelse(r2_show >= 0.5, sprintf("%.1f%%", r2_show), "")),
              position = position_stack(vjust = 0.5),
              size = 3.5, family = "heiti", color = "white") +
    scale_fill_manual(
      values = c("[a] 管理（独立）"="#D73027", "[b] Local（独立）"="#2CA25F",
                 "[c] 背景（独立）"="#4575B4", "[d] 未解释/共享"="#CCCCCC"),
      name = "方差分量") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title    = "Variance Partitioning (3 Groups): Mgmt vs Local vs Background",
         subtitle = paste0("Outcome: TRRI (all stations, 1–18)\n",
                           "[a] Mgmt = Investment (pa_built_10y) + CGI Governance\n",
                           "[b] Local = Precipitation; [c] Background = Population + Road Density + Building Footprint (2D)\n",
                           "Adj. R² unique fractions (negative shown as 0)"),
         x = NULL, y = "调整R²（%）") +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "varpart_stack.png"),
         p_vp_stack, width = 11, height = 5, dpi = 300)
  cat("-> varpart_stack.png\n")

  # 图2：三组独立解释力并排比较
  p_vp_compare <- vp_plot_df %>%
    pivot_longer(c(mgmt_only, local_only, bg_only),
                 names_to = "group", values_to = "r2") %>%
    mutate(
      r2_show   = pmax(r2, 0),
      grp_lbl   = factor(group,
                          levels = c("mgmt_only","local_only","bg_only"),
                          labels = c("管理","Local","背景")),
      r2_label  = sprintf("%.2f%%", r2)
    ) %>%
    ggplot(aes(x = grp_label, y = r2_show, fill = grp_lbl)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    geom_text(aes(label = r2_label),
              position = position_dodge(0.7), vjust = -0.4,
              size = 3.2, family = "heiti") +
    scale_fill_manual(
      values = c("管理"="#D73027","Local"="#2CA25F","背景"="#4575B4"),
      name = "变量组") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title    = "Variance Partitioning: Unique Fractions by Climate Zone",
         subtitle = paste0("Outcome: TRRI (all stations, 1–18)\n",
                           "Mgmt = Investment + CGI; Local = Precipitation; Background = Population + Road Density + Building Footprint (2D)\n",
                           "(negative = unique contribution below random, shown as 0)"),
         x = NULL, y = "独立调整R²（%）") +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "varpart_compare.png"),
         p_vp_compare, width = 11, height = 5, dpi = 300)
  cat("-> varpart_compare.png\n")
}


# =============================================================================
# L. 气候区 × CCM类型 分层回归（有序Logit：TRRI ~ 管理变量）
#    分层：koppen_group（A/B/C/D）× stype（4类）→ 最多16个子组
#    仅站点级；控制变量：经纬度（子组内不再含Köppen哑变量）
#    输出：子组系数表 + 气泡热图（效应大小/方向/显著性）
# =============================================================================

cat("\n=== L. 气候区×CCM类型 分层回归 ===\n")

koppen_groups <- sort(unique(anal_df$koppen_group))
# stype_levels 已在开头定义

strata_tbl <- map_dfr(koppen_groups, function(kg) {
  map_dfr(stype_levels, function(st) {
    df_sub <- anal_df %>%
      filter(koppen_group == kg, stype == st)
    grp_label <- sprintf("%s|%s", kg, st)

    # 控制：经纬度（无Köppen哑变量，子组内气候区已固定）
    ctrl_sub <- c("cgi_score_w", "longitude", "latitude")

    result <- run_olr_full(df_sub, "pa_built_10y_w", ctrl_sub, grp_label)
    if (!is.null(result)) {
      result %>% mutate(koppen_group = kg, stype = st, .before = 1)
    }
  })
})

if (nrow(strata_tbl) > 0) {
  cat("\n=== L. 分层回归：管理变量系数汇总 ===\n")
  print(strata_tbl %>%
          filter(variable %in% c("inv_w", "cgi_score_w")) %>%
          dplyr::select(koppen_group, stype, variable, n, coef, se, p_val, sig))

  write_csv(strata_tbl, file.path(OUT, "olr_strata_koppen_stype_coef.csv"))
  cat("-> olr_strata_koppen_stype_coef.csv\n")

  # 气泡热图：横轴=stype，纵轴=气候区，大小=|coef|，颜色=方向，透明度=显著性
  stype_lmap_short <- c(
    always_inhibit  = "全程抑制",
    inhibit_promote = "抑制→促进",
    promote_inhibit = "促进→抑制",
    always_promote  = "全程促进"
  )
  koppen_nm_short <- c(A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")

  bubble_df <- strata_tbl %>%
    filter(variable %in% c("inv_w", "cgi_score_w")) %>%
    mutate(
      stype_lbl  = factor(stype_lmap_short[stype], levels = stype_labels),
      koppen_lbl = factor(koppen_nm_short[koppen_group],
                          levels = koppen_nm_short[c("A","B","C","D")]),
      var_lbl    = ifelse(variable == "inv_w", "投资强度", "CGI治理指数"),
      sig_alpha  = case_when(
        sig %in% c("***","**","*") ~ 1.0,
        sig == "."                 ~ 0.65,
        TRUE                       ~ 0.30
      ),
      coef_dir   = ifelse(coef > 0, "正（促进TRRI）", "负（抑制TRRI）")
    ) %>%
    filter(!is.na(stype_lbl), !is.na(koppen_lbl))

  p_bubble <- ggplot(bubble_df,
                     aes(x = stype_lbl, y = koppen_lbl,
                         size = pmin(abs(coef), 2),
                         color = coef_dir, alpha = sig_alpha)) +
    geom_point(shape = 16) +
    geom_text(aes(label = sprintf("β=%.2f\n%s", coef, sig)),
              size = 2.6, family = "heiti", vjust = 2.5, alpha = 1) +
    scale_size_continuous(range = c(2, 12), name = "|β|（截断≤2）") +
    scale_color_manual(
      values = c("正（促进TRRI）"="#D73027","负（抑制TRRI）"="#4575B4"),
      name = "效应方向") +
    scale_alpha_identity() +
    facet_wrap(~var_lbl, nrow = 1) +
    labs(title    = "气候区 × CCM类型 分层有序Logit系数",
         subtitle = "气泡大小=|β|；颜色=方向；透明度：实=p<0.05, 半透=p<0.1, 淡=不显著\nn<15的子组已跳过",
         x = "CCM响应类型（stype）", y = "气候区") +
    theme_cn() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          panel.grid.major = element_line(color = "grey92"),
          legend.position = "right")

  ggsave(file.path(OUT, "strata_bubble_koppen_stype.png"),
         p_bubble, width = 14, height = 6, dpi = 300)
  cat("-> strata_bubble_koppen_stype.png\n")
} else {
  cat("  [警告] 分层回归无有效结果（各子组n均不足15）\n")
}


# =============================================================================
# M. 促进/抑制分组分析
#    促进组：stype %in% c("always_promote", "promote_inhibit")
#    抑制组：stype %in% c("always_inhibit", "inhibit_promote")
#    对每组分别做：
#      M1. 有序Logit（全国 + 分气候区）
#      M2. 三组方差分解（varpart）
# =============================================================================

cat("\n=== M. 促进/抑制分组分析 ===\n")

# TRRI 组内重新编码（统一为 1=最差, 9=最好）：
#   抑制组 (TRRI 1-9) : 1=全程抑制, 9=很快转促进  → 保持原值
#   促进组 (TRRI 10-18): 10=很快转抑制, 18=全程促进 → 减9，映射到 1-9
anal_df <- anal_df %>%
  mutate(
    response_group = ifelse(
      stype %in% c("always_promote", "promote_inhibit"), "促进", "抑制"
    ),
    TRRI_g = ifelse(response_group == "促进", TRRI - 9L, TRRI),
    TRRI_g_ord = factor(TRRI_g, levels = 1:9, ordered = TRUE)
  )

group_labels <- c("促进", "抑制")

# 导出站点级协变量表（供 15_hcsif 等新分析按 meteo_stat_id 复用，避免重算协变量）
saveRDS(
  anal_df %>%
    dplyr::select(meteo_stat_id, city_name, longitude, latitude,
                  koppen_class, koppen_group,
                  pa_built_10y, cgi_score, precip_mean,
                  soil_clay, soil_sand, pop_10y, road_density,
                  building_footprint, mean_height, building_vol_density,
                  pgdp_10y),
  file.path(OUT, "station_covariates.rds")
)
cat("-> station_covariates.rds（站点级协变量，供新分析复用）\n")

# M0. 原始数据关系图：投资/CGI vs TRRI_g（分促进/抑制组）
plot_df_m0 <- anal_df %>%
  filter(!is.na(pa_built_10y_w), !is.na(cgi_score_w), !is.na(TRRI_g)) %>%
  mutate(
    rg_label = factor(response_group,
                      levels = c("抑制", "促进"),
                      labels = c("抑制组（1=全程抑制, 9=很快转促进）",
                                 "促进组（1=很快转抑制, 9=全程促进）"))
  )

# 投资强度 vs TRRI_g
p_raw_inv <- plot_df_m0 %>%
  ggplot(aes(x = pa_built_10y_w, y = TRRI_g)) +
  geom_jitter(aes(color = response_group), height = 0.25, width = 0,
              alpha = 0.35, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                     guide = "none") +
  scale_y_continuous(breaks = 1:9) +
  facet_wrap(~rg_label, scales = "free_x") +
  labs(title    = "投资强度 vs 组内TRRI（原始数据）",
       subtitle = "黑线为OLS趋势线±95%CI；点有垂直抖动",
       x = "投资强度（pa_built_10y，截尾）", y = "TRRI_g（组内重编码）") +
  theme_cn()

ggsave(file.path(OUT, "raw_invest_vs_trri_g.png"),
       p_raw_inv, width = 12, height = 5, dpi = 300)
cat("-> raw_invest_vs_trri_g.png\n")

# CGI vs TRRI_g
p_raw_cgi <- plot_df_m0 %>%
  ggplot(aes(x = cgi_score_w, y = TRRI_g)) +
  geom_jitter(aes(color = response_group), height = 0.25, width = 0,
              alpha = 0.35, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                     guide = "none") +
  scale_y_continuous(breaks = 1:9) +
  facet_wrap(~rg_label, scales = "free_x") +
  labs(title    = "CGI治理指数 vs 组内TRRI（原始数据）",
       subtitle = "黑线为OLS趋势线±95%CI；点有垂直抖动",
       x = "CGI治理指数（截尾）", y = "TRRI_g（组内重编码）") +
  theme_cn()

ggsave(file.path(OUT, "raw_cgi_vs_trri_g.png"),
       p_raw_cgi, width = 12, height = 5, dpi = 300)
cat("-> raw_cgi_vs_trri_g.png\n")

# 箱线图：按TRRI_g分组，展示投资/CGI分布
p_box_inv <- plot_df_m0 %>%
  ggplot(aes(x = factor(TRRI_g), y = pa_built_10y_w, fill = response_group)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                    name = "响应类型") +
  facet_wrap(~rg_label, scales = "free") +
  labs(title = "各TRRI_g等级的投资强度分布",
       x = "TRRI_g（组内重编码）", y = "投资强度（pa_built_10y，截尾）") +
  theme_cn() + theme(legend.position = "none")

p_box_cgi <- plot_df_m0 %>%
  ggplot(aes(x = factor(TRRI_g), y = cgi_score_w, fill = response_group)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                    name = "响应类型") +
  facet_wrap(~rg_label, scales = "free") +
  labs(title = "各TRRI_g等级的CGI治理指数分布",
       x = "TRRI_g（组内重编码）", y = "CGI治理指数（截尾）") +
  theme_cn() + theme(legend.position = "none")

p_box_combined <- p_box_inv / p_box_cgi
ggsave(file.path(OUT, "raw_boxplot_trri_g.png"),
       p_box_combined, width = 12, height = 9, dpi = 300)
cat("-> raw_boxplot_trri_g.png\n")

# run_olr_full 使用 TRRI_ord 作为因变量，需临时替换为 TRRI_g_ord
# 在 M 段用专用包装函数，传入重编码因变量
run_olr_group <- function(df, inv_var, ctrl_vars, grp_label) {
  df <- df %>% mutate(TRRI = TRRI_g)  # 用组内重编码值替换，1-9统一方向
  run_olr_full(df, inv_var, ctrl_vars, grp_label)
}

# M1b 控制变量：Full Control (Footprint) — precip + pop + road + 建筑占地面积（2D）
all_vars_ctrl_all   <- c("cgi_score_w",
                         "precip_mean_w",
                         "pop_10y_w", "road_density_w", "building_footprint_w",
                         "longitude", "latitude",
                         "koppen_B", "koppen_C", "koppen_D")
all_vars_ctrl_inner <- c("cgi_score_w",
                         "precip_mean_w",
                         "pop_10y_w", "road_density_w", "building_footprint_w",
                         "longitude", "latitude")

# M1c 控制变量：Full Control (Volume) — precip + pop + road + 建筑体积密度（3D）
vol_ctrl_all   <- c("cgi_score_w",
                    "precip_mean_w",
                    "pop_10y_w", "road_density_w", "building_vol_density_w",
                    "longitude", "latitude",
                    "koppen_B", "koppen_C", "koppen_D")
vol_ctrl_inner <- c("cgi_score_w",
                    "precip_mean_w",
                    "pop_10y_w", "road_density_w", "building_vol_density_w",
                    "longitude", "latitude")

# M1. 有序Logit（基础版）：仅地理+Köppen控制，分促进/抑制组
olr_rg_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# M1b. 有序Logit（全变量版）：加入 local + bg 所有变量作为控制
olr_rg_full_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") all_vars_ctrl_all else all_vars_ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

if (nrow(olr_rg_tbl) > 0) {
  cat("\n--- M1. 有序Logit（分促进/抑制）---\n")
  print(olr_rg_tbl %>%
          filter(variable %in% c("inv_w", "cgi_score_w")) %>%
          dplyr::select(response_group, koppen, n, variable, coef, se, p_val, sig))

  write_csv(olr_rg_tbl, file.path(OUT, "olr_response_group_coef.csv"))
  cat("-> olr_response_group_coef.csv\n")

  # 森林图：促进 vs 抑制，分气候区（M1基础版）
  koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                  C = "C(温带)", D = "D(大陆)")
  m1_vars_display <- c("inv_w", "cgi_score_w",
                       "longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")
  m1_var_label_map <- c(
    inv_w       = "Investment",
    cgi_score_w = "CGI Governance",
    longitude   = "Longitude",
    latitude    = "Latitude",
    koppen_B    = "Köppen B (Arid)",
    koppen_C    = "Köppen C (Temperate)",
    koppen_D    = "Köppen D (Continental)"
  )
  m1_is_mgmt <- c(inv_w = TRUE, cgi_score_w = TRUE,
                  longitude = FALSE, latitude = FALSE,
                  koppen_B = FALSE, koppen_C = FALSE, koppen_D = FALSE)

  m1_plot_df <- olr_rg_tbl %>% filter(variable %in% m1_vars_display)
  m1_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(m1_plot_df$koppen))
  m1_grp_levels   <- koppen_nm2[m1_koppen_order]

  p_olr_rg <- m1_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = m1_grp_levels),
      var_label2 = factor(m1_var_label_map[variable],
                          levels = rev(m1_var_label_map[m1_vars_display])),
      is_mgmt    = m1_is_mgmt[variable],
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,   # filled diamond
        is_mgmt & !sig ~ 5L,    # open diamond
       !is_mgmt &  sig ~ 16L,   # filled circle
       !is_mgmt & !sig ~ 1L     # open circle
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "M1 Basic Model: Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Lon/Lat/Köppen (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group.png"),
         p_olr_rg, width = 14, height = 6, dpi = 300)
  cat("-> olr_response_group.png\n")
}

# M1b 输出（Full Control - Footprint）
mgmt_vars_display <- c("inv_w", "cgi_score_w")
all_vars_display  <- c("inv_w", "cgi_score_w",
                       "precip_mean_w",
                       "pop_10y_w", "road_density_w", "building_footprint_w",
                       "longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")
vol_vars_display  <- c("inv_w", "cgi_score_w",
                       "precip_mean_w",
                       "pop_10y_w", "road_density_w", "building_vol_density_w",
                       "longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

var_label_map <- c(
  inv_w                  = "Investment",
  cgi_score_w            = "CGI Governance",
  precip_mean_w          = "Precipitation",
  soil_clay_w            = "Soil Clay",
  soil_sand_w            = "Soil Sand",
  pop_10y_w              = "Population",
  road_density_w         = "Road Density",
  building_footprint_w   = "Building Footprint (2D)",
  building_vol_density_w = "Building Volume Density (3D)",
  longitude              = "Longitude",
  latitude               = "Latitude",
  koppen_B               = "Köppen B (Arid)",
  koppen_C               = "Köppen C (Temperate)",
  koppen_D               = "Köppen D (Continental)"
)

if (nrow(olr_rg_full_tbl) > 0) {
  cat("\n--- M1b. Full Control (Footprint): Promote vs Inhibit ---\n")
  print(olr_rg_full_tbl %>%
          filter(variable %in% mgmt_vars_display, koppen == "ALL") %>%
          dplyr::select(response_group, koppen, n, mcfadden,
                        variable, coef, se, p_val, sig))

  write_csv(olr_rg_full_tbl, file.path(OUT, "olr_response_group_fullvar_coef.csv"))
  cat("-> olr_response_group_fullvar_coef.csv\n")

  # Forest plot：全变量版，展示所有7个关键变量系数（全国ALL）
  koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                  C = "C(温带)", D = "D(大陆)")

  full_plot_df <- olr_rg_full_tbl %>% filter(variable %in% all_vars_display)
  full_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(full_plot_df$koppen))
  full_grp_levels   <- koppen_nm2[full_koppen_order]

  p_forest_full <- full_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = full_grp_levels),
      var_label2 = factor(var_label_map[variable],
                          levels = rev(var_label_map[all_vars_display])),
      is_mgmt    = variable %in% mgmt_vars_display,
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,
        is_mgmt & !sig ~ 5L,
       !is_mgmt &  sig ~ 16L,
       !is_mgmt & !sig ~ 1L
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "M1b Full Control (Footprint): Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Precip, Pop, Road, Footprint 2D, Lon/Lat/Köppen (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group_fullvar.png"),
         p_forest_full, width = 14, height = 8, dpi = 300)
  cat("-> olr_response_group_fullvar.png\n")
}

# M1c. Full Control (Volume)：与M1b相同结构，但用建筑体积密度（3D）代替占地面积（2D）
olr_rg_nopop_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") vol_ctrl_all else vol_ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

bg_vars_vol <- c("pop_10y_w", "road_density_w", "building_vol_density_w")  # M1c背景变量

vp_rg_nonpop_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    result <- run_vp3(df_g, "TRRI_g", mgmt_vars, local_vars, bg_vars_vol, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

if (nrow(olr_rg_nopop_tbl) > 0) {
  cat("\n--- M1c. Full Control (Volume): Promote vs Inhibit ---\n")
  print(olr_rg_nopop_tbl %>%
          filter(variable %in% c("inv_w", "cgi_score_w"), koppen == "ALL") %>%
          dplyr::select(response_group, koppen, n, mcfadden,
                        variable, coef, se, p_val, sig))

  write_csv(olr_rg_nopop_tbl, file.path(OUT, "olr_response_group_vol_coef.csv"))
  cat("-> olr_response_group_vol_coef.csv\n")

  koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                  C = "C(温带)", D = "D(大陆)")

  vol_plot_df <- olr_rg_nopop_tbl %>% filter(variable %in% vol_vars_display)
  vol_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(vol_plot_df$koppen))
  vol_grp_levels   <- koppen_nm2[vol_koppen_order]

  p_forest_nopop <- vol_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = vol_grp_levels),
      var_label2 = factor(var_label_map[variable],
                          levels = rev(var_label_map[vol_vars_display])),
      is_mgmt    = variable %in% mgmt_vars_display,
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,
        is_mgmt & !sig ~ 5L,
       !is_mgmt &  sig ~ 16L,
       !is_mgmt & !sig ~ 1L
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "M1c Full Control (Volume): Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Precip, Pop, Road, Vol.Density 3D, Lon/Lat/Köppen (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group_vol.png"),
         p_forest_nopop, width = 14, height = 8, dpi = 300)
  cat("-> olr_response_group_vol.png\n")
  
  # M1c varpart plot saved below in combined section
}

# M2. 方差分解：分促进/抑制组（M1b: footprint；已在上面M1c section中做了vol版本）
vp_rg_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    result <- run_vp3(df_g, "TRRI_g", mgmt_vars, local_vars, bg_vars, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

if (nrow(vp_rg_tbl) > 0) {
  cat("\n--- M2. Variance Partitioning (Promote / Inhibit, M1b footprint + M1c volume) ---\n")
  print(vp_rg_tbl %>%
          dplyr::select(response_group, koppen, n, mgmt_only, local_only,
                        bg_only, mgmt_local, mgmt_bg, local_bg, triple, unexplained))

  write_csv(vp_rg_tbl, file.path(OUT, "varpart_response_group.csv"))
  cat("-> varpart_response_group.csv\n")
  if (nrow(vp_rg_nonpop_tbl) > 0)
    write_csv(vp_rg_nonpop_tbl, file.path(OUT, "varpart_response_group_vol.csv"))

  # ── helper: build full-fraction stacked bar for one varpart table ──────────
  vp_comp_cols  <- c("mgmt_only","local_only","bg_only",
                     "mgmt_local","mgmt_bg","local_bg","triple","unexplained")
  vp_comp_labels <- c("Mgmt (unique)","Local (unique)","Background (unique)",
                      "Mgmt ∩ Local","Mgmt ∩ Background","Local ∩ Background",
                      "All Three","Unexplained")
  vp_comp_colors <- c(
    "Mgmt (unique)"       = "#D73027",
    "Local (unique)"      = "#2CA25F",
    "Background (unique)" = "#4575B4",
    "Mgmt ∩ Local"        = "#FC8D59",
    "Mgmt ∩ Background"   = "#984EA3",
    "Local ∩ Background"  = "#2CAAAA",
    "All Three"           = "#4D4D4D",
    "Unexplained"         = "#CCCCCC"
  )

  make_vp_plot <- function(vp_df, title_str, subtitle_str) {
    koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                    C = "C(温带)", D = "D(大陆)")
    ko  <- intersect(c("ALL","A","B","C","D"), unique(vp_df$koppen))
    lvl <- as.vector(outer(koppen_nm2[ko], group_labels,
                           function(k, g) paste0(k, "\n(", g, ")")))
    vp_df %>%
      mutate(grp_label = factor(
        paste0(koppen_nm2[koppen], "\n(", response_group, ")"), levels = lvl)) %>%
      filter(!is.na(grp_label)) %>%
      pivot_longer(all_of(vp_comp_cols), names_to = "component", values_to = "r2") %>%
      mutate(
        r2_show    = pmax(r2, 0),
        comp_label = factor(component, levels = vp_comp_cols, labels = vp_comp_labels)
      ) %>%
      ggplot(aes(x = grp_label, y = r2_show, fill = comp_label)) +
      geom_col(position = "stack", alpha = 0.88, width = 0.7) +
      geom_text(aes(label = ifelse(r2_show >= 0.3, sprintf("%.1f%%", r2_show), "")),
                position = position_stack(vjust = 0.5),
                size = 2.5, color = "white") +
      scale_fill_manual(values = vp_comp_colors, name = "Fraction") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(title = title_str, subtitle = subtitle_str,
           x = NULL, y = "Adj. R² (%)") +
      theme_cn() +
      theme(legend.position = "right", axis.text.x = element_text(size = 8))
  }

  p_vp_rg <- make_vp_plot(
    vp_rg_tbl,
    "Variance Partitioning – M1b Full Control (Footprint 2D): Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Mgmt = Investment + CGI; Local = Precipitation; ",
           "Background = Population + Road Density + Building Footprint (2D)\n",
           "Stacked fractions (negative shown as 0); shared fractions show collinear overlap")
  )

  ggsave(file.path(OUT, "varpart_response_group.png"),
         p_vp_rg, width = 14, height = 5, dpi = 300)
  cat("-> varpart_response_group.png\n")

  if (nrow(vp_rg_nonpop_tbl) > 0) {
    p_vp_rg_vol <- make_vp_plot(
      vp_rg_nonpop_tbl,
      "Variance Partitioning – M1c Full Control (Volume 3D): Promote vs Inhibit",
      paste0("Outcome: TRRI_g (within-group 1–9)\n",
             "Mgmt = Investment + CGI; Local = Precipitation; ",
             "Background = Population + Road Density + Building Volume Density (3D)\n",
             "Stacked fractions (negative shown as 0); shared fractions show collinear overlap")
    )
    ggsave(file.path(OUT, "varpart_response_group_vol.png"),
           p_vp_rg_vol, width = 14, height = 5, dpi = 300)
    cat("-> varpart_response_group_vol.png\n")

    # ── Combined: M1b (footprint) + M1c (volume) in one figure ───────────────
    combined_vp_df <- bind_rows(
      vp_rg_tbl       %>% mutate(model = "M1b: Full Control (Footprint 2D)"),
      vp_rg_nonpop_tbl %>% mutate(model = "M1c: Full Control (Volume 3D)")
    )
    koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                    C = "C(温带)", D = "D(大陆)")
    ko_all  <- intersect(c("ALL","A","B","C","D"), unique(combined_vp_df$koppen))
    lvl_all <- as.vector(outer(koppen_nm2[ko_all], group_labels,
                               function(k, g) paste0(k, "\n(", g, ")")))

    p_vp_combined <- combined_vp_df %>%
      mutate(
        grp_label = factor(
          paste0(koppen_nm2[koppen], "\n(", response_group, ")"), levels = lvl_all),
        model = factor(model, levels = c("M1b: Full Control (Footprint 2D)",
                                         "M1c: Full Control (Volume 3D)"))
      ) %>%
      filter(!is.na(grp_label)) %>%
      pivot_longer(all_of(vp_comp_cols), names_to = "component", values_to = "r2") %>%
      mutate(
        r2_show    = pmax(r2, 0),
        comp_label = factor(component, levels = vp_comp_cols, labels = vp_comp_labels)
      ) %>%
      ggplot(aes(x = grp_label, y = r2_show, fill = comp_label)) +
      geom_col(position = "stack", alpha = 0.88, width = 0.7) +
      geom_text(aes(label = ifelse(r2_show >= 0.5, sprintf("%.1f%%", r2_show), "")),
                position = position_stack(vjust = 0.5),
                size = 2.3, color = "white") +
      scale_fill_manual(values = vp_comp_colors, name = "Fraction") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      facet_wrap(~model, ncol = 1) +
      labs(
        title    = "Variance Partitioning: M1b vs M1c (All Fractions incl. Shared)",
        subtitle = paste0(
          "Outcome: TRRI_g (within-group 1–9)\n",
          "Mgmt = Investment + CGI; Local = Precipitation\n",
          "Background (M1b) = Population + Road + Building Footprint (2D)\n",
          "Background (M1c) = Population + Road + Building Volume Density (3D)\n",
          "Shared fractions reflect variance jointly explained by two or more variable groups"
        ),
        x = NULL, y = "Adj. R² (%)"
      ) +
      theme_cn() +
      theme(legend.position = "right", axis.text.x = element_text(size = 8))

    ggsave(file.path(OUT, "varpart_combined_models.png"),
           p_vp_combined, width = 16, height = 10, dpi = 300)
    cat("-> varpart_combined_models.png\n")
  }
}


# =============================================================================
# MX. 模型表现比较图（各OLR模型McFadden R² 对比）
# =============================================================================

cat("\n=== MX. 模型表现比较图 ===\n")

model_perf_list <- list()

# M1：Basic Model
if (exists("olr_rg_tbl") && nrow(olr_rg_tbl) > 0) {
  model_perf_list[["M1:\nBasic Model"]] <- olr_rg_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>%
    distinct() %>%
    mutate(model = "M1:\nBasic Model")
}

# M1b：Full Control (Footprint)
if (exists("olr_rg_full_tbl") && nrow(olr_rg_full_tbl) > 0) {
  model_perf_list[["M1b:\nFull Control\n(Footprint)"]] <- olr_rg_full_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>%
    distinct() %>%
    mutate(model = "M1b:\nFull Control\n(Footprint)")
}

# M1c：Full Control (Volume)
if (exists("olr_rg_nopop_tbl") && nrow(olr_rg_nopop_tbl) > 0) {
  model_perf_list[["M1c:\nFull Control\n(Volume)"]] <- olr_rg_nopop_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>%
    distinct() %>%
    mutate(model = "M1c:\nFull Control\n(Volume)")
}

if (length(model_perf_list) > 0) {
  koppen_nm_mx <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                    C = "C(温带)", D = "D(大陆)")
  model_levels <- names(model_perf_list)

  perf_df <- bind_rows(model_perf_list) %>%
    mutate(
      model     = factor(model, levels = model_levels),
      grp_label = factor(koppen_nm_mx[koppen],
                         levels = koppen_nm_mx[intersect(c("ALL","A","B","C","D"), koppen)])
    ) %>%
    filter(!is.na(grp_label))

  p_model_compare <- perf_df %>%
    ggplot(aes(x = grp_label, y = mcfadden, color = model, group = model)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.5) +
    geom_text(aes(label = sprintf("%.3f", mcfadden)),
              vjust = -0.7, size = 2.8, family = "heiti",
              position = position_dodge(0.1)) +
    scale_color_manual(
      values = c("M1:\nBasic Model"           = "#4575B4",
                 "M1b:\nFull Control\n(Footprint)" = "#D73027",
                 "M1c:\nFull Control\n(Volume)"    = "#1A9641"),
      name = "Model"
    ) +
    facet_wrap(~response_group, ncol = 2) +
    labs(
      title    = "OLR Model Comparison: McFadden Pseudo-R²",
      subtitle = paste0(
        "Outcome: TRRI_g (within-group 1–9)\n",
        "M1 (Basic Model): Investment + CGI + Lon/Lat/Köppen\n",
        "M1b (Full Control, Footprint): M1 + Precip + Pop + Road + Building Footprint (2D)\n",
        "M1c (Full Control, Volume): M1 + Precip + Pop + Road + Building Volume Density (3D)"
      ),
      x = NULL, y = "McFadden Pseudo-R²"
    ) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "model_compare_mcfadden.png"),
         p_model_compare, width = 12, height = 6, dpi = 300)
  cat("-> model_compare_mcfadden.png\n")

  # 同时输出汇总表
  perf_summary <- perf_df %>%
    dplyr::select(model, response_group, koppen, n, mcfadden) %>%
    arrange(model, response_group, koppen)
  write_csv(perf_summary, file.path(OUT, "model_compare_mcfadden.csv"))
  cat("-> model_compare_mcfadden.csv\n")
  print(perf_summary)
}

# =============================================================================
# MN. No-Geography Models（去掉经度纬度的版本）
# =============================================================================

cat("\n=== MN. No-Geography Models (without Lon/Lat) ===\n")

# MN1_nogeo 控制变量
m1ng_ctrl_all   <- c("cgi_score_w", "koppen_B", "koppen_C", "koppen_D")
m1ng_ctrl_inner <- c("cgi_score_w")

# MN1b_nogeo 控制变量
allng_ctrl_all   <- c("cgi_score_w", "precip_mean_w",
                      "pop_10y_w", "road_density_w", "building_footprint_w",
                      "koppen_B", "koppen_C", "koppen_D")
allng_ctrl_inner <- c("cgi_score_w", "precip_mean_w",
                      "pop_10y_w", "road_density_w", "building_footprint_w")

# MN1c_nogeo 控制变量
volng_ctrl_all   <- c("cgi_score_w", "precip_mean_w",
                      "pop_10y_w", "road_density_w", "building_vol_density_w",
                      "koppen_B", "koppen_C", "koppen_D")
volng_ctrl_inner <- c("cgi_score_w", "precip_mean_w",
                      "pop_10y_w", "road_density_w", "building_vol_density_w")

# ── MN1. OLR (Basic, No Geo) ─────────────────────────────────────────────────
olr_rg_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") m1ng_ctrl_all else m1ng_ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── MN1b. OLR (Full Control Footprint, No Geo) ───────────────────────────────
olr_rg_full_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") allng_ctrl_all else allng_ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── MN1c. OLR (Full Control Volume, No Geo) ──────────────────────────────────
olr_rg_vol_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") volng_ctrl_all else volng_ctrl_inner
    result <- run_olr_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── Varpart: MN1b (Footprint, No Geo) ────────────────────────────────────────
vp_rg_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    result <- run_vp3(df_g, "TRRI_g", mgmt_vars, local_vars, bg_vars, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── Varpart: MN1c (Volume, No Geo) ───────────────────────────────────────────
vp_rg_vol_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    result <- run_vp3(df_g, "TRRI_g", mgmt_vars, local_vars, bg_vars_vol, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                C = "C(温带)", D = "D(大陆)")

# Display vars for nogeo models
m1ng_vars_display  <- c("inv_w", "cgi_score_w",
                         "koppen_B", "koppen_C", "koppen_D")
allng_vars_display <- c("inv_w", "cgi_score_w",
                         "precip_mean_w",
                         "pop_10y_w", "road_density_w", "building_footprint_w",
                         "koppen_B", "koppen_C", "koppen_D")
volng_vars_display <- c("inv_w", "cgi_score_w",
                         "precip_mean_w",
                         "pop_10y_w", "road_density_w", "building_vol_density_w",
                         "koppen_B", "koppen_C", "koppen_D")

# ── Forest plot: MN1 (Basic, No Geo) ─────────────────────────────────────────
if (nrow(olr_rg_ng_tbl) > 0) {
  write_csv(olr_rg_ng_tbl, file.path(OUT, "olr_response_group_nogeo_coef.csv"))
  cat("-> olr_response_group_nogeo_coef.csv\n")

  ng_plot_df      <- olr_rg_ng_tbl %>% filter(variable %in% m1ng_vars_display)
  ng_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(ng_plot_df$koppen))
  ng_grp_levels   <- koppen_nm2[ng_koppen_order]

  m1ng_var_label_map <- c(
    inv_w       = "Investment",
    cgi_score_w = "CGI Governance",
    koppen_B    = "Köppen B (Arid)",
    koppen_C    = "Köppen C (Temperate)",
    koppen_D    = "Köppen D (Continental)"
  )
  m1ng_is_mgmt <- c(inv_w = TRUE, cgi_score_w = TRUE,
                    koppen_B = FALSE, koppen_C = FALSE, koppen_D = FALSE)

  p_olr_ng <- ng_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = ng_grp_levels),
      var_label2 = factor(m1ng_var_label_map[variable],
                          levels = rev(m1ng_var_label_map[m1ng_vars_display])),
      is_mgmt    = m1ng_is_mgmt[variable],
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,
        is_mgmt & !sig ~ 5L,
       !is_mgmt &  sig ~ 16L,
       !is_mgmt & !sig ~ 1L
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "MN1 Basic Model (No Geo): Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Köppen only — NO Lon/Lat (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group_nogeo.png"),
         p_olr_ng, width = 14, height = 6, dpi = 300)
  cat("-> olr_response_group_nogeo.png\n")
}

# ── Forest plot: MN1b (Footprint, No Geo) ────────────────────────────────────
if (nrow(olr_rg_full_ng_tbl) > 0) {
  write_csv(olr_rg_full_ng_tbl, file.path(OUT, "olr_response_group_fullvar_nogeo_coef.csv"))
  cat("-> olr_response_group_fullvar_nogeo_coef.csv\n")

  fullng_plot_df      <- olr_rg_full_ng_tbl %>% filter(variable %in% allng_vars_display)
  fullng_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(fullng_plot_df$koppen))
  fullng_grp_levels   <- koppen_nm2[fullng_koppen_order]

  p_forest_full_ng <- fullng_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = fullng_grp_levels),
      var_label2 = factor(var_label_map[variable],
                          levels = rev(var_label_map[allng_vars_display])),
      is_mgmt    = variable %in% mgmt_vars_display,
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,
        is_mgmt & !sig ~ 5L,
       !is_mgmt &  sig ~ 16L,
       !is_mgmt & !sig ~ 1L
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "MN1b Full Control Footprint (No Geo): Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Precip, Pop, Road, Footprint 2D, Köppen — NO Lon/Lat (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group_fullvar_nogeo.png"),
         p_forest_full_ng, width = 14, height = 8, dpi = 300)
  cat("-> olr_response_group_fullvar_nogeo.png\n")
}

# ── Forest plot: MN1c (Volume, No Geo) ───────────────────────────────────────
if (nrow(olr_rg_vol_ng_tbl) > 0) {
  write_csv(olr_rg_vol_ng_tbl, file.path(OUT, "olr_response_group_vol_nogeo_coef.csv"))
  cat("-> olr_response_group_vol_nogeo_coef.csv\n")

  volng_plot_df      <- olr_rg_vol_ng_tbl %>% filter(variable %in% volng_vars_display)
  volng_koppen_order <- intersect(c("ALL","A","B","C","D"), unique(volng_plot_df$koppen))
  volng_grp_levels   <- koppen_nm2[volng_koppen_order]

  p_forest_vol_ng <- volng_plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = volng_grp_levels),
      var_label2 = factor(var_label_map[variable],
                          levels = rev(var_label_map[volng_vars_display])),
      is_mgmt    = variable %in% mgmt_vars_display,
      sig        = p_val < 0.05,
      pt_shape   = case_when(
        is_mgmt &  sig ~ 18L,
        is_mgmt & !sig ~ 5L,
       !is_mgmt &  sig ~ 16L,
       !is_mgmt & !sig ~ 1L
      )
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = coef, y = var_label2, color = response_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(0.6)) +
    geom_point(aes(shape = pt_shape), size = 2.5,
               position = position_dodge(0.6)) +
    scale_color_manual(values = c("促进" = "#E6550D", "抑制" = "#3182BD"),
                       name = "Response") +
    scale_shape_identity(guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = "MN1c Full Control Volume (No Geo): Promote vs Inhibit Groups",
         subtitle = paste0("Outcome: TRRI_g (within-group 1–9, 1=worst, 9=best)\n",
                           "Management: Investment, CGI (◆ filled=sig, ◇ open=ns)\n",
                           "Controls: Precip, Pop, Road, Vol.Density 3D, Köppen — NO Lon/Lat (● filled=sig, ○ open=ns); error bars = 95% CI"),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "olr_response_group_vol_nogeo.png"),
         p_forest_vol_ng, width = 14, height = 8, dpi = 300)
  cat("-> olr_response_group_vol_nogeo.png\n")
}

# ── Varpart plots: nogeo combined ────────────────────────────────────────────
if (nrow(vp_rg_ng_tbl) > 0 && nrow(vp_rg_vol_ng_tbl) > 0) {
  write_csv(vp_rg_ng_tbl,     file.path(OUT, "varpart_response_group_nogeo.csv"))
  write_csv(vp_rg_vol_ng_tbl, file.path(OUT, "varpart_response_group_vol_nogeo.csv"))

  p_vp_ng     <- make_vp_plot(
    vp_rg_ng_tbl,
    "Variance Partitioning – MN1b Full Control Footprint (No Geo): Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Mgmt = Investment + CGI; Local = Precipitation; ",
           "Background = Population + Road Density + Building Footprint (2D)\n",
           "No Lon/Lat control; stacked fractions (negative shown as 0)")
  )
  p_vp_vol_ng <- make_vp_plot(
    vp_rg_vol_ng_tbl,
    "Variance Partitioning – MN1c Full Control Volume (No Geo): Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Mgmt = Investment + CGI; Local = Precipitation; ",
           "Background = Population + Road Density + Building Volume Density (3D)\n",
           "No Lon/Lat control; stacked fractions (negative shown as 0)")
  )

  # Combined nogeo varpart
  combined_vp_ng_df <- bind_rows(
    vp_rg_ng_tbl     %>% mutate(model = "MN1b: Full Control (Footprint 2D, No Geo)"),
    vp_rg_vol_ng_tbl %>% mutate(model = "MN1c: Full Control (Volume 3D, No Geo)")
  )
  ko_ng  <- intersect(c("ALL","A","B","C","D"), unique(combined_vp_ng_df$koppen))
  lvl_ng <- as.vector(outer(koppen_nm2[ko_ng], group_labels,
                             function(k, g) paste0(k, "\n(", g, ")")))

  p_vp_ng_combined <- combined_vp_ng_df %>%
    mutate(
      grp_label = factor(
        paste0(koppen_nm2[koppen], "\n(", response_group, ")"), levels = lvl_ng),
      model = factor(model, levels = c("MN1b: Full Control (Footprint 2D, No Geo)",
                                       "MN1c: Full Control (Volume 3D, No Geo)"))
    ) %>%
    filter(!is.na(grp_label)) %>%
    pivot_longer(all_of(vp_comp_cols), names_to = "component", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(component, levels = vp_comp_cols, labels = vp_comp_labels)
    ) %>%
    ggplot(aes(x = grp_label, y = r2_show, fill = comp_label)) +
    geom_col(position = "stack", alpha = 0.88, width = 0.7) +
    geom_text(aes(label = ifelse(r2_show >= 0.5, sprintf("%.1f%%", r2_show), "")),
              position = position_stack(vjust = 0.5),
              size = 2.3, color = "white") +
    scale_fill_manual(values = vp_comp_colors, name = "Fraction") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(~model, ncol = 1) +
    labs(
      title    = "Variance Partitioning: MN1b vs MN1c (No Geo, All Fractions incl. Shared)",
      subtitle = paste0(
        "Outcome: TRRI_g (within-group 1–9); NO Lon/Lat control\n",
        "Mgmt = Investment + CGI; Local = Precipitation\n",
        "Background (MN1b) = Population + Road + Building Footprint (2D)\n",
        "Background (MN1c) = Population + Road + Building Volume Density (3D)\n",
        "Shared fractions reflect variance jointly explained by two or more variable groups"
      ),
      x = NULL, y = "Adj. R² (%)"
    ) +
    theme_cn() +
    theme(legend.position = "right", axis.text.x = element_text(size = 8))

  ggsave(file.path(OUT, "varpart_nogeo_combined.png"),
         p_vp_ng_combined, width = 16, height = 10, dpi = 300)
  cat("-> varpart_nogeo_combined.png\n")
}

# ── McFadden comparison: all 6 models ────────────────────────────────────────
model_perf_all_list <- list()

if (exists("olr_rg_tbl") && nrow(olr_rg_tbl) > 0)
  model_perf_all_list[["M1:\nBasic\n(Geo)"]] <- olr_rg_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "M1:\nBasic\n(Geo)")

if (exists("olr_rg_ng_tbl") && nrow(olr_rg_ng_tbl) > 0)
  model_perf_all_list[["MN1:\nBasic\n(No Geo)"]] <- olr_rg_ng_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "MN1:\nBasic\n(No Geo)")

if (exists("olr_rg_full_tbl") && nrow(olr_rg_full_tbl) > 0)
  model_perf_all_list[["M1b:\nFull Ctrl\nFootprint (Geo)"]] <- olr_rg_full_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "M1b:\nFull Ctrl\nFootprint (Geo)")

if (exists("olr_rg_full_ng_tbl") && nrow(olr_rg_full_ng_tbl) > 0)
  model_perf_all_list[["MN1b:\nFull Ctrl\nFootprint (No Geo)"]] <- olr_rg_full_ng_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "MN1b:\nFull Ctrl\nFootprint (No Geo)")

if (exists("olr_rg_nopop_tbl") && nrow(olr_rg_nopop_tbl) > 0)
  model_perf_all_list[["M1c:\nFull Ctrl\nVolume (Geo)"]] <- olr_rg_nopop_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "M1c:\nFull Ctrl\nVolume (Geo)")

if (exists("olr_rg_vol_ng_tbl") && nrow(olr_rg_vol_ng_tbl) > 0)
  model_perf_all_list[["MN1c:\nFull Ctrl\nVolume (No Geo)"]] <- olr_rg_vol_ng_tbl %>%
    dplyr::select(response_group, koppen, n, mcfadden) %>% distinct() %>%
    mutate(model = "MN1c:\nFull Ctrl\nVolume (No Geo)")

if (length(model_perf_all_list) > 0) {
  model_all_levels <- names(model_perf_all_list)
  model_all_colors <- c(
    "M1:\nBasic\n(Geo)"                  = "#4575B4",
    "MN1:\nBasic\n(No Geo)"              = "#74ADD1",
    "M1b:\nFull Ctrl\nFootprint (Geo)"   = "#D73027",
    "MN1b:\nFull Ctrl\nFootprint (No Geo)" = "#F46D43",
    "M1c:\nFull Ctrl\nVolume (Geo)"      = "#1A9641",
    "MN1c:\nFull Ctrl\nVolume (No Geo)"  = "#74C476"
  )

  perf_all_df <- bind_rows(model_perf_all_list) %>%
    mutate(
      model     = factor(model, levels = model_all_levels),
      grp_label = factor(koppen_nm2[koppen],
                         levels = koppen_nm2[intersect(c("ALL","A","B","C","D"), koppen)])
    ) %>%
    filter(!is.na(grp_label))

  p_model_compare_all <- perf_all_df %>%
    ggplot(aes(x = grp_label, y = mcfadden, color = model, group = model)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.5) +
    geom_text(aes(label = sprintf("%.3f", mcfadden)),
              vjust = -0.7, size = 2.5, family = "heiti",
              position = position_dodge(0.15)) +
    scale_color_manual(values = model_all_colors, name = "Model") +
    facet_wrap(~response_group, ncol = 2) +
    labs(
      title    = "OLR Model Comparison: McFadden Pseudo-R² (Geo vs No-Geo)",
      subtitle = paste0(
        "Outcome: TRRI_g (within-group 1–9)\n",
        "Dark colors = with Lon/Lat; Light colors = without Lon/Lat\n",
        "M1/MN1: Basic; M1b/MN1b: + Footprint; M1c/MN1c: + Volume Density"
      ),
      x = NULL, y = "McFadden Pseudo-R²"
    ) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "model_compare_mcfadden_all.png"),
         p_model_compare_all, width = 14, height = 6, dpi = 300)
  cat("-> model_compare_mcfadden_all.png\n")

  write_csv(perf_all_df %>%
              dplyr::select(model, response_group, koppen, n, mcfadden) %>%
              arrange(model, response_group, koppen),
            file.path(OUT, "model_compare_mcfadden_all.csv"))
  cat("-> model_compare_mcfadden_all.csv\n")
}

# =============================================================================
# MR. 随机森林分析（补充有序Logit，非参数框架下验证各变量的预测力）
# =============================================================================

cat("\n=== MR. Random Forest Analysis ===\n")

if (!requireNamespace("ranger", quietly = TRUE))
  install.packages("ranger", repos = "https://cloud.r-project.org")
library(ranger)

# ── Helper: run RF for one df/variable set ────────────────────────────────────
run_rf <- function(df, inv_var, ctrl_vars, grp_label, num_trees = 500) {
  all_vars <- c(inv_var, ctrl_vars)
  df2 <- df %>%
    dplyr::select(TRRI_g, all_of(all_vars)) %>%
    drop_na()
  n <- nrow(df2)
  if (n < 20) {
    message(sprintf("  [RF skip] %s: n=%d < 20", grp_label, n))
    return(NULL)
  }
  df2$TRRI_g <- as.numeric(df2$TRRI_g)
  tryCatch({
    fit <- ranger(
      formula   = TRRI_g ~ .,
      data      = df2,
      num.trees = num_trees,
      importance = "permutation",
      seed      = 42
    )
    imp <- importance(fit)
    # Direction via Spearman correlation
    cors <- sapply(names(imp), function(v) {
      cor(df2[[v]], df2$TRRI_g, method = "spearman", use = "complete.obs")
    })
    tibble(
      koppen    = grp_label,
      n         = n,
      oob_r2    = fit$r.squared,
      variable  = names(imp),
      importance = as.numeric(imp),
      direction  = sign(cors)
    )
  }, error = function(e) {
    message(sprintf("  [RF error] %s: %s", grp_label, conditionMessage(e)))
    NULL
  })
}

run_rf_group <- function(df, inv_var, ctrl_vars, grp_label) {
  df <- df %>% mutate(TRRI_g = TRRI_g)
  run_rf(df, inv_var, ctrl_vars, grp_label)
}

# ── MR1. RF Basic Model ───────────────────────────────────────────────────────
rf_rg_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── MR1b. RF Full Control (Footprint) ────────────────────────────────────────
rf_rg_full_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") all_vars_ctrl_all else all_vars_ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── MR1c. RF Full Control (Volume) ───────────────────────────────────────────
rf_rg_vol_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") vol_ctrl_all else vol_ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── MR1_nogeo / MR1b_nogeo / MR1c_nogeo ──────────────────────────────────────
rf_rg_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") m1ng_ctrl_all else m1ng_ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

rf_rg_full_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") allng_ctrl_all else allng_ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

rf_rg_vol_ng_tbl <- map_dfr(group_labels, function(rg) {
  df_rg <- filter(anal_df, response_group == rg)
  map_dfr(c("ALL", sort(unique(df_rg$koppen_group))), function(grp) {
    df_g   <- if (grp == "ALL") df_rg else filter(df_rg, koppen_group == grp)
    ctrl_v <- if (grp == "ALL") volng_ctrl_all else volng_ctrl_inner
    result <- run_rf_group(df_g, "pa_built_10y_w", ctrl_v, grp)
    if (!is.null(result)) result %>% mutate(response_group = rg, .before = 1)
  })
})

# ── Helper: RF importance bar plot ───────────────────────────────────────────
make_rf_plot <- function(rf_tbl, vars_display, title_str, subtitle_str,
                          mgmt_vars = c("inv_w", "cgi_score_w"),
                          label_map = var_label_map) {
  koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                  C = "C(温带)", D = "D(大陆)")

  plot_df     <- rf_tbl %>% filter(variable %in% vars_display)
  ko_order    <- intersect(c("ALL","A","B","C","D"), unique(plot_df$koppen))
  grp_levels  <- koppen_nm2[ko_order]

  plot_df %>%
    mutate(
      grp_label  = factor(koppen_nm2[koppen], levels = grp_levels),
      var_label2 = factor(label_map[variable],
                          levels = rev(label_map[vars_display])),
      is_mgmt    = variable %in% mgmt_vars,
      dir_label  = ifelse(direction >= 0, "Positive", "Negative"),
      # importance can be negative (permutation); clip at 0 for display length,
      # but keep value for color
      imp_show   = importance,
      pt_shape   = ifelse(is_mgmt, 18L, 16L)
    ) %>%
    filter(!is.na(grp_label), !is.na(var_label2)) %>%
    ggplot(aes(x = imp_show, y = var_label2,
               fill = interaction(dir_label, response_group),
               alpha = is_mgmt)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(
      values = c(
        "Positive.促进" = "#E6550D", "Negative.促进" = "#FDAE6B",
        "Positive.抑制" = "#3182BD", "Negative.抑制" = "#9ECAE1"
      ),
      labels = c(
        "Positive.促进" = "Promote (+)",
        "Negative.促进" = "Promote (−)",
        "Positive.抑制" = "Inhibit (+)",
        "Negative.抑制" = "Inhibit (−)"
      ),
      name = "Group & Direction"
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.65),
                       guide = "none") +
    facet_wrap(~grp_label, nrow = 1) +
    labs(title    = title_str,
         subtitle = subtitle_str,
         x = "Permutation Importance (% ↑ MSE)", y = NULL) +
    theme_cn() +
    theme(legend.position = "right")
}

# ── Plots ─────────────────────────────────────────────────────────────────────
# Extend var_label_map for M1 basic (which uses m1_var_label_map)
rf_label_map <- c(
  inv_w                  = "Investment",
  cgi_score_w            = "CGI Governance",
  precip_mean_w          = "Precipitation",
  pop_10y_w              = "Population",
  road_density_w         = "Road Density",
  building_footprint_w   = "Building Footprint (2D)",
  building_vol_density_w = "Building Volume Density (3D)",
  longitude              = "Longitude",
  latitude               = "Latitude",
  koppen_B               = "Köppen B (Arid)",
  koppen_C               = "Köppen C (Temperate)",
  koppen_D               = "Köppen D (Continental)"
)

if (nrow(rf_rg_tbl) > 0) {
  write_csv(rf_rg_tbl, file.path(OUT, "rf_response_group_coef.csv"))
  p_rf1 <- make_rf_plot(
    rf_rg_tbl, m1_vars_display,
    "MR1 Basic Model (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Importance: permutation (% increase in OOB MSE); ",
           "Direction: sign of Spearman correlation\n",
           "Full bars = management variables; translucent = controls; color saturation = direction"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group.png"),
         p_rf1, width = 14, height = 6, dpi = 300)
  cat("-> rf_response_group.png\n")
}

if (nrow(rf_rg_full_tbl) > 0) {
  write_csv(rf_rg_full_tbl, file.path(OUT, "rf_response_group_fullvar_coef.csv"))
  p_rf1b <- make_rf_plot(
    rf_rg_full_tbl, all_vars_display,
    "MR1b Full Control Footprint (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Importance: permutation (% increase in OOB MSE); ",
           "Direction: sign of Spearman correlation\n",
           "Full bars = management variables; translucent = controls"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group_fullvar.png"),
         p_rf1b, width = 14, height = 8, dpi = 300)
  cat("-> rf_response_group_fullvar.png\n")
}

if (nrow(rf_rg_vol_tbl) > 0) {
  write_csv(rf_rg_vol_tbl, file.path(OUT, "rf_response_group_vol_coef.csv"))
  p_rf1c <- make_rf_plot(
    rf_rg_vol_tbl, vol_vars_display,
    "MR1c Full Control Volume (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9)\n",
           "Importance: permutation (% increase in OOB MSE); ",
           "Direction: sign of Spearman correlation\n",
           "Full bars = management variables; translucent = controls"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group_vol.png"),
         p_rf1c, width = 14, height = 8, dpi = 300)
  cat("-> rf_response_group_vol.png\n")
}

if (nrow(rf_rg_ng_tbl) > 0) {
  write_csv(rf_rg_ng_tbl, file.path(OUT, "rf_response_group_nogeo_coef.csv"))
  p_rf_ng1 <- make_rf_plot(
    rf_rg_ng_tbl, m1ng_vars_display,
    "MR-N1 Basic Model No Geo (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9); NO Lon/Lat\n",
           "Importance: permutation (% increase in OOB MSE); Direction: sign of Spearman correlation"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group_nogeo.png"),
         p_rf_ng1, width = 14, height = 6, dpi = 300)
  cat("-> rf_response_group_nogeo.png\n")
}

if (nrow(rf_rg_full_ng_tbl) > 0) {
  write_csv(rf_rg_full_ng_tbl, file.path(OUT, "rf_response_group_fullvar_nogeo_coef.csv"))
  p_rf_ng1b <- make_rf_plot(
    rf_rg_full_ng_tbl, allng_vars_display,
    "MR-N1b Full Footprint No Geo (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9); NO Lon/Lat\n",
           "Importance: permutation (% increase in OOB MSE); Direction: sign of Spearman correlation"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group_fullvar_nogeo.png"),
         p_rf_ng1b, width = 14, height = 8, dpi = 300)
  cat("-> rf_response_group_fullvar_nogeo.png\n")
}

if (nrow(rf_rg_vol_ng_tbl) > 0) {
  write_csv(rf_rg_vol_ng_tbl, file.path(OUT, "rf_response_group_vol_nogeo_coef.csv"))
  p_rf_ng1c <- make_rf_plot(
    rf_rg_vol_ng_tbl, volng_vars_display,
    "MR-N1c Full Volume No Geo (RF): Variable Importance – Promote vs Inhibit",
    paste0("Outcome: TRRI_g (within-group 1–9); NO Lon/Lat\n",
           "Importance: permutation (% increase in OOB MSE); Direction: sign of Spearman correlation"),
    label_map = rf_label_map
  )
  ggsave(file.path(OUT, "rf_response_group_vol_nogeo.png"),
         p_rf_ng1c, width = 14, height = 8, dpi = 300)
  cat("-> rf_response_group_vol_nogeo.png\n")
}

# ── OOB R² comparison across all RF models ───────────────────────────────────
rf_all_perf <- bind_rows(
  if (nrow(rf_rg_tbl)        > 0) rf_rg_tbl        %>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR1:\nBasic (Geo)"),
  if (nrow(rf_rg_ng_tbl)     > 0) rf_rg_ng_tbl     %>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR-N1:\nBasic (No Geo)"),
  if (nrow(rf_rg_full_tbl)   > 0) rf_rg_full_tbl   %>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR1b:\nFootprint (Geo)"),
  if (nrow(rf_rg_full_ng_tbl)> 0) rf_rg_full_ng_tbl%>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR-N1b:\nFootprint (No Geo)"),
  if (nrow(rf_rg_vol_tbl)    > 0) rf_rg_vol_tbl    %>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR1c:\nVolume (Geo)"),
  if (nrow(rf_rg_vol_ng_tbl) > 0) rf_rg_vol_ng_tbl %>% dplyr::select(response_group, koppen, n, oob_r2) %>% distinct() %>% mutate(model = "MR-N1c:\nVolume (No Geo)")
)

if (nrow(rf_all_perf) > 0) {
  koppen_nm2 <- c(ALL = "全部", A = "A(热带)", B = "B(干旱)",
                  C = "C(温带)", D = "D(大陆)")
  rf_model_levels  <- unique(rf_all_perf$model)
  rf_model_colors  <- c(
    "MR1:\nBasic (Geo)"          = "#4575B4",
    "MR-N1:\nBasic (No Geo)"     = "#74ADD1",
    "MR1b:\nFootprint (Geo)"     = "#D73027",
    "MR-N1b:\nFootprint (No Geo)"= "#F46D43",
    "MR1c:\nVolume (Geo)"        = "#1A9641",
    "MR-N1c:\nVolume (No Geo)"   = "#74C476"
  )

  p_rf_compare <- rf_all_perf %>%
    mutate(
      model     = factor(model, levels = rf_model_levels),
      grp_label = factor(koppen_nm2[koppen],
                         levels = koppen_nm2[intersect(c("ALL","A","B","C","D"), koppen)])
    ) %>%
    filter(!is.na(grp_label)) %>%
    ggplot(aes(x = grp_label, y = oob_r2, color = model, group = model)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.5) +
    geom_text(aes(label = sprintf("%.3f", oob_r2)),
              vjust = -0.7, size = 2.5, family = "heiti",
              position = position_dodge(0.15)) +
    scale_color_manual(values = rf_model_colors, name = "Model") +
    facet_wrap(~response_group, ncol = 2) +
    labs(
      title    = "RF Model Comparison: OOB R² (Geo vs No-Geo)",
      subtitle = paste0(
        "Outcome: TRRI_g (within-group 1–9)\n",
        "Dark = with Lon/Lat; Light = without Lon/Lat; ",
        "OOB R² from ranger permutation forest (500 trees)"
      ),
      x = NULL, y = "OOB R²"
    ) +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "rf_model_compare_oobr2.png"),
         p_rf_compare, width = 14, height = 6, dpi = 300)
  cat("-> rf_model_compare_oobr2.png\n")

  write_csv(rf_all_perf %>% arrange(model, response_group, koppen),
            file.path(OUT, "rf_model_compare_oobr2.csv"))
  cat("-> rf_model_compare_oobr2.csv\n")
}

# =============================================================================
# 完成
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("05_1_analysis.R 运行完成\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("输出目录  : %s\n", OUT))
cat(sprintf("CCM站点数 : %d\n", length(unique(ccm_results$meteo_stat_id))))
cat(sprintf("分析站点数: %d\n", nrow(anal_df)))
cat(sprintf("分析城市数: %d（仅用于G3热图）\n", nrow(city_agg)))
cat("\n输出文件:\n")
for (f in sort(list.files(OUT, pattern = "\\.(csv|png)$")))
  cat(sprintf("  %s\n", f))
