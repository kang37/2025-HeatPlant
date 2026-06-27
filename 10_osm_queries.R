# 10_osm_queries.R
# 逐站点查询 Overpass API（httr + jsonlite，绕开 osmdata 重试逻辑）：
#   - 道路密度（10km缓冲区，km/km²）
#   - 建筑密度（500m缓冲区，建筑占地/缓冲区面积）
# 结果缓存到 data_raw/road_density_station.rds，供 05_1_analysis.R 读取

library(tidyverse)
library(sf)
library(httr)
library(jsonlite)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")

OVERPASS_URL <- "https://maps.mail.ru/osm/tools/overpass/api/interpreter"
ROAD_CACHE   <- "data_raw/road_density_station.rds"

# 直接 HTTP POST 查询，返回原始 JSON 字符串
overpass_raw <- function(ql, timeout_s = 40) {
  resp <- tryCatch(
    httr::POST(OVERPASS_URL, body = list(data = ql), encode = "form",
               httr::timeout(timeout_s)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
  httr::content(resp, as = "text", encoding = "UTF-8")
}

# 解析 Overpass JSON (out geom) → sf LINESTRING（用于道路）
parse_lines <- function(json_text) {
  if (is.null(json_text)) return(NULL)
  dat <- tryCatch(jsonlite::fromJSON(json_text, simplifyVector = FALSE),
                  error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  els <- dat$elements
  ways <- Filter(function(e) e$type == "way" && !is.null(e$geometry), els)
  if (length(ways) == 0) return(NULL)
  geoms <- lapply(ways, function(w) {
    coords <- tryCatch(
      do.call(rbind, lapply(w$geometry, function(g) c(g$lon, g$lat))),
      error = function(e) NULL
    )
    if (is.null(coords) || nrow(coords) < 2) return(NULL)
    sf::st_linestring(coords)
  })
  geoms <- Filter(Negate(is.null), geoms)
  if (length(geoms) == 0) return(NULL)
  sf::st_sfc(geoms, crs = 4326)
}

# 解析 Overpass JSON (out geom) → sf POLYGON（用于建筑）
parse_polygons <- function(json_text) {
  if (is.null(json_text)) return(NULL)
  dat <- tryCatch(jsonlite::fromJSON(json_text, simplifyVector = FALSE),
                  error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  els <- dat$elements
  ways <- Filter(function(e) e$type == "way" && !is.null(e$geometry), els)
  if (length(ways) == 0) return(NULL)
  geoms <- lapply(ways, function(w) {
    coords <- tryCatch(
      do.call(rbind, lapply(w$geometry, function(g) c(g$lon, g$lat))),
      error = function(e) NULL
    )
    if (is.null(coords) || nrow(coords) < 3) return(NULL)
    # 确保首尾闭合
    if (!all(coords[1, ] == coords[nrow(coords), ]))
      coords <- rbind(coords, coords[1, ])
    tryCatch(sf::st_polygon(list(coords)), error = function(e) NULL)
  })
  geoms <- Filter(Negate(is.null), geoms)
  if (length(geoms) == 0) return(NULL)
  sf::st_sfc(geoms, crs = 4326)
}

# 读取站点坐标
station_coords <- read_csv("data_proc/output_10y_built_up_05_01/trri_station.csv",
                           show_col_types = FALSE) %>%
  dplyr::select(meteo_stat_id, longitude, latitude) %>%
  filter(!is.na(longitude), !is.na(latitude))

cat(sprintf("共 %d 个站点\n", nrow(station_coords)))

main_road_types <- paste(
  c("motorway","trunk","primary","secondary","tertiary",
    "unclassified","residential","motorway_link","trunk_link",
    "primary_link","secondary_link","tertiary_link"),
  collapse = "|"
)

road_buf_m     <- 10000
road_buf_area  <- pi * 10^2     # ≈ 314.16 km²
build_buf_m    <- 500
build_buf_area <- pi * 0.5^2   # ≈ 0.785 km²

# 断点续传
if (file.exists(ROAD_CACHE)) {
  done <- readRDS(ROAD_CACHE)
  if (!"building_density" %in% names(done)) done$building_density <- NA_real_
  cat(sprintf("读取已有缓存：%d 个站点已完成\n", nrow(done)))
  remaining <- station_coords %>% filter(!meteo_stat_id %in% done$meteo_stat_id)
} else {
  done      <- tibble()
  remaining <- station_coords
}

cat(sprintf("剩余待查询：%d 个站点\n", nrow(remaining)))

results <- vector("list", nrow(remaining))

for (i in seq_len(nrow(remaining))) {
  row <- remaining[i, ]
  lon <- row$longitude; lat <- row$latitude
  cat(sprintf("[%d/%d] 站点 %s\n", i, nrow(remaining), row$meteo_stat_id))

  zone    <- (floor((lon + 180) / 6) %% 60) + 1
  utm_crs <- paste0("+proj=utm +zone=", zone, " +ellps=WGS84")
  center_sf <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326) %>%
    sf::st_transform(utm_crs)

  # --- 道路密度（10km）---
  d <- 0.09
  ql_road <- sprintf(
    '[out:json][timeout:35][bbox:%f,%f,%f,%f];(way["highway"~"^(%s)$"];);out geom;',
    lat - d, lon - d, lat + d, lon + d, main_road_types
  )
  road_density <- tryCatch({
    lines <- parse_lines(overpass_raw(ql_road))
    if (is.null(lines)) { 0 } else {
      buf        <- sf::st_buffer(center_sf, dist = road_buf_m)
      lines_proj <- sf::st_transform(lines, utm_crs)
      clipped    <- suppressWarnings(sf::st_intersection(lines_proj, buf))
      if (length(clipped) == 0) 0
      else as.numeric(sum(sf::st_length(clipped), na.rm = TRUE) / 1000) / road_buf_area
    }
  }, error = function(e) { cat(sprintf("  [道路失败: %s]\n", e$message)); NA_real_ })

  # --- 建筑密度（500m）---
  d2 <- 0.005
  ql_build <- sprintf(
    '[out:json][timeout:35][bbox:%f,%f,%f,%f];(way["building"];);out geom;',
    lat - d2, lon - d2, lat + d2, lon + d2
  )
  building_density <- tryCatch({
    polys <- parse_polygons(overpass_raw(ql_build))
    if (is.null(polys)) { 0 } else {
      buf        <- sf::st_buffer(center_sf, dist = build_buf_m)
      polys_proj <- sf::st_transform(polys, utm_crs)
      clipped    <- suppressWarnings(sf::st_intersection(polys_proj, buf))
      if (length(clipped) == 0) 0
      else as.numeric(sum(sf::st_area(clipped), na.rm = TRUE) / 1e6) / build_buf_area
    }
  }, error = function(e) { cat(sprintf("  [建筑失败: %s]\n", e$message)); NA_real_ })

  cat(sprintf("  道路: %.4f km/km²  建筑: %.4f\n",
              ifelse(is.na(road_density), NaN, road_density),
              ifelse(is.na(building_density), NaN, building_density)))

  results[[i]] <- tibble(
    meteo_stat_id    = row$meteo_stat_id,
    road_density     = road_density,
    building_density = building_density
  )

  Sys.sleep(1)

  if (i %% 10 == 0 || i == nrow(remaining)) {
    current <- bind_rows(done, bind_rows(results[1:i]))
    saveRDS(current, ROAD_CACHE)
    cat(sprintf("  -> 进度已保存（共 %d 站点）\n", nrow(current)))
  }
}

road_df <- bind_rows(done, bind_rows(results))
saveRDS(road_df, ROAD_CACHE)
cat(sprintf("\n完成！道路: %d / %d 有值，建筑: %d / %d 有值\n",
            sum(!is.na(road_df$road_density)),     nrow(road_df),
            sum(!is.na(road_df$building_density)), nrow(road_df)))
