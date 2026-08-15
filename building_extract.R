# ============================================================
# 计算空气质量站点 30 m 缓冲区内建筑物占地面积
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(lwgeom)
  library(foreach)
  library(doParallel)
})

sf_use_s2(TRUE)

# ---------------------------- 路径与参数 ----------------------------

province_shp_path <- paste0(
  "/share/home/wangl/buiding/",
  "中国标准行政区划数据GS（2024）0650号/",
  "shp格式/中国_省.shp"
)

station_csv_path <- "/share/home/wangl/buiding/station_regression.csv"
building_dir <- "/share/home/wangl/china_new"
result_dir <- file.path(building_dir, "station_results_30m")
final_output <- "/share/home/wangl/China_stations_buildings_with_coords30.csv"

buffer_m <- 30
workers_requested <- 12

# TRUE：全部重算；FALSE：已有省份结果则跳过
overwrite_existing <- TRUE

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(province_shp_path),
  file.exists(station_csv_path),
  dir.exists(building_dir)
)

# ---------------------------- 省界数据 ----------------------------

message("读取省界：", province_shp_path)

china_provinces <- st_read(province_shp_path, quiet = TRUE)

if (!inherits(china_provinces, "sf")) stop("省界未读取为空间对象。")
if (is.na(st_crs(china_provinces))) stop("省界缺少 CRS。")
if (!"name" %in% names(china_provinces)) {
  stop("省界中没有 name 字段。现有字段：", paste(names(china_provinces), collapse = ", "))
}

china_provinces <- st_zm(china_provinces, drop = TRUE, what = "ZM")

bad_province_geometry <- !st_is_valid(china_provinces)
if (any(bad_province_geometry)) {
  message("发现 ", sum(bad_province_geometry), " 个无效省界几何，正在修复。")
  china_provinces <- st_make_valid(china_provinces)
}

china_provinces <- st_transform(china_provinces, 4326)

province_name_map <- c(
  "河南" = "Henan", "浙江" = "Zhejiang", "海南" = "Hainan",
  "台湾" = "Taiwan", "甘肃" = "Gansu", "湖北" = "Hubei",
  "天津" = "Tianjin", "内蒙古" = "NeiMongol", "广西" = "Guangxi",
  "江苏" = "Jiangsu", "山东" = "Shandong", "北京" = "Beijing",
  "西藏" = "Tibet", "青海" = "Qinghai", "重庆" = "Chongqing",
  "澳门" = "Macau", "吉林" = "Jilin", "福建" = "Fujian",
  "广东" = "Guangdong", "上海" = "Shanghai", "江西" = "Jiangxi",
  "香港" = "HongKong", "黑龙江" = "Heilongjiang", "四川" = "Sichuan",
  "云南" = "Yunnan", "陕西" = "Shaanxi", "贵州" = "Guizhou",
  "安徽" = "Anhui", "山西" = "Shanxi", "河北" = "Hebei",
  "辽宁" = "Liaoning", "宁夏" = "Ningxia", "湖南" = "Hunan",
  "新疆" = "Xinjiang"
)

china_provinces <- china_provinces %>%
  mutate(
    clean_name = str_remove(
      as.character(name),
      "(壮族自治区|回族自治区|维吾尔自治区|特别行政区|自治区|省|市)$"
    ),
    english_name = unname(province_name_map[clean_name])
  )

unmapped_provinces <- china_provinces %>%
  st_drop_geometry() %>%
  filter(is.na(english_name) | english_name == "") %>%
  distinct(name, clean_name)

if (nrow(unmapped_provinces) > 0) {
  warning("以下省份未匹配英文名：\n", paste(unmapped_provinces$name, collapse = "\n"))
}

# ---------------------------- 建筑物 GPKG 清单 ----------------------------

gpkg_files <- list.files(
  building_dir,
  pattern = "_buildings_CRS_.*\\.gpkg$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(gpkg_files) == 0) {
  stop("没有找到建筑物 GPKG。目录：", building_dir)
}

message("找到 ", length(gpkg_files), " 个建筑物 GPKG。")

building_files <- tibble(
  gpkg_file = gpkg_files,
  file_name = basename(gpkg_files),
  file_size_GB = file.info(gpkg_files)$size / 1024^3
) %>%
  mutate(
    Province = str_match(
      file_name,
      regex("^(.+?)_buildings_CRS_(.+)\\.gpkg$", ignore_case = TRUE)
    )[, 2],
    CRS_Group = str_match(
      file_name,
      regex("^(.+?)_buildings_CRS_(.+)\\.gpkg$", ignore_case = TRUE)
    )[, 3]
  ) %>%
  filter(!is.na(Province), !is.na(CRS_Group))

crs_valid <- function(x) {
  inherits(x, "crs") && tryCatch(!is.na(x), error = function(e) FALSE)
}

inspect_gpkg <- function(gpkg_file, file_name, file_size_GB, Province, CRS_Group) {
  tryCatch({
    info <- st_layers(gpkg_file)
    if (length(info$name) == 0) stop("GPKG 中没有图层。")

    layer_names <- info$name
    name_match <- which(str_detect(str_to_lower(layer_names), "building"))

    crs_match <- which(vapply(seq_along(layer_names), function(i) {
      tryCatch(crs_valid(info$crs[[i]]), error = function(e) FALSE)
    }, logical(1)))

    selected_index <- unique(c(name_match, crs_match, seq_along(layer_names)))[1]
    selected_layer <- layer_names[selected_index]
    selected_crs <- tryCatch(info$crs[[selected_index]], error = function(e) NULL)

    if (!crs_valid(selected_crs)) {
      epsg_from_name <- suppressWarnings(
        as.integer(str_extract(as.character(CRS_Group), "\\d+"))
      )
      if (!is.na(epsg_from_name)) selected_crs <- st_crs(epsg_from_name)
    }

    if (!crs_valid(selected_crs)) stop("无法确定图层 CRS。")

    geom_text <- tryCatch(
      paste(as.character(info$geomtype[[selected_index]]), collapse = ", "),
      error = function(e) ""
    )

    tibble(
      gpkg_file = gpkg_file,
      file_name = file_name,
      file_size_GB = file_size_GB,
      Province = Province,
      CRS_Group = as.character(CRS_Group),
      layer_name = selected_layer,
      geometry_type = geom_text,
      crs_wkt = selected_crs$wkt,
      epsg = selected_crs$epsg,
      valid_file = TRUE,
      inspection_error = NA_character_
    )
  }, error = function(e) {
    tibble(
      gpkg_file = gpkg_file,
      file_name = file_name,
      file_size_GB = file_size_GB,
      Province = Province,
      CRS_Group = as.character(CRS_Group),
      layer_name = NA_character_,
      geometry_type = NA_character_,
      crs_wkt = NA_character_,
      epsg = NA_integer_,
      valid_file = FALSE,
      inspection_error = conditionMessage(e)
    )
  })
}

message("读取建筑物 GPKG 图层元数据。")

building_inventory <- pmap_dfr(
  building_files,
  inspect_gpkg
)

write_csv(
  building_inventory,
  file.path(result_dir, "Building_GPKG_Inventory.csv")
)

invalid_gpkg <- building_inventory %>% filter(!valid_file)
if (nrow(invalid_gpkg) > 0) {
  write_csv(invalid_gpkg, file.path(result_dir, "Invalid_Building_GPKG_Files.csv"))
  warning(nrow(invalid_gpkg), " 个 GPKG 检查失败。")
}

building_inventory <- building_inventory %>% filter(valid_file)
if (nrow(building_inventory) == 0) stop("没有可使用的建筑物 GPKG。")

message("有效建筑物 GPKG：", nrow(building_inventory))
print(
  building_inventory %>%
    select(Province, file_name, layer_name, geometry_type, epsg, file_size_GB),
  n = Inf
)

# ---------------------------- 站点数据 ----------------------------

station_data <- read_csv(
  station_csv_path,
  show_col_types = FALSE,
  name_repair = "unique"
)

message("站点文件字段：", paste(names(station_data), collapse = ", "))

required_station_columns <- c("meteo_stat_id", "longitude", "latitude")
missing_station_columns <- setdiff(required_station_columns, names(station_data))

if (length(missing_station_columns) > 0) {
  stop(
    "站点文件缺少字段：",
    paste(missing_station_columns, collapse = ", "),
    "\n实际字段：",
    paste(names(station_data), collapse = ", ")
  )
}

# 删除无名索引列，如 ...1
station_data <- station_data %>%
  select(-matches("^\\.\\.\\.[0-9]+$")) %>%
  rename(
    AirQualityStation = meteo_stat_id,
    Longitude = longitude,
    Latitude = latitude
  ) %>%
  mutate(
    station_row_id = row_number(),
    AirQualityStation = str_trim(as.character(AirQualityStation)),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude = suppressWarnings(as.numeric(Latitude))
  )

invalid_stations <- station_data %>%
  filter(
    is.na(AirQualityStation) | AirQualityStation == "" |
      is.na(Longitude) | is.na(Latitude) |
      !between(Longitude, -180, 180) |
      !between(Latitude, -90, 90)
  )

if (nrow(invalid_stations) > 0) {
  write_csv(
    invalid_stations,
    file.path(result_dir, "Stations_with_Invalid_Coordinates.csv")
  )
  warning(nrow(invalid_stations), " 个站点编号或坐标无效，已排除。")
}

station_data <- station_data %>%
  filter(
    !is.na(AirQualityStation),
    AirQualityStation != "",
    !is.na(Longitude),
    !is.na(Latitude),
    between(Longitude, -180, 180),
    between(Latitude, -90, 90)
  )

if (nrow(station_data) == 0) stop("清理后没有有效站点。")
message("有效站点数：", nrow(station_data))

duplicate_ids <- station_data %>%
  count(AirQualityStation, name = "record_number") %>%
  filter(record_number > 1)

if (nrow(duplicate_ids) > 0) {
  write_csv(duplicate_ids, file.path(result_dir, "Duplicate_Station_IDs.csv"))
}

stations_sf <- st_as_sf(
  station_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# ---------------------------- 匹配站点所属省份 ----------------------------

stations_with_province <- st_join(
  stations_sf,
  china_provinces %>% select(clean_name, english_name),
  join = st_intersects,
  left = TRUE
) %>%
  arrange(station_row_id, is.na(english_name)) %>%
  group_by(station_row_id) %>%
  slice(1) %>%
  ungroup()

unmatched_index <- which(
  is.na(stations_with_province$english_name) |
    stations_with_province$english_name == ""
)

if (length(unmatched_index) > 0) {
  message("有 ", length(unmatched_index), " 个站点未落入省界，匹配最近省份。")

  nearest_index <- st_nearest_feature(
    stations_with_province[unmatched_index, ],
    china_provinces
  )

  stations_with_province$clean_name[unmatched_index] <-
    china_provinces$clean_name[nearest_index]

  stations_with_province$english_name[unmatched_index] <-
    china_provinces$english_name[nearest_index]
}

still_unmatched <- stations_with_province %>%
  filter(is.na(english_name) | english_name == "")

if (nrow(still_unmatched) > 0) {
  write_csv(
    st_drop_geometry(still_unmatched),
    file.path(result_dir, "Stations_without_Province.csv")
  )
  warning(nrow(still_unmatched), " 个站点仍无法匹配省份，已排除。")
}

stations_with_province <- stations_with_province %>%
  filter(!is.na(english_name), english_name != "")

write_csv(
  st_drop_geometry(stations_with_province),
  file.path(result_dir, "Stations_with_Provinces.csv")
)

province_list <- sort(unique(stations_with_province$english_name))
building_provinces <- sort(unique(building_inventory$Province))

provinces_without_gpkg <- setdiff(province_list, building_provinces)

if (length(provinces_without_gpkg) > 0) {
  write_csv(
    tibble(Province = provinces_without_gpkg),
    file.path(result_dir, "Provinces_without_Building_GPKG.csv")
  )
  warning(
    "以下省份没有建筑物 GPKG，面积将保存为 NA：\n",
    paste(provinces_without_gpkg, collapse = "\n")
  )
}

# ---------------------------- 空间计算函数 ----------------------------

get_local_utm_epsg <- function(longitude, latitude) {
  zone <- floor((longitude + 180) / 6) + 1
  zone <- max(1, min(60, zone))
  if (latitude >= 0) 32600 + zone else 32700 + zone
}

clean_polygon_sf <- function(x) {
  if (is.null(x) || !inherits(x, "sf") || nrow(x) == 0) return(x)

  x <- st_zm(x, drop = TRUE, what = "ZM")

  if (any(!st_is_valid(x))) {
    x <- st_make_valid(x)
  }

  geom_type <- as.character(st_geometry_type(x))

  if (any(geom_type == "GEOMETRYCOLLECTION")) {
    x <- suppressWarnings(st_collection_extract(x, "POLYGON"))
  }

  geom_type <- as.character(st_geometry_type(x))
  x[geom_type %in% c("POLYGON", "MULTIPOLYGON"), , drop = FALSE]
}

calculate_area_in_one_gpkg <- function(
    station_row,
    gpkg_file,
    layer_name,
    layer_crs_wkt,
    buffer_m
) {
  longitude <- as.numeric(station_row$Longitude[[1]])
  latitude <- as.numeric(station_row$Latitude[[1]])
  local_utm <- get_local_utm_epsg(longitude, latitude)

  station_utm <- st_transform(station_row, local_utm)
  buffer_utm <- st_buffer(st_geometry(station_utm), dist = buffer_m)
  buffer_utm <- st_make_valid(buffer_utm)

  building_crs <- st_crs(layer_crs_wkt)
  if (is.na(building_crs)) stop("建筑物图层 CRS 无效。")

  filter_wkt <- st_as_text(
    st_transform(buffer_utm, building_crs)
  )[[1]]

  nearby_buildings <- st_read(
    gpkg_file,
    layer = layer_name,
    wkt_filter = filter_wkt,
    quiet = TRUE
  )

  if (!inherits(nearby_buildings, "sf") || nrow(nearby_buildings) == 0) {
    return(0)
  }

  nearby_buildings <- clean_polygon_sf(nearby_buildings)
  if (is.null(nearby_buildings) || nrow(nearby_buildings) == 0) return(0)

  nearby_buildings <- st_transform(nearby_buildings, local_utm)
  nearby_buildings <- clean_polygon_sf(nearby_buildings)
  if (is.null(nearby_buildings) || nrow(nearby_buildings) == 0) return(0)

  intersection_geometry <- suppressWarnings(
    st_intersection(
      st_geometry(nearby_buildings),
      buffer_utm
    )
  )

  if (length(intersection_geometry) == 0) return(0)

  intersection_geometry <- intersection_geometry[
    !st_is_empty(intersection_geometry)
  ]

  if (length(intersection_geometry) == 0) return(0)

  # 合并后计算面积，避免重叠建筑重复计数。
  union_geometry <- suppressWarnings(st_union(intersection_geometry))
  area <- as.numeric(st_area(union_geometry))

  if (length(area) == 0 || is.na(area)) 0 else area
}

process_one_station <- function(station_row, province_files, buffer_m) {
  results <- vector("list", nrow(province_files))

  for (i in seq_len(nrow(province_files))) {
    current_file <- province_files[i, ]

    results[[i]] <- tryCatch({
      area <- calculate_area_in_one_gpkg(
        station_row = station_row,
        gpkg_file = current_file$gpkg_file[[1]],
        layer_name = current_file$layer_name[[1]],
        layer_crs_wkt = current_file$crs_wkt[[1]],
        buffer_m = buffer_m
      )

      tibble(
        buildingFootprint = area,
        CRS_Group = as.character(current_file$CRS_Group[[1]]),
        Source_GPKG = current_file$file_name[[1]],
        Source_Layer = current_file$layer_name[[1]],
        Error = NA_character_
      )
    }, error = function(e) {
      tibble(
        buildingFootprint = NA_real_,
        CRS_Group = as.character(current_file$CRS_Group[[1]]),
        Source_GPKG = current_file$file_name[[1]],
        Source_Layer = current_file$layer_name[[1]],
        Error = conditionMessage(e)
      )
    })
  }

  results <- bind_rows(results)
  valid_results <- results %>% filter(!is.na(buildingFootprint))

  if (nrow(valid_results) > 0) {
    # 原始逻辑：多个 CRS 文件中选择面积最大的结果。
    best <- valid_results %>%
      slice_max(buildingFootprint, n = 1, with_ties = FALSE)
  } else {
    best <- tibble(
      buildingFootprint = NA_real_,
      CRS_Group = NA_character_,
      Source_GPKG = NA_character_,
      Source_Layer = NA_character_
    )
  }

  errors <- results %>% filter(!is.na(Error))
  combined_error <- if (nrow(errors) == 0) {
    NA_character_
  } else {
    paste(
      paste0(errors$Source_GPKG, ": ", errors$Error),
      collapse = " | "
    )
  }

  tibble(
    station_row_id = station_row$station_row_id[[1]],
    AirQualityStation = as.character(station_row$AirQualityStation[[1]]),
    Longitude = station_row$Longitude[[1]],
    Latitude = station_row$Latitude[[1]],
    Province = as.character(station_row$english_name[[1]]),
    CRS_Group = best$CRS_Group[[1]],
    buildingFootprint = best$buildingFootprint[[1]],
    Source_GPKG = best$Source_GPKG[[1]],
    Source_Layer = best$Source_Layer[[1]],
    Error = combined_error
  )
}

process_one_province <- function(province_name) {
  output_file <- file.path(
    result_dir,
    paste0(province_name, "_stations_buildings_30m.csv")
  )

  error_file <- file.path(
    result_dir,
    paste0(province_name, "_station_errors_30m.csv")
  )

  if (file.exists(output_file) && !overwrite_existing) {
    return(paste0(province_name, ": skipped"))
  }

  if (overwrite_existing && file.exists(error_file)) {
    file.remove(error_file)
  }

  province_stations <- stations_with_province %>%
    filter(english_name == province_name)

  province_files <- building_inventory %>%
    filter(Province == province_name)

  if (nrow(province_files) == 0) {
    output <- province_stations %>%
      st_drop_geometry() %>%
      transmute(
        station_row_id,
        AirQualityStation = as.character(AirQualityStation),
        Longitude,
        Latitude,
        Province = province_name,
        CRS_Group = NA_character_,
        buildingFootprint = NA_real_,
        Source_GPKG = NA_character_,
        Source_Layer = NA_character_
      )

    write_csv(output, output_file)

    return(
      paste0(
        province_name,
        ": no building gpkg; ",
        nrow(output),
        " stations saved as NA"
      )
    )
  }

  output <- bind_rows(
    lapply(seq_len(nrow(province_stations)), function(i) {
      process_one_station(
        province_stations[i, ],
        province_files,
        buffer_m
      )
    })
  )

  write_csv(output %>% select(-Error), output_file)

  station_errors <- output %>% filter(!is.na(Error))
  if (nrow(station_errors) > 0) {
    write_csv(station_errors, error_file)
  }

  paste0(
    province_name,
    ": completed; ",
    nrow(output),
    " stations; ",
    sum(!is.na(output$buildingFootprint)),
    " successful"
  )
}

# ---------------------------- 运行前测试 ----------------------------

common_provinces <- intersect(province_list, building_provinces)
if (length(common_provinces) == 0) {
  stop("站点所在省份与建筑物文件省份没有交集。")
}

test_province <- common_provinces[[1]]
test_station <- stations_with_province %>%
  filter(english_name == test_province) %>%
  slice(1)

test_file <- building_inventory %>%
  filter(Province == test_province) %>%
  slice(1)

message(
  "运行前测试：",
  test_province,
  "，站点 ",
  test_station$AirQualityStation[[1]]
)

test_area <- tryCatch(
  calculate_area_in_one_gpkg(
    station_row = test_station,
    gpkg_file = test_file$gpkg_file[[1]],
    layer_name = test_file$layer_name[[1]],
    layer_crs_wkt = test_file$crs_wkt[[1]],
    buffer_m = buffer_m
  ),
  error = function(e) {
    stop("运行前测试失败：", conditionMessage(e))
  }
)

message("运行前测试成功，建筑面积：", round(test_area, 3), " m²。")

# ---------------------------- 并行处理 ----------------------------

available_cores <- parallel::detectCores(logical = TRUE)

worker_number <- min(
  workers_requested,
  length(province_list),
  max(1, available_cores - 1)
)

message("处理省份数：", length(province_list))
message("并行进程数：", worker_number)

cl <- makeCluster(worker_number)
registerDoParallel(cl)

province_status <- tryCatch({
  foreach(
    province_name = province_list,
    .packages = c("dplyr", "readr", "stringr", "tibble", "sf", "lwgeom"),
    .export = c(
      "get_local_utm_epsg",
      "clean_polygon_sf",
      "calculate_area_in_one_gpkg",
      "process_one_station",
      "process_one_province",
      "stations_with_province",
      "building_inventory",
      "result_dir",
      "buffer_m",
      "overwrite_existing"
    ),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch(
      process_one_province(province_name),
      error = function(e) {
        paste0(province_name, ": ERROR: ", conditionMessage(e))
      }
    )
  }
}, finally = {
  stopCluster(cl)
  registerDoSEQ()
})

province_status_table <- tibble(
  Province = province_list,
  Status = vapply(province_status, function(x) {
    if (inherits(x, "error")) {
      paste0("ERROR: ", conditionMessage(x))
    } else if (is.null(x) || length(x) == 0) {
      NA_character_
    } else {
      as.character(x)[[1]]
    }
  }, character(1))
)

print(province_status_table, n = Inf)

write_csv(
  province_status_table,
  file.path(result_dir, "Province_Processing_Status.csv")
)

# ---------------------------- 合并错误日志 ----------------------------

station_error_files <- list.files(
  result_dir,
  pattern = "_station_errors_30m\\.csv$",
  full.names = TRUE
)

all_error_output <- file.path(
  result_dir,
  "All_Station_Processing_Errors.csv"
)

if (length(station_error_files) > 0) {
  all_station_errors <- map_dfr(
    station_error_files,
    ~ read_csv(.x, show_col_types = FALSE)
  )

  write_csv(all_station_errors, all_error_output)
  message("站点级错误记录数：", nrow(all_station_errors))
} else {
  if (file.exists(all_error_output)) file.remove(all_error_output)
  message("没有站点级计算错误。")
}

# ---------------------------- 合并全国结果 ----------------------------

province_output_files <- list.files(
  result_dir,
  pattern = "_stations_buildings_30m\\.csv$",
  full.names = TRUE
)

if (length(province_output_files) == 0) {
  stop("没有生成任何省份结果。")
}

final_result <- map_dfr(
  province_output_files,
  ~ read_csv(.x, show_col_types = FALSE)
) %>%
  arrange(station_row_id) %>%
  distinct(station_row_id, .keep_all = TRUE)

write_csv(final_result, final_output)

expected_n <- nrow(stations_with_province)
actual_n <- nrow(final_result)

if (actual_n != expected_n) {
  warning(
    "最终结果数量与有效站点数不一致：",
    "有效站点数=", expected_n,
    "；最终结果数=", actual_n
  )
}

message("================================================")
message("全部完成。")
message("省级结果目录：", result_dir)
message("全国结果：", final_output)
message("有效站点数：", expected_n)
message("最终结果数：", actual_n)
message("成功计算：", sum(!is.na(final_result$buildingFootprint)))
message("建筑面积 > 0：", sum(final_result$buildingFootprint > 0, na.rm = TRUE))
message("建筑面积 = 0：", sum(final_result$buildingFootprint == 0, na.rm = TRUE))
message("建筑面积为 NA：", sum(is.na(final_result$buildingFootprint)))
message("================================================")