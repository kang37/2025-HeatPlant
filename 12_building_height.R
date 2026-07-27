# 12_building_height.R
# 从CMAB中国建筑高度数据集提取各气象站点30m缓冲区内的建筑面积加权平均高度
# 来源：CMAB - A First National-Scale Multi-Attribute Building Dataset
#       (Lu et al. 2024, Scientific Data)
# 缓冲区：30m（与building_footprint保持一致，来源：building_extract.R / 王琳）
# 输出：data_raw/building_height_station.rds
#   字段：meteo_stat_id, mean_height（面积加权平均高度，m）, n_buildings

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")

ZIP_DIR   <- "data_raw/中国建筑高度/分省数据"
CACHE_OUT <- "data_raw/building_height_station.rds"
BUFFER_M  <- 30
TMP_BASE  <- "/tmp/cmab_extract"

# --------------------------------------------------------------------------
# 阶段1：用Python将所有省份zip解压到ASCII路径（处理中文目录名）
# --------------------------------------------------------------------------

py_extract_script <- tempfile(fileext = ".py")
writeLines(sprintf('
import zipfile, os, shutil

zip_dir  = "%s"
tmp_base = "%s"

zips = [f for f in os.listdir(zip_dir) if f.endswith(".zip")]
print(f"共 {len(zips)} 个省份zip")

for zname in sorted(zips):
    prov = zname.replace(".zip", "")
    out  = os.path.join(tmp_base, prov)
    os.makedirs(out, exist_ok=True)
    if len([f for f in os.listdir(out) if f.endswith(".shp")]) > 0:
        print(f"  [{prov}] 已解压，跳过")
        continue
    try:
        with zipfile.ZipFile(os.path.join(zip_dir, zname)) as zf:
            zf.extractall(out)
        # 平铺：将中文子目录内文件复制到out根目录
        for root, dirs, files in os.walk(out):
            if root == out:
                continue
            for f in files:
                shutil.copy2(os.path.join(root, f), os.path.join(out, f))
        n_shp = len([f for f in os.listdir(out) if f.endswith(".shp")])
        print(f"  [{prov}] 解压完成，{n_shp} 个shp")
    except Exception as e:
        print(f"  [{prov}] 解压失败: {e}")
', ZIP_DIR, TMP_BASE), py_extract_script)

cat("=== 阶段1：解压省份zip ===\n")
system2("python3", py_extract_script)

# --------------------------------------------------------------------------
# 阶段2：逐站点提取建筑高度
# --------------------------------------------------------------------------

cat("\n=== 阶段2：提取建筑高度 ===\n")

station_coords <- read_csv(
  "data_proc/output_10y_built_up_05_01/trri_station.csv",
  show_col_types = FALSE
) %>%
  dplyr::select(meteo_stat_id, longitude, latitude) %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(meteo_stat_id = as.character(meteo_stat_id))

cat(sprintf("站点总数: %d\n", nrow(station_coords)))

# 断点续传
if (file.exists(CACHE_OUT)) {
  done <- readRDS(CACHE_OUT) %>%
    mutate(meteo_stat_id = as.character(meteo_stat_id))
  cat(sprintf("已有缓存: %d 个站点\n", nrow(done)))
  remaining <- station_coords %>%
    filter(!meteo_stat_id %in% done$meteo_stat_id)
} else {
  done      <- tibble()
  remaining <- station_coords
}

cat(sprintf("待处理: %d 个站点\n", nrow(remaining)))
if (nrow(remaining) == 0) {
  cat("全部完成。\n")
  quit(save = "no")
}

# 站点缓冲区（EPSG:3857，与CMAB一致）
stations_sf <- remaining %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(3857)
stations_buf  <- st_buffer(stations_sf, dist = BUFFER_M)
station_cents <- st_coordinates(st_centroid(stations_buf))

# 结果初始化
results <- tibble(
  meteo_stat_id = remaining$meteo_stat_id,
  mean_height   = NA_real_,
  n_buildings   = NA_integer_
)

prov_dirs <- list.dirs(TMP_BASE, recursive = FALSE)
cat(sprintf("找到 %d 个省份目录\n\n", length(prov_dirs)))

for (pdir in sort(prov_dirs)) {
  prov_name <- basename(pdir)
  shp_files <- list.files(pdir, pattern = "[.]shp$", full.names = TRUE)
  if (length(shp_files) == 0) next

  # 读取省份所有建筑
  prov_sf <- tryCatch({
    maps <- lapply(shp_files, function(f) tryCatch(
      sf::st_read(f, quiet = TRUE) %>%
        dplyr::select(Height, geometry) %>%
        mutate(Height = as.numeric(Height)),
      error = function(e) NULL
    ))
    maps <- Filter(Negate(is.null), maps)
    if (length(maps) == 0) return(NULL)
    do.call(rbind, maps)
  }, error = function(e) NULL)

  if (is.null(prov_sf) || nrow(prov_sf) == 0) {
    cat(sprintf("  [%s] 读取失败\n", prov_name)); next
  }

  # 确保CRS为3857
  if (!isTRUE(sf::st_crs(prov_sf)$epsg == 3857))
    prov_sf <- sf::st_transform(prov_sf, 3857)

  # 粗筛：站点中心在省份bbox ± 50km范围内
  bb <- sf::st_bbox(prov_sf)
  idx_in <- which(
    station_cents[, 1] >= bb["xmin"] - 50000 &
    station_cents[, 1] <= bb["xmax"] + 50000 &
    station_cents[, 2] >= bb["ymin"] - 50000 &
    station_cents[, 2] <= bb["ymax"] + 50000
  )

  if (length(idx_in) == 0) {
    cat(sprintf("  [%s] 无站点，跳过\n", prov_name)); next
  }

  n_updated <- 0
  for (idx in idx_in) {
    sid  <- remaining$meteo_stat_id[idx]
    buf  <- stations_buf[idx, ]
    inter <- tryCatch(
      suppressWarnings(sf::st_intersection(prov_sf, buf)),
      error = function(e) NULL
    )
    if (is.null(inter) || nrow(inter) == 0) {
      # 缓冲区内无建筑 → 高度记为0（与footprint=0一致）
      if (is.na(results$mean_height[results$meteo_stat_id == sid])) {
        results$mean_height[results$meteo_stat_id == sid]  <- 0
        results$n_buildings[results$meteo_stat_id == sid]  <- 0L
      }
      next
    }
    areas  <- as.numeric(sf::st_area(inter))
    heights <- as.numeric(inter$Height)
    valid   <- !is.na(heights) & heights > 0 & areas > 0
    if (!any(valid)) {
      results$mean_height[results$meteo_stat_id == sid]  <- 0
      results$n_buildings[results$meteo_stat_id == sid]  <- 0L
      next
    }
    results$mean_height[results$meteo_stat_id == sid] <-
      sum(areas[valid] * heights[valid]) / sum(areas[valid])
    results$n_buildings[results$meteo_stat_id == sid] <- sum(valid)
    n_updated <- n_updated + 1
  }

  cat(sprintf("  [%s] %d shp, %d 候选站点, %d 有建筑\n",
              prov_name, length(shp_files), length(idx_in), n_updated))

  # 每省保存进度
  current <- bind_rows(done, results %>% filter(!is.na(mean_height)))
  saveRDS(current, CACHE_OUT)
}

# 最终保存
final <- bind_rows(done, results %>% mutate(
  mean_height  = replace_na(mean_height, 0),
  n_buildings  = replace_na(n_buildings, 0L)
))
saveRDS(final, CACHE_OUT)

cat(sprintf(
  "\n=== 完成 ===\nmean_height > 0: %d / %d 站点\n高度范围: %.1f ~ %.1f m\n",
  sum(final$mean_height > 0, na.rm = TRUE), nrow(final),
  min(final$mean_height[final$mean_height > 0], na.rm = TRUE),
  max(final$mean_height, na.rm = TRUE)
))
