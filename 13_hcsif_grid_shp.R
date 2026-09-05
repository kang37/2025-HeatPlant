#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_hcsif_grid_shp.R
#
# 输出"各站点 SIF 取值所用网格"的矢量文件，并附带缺失统计。
#
#   hcsif_grid_point.shp    每站 1 个多边形 —— 站点落入的那个 500 m 像元(点值来源)
#   hcsif_grid_buf750.shp   每站若干个多边形 —— 750 m 缓冲区覆盖的像元
#   hcsif_grid_buf1000.shp  同上，1000 m 缓冲区
#   (BUFFER_M 里每加一个半径就多输出一个 hcsif_grid_buf<r>.shp)
#
# 用法:
#   Rscript 13_hcsif_grid_shp.R [out_dir]
#
# 网格几何优先从现存的 HCSIF 栅格读取；若原始栅格已被流式流程删除，
# 则回退到实测确认的网格参数(见下方 GRID)。两条路径给出的结果一致。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
RAW_DIR  <- file.path(PROJ, "data_raw/hcsif/tmp")
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v2")
BUFFER_M <- c(750, 1000, 2000, 3000)   # 每个半径生成一个 hcsif_grid_buf<r>.shp

# 实测自 2000196.tif / .tfw 的网格参数(EPSG:4326)
# .tfw 给出左上角像元中心，栅格边界需再外扩半个像元。
GRID <- list(
  res      = 0.0044915764,
  ctr_x    = 73.5024023349,   # 左上角像元中心 x
  ctr_y    = 53.5598030274,   # 左上角像元中心 y
  nrow     = 10518L,
  ncol     = 13712L,
  crs      = "EPSG:4326"
)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# --- 网格模板 --------------------------------------------------------------

build_template <- function() {
  tifs <- list.files(RAW_DIR, pattern = "\\.tif$", full.names = TRUE)
  if (length(tifs)) {
    r <- terra::rast(tifs[1])
    log_msg("网格来源: 实际栅格 ", basename(tifs[1]))
    return(terra::rast(terra::ext(r), resolution = terra::res(r), crs = terra::crs(r)))
  }
  half <- GRID$res / 2
  xmin <- GRID$ctr_x - half
  ymax <- GRID$ctr_y + half
  e <- terra::ext(xmin, xmin + GRID$ncol * GRID$res,
                  ymax - GRID$nrow * GRID$res, ymax)
  log_msg("网格来源: 内置参数(原始栅格已清理)")
  terra::rast(e, resolution = GRID$res, crs = GRID$crs)
}

# 由像元号构造方格多边形
cells_to_polys <- function(tmpl, cells, ids) {
  xy   <- terra::xyFromCell(tmpl, cells)
  half <- terra::res(tmpl) / 2
  n    <- length(cells)
  geom <- cbind(
    object = rep(seq_len(n), each = 4L),
    part   = 1L,
    x = as.vector(t(cbind(xy[, 1] - half[1], xy[, 1] + half[1],
                          xy[, 1] + half[1], xy[, 1] - half[1]))),
    y = as.vector(t(cbind(xy[, 2] + half[2], xy[, 2] + half[2],
                          xy[, 2] - half[2], xy[, 2] - half[2]))),
    hole = 0L
  )
  v <- terra::vect(geom, type = "polygons", crs = terra::crs(tmpl))
  v$stat_id <- ids
  v$cell    <- cells
  v$cx      <- round(xy[, 1], 7)
  v$cy      <- round(xy[, 2], 7)
  v
}

# --- 站点缺失统计 ----------------------------------------------------------

station_stats <- function(st) {
  fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
  if (!length(fs)) {
    log_msg("警告: 尚无提取结果，shapefile 只含几何与站点标识")
    return(data.table(meteo_stat = st$meteo_stat, n_obs = 0L))
  }
  log_msg("汇总 ", length(fs), " 个年度文件: ",
          paste(sub("hcsif_station_", "", tools::file_path_sans_ext(basename(fs))), collapse = ", "))
  d <- rbindlist(lapply(fs, fread), fill = TRUE)

  s <- d[, .(
    n_obs      = .N,
    n_valid    = sum(!is.na(SIF)),
    n_miss     = sum(is.na(SIF)),
    n_masked   = if ("pt_masked" %in% names(d)) sum(pt_masked %in% c(TRUE, "TRUE")) else NA_integer_,
    sif_mean   = round(mean(SIF, na.rm = TRUE), 5),
    buf_mean   = round(mean(SIF_buf750, na.rm = TRUE), 5),
    nvalid_buf = sum(!is.na(SIF_buf750))
  ), by = meteo_stat]
  s[, miss_pct := round(100 * n_miss / n_obs, 2)]
  s[is.nan(sif_mean), sif_mean := NA_real_]
  s[is.nan(buf_mean), buf_mean := NA_real_]
  s
}

# --- 主流程 ----------------------------------------------------------------

main <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  st <- fread(file.path(PROJ, "data_raw/hcsif/stations_924.csv"))
  log_msg("站点数: ", nrow(st))

  tmpl <- build_template()
  log_msg("分辨率 ", signif(terra::res(tmpl)[1], 9),
          " | 范围 ", paste(signif(as.vector(terra::ext(tmpl)), 8), collapse = " "))

  stats <- station_stats(st)

  # --- 1. 点值像元 ---
  cells <- terra::cellFromXY(tmpl, cbind(st$longitude, st$latitude))
  if (any(is.na(cells))) stop("有站点落在栅格范围之外: ", sum(is.na(cells)))

  vp <- cells_to_polys(tmpl, cells, st$meteo_stat)
  vp$lon <- st$longitude
  vp$lat <- st$latitude
  att <- merge(data.table(meteo_stat = st$meteo_stat), stats, by = "meteo_stat",
               all.x = TRUE, sort = FALSE)
  for (cl in setdiff(names(att), "meteo_stat")) vp[[cl]] <- att[[cl]]

  fp1 <- file.path(out_dir, "hcsif_grid_point.shp")
  terra::writeVector(vp, fp1, overwrite = TRUE)
  log_msg("写出 ", basename(fp1), "  (", nrow(vp), " 个像元多边形)")

  # 同一像元被多个站点共用的情况(站点相距 < 500 m)
  dup <- sum(duplicated(cells))
  if (dup > 0) log_msg("注意: ", dup, " 个站点与其他站点共用同一像元")

  # --- 2. 各级缓冲区像元 ---
  pts <- terra::vect(as.data.frame(st), geom = c("longitude", "latitude"), crs = "EPSG:4326")
  for (bm in BUFFER_M) {
    buf <- terra::buffer(pts, width = bm)
    cb  <- terra::cells(tmpl, buf)          # 返回 ID(缓冲区序号) 与 cell
    cb  <- as.data.table(cb)
    setnames(cb, 1:2, c("ID", "cell"))

    vb  <- cells_to_polys(tmpl, cb$cell, st$meteo_stat[cb$ID])
    fpb <- file.path(out_dir, sprintf("hcsif_grid_buf%d.shp", bm))
    terra::writeVector(vb, fpb, overwrite = TRUE)
    log_msg("写出 ", basename(fpb), "  (", nrow(vb), " 个像元，平均每站 ",
            round(nrow(vb) / nrow(st), 1), " 个)")
  }

  # --- 3. 缺失概览 ---
  if ("n_valid" %in% names(stats)) {
    cat("\n=== SIF 缺失概览 ===\n")
    cat("记录总数        :", sum(stats$n_obs), "\n")
    cat("有效点值        :", sum(stats$n_valid),
        sprintf("(%.1f%%)\n", 100 * sum(stats$n_valid) / sum(stats$n_obs)))
    cat("缺失点值        :", sum(stats$n_miss),
        sprintf("(%.1f%%)\n", 100 * sum(stats$n_miss) / sum(stats$n_obs)))
    cat("其中被掩膜(0值) :", sum(stats$n_masked, na.rm = TRUE), "\n")
    cat("缓冲均值有效    :", sum(stats$nvalid_buf),
        sprintf("(%.1f%%)\n", 100 * sum(stats$nvalid_buf) / sum(stats$n_obs)))
    cat("\n按站点看:\n")
    cat("  全部时相都有值  :", sum(stats$n_miss == 0), "站\n")
    cat("  部分缺失        :", sum(stats$n_miss > 0 & stats$n_valid > 0), "站\n")
    cat("  全部缺失(完全无值):", sum(stats$n_valid == 0), "站\n")
    cat("\n缺失率分布:\n")
    print(summary(stats$miss_pct))
    worst <- stats[order(-miss_pct)][1:min(10, .N),
                   .(meteo_stat, n_obs, n_miss, miss_pct, nvalid_buf)]
    cat("\n缺失最严重的站点:\n"); print(worst)
  }
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
main(if (length(args)) args[1] else file.path(PROJ, "data_raw/hcsif/grid"))
