#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 16_ntl_1km.R — 夜间灯光(SVNL / NPP-VIIRS-like)按 1000 m 缓冲区提取
#
# 数据: Chen et al., A history reconstructed time series of annual global
#       NPP-VIIRS-like nighttime light data (SRUNet)
#       Figshare DOI 10.6084/m9.figshare.22262545
#       1992-2023, 15 角秒(约 500 m), WGS84, GeoTIFF
#       1992-2011 由 DMSP-OLS 经超分辨率 U-Net 重建; 2012-2023 为 NPP-VIIRS V2
#
# 分辨率 500 m 与 HCSIF 一致，可直接落到既有的取值网格上。
#
# 输出两级(与 15_covariates_1km.R 一致):
#   ntl_cell_<年>.csv     每个 HCSIF 500 m 像元一行
#   ntl_station_<年>.csv  每个站点一行(按像元数加权)
#
# 用法: Rscript 16_ntl_1km.R <年份> <tif路径>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
GRID_SHP <- file.path(PROJ, "data_raw/hcsif/grid/hcsif_grid_buf1000.shp")
OUT_DIR  <- file.path(PROJ, "data_raw/covariates_1km")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript 16_ntl_1km.R <年份> <tif路径>")
year <- as.integer(args[1]); tif <- args[2]
if (!file.exists(tif)) stop("找不到栅格: ", tif)

v <- terra::vect(GRID_SHP)
r <- terra::rast(tif)
log_msg("NTL ", year, ": ", paste(dim(r)[1:2], collapse = " x "),
        " | 分辨率 ", signif(terra::res(r)[1], 6),
        " | CRS ", terra::crs(r, describe = TRUE)$code)

if (!terra::same.crs(v, terra::crs(r))) v <- terra::project(v, terra::crs(r))

# 灯光是连续值，按像元面积均值聚合即可(不像土地覆盖要算类别占比)
ex <- data.table::as.data.table(terra::extract(r, v))
data.table::setnames(ex, 1:2, c("ID", "ntl"))

agg <- ex[, .(ntl_mean = mean(ntl, na.rm = TRUE),
              ntl_max  = max(ntl, na.rm = TRUE),
              n_px     = sum(!is.na(ntl))), by = ID]
agg <- agg[data.table::data.table(ID = seq_len(nrow(v))), on = "ID"]
agg[is.nan(ntl_mean), ntl_mean := NA_real_]
agg[is.infinite(ntl_max), ntl_max := NA_real_]
agg[is.na(n_px), n_px := 0L]

cell <- data.table(stat_id = v$stat_id, cell = v$cell,
                   ntl_mean = agg$ntl_mean, ntl_max = agg$ntl_max,
                   n_px = agg$n_px, year = year)
fwrite(cell, file.path(OUT_DIR, sprintf("ntl_cell_%d.csv", year)))
log_msg("网格级写出: ntl_cell_", year, ".csv (", nrow(cell), " 行)")

st <- cell[, .(ntl_mean = sum(ntl_mean * n_px, na.rm = TRUE) / sum(n_px),
               ntl_max  = max(ntl_max, na.rm = TRUE),
               n_cell   = .N, n_px = sum(n_px)), by = stat_id]
st[is.nan(ntl_mean), ntl_mean := NA_real_]
st[is.infinite(ntl_max), ntl_max := NA_real_]
st[, year := year]

if (nrow(st) != 924) log_msg("!! 站点数 ", nrow(st), " != 924")
fwrite(st, file.path(OUT_DIR, sprintf("ntl_station_%d.csv", year)))
log_msg("站点级写出: ntl_station_", year, ".csv (", nrow(st), " 站)")
log_msg(sprintf("ntl_mean 均值 %.2f | 中位 %.2f | 范围 [%.2f, %.2f]",
                mean(st$ntl_mean, na.rm = TRUE), median(st$ntl_mean, na.rm = TRUE),
                min(st$ntl_mean, na.rm = TRUE), max(st$ntl_mean, na.rm = TRUE)))
