#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 17_glc_extract.R — 从单个 GLC_FCS30D 瓦片提取站点 1000 m 缓冲区的类别占比
#
# 瓦片是 5°x5°、30 m、23 波段(2000-2022 逐年)的 GeoTIFF，由
# 17_glc_fetch_tile.py 从远程 zip 中按字节范围取回。
#
# 输出按瓦片分片写出，最后由 17_glc_run.sh 合并;
# 落在多个瓦片交界处的站点会在各瓦片各得一部分像元，合并时按像元数加权。
#
# 用法: Rscript 17_glc_extract.R <tile.tif> <tile_id> <out_dir>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(terra); library(data.table) })

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
GRID_SHP <- file.path(PROJ, "data_raw/hcsif/grid/hcsif_grid_buf1000.shp")
YEARS    <- 2000:2022

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("用法: Rscript 17_glc_extract.R <tile.tif> <tile_id> <out_dir>")
tif <- a[1]; tile_id <- a[2]; out_dir <- a[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

r <- rast(tif)
if (nlyr(r) != length(YEARS))
  log_msg("注意: 波段数 ", nlyr(r), " 与预期 ", length(YEARS), " 不符")

v <- vect(GRID_SHP)
# 只保留与本瓦片相交的网格，避免全量求交
keep <- terra::is.related(v, as.polygons(ext(r), crs = crs(r)), "intersects")
v <- v[keep]
if (!nrow(v)) { log_msg(tile_id, ": 无网格落入，跳过"); quit(status = 0) }
log_msg(tile_id, ": ", nrow(v), " 个网格, ", uniqueN(v$stat_id), " 个站点")

res_l <- vector("list", nlyr(r))
for (b in seq_len(nlyr(r))) {
  ex <- as.data.table(terra::extract(r[[b]], v))
  setnames(ex, 1:2, c("ID", "cls"))
  ex <- ex[!is.na(cls) & cls > 0]
  if (!nrow(ex)) next
  tab <- ex[, .N, by = .(ID, cls)]
  tab[, `:=`(stat_id = v$stat_id[ID], cell = v$cell[ID], year = YEARS[b])]
  res_l[[b]] <- tab[, .(stat_id, cell, year, cls, n = N)]
}
out <- rbindlist(res_l[!vapply(res_l, is.null, logical(1))])
if (!nrow(out)) { log_msg(tile_id, ": 无有效像元"); quit(status = 0) }

fp <- file.path(out_dir, paste0("glc_long_", tile_id, ".csv"))
fwrite(out, fp)
log_msg("写出 ", basename(fp), " (", nrow(out), " 行, ",
        uniqueN(out$cls), " 个类别, ", uniqueN(out$year), " 年)")
