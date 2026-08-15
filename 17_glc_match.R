#!/usr/bin/env Rscript
# 从各包的瓦片清单中，挑出覆盖 924 个站点 1000 m 缓冲区的 Annual 瓦片。
#
# 包名按 10 度经度命名(E80-E85.zip 含 E80 与 E85 两列瓦片)；
# 包内同时有 GLC_FCS30D_20002022_*_Annual.tif(逐年，本项目所需)
# 与 GLC_FCS30_19852000_5yearsMap_*(5 年一期，忽略)。
#
# 瓦片命名约定(2026-08-14 实测更正):
#   瓦片是 5x5 度。Exxx 是**西**边界(覆盖 xxx .. xxx+5)，
#   Nyy 是**北**边界(覆盖 yy-5 .. yy)——不是南边界。
#   证据: 已解出的 glc_long_E105N25.csv 里站点纬度 21.48-24.98，
#         glc_long_E115N40.csv 里 35.07-40.00，均落在 [yy-5, yy]。
#   早先按南边界匹配，导致每列都选错一行、下载了 10 个空瓦片(各约 400 MB)。
#
# 另: 按缓冲区网格的外接矩形而非站点坐标匹配，跨瓦片边界的缓冲区会同时选中
#     两个瓦片，各自提取一部分像元，合并时按像元数加权即可。

suppressPackageStartupMessages({ library(data.table); library(terra) })
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")

GRID  <- "data_raw/hcsif/grid/hcsif_grid_buf1000.shp"
PARTS <- "data_raw/glc/parts"
TILE  <- 5L

# --- 1. 可用瓦片清单 -------------------------------------------------------
fs <- list.files("data_raw/glc/lists", pattern = "^E[0-9]+[.]txt$", full.names = TRUE)
if (!length(fs)) stop("没有瓦片清单")
tl <- rbindlist(lapply(fs, function(f) {
  d <- fread(f, header = FALSE, sep = "\t", col.names = c("name", "csize"),
             showProgress = FALSE)
  d[, band := sub("[.]txt$", "", basename(f))]
  d[grepl("^GLC_FCS30D_20002022_E[0-9]+N[0-9]+_Annual[.]tif$", name)]
}))
if (!nrow(tl)) stop("清单里没有北半球 Annual 瓦片")
tl[, `:=`(lon_w = as.integer(sub("^GLC_FCS30D_20002022_E([0-9]+)N.*", "\\1", name)),
          lat_n = as.integer(sub("^GLC_FCS30D_20002022_E[0-9]+N([0-9]+)_.*", "\\1", name)))]
tl[, tid := sub(".*_(E[0-9]+N[0-9]+)_Annual[.]tif$", "\\1", name)]

# --- 2. 缓冲区网格的外接矩形 -> 覆盖到哪些瓦片 -----------------------------
v  <- vect(GRID)
g  <- as.data.table(geom(v))          # 顶点坐标，geom 列是网格序号
bb <- g[, .(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)), by = geom]
bb[, `:=`(stat_id = v$stat_id[geom], cell = v$cell[geom])]

need <- bb[, {
  lw <- seq(floor(xmin / TILE) * TILE, floor(xmax / TILE) * TILE, by = TILE)
  ln <- seq(ceiling(ymin / TILE) * TILE, ceiling(ymax / TILE) * TILE, by = TILE)
  CJ(lon_w = as.integer(lw), lat_n = as.integer(ln))
}, by = .(geom, stat_id)]

need <- need[, .(n_cell = .N, n_station = uniqueN(stat_id)), by = .(lon_w, lat_n)]

m <- merge(need, tl[, .(band, name, csize, tid, lon_w, lat_n)],
           by = c("lon_w", "lat_n"), all.x = TRUE)
if (m[is.na(name), .N]) {
  warning("以下瓦片站点需要但清单里没有:\n",
          paste(m[is.na(name), sprintf("  E%dN%d (%d 站)", lon_w, lat_n, n_station)],
                collapse = "\n"))
  m <- m[!is.na(name)]
}

# --- 3. 已完成的分片直接跳过 -----------------------------------------------
done <- sub("^glc_long_(.*)[.]csv$", "\\1",
            list.files(PARTS, "^glc_long_.*[.]csv$"))
m[, done := tid %in% done]
setorder(m, -csize)

fwrite(m[, .(band, name, csize, n_station, tid, done)], "data_raw/glc/needed_tiles.csv")

cat(sprintf("需要瓦片 %d 个 | 已完成 %d | 待取 %d\n",
            nrow(m), sum(m$done), sum(!m$done)))
cat(sprintf("待取压缩合计 %.1f GB\n", sum(as.numeric(m[!(done)]$csize)) / 1e9))
cat("覆盖站点(去重):", uniqueN(need$lon_w) * 0 + nrow(bb[, .N, by = stat_id]), "\n")
cat("\n待取瓦片:\n")
print(m[!(done), .(tid, MB = round(csize / 1048576), n_station)])
