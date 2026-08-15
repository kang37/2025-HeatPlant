#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 19_era5_extract.R — 把 ERA5-Land 下行短波辐射提取到各气象站
#
# 输入: data_raw/era5_land/era5land_rsds_<year>.nc
#       每年 153 层 = 5/1 - 9/30 逐日，0.1 度，单位 J m-2
#
# 三个要点:
#
# 1) 单位: ssrd 是从当日 00 UTC 起的累积量，取 23:00 时次即当日总量(J m-2)。
#    除以 86400 得日均通量 W m-2。中国夏季日均应在 150-300 W m-2。
#
# 2) 时区: ERA5-Land 按 UTC 日累积，中国是 UTC+8，故这个"日总量"实际覆盖
#    当地 08:00 当日 - 07:00 次日。对逐日分析会有约 8 小时相位偏移；但本项目
#    的 SIF 是 8 天合成，8 天窗口内这个偏移只影响首尾两端，量级约 4%，可接受。
#    脚本同时输出逐日值和 8 天合成均值，后者才是与 SIF 配对用的。
#
# 3) 空间尺度: ERA5-Land 约 9 km，远粗于 1000 m 缓冲区。同一站点缓冲区内所有
#    网格落在同一或相邻像元上，逐网格提取无信息增益(脚本会实测并报告站内离散度)。
#    故按站点提取，用双线性插值取站点坐标处的值。
#
# 用法: Rscript 19_era5_extract.R [起始年] [结束年]
#       Rscript 19_era5_extract.R 2020 2022
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(terra); library(data.table) })

PROJ    <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
NC_DIR  <- file.path(PROJ, "data_raw/era5_land")
OUT_DIR <- file.path(PROJ, "data_raw/covariates_1km")
ST_CSV  <- file.path(PROJ, "data_raw/hcsif/stations_924.csv")
GRID    <- file.path(PROJ, "data_raw/hcsif/grid/hcsif_grid_buf1000.shp")

SIF_DOY  <- seq(124, 268, by = 8)  # 与 13_hcsif_extract.R 输出的时相网格一致
J_TO_W   <- 86400                  # J m-2 (日累积) -> W m-2 (日均通量)
PLAUS    <- c(20, 400)             # 5-9 月日均辐射的合理区间，仅作量级兜底
COAST_R  <- 25000                  # 落在 ERA5-Land 海洋掩膜上的站，取此半径内陆地像元均值

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

a  <- commandArgs(trailingOnly = TRUE)
y0 <- if (length(a) >= 1) as.integer(a[1]) else 2020
y1 <- if (length(a) >= 2) as.integer(a[2]) else 2022

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

st <- fread(ST_CSV)
setnames(st, "meteo_stat", "stat_id")
pts <- vect(st, geom = c("longitude", "latitude"), crs = "EPSG:4326")
log_msg("站点 ", nrow(st), " 个")

# --- 层名里的 valid_time 是 epoch 秒，据此还原日期 -------------------------
layer_dates <- function(r) {
  ep <- as.numeric(sub(".*valid_time=", "", names(r)))
  if (anyNA(ep)) stop("层名里解析不出 valid_time")
  as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")
}

# --- 站内离散度自检: ERA5 这么粗，缓冲区内各网格是否真的等值 ---------------
check_within_station <- function(r) {
  if (!file.exists(GRID)) { log_msg("无缓冲区 shp，跳过站内离散度自检"); return(invisible()) }
  v  <- vect(GRID)
  cn <- centroids(v)
  ex <- as.data.table(terra::extract(r[[1]], cn, method = "bilinear"))
  ex[, stat_id := cn$stat_id]
  setnames(ex, 2, "val")
  s <- ex[!is.na(val), .(rng = max(val) - min(val), mu = mean(val)), by = stat_id]
  log_msg(sprintf("站内离散度自检(第1天): 缓冲区内极差/均值 中位数 %.4f%%, 最大 %.4f%%",
                  100 * median(s$rng / s$mu), 100 * max(s$rng / s$mu)))
}

first <- TRUE
for (y in y0:y1) {
  fp <- file.path(NC_DIR, sprintf("era5land_rsds_%d.nc", y))
  if (!file.exists(fp)) { log_msg("缺文件，跳过: ", basename(fp)); next }

  r  <- rast(fp)
  ts <- layer_dates(r)
  log_msg(y, ": ", nlyr(r), " 层, ", format(min(ts), "%m-%d %H:%M"), " UTC 至 ",
          format(max(ts), "%m-%d %H:%M"), " UTC")

  if (first) { check_within_station(r); first <- FALSE }

  # 双线性插值到站点坐标; 转置成 长表(站 x 日)
  ex <- as.data.table(terra::extract(r, pts, method = "bilinear", ID = FALSE))
  ex[, stat_id := st$stat_id]

  # ERA5-Land 只覆盖陆地。个别近海站(如 54623 渤海湾)整站落在海洋掩膜上，
  # 双线性拿不到值，改取周边 25 km 内陆地像元的均值。
  sea <- which(rowSums(!is.na(ex[, -"stat_id"])) == 0)
  if (length(sea)) {
    bufv <- buffer(pts[sea], COAST_R)
    fb <- as.data.table(terra::extract(r, bufv, fun = mean, na.rm = TRUE, ID = FALSE))
    for (j in seq_len(ncol(fb))) set(ex, i = sea, j = j, value = fb[[j]])
    still <- sum(rowSums(!is.na(ex[sea, -"stat_id", with = FALSE])) == 0)
    log_msg("  海洋掩膜补值: ", length(sea), " 站 (", paste(st$stat_id[sea], collapse = ","),
            "), 补后仍缺 ", still, " 站")
  }
  d <- melt(ex, id.vars = "stat_id", variable.name = "lyr", value.name = "ssrd_J")
  d[, date := as.Date(ts[as.integer(lyr)])]          # UTC 累积窗口所属的那一天
  d[, lyr := NULL]
  d[, rsds := ssrd_J / J_TO_W]
  d[, ssrd_J := NULL]
  d[, doy := as.integer(format(date, "%j"))]
  d[, year := y]

  bad <- d[!is.na(rsds) & (rsds < PLAUS[1] | rsds > PLAUS[2]), .N]
  log_msg(sprintf("  日均辐射 %.1f - %.1f W m-2 (中位 %.1f), 越界 %d 个, 缺测 %d 个",
                  min(d$rsds, na.rm = TRUE), max(d$rsds, na.rm = TRUE),
                  median(d$rsds, na.rm = TRUE), bad, d[is.na(rsds), .N]))

  fwrite(d[, .(stat_id, year, doy, date, rsds)],
         file.path(OUT_DIR, sprintf("era5_rsds_daily_%d.csv", y)))

  # --- 聚合到 SIF 的 8 天时相: 时相 d 覆盖 DOY d..d+7 -----------------------
  # findInterval 对早于第一个时相的日子返回 0，直接下标会丢元素，故先算下标再置 NA
  k <- findInterval(d$doy, SIF_DOY)
  d[, sif_doy := ifelse(k == 0L | doy > max(SIF_DOY) + 7L, NA_integer_, SIF_DOY[pmax(k, 1L)])]
  agg <- d[!is.na(sif_doy), .(rsds_mean = mean(rsds, na.rm = TRUE),
                              rsds_max  = max(rsds,  na.rm = TRUE),
                              n_day     = sum(!is.na(rsds))),
           by = .(stat_id, year, doy = sif_doy)]
  agg[is.nan(rsds_mean), rsds_mean := NA_real_]      # 全 NA 时 mean 返回 NaN
  agg[is.infinite(rsds_max), rsds_max := NA_real_]
  fwrite(agg, file.path(OUT_DIR, sprintf("era5_rsds_sif8d_%d.csv", y)))
  log_msg("  写出 ", nrow(d), " 行逐日 / ", nrow(agg), " 行 8 天合成 (n_day 应为 8: ",
          paste(sort(unique(agg$n_day)), collapse = ","), ")")
}
log_msg("完成")
