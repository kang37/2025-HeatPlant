#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 21_meteo_agg.R — 从台站原始气象记录聚合出站点级的气温/湿度/降水/云量
#
# 原始文件: data_raw/meteo_data_1961-2023/<站号>.txt
#   第 1 行  站号,经度,纬度
#   第 2 行  date,tmax,tmin,tavg,RH,CF,precip
#   缺测标记 999999.0
#
# 与 SIF 窗口一致，只取 2000-2022 年 5-9 月。
# 输出站点级多年均值，同时给出各要素的有效天数占比，便于判断可用性。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_raw/covariates_1km"
MISS <- 999999.0
# 物理量程: 中国 5-9 月的合理范围, 超界即判缺测
RANGE <- list(tavg = c(-60, 60), tmax = c(-60, 60), RH = c(0, 100),
              CF = c(0, 100), precip = c(0, 2000))
Y0 <- 2000L; Y1 <- 2022L
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

st <- fread("data_raw/hcsif/stations_924.csv")
setnames(st, "meteo_stat", "stat_id")

read_one <- function(sid) {
  fp <- file.path("data_raw/meteo_data_1961-2023", paste0(sid, ".txt"))
  if (!file.exists(fp)) return(NULL)
  d <- fread(fp, skip = 1, showProgress = FALSE)      # 跳过站号那行, 第2行即表头
  if (!"date" %in% names(d)) return(NULL)
  d[, date := as.IDate(date)]
  d <- d[year(date) >= Y0 & year(date) <= Y1 & month(date) %in% 5:9]
  if (!nrow(d)) return(NULL)
  # 只挡 999999 不够: 个别站(如 59358)存在 250022.6 / 750006.6 这类取值，
  # 是缺测码混入日均计算的产物。改用物理量程过滤，顺带挡掉负降水等脏值。
  for (v in names(RANGE))
    if (v %in% names(d))
      d[!is.finite(get(v)) | get(v) < RANGE[[v]][1] | get(v) > RANGE[[v]][2],
        (v) := NA_real_]
  d[, .(stat_id = sid,
        tavg   = mean(tavg,   na.rm = TRUE), tavg_n   = mean(!is.na(tavg)),
        tmax   = mean(tmax,   na.rm = TRUE), tmax_n   = mean(!is.na(tmax)),
        rh     = mean(RH,     na.rm = TRUE), rh_n     = mean(!is.na(RH)),
        cloud  = mean(CF,     na.rm = TRUE), cloud_n  = mean(!is.na(CF)),
        precip = mean(precip, na.rm = TRUE), precip_n = mean(!is.na(precip)),
        n_day  = .N)]
}

res <- rbindlist(lapply(st$stat_id, read_one), fill = TRUE)
log_msg("有原始记录的站: ", nrow(res), " / ", nrow(st))

for (v in c("tavg", "tmax", "rh", "cloud", "precip"))
  set(res, i = which(is.nan(res[[v]])), j = v, value = NA_real_)

log_msg("各要素有效天数占比(中位) 与 站点缺失率:")
for (v in c("tavg", "tmax", "rh", "cloud", "precip"))
  log_msg(sprintf("  %-7s 有效天占比中位 %.3f | 该列全缺的站 %d",
                  v, median(res[[paste0(v, "_n")]], na.rm = TRUE), sum(is.na(res[[v]]))))

# 有效天数过少的站, 均值不可靠, 置为缺失
for (v in c("tavg", "tmax", "rh", "cloud", "precip"))
  set(res, i = which(res[[paste0(v, "_n")]] < 0.5), j = v, value = NA_real_)

fwrite(res, file.path(OUT, "meteo_station_agg.csv"))
log_msg("写出 meteo_station_agg.csv")
log_msg("量级检查: 气温 ", paste(round(range(res$tavg, na.rm = TRUE), 1), collapse = "-"),
        " C | 湿度 ", paste(round(range(res$rh, na.rm = TRUE), 1), collapse = "-"), " %")
