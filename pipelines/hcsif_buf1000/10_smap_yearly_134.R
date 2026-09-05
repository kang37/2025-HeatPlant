#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_smap_yearly_134.R — 134站二变量S-map系数按年拆分
#
# 06_smap_bivar_134.R 把每站整个 2000-2022 的S-map系数序列汇总成一个方向
# (median_coef_2v/frac_positive_2v pooling全部年份)。用户想看这个方向是不是
# 每年都一样，还是像 smap_timevarying_coef.R 发现的单变量系数那样逐年翻转。
#
# 方法: 完全复用 06 的状态空间/E/theta/optimal_tp/随机种子(结果应与06逐点一致，
# 只是这里额外把每个S-map系数对应的日历年份保留下来，06为了省事在生成pooled
# 统计量后就丢弃了这层信息)，按 (meteo_stat, year) 汇总 median/frac_positive/n。
# "Overall"列直接复用06已经算好的 median_coef_2v/frac_positive_2v，不重新算，
# 保证与已发布结果完全一致。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
BUF_R    <- 1000L
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v3")
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")
SIF_COL  <- sprintf("SIF_buf%d", BUF_R)
OUT_DIR  <- file.path(PROJ, "data_proc/smap_bivar_134")

MIN_PTS      <- 40L   # 与06一致: 全序列最少点数门槛
MIN_PTS_YEAR <- 5L    # 新增: 单年最少点数门槛(生长季约19个8天窗口, 5是约1/4)
MIN_DAYS_WIN <- 5L
SEED_BASE    <- 20260905L
THETA        <- 2

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ===========================================================================
# 1. 134 站名单 + optimal tp + 已发布的 Overall 方向(不重算, 直接复用06结果)
# ===========================================================================
st134 <- readRDS(file.path(OUT_DIR, "smap_bivar_134.rds"))
st134 <- st134[!is.na(direction_2v)]
log_msg("134 站名单读取完成, 实际有方向结果: ", nrow(st134), " 站")

# ===========================================================================
# 2. SIF + VPD 合并, 归一化(与06一致)
# ===========================================================================
fs  <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[meteo_stat %in% st134$meteo_stat, c("meteo_stat", "year", "doy", "date", SIF_COL), with = FALSE]

vpd <- readRDS(CACHE)[n_days >= MIN_DAYS_WIN & meteo_stat %in% st134$meteo_stat]

d <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
d[, idx8 := as.integer(round(as.numeric(
  as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]

znorm <- function(x) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
for (v in c(SIF_COL, "vpd_mean")) d[, paste0(v, "_dt") := znorm(get(v)), by = meteo_stat]
SIF_V <- paste0(SIF_COL, "_dt"); VPD_V <- "vpd_mean_dt"

log_msg("SIF-VPD 匹配后: ", nrow(d), " 条, ", uniqueN(d$meteo_stat), " 站")

# ===========================================================================
# 3. 逐站二变量 S-map, 保留 (time -> idx8 -> year) 映射, 输出逐年系数
# ===========================================================================
bivar_smap_yearly <- function(sid, tp_x, theta = THETA) {
  dd <- d[meteo_stat == sid]
  dd[, x_lag := shift(get(VPD_V), n = tp_x, type = "lag"), by = year]
  dd <- dd[!is.na(idx8)]
  year_of_idx8 <- unique(dd[, .(idx8, year)])          # idx8 在本项目场景下与年份一一对应

  tmp  <- dd[, .(idx8, y = get(SIF_V), x = x_lag)]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)
  n_valid <- full[!is.na(y) & !is.na(x), .N]
  if (n_valid < MIN_PTS) return(NULL)

  full[, time := .I]
  full <- merge(full, year_of_idx8, by = "idx8", all.x = TRUE)   # 补回年份, 空档行year=NA
  setorder(full, time)

  df <- data.frame(time = full$time, sif = full$y, vpd = full$x)
  n <- nrow(df)
  set.seed(SEED_BASE + as.integer(sid) + tp_x * 1000L)

  sm <- tryCatch(rEDM::SMap(dataFrame = df, E = 2, theta = theta,
                             lib = paste("1", n), pred = paste("1", n),
                             columns = c("sif", "vpd"), target = "sif", embedded = TRUE),
                 error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  co <- as.data.table(sm$coefficients)
  cc <- which(grepl("vpd", colnames(co), ignore.case = TRUE))[1]
  if (is.na(cc)) return(NULL)

  # rEDM 2.0.2 的系数表列名是 "Time"(大写), 且行序本来就与 df 的行序(1:n)一致,
  # 直接用行号做 time 更稳妥(避免大小写/列名不存在导致静默取到 NULL)
  out <- data.table(time = seq_len(nrow(co)), coef = co[[cc]])
  out <- merge(out, full[, .(time, year)], by = "time")
  out <- out[!is.na(coef) & !is.nan(coef) & !is.na(year)]
  if (!nrow(out)) return(NULL)
  out[, .(meteo_stat = sid, time, year, coef)]
}

log_msg("开始逐站二变量 S-map(保留年份), 共 ", nrow(st134), " 站 ...")
coef_all <- rbindlist(lapply(seq_len(nrow(st134)), function(i) {
  r <- st134[i]
  out <- tryCatch(bivar_smap_yearly(r$meteo_stat, r$optimal_tp), error = function(e) NULL)
  if (i %% 20 == 0) log_msg("  ", i, "/", nrow(st134))
  out
}), fill = TRUE)
log_msg("完成, 系数点总数: ", nrow(coef_all), ", 覆盖站数: ", uniqueN(coef_all$meteo_stat))

# ===========================================================================
# 4. 按 (站, 年) 汇总: median/mean/frac_positive/n
# ===========================================================================
yearly <- coef_all[, .(n_pts = .N,
                        median_coef = median(coef),
                        mean_coef   = mean(coef),
                        frac_positive = mean(coef > 0)),
                    by = .(meteo_stat, year)]
yearly[, valid := n_pts >= MIN_PTS_YEAR]
log_msg("站x年 格子总数: ", nrow(yearly), "; 达到 n_pts>=", MIN_PTS_YEAR, " 门槛: ", sum(yearly$valid))

# 每站: 有多少年"方向"与Overall(06主线)相反, 作为一个直接的稳健性数字
dir_from_val <- function(v, thresh = 0.75) fifelse(v >= thresh, "Promote", fifelse(v <= 1 - thresh, "Inhibit", "Ambiguous"))
yearly[, direction_year := dir_from_val(frac_positive)]
chk <- merge(yearly[valid == TRUE], st134[, .(meteo_stat, direction_2v)], by = "meteo_stat")
chk[, opposite := (direction_2v == "Promote" & direction_year == "Inhibit") |
                   (direction_2v == "Inhibit" & direction_year == "Promote")]
n_st_with_opposite <- chk[, .(any_opp = any(opposite)), by = meteo_stat][any_opp == TRUE, .N]
log_msg("在 Promote/Inhibit 站里, 至少1年方向与Overall直接相反的站数: ", n_st_with_opposite,
        " / ", uniqueN(chk[direction_2v != "Ambiguous", meteo_stat]))

saveRDS(list(coef_all = coef_all, yearly = yearly), file.path(OUT_DIR, "smap_yearly_134.rds"))
fwrite(yearly, file.path(OUT_DIR, "smap_yearly_134.csv"))
log_msg("已写出 ", file.path(OUT_DIR, "smap_yearly_134.{rds,csv}"))
