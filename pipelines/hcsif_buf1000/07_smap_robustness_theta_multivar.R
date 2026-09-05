#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 07_smap_robustness_theta_multivar.R
#
# 06_smap_bivar_134.R 的两项稳健性检查，抽样约 30 站(方向分层抽样，
# 而非随机——重点是压力测试各方向类别是否稳定):
#
#   (A) theta: 06 脚本固定 theta=2, 这里对每个抽样站在 optimal tp 上扫描
#       theta 网格(0,0.5,1,2,3,4,6,8), 用样本内预测技巧(cor(Obs,Pred))选
#       theta_cv=argmax, 比较 theta=2 vs theta_cv 的方向判断是否一致。
#
#   (B) 状态空间维度: 二变量 [sif,vpd_lag] vs 三变量 [sif,vpd_lag,tavg_lag]
#       vs 四变量 [sif,vpd_lag,tavg_lag,rsds_lag](温度=站点气象日值8天均,
#       辐射=ERA5-Land era5_rsds_sif8d_<year>.csv, 均按与vpd相同的optimal tp
#       滞后), 比较系数方向是否随维度增加而改变。
#
# 只是稳健性抽查，不替换 06 的主线结果；若这里发现大比例翻转，需要回头改主线。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
BUF_R    <- 1000L
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v3")
METEO    <- file.path(PROJ, "data_raw/meteo_data_1961-2023")
ERA5_DIR <- file.path(PROJ, "data_raw/covariates_1km")
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")
SIF_COL  <- sprintf("SIF_buf%d", BUF_R)
OUT_DIR  <- file.path(PROJ, "data_proc/smap_bivar_134")

MIN_PTS      <- 40L
MIN_DAYS_WIN <- 5L
SEED_BASE    <- 20260905L
THETA_GRID   <- c(0, 0.5, 1, 2, 3, 4, 6, 8)
N_PER_GROUP  <- 10L

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ===========================================================================
# 1. 分层抽样 ~30 站(按 06 的方向分组)
# ===========================================================================
st134 <- readRDS(file.path(OUT_DIR, "smap_bivar_134.rds"))
set.seed(20260905L)
sub <- st134[!is.na(direction_2v), .SD[sample(.N, min(.N, N_PER_GROUP))], by = direction_2v]
sub <- sub[, .(meteo_stat, optimal_tp, direction_2v, mean_coef_2v, median_coef_2v, frac_positive_2v)]
log_msg("抽样: ", nrow(sub), " 站 (", paste(sub[, .N, by = direction_2v][, paste0(direction_2v, "=", N)], collapse = ", "), ")")

# ===========================================================================
# 2. 载入 SIF/VPD/温度/辐射, 归一化
# ===========================================================================
ids <- sub$meteo_stat
fs  <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[meteo_stat %in% ids, c("meteo_stat", "year", "doy", "date", SIF_COL), with = FALSE]

vpd <- readRDS(CACHE)[n_days >= MIN_DAYS_WIN & meteo_stat %in% ids]

# --- 温度: 站点气象日值 -> 8天窗口均值(与 vpd_8day_cache 同一套窗口规则) ----
build_tavg <- function(ids) {
  fl <- file.path(METEO, paste0(ids, ".txt"))
  fl <- fl[file.exists(fl)]
  win_doys <- sort(unique(sif$doy))
  one <- function(p) {
    sid <- as.integer(sub("\\.txt$", "", basename(p)))
    dd <- tryCatch(fread(p, skip = 1, showProgress = FALSE), error = function(e) NULL)
    if (is.null(dd) || !all(c("date", "tavg") %in% names(dd))) return(NULL)
    dd <- dd[, .(date, tavg)]
    dd[tavg >= 999990, tavg := NA_real_]
    dd[, date := as.Date(date)]
    dd[, doy := as.integer(format(date, "%j"))]
    dd[, win := 60L + 8L * ((doy - 60L) %/% 8L)]
    dd <- dd[win %in% win_doys & !is.na(tavg)]
    if (!nrow(dd)) return(NULL)
    dd[, .(meteo_stat = sid, tavg_mean = mean(tavg), n_days_t = .N),
       by = .(year = as.integer(format(date, "%Y")), doy = win)]
  }
  rbindlist(lapply(fl, one), fill = TRUE)
}
tavg <- build_tavg(ids)
log_msg("温度8天序列: ", nrow(tavg), " 条, ", uniqueN(tavg$meteo_stat), " 站")

# --- 辐射: ERA5-Land era5_rsds_sif8d_<year>.csv, 已是与 SIF 同网格的8天合成 --
yrs_needed <- 2000:2022
rsds_files <- file.path(ERA5_DIR, sprintf("era5_rsds_sif8d_%d.csv", yrs_needed))
rsds_files <- rsds_files[file.exists(rsds_files)]
rsds <- rbindlist(lapply(rsds_files, fread), fill = TRUE)
rsds <- rsds[stat_id %in% ids, .(meteo_stat = stat_id, year, doy, rsds_mean)]
log_msg("辐射8天序列: ", nrow(rsds), " 条, ", uniqueN(rsds$meteo_stat), " 站")

d <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
d <- merge(d, tavg, by = c("meteo_stat", "year", "doy"), all.x = TRUE)
d <- merge(d, rsds, by = c("meteo_stat", "year", "doy"), all.x = TRUE)
setorder(d, meteo_stat, year, doy)
d[, idx8 := as.integer(round(as.numeric(
  as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]

znorm <- function(x) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
for (v in c(SIF_COL, "vpd_mean", "tavg_mean", "rsds_mean")) {
  if (v %in% names(d)) d[, paste0(v, "_dt") := znorm(get(v)), by = meteo_stat]
}
SIF_V <- paste0(SIF_COL, "_dt"); VPD_V <- "vpd_mean_dt"
TAVG_V <- "tavg_mean_dt"; RSDS_V <- "rsds_mean_dt"
log_msg("四要素匹配后: ", nrow(d), " 条 (缺 tavg: ", d[is.na(get(TAVG_V)), .N],
        ", 缺 rsds: ", d[is.na(get(RSDS_V)), .N], ")")

# ===========================================================================
# 3. 通用: 构造某站某 tp 的状态空间数据框
# ===========================================================================
build_df <- function(sid, tp_x, extra_vars = character(0)) {
  extra_lag <- if (length(extra_vars)) paste0(extra_vars, "_lag") else character(0)
  dd <- d[meteo_stat == sid]
  dd[, x_lag := shift(get(VPD_V), n = tp_x, type = "lag"), by = year]
  for (v in extra_vars) dd[, paste0(v, "_lag") := shift(get(v), n = tp_x, type = "lag"), by = year]
  dd <- dd[!is.na(idx8)]
  tmp <- dd[, .(idx8, y = get(SIF_V), x = x_lag)]
  for (v in extra_vars) tmp[[paste0(v, "_lag")]] <- dd[[paste0(v, "_lag")]]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)
  keycols <- c("y", "x", extra_lag)
  n_valid <- sum(complete.cases(as.data.frame(full)[, keycols, drop = FALSE]))
  df <- data.frame(time = seq_len(nrow(full)), sif = full$y, vpd = full$x)
  for (v in extra_vars) df[[v]] <- full[[paste0(v, "_lag")]]
  list(df = df, n_valid = n_valid)
}

run_smap <- function(df, cols, theta, target = "sif") {
  n <- nrow(df)
  sm <- tryCatch(rEDM::SMap(dataFrame = df, E = length(cols), theta = theta,
                             lib = paste("1", n), pred = paste("1", n),
                             columns = cols, target = target, embedded = TRUE),
                 error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  co <- as.data.table(sm$coefficients)
  cc <- which(grepl("vpd", colnames(co), ignore.case = TRUE))[1]
  if (is.na(cc)) return(NULL)
  cv <- co[[cc]]; cv <- cv[!is.na(cv) & !is.nan(cv)]
  pr <- as.data.table(sm$predictions)
  pr_ok <- pr[!is.na(Observations) & !is.na(Predictions)]
  rho_pred <- if (nrow(pr_ok) >= 5) cor(pr_ok$Observations, pr_ok$Predictions) else NA_real_
  list(coef = cv, rho_pred = rho_pred)
}
direction_call <- function(cv, thresh = 0.75) {
  if (!length(cv)) return(NA_character_)
  fp <- mean(cv > 0)
  if (fp >= thresh) "Promote" else if (fp <= (1 - thresh)) "Inhibit" else "Ambiguous"
}

# ===========================================================================
# 4A. theta 网格搜索(二变量状态空间, 不同 theta 比预测技巧)
# ===========================================================================
log_msg("=== (A) theta 稳健性 ===")
theta_res <- rbindlist(lapply(seq_len(nrow(sub)), function(i) {
  r <- sub[i]
  bd <- build_df(r$meteo_stat, r$optimal_tp)
  if (bd$n_valid < MIN_PTS) return(NULL)
  scan <- rbindlist(lapply(THETA_GRID, function(th) {
    set.seed(SEED_BASE + as.integer(r$meteo_stat) + r$optimal_tp * 1000L + round(th * 10))
    out <- run_smap(bd$df, c("sif", "vpd"), th)
    if (is.null(out)) return(NULL)
    data.table(theta = th, rho_pred = out$rho_pred, median_coef = median(out$coef),
               direction = direction_call(out$coef))
  }))
  if (!nrow(scan)) return(NULL)
  best <- scan[which.max(rho_pred)]
  base <- scan[theta == 2]
  data.table(meteo_stat = r$meteo_stat, direction_theta2 = base$direction,
             theta_cv = best$theta, rho_pred_theta2 = base$rho_pred,
             rho_pred_cv = best$rho_pred, direction_cv = best$direction,
             flipped = !identical(base$direction, best$direction))
}))
print(theta_res)
n_flip <- sum(theta_res$flipped, na.rm = TRUE)
log_msg("theta=2 vs theta_cv 方向不一致: ", n_flip, "/", nrow(theta_res), " 站")

# ===========================================================================
# 4B. 多变量稳健性: 二/三/四变量
# ===========================================================================
log_msg("=== (B) 多变量稳健性 ===")
multivar_res <- rbindlist(lapply(seq_len(nrow(sub)), function(i) {
  r <- sub[i]
  out2 <- { bd <- build_df(r$meteo_stat, r$optimal_tp); if (bd$n_valid < MIN_PTS) NULL else
            { set.seed(SEED_BASE + as.integer(r$meteo_stat) + r$optimal_tp * 1000L)
              run_smap(bd$df, c("sif", "vpd"), THETA <- 2) } }
  out3 <- { bd <- build_df(r$meteo_stat, r$optimal_tp, extra_vars = TAVG_V); if (bd$n_valid < MIN_PTS) NULL else
            { set.seed(SEED_BASE + as.integer(r$meteo_stat) + r$optimal_tp * 1000L + 1L)
              run_smap(bd$df, c("sif", "vpd", TAVG_V), 2) } }
  out4 <- { bd <- build_df(r$meteo_stat, r$optimal_tp, extra_vars = c(TAVG_V, RSDS_V)); if (bd$n_valid < MIN_PTS) NULL else
            { set.seed(SEED_BASE + as.integer(r$meteo_stat) + r$optimal_tp * 1000L + 2L)
              run_smap(bd$df, c("sif", "vpd", TAVG_V, RSDS_V), 2) } }
  data.table(meteo_stat = r$meteo_stat, direction_baseline = r$direction_2v,
             direction_2v_bis = direction_call(out2$coef),
             direction_3v = direction_call(out3$coef),
             direction_4v = direction_call(out4$coef),
             n3 = if (is.null(out3)) 0L else length(out3$coef),
             n4 = if (is.null(out4)) 0L else length(out4$coef))
}))
print(multivar_res)
log_msg("2v(重算) vs 06主线 方向不一致: ",
        sum(multivar_res$direction_baseline != multivar_res$direction_2v_bis, na.rm = TRUE), "/", nrow(multivar_res))
log_msg("2v vs 3v 方向不一致: ",
        sum(multivar_res$direction_2v_bis != multivar_res$direction_3v, na.rm = TRUE), "/", nrow(multivar_res))
log_msg("2v vs 4v 方向不一致: ",
        sum(multivar_res$direction_2v_bis != multivar_res$direction_4v, na.rm = TRUE), "/",
        multivar_res[!is.na(direction_4v), .N])

saveRDS(list(theta = theta_res, multivar = multivar_res), file.path(OUT_DIR, "robustness_check.rds"))
fwrite(theta_res, file.path(OUT_DIR, "robustness_theta.csv"))
fwrite(multivar_res, file.path(OUT_DIR, "robustness_multivar.csv"))
log_msg("已写出稳健性检查结果到 ", OUT_DIR)
