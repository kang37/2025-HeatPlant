#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 06_smap_bivar_134.R
#
# 对"CCM确认因果耦合"的站(2026-09-07起改用 12_ccm_causal_confirmation.R 统一
# 判据的产出，不再是旧的134站独立判据；脚本名沿用历史命名，实际站数以
# stations_confirmed.csv 为准)重新用二变量 S-map 判断因果方向，修正
# 00_cross_module_issues.md #1 指出的单变量 S-map 偏差(单变量只嵌入 VPD 自己，
# 未把 SIF 自身状态纳入状态空间；smap_2var_vs_1var.R 已验证二变量 [SIF,VPD]
# 与单变量符号一致率仅 63%)。
#
# 站点范围 + optimal tp: 直接读 12_ccm_causal_confirmation.R 的产出
# (data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv)，
# 不在本脚本里重复判据逻辑——判据改动只改12号脚本，这里自动跟着变。
#
# 方法: 状态空间 [SIF(t), VPD_lag(t)] (embedded=TRUE, E=2, theta=2 先固定，
# 稳健性见 07_smap_robustness_theta_multivar.R)，在每站的 optimal tp 上算，
# 读時变 ∂SIF/∂VPD 系数序列，用 median_coef + 主导符号占比(≥75%阈值，沿用
# smap_timevarying_coef.R 的规则)判方向: Promote / Inhibit / Ambiguous。
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
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_PTS      <- 40L
MIN_DAYS_WIN <- 5L
SEED_BASE    <- 20260905L
THETA        <- 2

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ===========================================================================
# 1. 站名单 + optimal tp: 直接读 12_ccm_causal_confirmation.R 的统一判据产出
# ===========================================================================
CONFIRMED_F <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
stopifnot(file.exists(CONFIRMED_F))
st134 <- fread(CONFIRMED_F)[, .(meteo_stat, optimal_tp, optimal_rho, E_uni = E, nsig_0to8, nsig_all, tp_sig_list)]
stopifnot(nrow(st134) > 0)
log_msg("读取统一因果确认站名单 n=", nrow(st134), "(脚本/文件名沿用历史 134 命名，实际站数以此为准)，",
        "optimal_tp 范围: ", min(st134$optimal_tp), "..", max(st134$optimal_tp))

# ===========================================================================
# 2. SIF + VPD 合并, 归一化(与主线 znorm 一致)
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
# 3. 逐站二变量 S-map: 状态空间 [sif, vpd_lag], E=2, theta=2, embedded=TRUE
# ===========================================================================
bivar_smap <- function(sid, tp_x, theta = THETA) {
  dd <- d[meteo_stat == sid]
  dd[, x_lag := shift(get(VPD_V), n = tp_x, type = "lag"), by = year]
  dd <- dd[!is.na(idx8)]
  tmp  <- dd[, .(idx8, y = get(SIF_V), x = x_lag)]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)
  n_valid <- full[!is.na(y) & !is.na(x), .N]
  if (n_valid < MIN_PTS) return(NULL)

  df <- data.frame(time = seq_len(nrow(full)), sif = full$y, vpd = full$x)
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
  cv <- co[[cc]]
  ok <- !is.na(cv) & !is.nan(cv)
  cv <- cv[ok]
  if (!length(cv)) return(NULL)

  # theta CV时同时需要预测技巧: 用 predictions 里 Observations/Predictions 的相关
  pr <- as.data.table(sm$predictions)
  pr_ok <- pr[!is.na(Observations) & !is.na(Predictions)]
  rho_pred <- if (nrow(pr_ok) >= 5) cor(pr_ok$Observations, pr_ok$Predictions) else NA_real_

  list(coef = cv, rho_pred = rho_pred, n_valid = n_valid)
}

direction_call <- function(cv, thresh = 0.75) {
  fp <- mean(cv > 0)
  if (fp >= thresh) "Promote" else if (fp <= (1 - thresh)) "Inhibit" else "Ambiguous"
}

log_msg("开始逐站二变量 S-map, ", nrow(st134), " 站 x optimal tp ...")
res <- rbindlist(lapply(seq_len(nrow(st134)), function(i) {
  r <- st134[i]
  out <- bivar_smap(r$meteo_stat, r$optimal_tp)
  if (is.null(out)) {
    return(data.table(meteo_stat = r$meteo_stat, n_coef_pts = 0L,
                       mean_coef_2v = NA_real_, median_coef_2v = NA_real_,
                       frac_positive_2v = NA_real_, direction_2v = NA_character_,
                       rho_pred_2v = NA_real_))
  }
  data.table(meteo_stat = r$meteo_stat, n_coef_pts = length(out$coef),
             mean_coef_2v = mean(out$coef), median_coef_2v = median(out$coef),
             frac_positive_2v = mean(out$coef > 0),
             direction_2v = direction_call(out$coef),
             rho_pred_2v = out$rho_pred)
}), fill = TRUE)

st134 <- merge(st134, res, by = "meteo_stat", all.x = TRUE)
setorder(st134, meteo_stat)

log_msg("完成: ", st134[!is.na(direction_2v), .N], "/", nrow(st134), " 站算出二变量方向")
print(st134[, .N, by = direction_2v])

saveRDS(st134, file.path(OUT_DIR, "smap_bivar_134.rds"))
fwrite(st134, file.path(OUT_DIR, "smap_bivar_134.csv"))
log_msg("已写出 ", file.path(OUT_DIR, "smap_bivar_134.{rds,csv}"))
