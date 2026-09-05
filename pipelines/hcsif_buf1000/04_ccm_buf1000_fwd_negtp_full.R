#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04_ccm_buf1000_fwd_negtp_full.R
#
# 全量重跑 fwd(VPD->SIF) 方向的 CCM，把 tp 扫描范围从主线的 0..8 扩到 -8..8。
#
# 背景(见 docs/02_ccm.md "确认 VPD→SIF 方向判据" 相关条目)：判断"方向是否确认"
# 需要两步 —— ①存在性(收敛+p_surr<0.05) ②optimal tp(显著tp中rho最高者)必须≥0。
# 第②步要求知道"全部候选滞后里 rho 最高的到底是不是负的"，主线 01_ccm_buf1000.R
# 从未测过负 tp，这一步在全量上做不了。本脚本补上这个缺口。
#
# 不测反向(SIF->VPD)——2026-09-04 已确认不需要(见 docs/02_ccm.md)。
# 不做多重比较的形式化校正——已查证 Ye 2015/Ushio 2018/Guo 2024 三篇参考文献都
# 没有对"扫描多个tp"做形式化校正，这个领域本身没有公认解法(见 docs/02_ccm.md)。
#
# 方法(E选取/CCM/S-map/替代检验/PREPROC=normalize)与 01_ccm_buf1000.R、
# 02_ccm_buf1000_dirscan.R 完全一致，只是站点范围改回全量、方向只测 fwd、
# tp 范围扩到 -8..8。用的是 rEDM 2.0.2(新API，见"环境坑")。
#
# 用法: Rscript 04_ccm_buf1000_fwd_negtp_full.R
#       HCSIF_TEST=1 Rscript ... 冒烟测试(小站点/tp子集)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
BUF_R    <- 1000L
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v3")
OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_fwd_negtp_full")
SIF_COL  <- sprintf("SIF_buf%d", BUF_R)
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")

TP_SEQ        <- -8:8
MIN_PTS       <- 40
MIN_DAYS_WIN  <- 5
N_CORES       <- max(1L, parallel::detectCores() - 2L)
PREPROC       <- "normalize"          # 与 norm_surr 主线口径一致
SEED_BASE     <- 20260904L
N_SURR        <- as.integer(Sys.getenv("HCSIF_N_SURR", "199"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ===========================================================================
# 1. HCSIF 全量站点序列 + VPD(复用缓存，与主线完全一致)
# ===========================================================================
fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有 HCSIF 提取结果")
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[, c("meteo_stat", "year", "doy", "date", SIF_COL), with = FALSE]
log_msg("HCSIF: ", nrow(sif), " 条, ", uniqueN(sif$meteo_stat), " 站")

if (!file.exists(CACHE)) stop("缺少 VPD 缓存: ", CACHE, " —— 先跑一次 01_ccm_buf1000.R 生成")
vpd <- readRDS(CACHE)
vpd <- vpd[n_days >= MIN_DAYS_WIN]
log_msg("VPD 8天序列(缓存): ", nrow(vpd), " 条")

# ===========================================================================
# 2. 合并 + 归一化预处理(与主线 znorm 完全一致)
# ===========================================================================
d <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
log_msg("SIF-VPD 匹配后: ", nrow(d), " 条, ", uniqueN(d$meteo_stat), " 站")

znorm <- function(x) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
d[, idx8 := as.integer(round(as.numeric(
  as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]
for (v in c(SIF_COL, "vpd_mean")) d[, paste0(v, "_dt") := znorm(get(v)), by = meteo_stat]

SIF_V <- paste0(SIF_COL, "_dt"); VPD_V <- "vpd_mean_dt"
STATIONS <- sort(unique(d$meteo_stat))

# 排除疑似地表覆盖断点站(同 01_ccm_buf1000.R；见 00_landcover_break_buf1000.R)
LC_BREAK_F <- file.path(PROJ, "data_proc/output_hcsif_buf1000/landcover_break_stations.csv")
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F)
  excl <- lcb[has_break == TRUE]$stat_id
  n0 <- length(STATIONS)
  STATIONS <- setdiff(STATIONS, excl)
  log_msg("地表覆盖断点过滤: 候选 ", length(excl), " 站, 命中 ", n0 - length(STATIONS),
          " 站, 排除后剩 ", length(STATIONS), " 站")
}

if (Sys.getenv("HCSIF_TEST", "0") == "1") {
  STATIONS <- head(STATIONS, 6)
  TP_SEQ   <- c(-2, 0, 2)
  OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_fwd_negtp_full_TEST")
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
}
log_msg("站点数: ", length(STATIONS), " | tp 范围: ", min(TP_SEQ), "..", max(TP_SEQ))

# ===========================================================================
# 3. CCM: fwd 方向(y=SIF, x=VPD), signed tp: >=0 为"x 滞后 tp 步"(过去驱动现在)，
#    <0 为"x 领先 |tp| 步"(用来判断 optimal tp 是否落在时间顺序不成立的负区间)
# ===========================================================================
lag_or_lead <- function(x, tp) {
  if (tp >= 0) shift(x, n = tp,  type = "lag")
  else         shift(x, n = -tp, type = "lead")
}

run_ccm <- function(sid, tp_x) {
  # 单站超时保护，见 02_ccm_buf1000_dirscan.R 的教训(numProcess 与 mclapply
  # fork 混用曾观察到整批停滞)：全部显式 numProcess=1 + 90s 硬超时兜底。
  setTimeLimit(elapsed = 90, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  dd <- d[meteo_stat == sid]
  dd[, x_shift := lag_or_lead(get(VPD_V), tp_x), by = year]

  dd <- dd[!is.na(idx8)]
  if (!nrow(dd)) return(NULL)
  tmp  <- dd[, .(idx8, y = get(SIF_V), x = x_shift)]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)

  n_valid_pairs <- full[!is.na(y) & !is.na(x), .N]
  if (n_valid_pairs < MIN_PTS) return(NULL)

  df <- data.frame(time = seq_len(nrow(full)), resp = full$y, driv = full$x)
  n <- nrow(df)

  set.seed(SEED_BASE + as.integer(sid) + tp_x * 1000L)

  tryCatch({
    ed <- rEDM::EmbedDimension(dataFrame = df, columns = "resp", target = "resp",
                               lib = paste("1", n), pred = paste("1", n),
                               maxE = max(2, min(8, floor(n_valid_pairs / 10))),
                               numProcess = 1, showPlot = FALSE)
    best_E <- ed$E[which.max(ed$rho)]

    lib_max <- min(n - best_E, n_valid_pairs)
    cm <- rEDM::CCM(dataFrame = df, E = best_E, Tp = 0,
                    columns = "driv", target = "resp",
                    libSizes = paste(best_E + 2, lib_max,
                                     max(2, floor((lib_max - best_E - 2) / 10))),
                    sample = 50, seed = SEED_BASE + as.integer(sid) + tp_x * 1000L,
                    numProcess = 1, showPlot = FALSE)
    cs <- as.data.table(cm)[, .(rho_m = mean(`driv:resp`, na.rm = TRUE)), by = LibSize]
    setorder(cs, LibSize)

    sm <- rEDM::SMap(dataFrame = df, E = best_E, theta = 2,
                     lib = paste("1", n), pred = paste("1", n),
                     columns = "driv", target = "resp", embedded = FALSE)
    cc <- which(grepl("driv", colnames(sm$coefficients), ignore.case = TRUE))[1]
    if (is.na(cc)) cc <- 2L

    xmap_lmax <- function(driv_vec, surr_seed) {
      dfa <- data.frame(time = seq_len(n), resp = full$y, driv = driv_vec)
      cm2 <- tryCatch(rEDM::CCM(dataFrame = dfa, E = best_E, Tp = 0,
                       columns = "driv", target = "resp",
                       libSizes = paste(lib_max, lib_max, 1),
                       sample = 1, seed = surr_seed,
                       numProcess = 1, showPlot = FALSE),
                     error = function(e) NULL)
      if (is.null(cm2)) return(NA_real_)
      mean(as.data.table(cm2)[["driv:resp"]], na.rm = TRUE)
    }
    p_surr <- NA_real_; rho_null95 <- NA_real_
    base_seed <- SEED_BASE + as.integer(sid) + tp_x * 1000L
    if (N_SURR > 0) {
    rho_obs_s <- xmap_lmax(full$x, base_seed)

    doys <- sort(unique(dd$doy)); yrs <- sort(unique(dd$year))
    if (length(yrs) >= 3 && is.finite(rho_obs_s)) {
      di <- match(dd$doy, doys); yi <- match(dd$year, yrs)
      M  <- matrix(NA_real_, length(doys), length(yrs))
      M[cbind(di, yi)] <- dd$x_shift
      pos <- match(dd$idx8, full$idx8)
      rho_null <- rep(NA_real_, N_SURR)
      for (s in seq_len(N_SURR)) {
        src <- sample(length(yrs))
        driv_full <- rep(NA_real_, n)
        driv_full[pos] <- M[cbind(di, src[yi])]
        rho_null[s] <- xmap_lmax(driv_full, base_seed + s)
      }
      nv <- sum(!is.na(rho_null))
      if (nv > 0) {
        p_surr     <- (1 + sum(rho_null >= rho_obs_s, na.rm = TRUE)) / (1 + nv)
        rho_null95 <- as.numeric(quantile(rho_null, 0.95, na.rm = TRUE))
      }
    }
    }

    data.table(meteo_stat = sid, tp = tp_x,
               E = best_E, n_obs = n_valid_pairs, n_grid = n,
               rho       = cs$rho_m[nrow(cs)],
               rho_min   = cs$rho_m[1],
               trend     = cor(cs$LibSize, cs$rho_m),
               mean_coef = mean(sm$coefficients[, cc], na.rm = TRUE),
               p_surr    = p_surr,
               rho_null95 = rho_null95)
  }, error = function(e) NULL)
}

log_msg("CCM 规模: ", length(STATIONS), " 站 x ", length(TP_SEQ), " 个 tp = ",
        length(STATIONS) * length(TP_SEQ), " 次")

safe_io <- function(expr, what, tries = 6) {
  for (k in seq_len(tries)) {
    r <- tryCatch(force(expr), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    log_msg("IO失败[", what, "]: ", conditionMessage(r), " —— 重试 ", k, "/", tries)
    Sys.sleep(3)
  }
  stop("多次重试仍失败: ", what)
}

PARTS_DIR <- file.path(OUT_DIR, "parts")
dir.create(PARTS_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
for (i in seq_along(TP_SEQ)) {
  tp_i <- TP_SEQ[i]
  part_f <- file.path(PARTS_DIR, sprintf("part_%02d.rds", i))
  if (file.exists(part_f)) { log_msg(sprintf("  [%2d/%2d] tp=%+d 已有检查点，跳过", i, length(TP_SEQ), tp_i)); next }
  r <- parallel::mclapply(STATIONS, function(s) run_ccm(s, tp_i), mc.cores = N_CORES)
  part <- rbindlist(r[!vapply(r, is.null, logical(1))])
  safe_io(saveRDS(part, part_f), paste0("checkpoint part_", i))
  log_msg(sprintf("  [%2d/%2d] tp=%+d : %d/%d 站成功 → %s (累计 %.1f 分钟)",
                  i, length(TP_SEQ), tp_i, nrow(part), length(STATIONS),
                  basename(part_f), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

parts <- sort(list.files(PARTS_DIR, pattern = "^part_[0-9]+\\.rds$", full.names = TRUE))
res <- rbindlist(lapply(parts, readRDS), fill = TRUE)
if (!nrow(res)) stop("所有组合均无成功结果，检查上面每行的失败日志")

stamp <- format(Sys.time(), "%Y%m%d_%H%M")
safe_io(saveRDS(res, file.path(OUT_DIR, sprintf("fwd_negtp_full_%s.rds", stamp))), "save rds")
safe_io(fwrite(res,  file.path(OUT_DIR, sprintf("fwd_negtp_full_%s.csv", stamp))), "save csv")
log_msg("写出: ", OUT_DIR, "/fwd_negtp_full_", stamp, ".{rds,csv}  (", nrow(res), " 行)")

cat("\n=== 判据核验(p_surr<0.05 且 Δρ>0；第二步要求 optimal tp>=0) ===\n")
res[, sig := !is.na(p_surr) & p_surr < 0.05 & (rho - rho_min) > 0]
step1 <- res[sig == TRUE, uniqueN(meteo_stat)]
opt <- res[sig == TRUE, .SD[which.max(rho)], by = meteo_stat][, .(meteo_stat, opt_tp = tp)]
step2 <- opt[opt_tp >= 0]
cat("站点总数:", uniqueN(res$meteo_stat), "\n")
cat("第一步(至少1个tp: p_surr<0.05且Δρ>0)通过:", step1, "站\n")
cat("第二步(+optimal tp>=0)通过:", nrow(step2), "站\n")
