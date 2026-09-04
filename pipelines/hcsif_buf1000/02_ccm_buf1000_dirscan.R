#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02_ccm_buf1000_dirscan.R
#
# 诊断性检验：反向 CCM(SIF -> VPD) + 负 tp 扫描，用于排查"强耦合/同步导致的
# 双向伪因果"(Ye et al. 2015, Sci Rep)。只在分层抽样的站点上跑，不是全量重跑。
#
# 背景: 现有主线 01_ccm_buf1000.R 只测 VPD->SIF 一个方向、tp 只扫正值(0..8,
# 代表"VPD 领先 SIF")。这套设计无法排除"两变量强耦合到接近同步，导致双向
# CCM 都显著收敛"这种伪影(Ye et al. 2015 的核心论点)。
#
# 分层抽样(基于宽松判据 p_surr<0.1 且 rho-rho_min>0 下的显著 tp 个数 nsig,
# 数据源 data_proc/ccm_hcsif_buf1000_norm_surr):
#   Tier A: nsig>=4 的全部 8 站(同步伪影风险最高，全测)
#   Tier B: nsig in 2:3 的 115 站中随机抽 25 站
#   Tier C: nsig==1  的 260 站中随机抽 20 站(基线对照:若这层也测出同步伪影迹象,
#           说明"多tp显著=高风险"这个分层假设本身站不住脚，需要全测)
#
# 每站测两个方向 x tp=-8..8(共 17 个 tp，含原有 0..8 的重复以保证自洽性和
# 直接可比):
#   fwd: y=SIF_dt, x=VPD_dt   (与主线一致，但补上负 tp)
#   rev: y=VPD_dt, x=SIF_dt   (主线从未测过的反向)
#
# 预处理固定用 normalize(与产出 p_surr 的主线结果 norm_surr 口径一致，可比)。
# 方法(E选取/CCM/S-map/替代检验)与 01_ccm_buf1000.R 完全一致，唯一区别是
# 站点范围、tp 范围(含负值)、可测方向(fwd/rev)。
#
# 用法: Rscript 02_ccm_buf1000_dirscan.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
BUF_R    <- 1000L
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v3")
OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_dirscan")
SIF_COL  <- sprintf("SIF_buf%d", BUF_R)
METEO    <- file.path(PROJ, "data_raw/meteo_data_1961-2023")
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")

VPD_THRESHOLD <- 2.0
TP_SEQ        <- -8:8
MIN_PTS       <- 40
MIN_DAYS_WIN  <- 5
N_CORES       <- max(1L, parallel::detectCores() - 2L)
PREPROC       <- "normalize"          # 与 norm_surr 主线口径一致
SEED_BASE     <- 20260904L
N_SURR        <- as.integer(Sys.getenv("HCSIF_N_SURR", "199"))

# --- 三层抽样站点(2026-09-04 基于 ccm_hcsif_buf1000_norm_surr 的 nsig 分层，见对话记录) ---
TIER_A <- c(54503,57514,57517,57520,57612,57978,58337,58502)          # nsig>=4，全测
TIER_B <- c(50953,51431,51814,51829,53885,54092,54324,54398,54419,54568,
            54818,54830,56666,57178,57206,57273,57431,57814,58362,58467,
            59022,59134,59297,59470,59632)                             # nsig 2-3，抽25/115
TIER_C <- c(53582,53986,54076,54453,54514,54611,54702,54736,56187,56396,
            57171,57710,58158,58222,58566,58642,59081,59238,59452,59493) # nsig==1，抽20/260
STATIONS <- sort(unique(c(TIER_A, TIER_B, TIER_C)))

# 冒烟测试开关: HCSIF_DIRSCAN_TEST=1 时只用极小站点/tp范围快速验证脚本可跑通
if (Sys.getenv("HCSIF_DIRSCAN_TEST", "0") == "1") {
  STATIONS <- head(STATIONS, 2)
  TP_SEQ   <- c(-2, 0, 2)
  OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_dirscan_TEST")
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
log_msg("分层抽样站点: ", length(STATIONS), " 站 (TierA=", length(TIER_A),
        " TierB=", length(TIER_B), " TierC=", length(TIER_C), ")")

# ===========================================================================
# 1. HCSIF 站点序列 + VPD(复用缓存，与主线完全一致)
# ===========================================================================
fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有 HCSIF 提取结果")
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[meteo_stat %in% STATIONS, c("meteo_stat", "year", "doy", "date", SIF_COL), with = FALSE]
log_msg("HCSIF: ", nrow(sif), " 条, ", uniqueN(sif$meteo_stat), " 站")

if (!file.exists(CACHE)) stop("缺少 VPD 缓存: ", CACHE, " —— 先跑一次 01_ccm_buf1000.R 生成")
vpd <- readRDS(CACHE)
vpd <- vpd[meteo_stat %in% STATIONS & n_days >= MIN_DAYS_WIN]
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

# ===========================================================================
# 3. CCM(方向可切换: fwd = y:SIF,x:VPD ; rev = y:VPD,x:SIF)
#    signed tp: >=0 为"x 滞后 tp 步"(过去驱动现在)，<0 为"x 领先 |tp| 步"
#    (未来/领先驱动现在——用来对称地检验反向时间关系，Ye et al. 2015 的思路)
# ===========================================================================
lag_or_lead <- function(x, tp) {
  if (tp >= 0) shift(x, n = tp,  type = "lag")
  else         shift(x, n = -tp, type = "lead")
}

run_ccm <- function(sid, y_col, x_col, tp_x) {
  # 单站超时保护：rEDM 2.0.2 的内部线程池(numProcess)与外层 mclapply(fork)
  # 混用曾观察到整批卡死(CPU 时间不涨)，此处加硬超时 + 全部显式 numProcess=1
  # 双重保险，超时的站直接算失败(返回 NULL)，不拖累整批。
  setTimeLimit(elapsed = 90, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  dd <- d[meteo_stat == sid]
  dd[, x_shift := lag_or_lead(get(x_col), tp_x), by = year]

  dd <- dd[!is.na(idx8)]
  if (!nrow(dd)) return(NULL)
  tmp  <- dd[, .(idx8, y = get(y_col), x = x_shift)]
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

    data.table(meteo_stat = sid, y_var = y_col, x_var = x_col, tp = tp_x,
               E = best_E, n_obs = n_valid_pairs, n_grid = n,
               rho       = cs$rho_m[nrow(cs)],
               rho_min   = cs$rho_m[1],
               trend     = cor(cs$LibSize, cs$rho_m),
               mean_coef = mean(sm$coefficients[, cc], na.rm = TRUE),
               p_surr    = p_surr,
               rho_null95 = rho_null95)
  }, error = function(e) NULL)
}

grid <- rbindlist(list(
  data.table(y_col = SIF_V, x_col = VPD_V, tp = TP_SEQ, direction = "fwd"),  # VPD -> SIF
  data.table(y_col = VPD_V, x_col = SIF_V, tp = TP_SEQ, direction = "rev")   # SIF -> VPD
))
log_msg("CCM 规模: ", length(STATIONS), " 站 x ", nrow(grid), " 组合 = ",
        length(STATIONS) * nrow(grid), " 次")

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
for (i in seq_len(nrow(grid))) {
  part_f <- file.path(PARTS_DIR, sprintf("part_%02d.rds", i))
  if (file.exists(part_f)) { log_msg(sprintf("  [%2d/%2d] 已有检查点，跳过", i, nrow(grid))); next }
  g <- grid[i]
  r <- parallel::mclapply(STATIONS, function(s) run_ccm(s, g$y_col, g$x_col, g$tp),
                          mc.cores = N_CORES)
  part <- rbindlist(r[!vapply(r, is.null, logical(1))])
  if (nrow(part)) part[, direction := g$direction]
  safe_io(saveRDS(part, part_f), paste0("checkpoint part_", i))
  log_msg(sprintf("  [%2d/%2d] %s(%s) <- %s tp=%+d : %d/%d 站成功 → %s (累计 %.1f 分钟)",
                  i, nrow(grid), g$direction, g$y_col, g$x_col, g$tp,
                  nrow(part), length(STATIONS), basename(part_f),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

parts <- sort(list.files(PARTS_DIR, pattern = "^part_[0-9]+\\.rds$", full.names = TRUE))
res <- rbindlist(lapply(parts, readRDS), fill = TRUE)
if (!nrow(res)) stop("所有组合均无成功结果，检查上面每行的失败日志")

# 附加分层标签，方便后续按层汇总
res[, tier := fifelse(meteo_stat %in% TIER_A, "A_multi_tp",
              fifelse(meteo_stat %in% TIER_B, "B_2to3_tp", "C_single_tp"))]

stamp <- format(Sys.time(), "%Y%m%d_%H%M")
safe_io(saveRDS(res, file.path(OUT_DIR, sprintf("dirscan_buf1000_%s.rds", stamp))), "save rds")
safe_io(fwrite(res,  file.path(OUT_DIR, sprintf("dirscan_buf1000_%s.csv", stamp))), "save csv")
log_msg("写出: ", OUT_DIR, "/dirscan_buf1000_", stamp, ".{rds,csv}  (", nrow(res), " 行)")

cat("\n=== 各层/各方向 显著 tp 覆盖情况(宽松判据 p_surr<0.1 且 rho-rho_min>0) ===\n")
res[, sig := !is.na(p_surr) & p_surr < 0.1 & (rho - rho_min) > 0]
print(res[, .(n_station = uniqueN(meteo_stat),
              n_sig_combo = sum(sig, na.rm = TRUE),
              pct_sig = round(100 * mean(sig, na.rm = TRUE), 1)),
          by = .(tier, direction)][order(tier, direction)])
