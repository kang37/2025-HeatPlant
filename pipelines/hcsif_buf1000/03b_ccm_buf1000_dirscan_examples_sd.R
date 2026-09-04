#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03b_ccm_buf1000_dirscan_examples_sd.R
#
# 02_ccm_buf1000_dirscan.R 的精简变体：只在画图用的 9 个示例站上重跑，
# 唯一区别是额外保存 rho_sd —— 50 次随机子抽样(sample=50)在最大库长处的
# 标准差，用来给 fig_dirscan_examples.png 画误差棒(参考 Ye et al. 2015 的画法；
# 主线 02 脚本没存这个量，403_fig 目前没有误差棒)。
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
OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_dirscan_examples")
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

# --- 9 个示例站(2026-09-04，从 dirscan_buf1000_20260904_0046 结果里按类型挑选) ---
# 每类 3 站，用于 fig_dirscan_examples.png：
EX_CLEAN  <- c(57612, 57517, 58467)   # 清晰方向型: tp=0 处仅 fwd 显著
EX_AMBIG  <- c(59238, 57273, 54818)   # 双向可疑型: tp=0 处 fwd 与 rev 都显著
EX_NULL   <- c(57171, 59493, 57206)   # 无信号对照: nsig=0 或 1
STATIONS  <- sort(unique(c(EX_CLEAN, EX_AMBIG, EX_NULL)))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
log_msg("示例站: ", length(STATIONS), " 站 (clean=", length(EX_CLEAN),
        " ambig=", length(EX_AMBIG), " null=", length(EX_NULL), ")")

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
    # includeData=TRUE 才会附带 `driv:resp_var` 逐库长的抽样方差列——
    # rEDM 2.0.2 的 CCM() 每个库长只回一行(已经是50次子抽样的均值)，不像旧版
    # 那样能拿到50次各自的原始值自己算sd，必须显式要它把方差带出来。
    cm <- rEDM::CCM(dataFrame = df, E = best_E, Tp = 0,
                    columns = "driv", target = "resp",
                    libSizes = paste(best_E + 2, lib_max,
                                     max(2, floor((lib_max - best_E - 2) / 10))),
                    sample = 50, seed = SEED_BASE + as.integer(sid) + tp_x * 1000L,
                    numProcess = 1, includeData = TRUE, showPlot = FALSE)
    cs <- as.data.table(cm)[, .(rho_m = mean(`driv:resp`, na.rm = TRUE),
                                 rho_var = mean(`driv:resp_var`, na.rm = TRUE)), by = LibSize]
    setorder(cs, LibSize)

    # --- Ye et al. 2015 风格误差棒：整条 tp 曲线固定用同一个库长(60% lib_max)，
    # 点估计和误差棒来自同一批 50 次抽样，不是像上面那样"点估计用最大库长、
    # sd 从别的库长借"的拼凑做法。只用于画图，不改主线的 rho/rho_min/p_surr 定义
    # 和显著性判据(sig 仍然基于上面 cs 的最大库长)。
    lib_ye <- max(best_E + 2, round(0.6 * lib_max))
    cm_ye <- rEDM::CCM(dataFrame = df, E = best_E, Tp = 0,
                       columns = "driv", target = "resp",
                       libSizes = paste(lib_ye, lib_ye, 1),
                       sample = 50, seed = SEED_BASE + as.integer(sid) + tp_x * 1000L,
                       numProcess = 1, includeData = TRUE, showPlot = FALSE)
    cs_ye <- as.data.table(cm_ye)

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
               # 库长=全部有效点数时子抽样没有随机性(只有一种取法)，方差恒为0。
               # 用方差>0的最大库长处的 sd 作为局部抽样不确定性的近似，只用于
               # 画误差棒，不影响 rho 点估计本身(仍取自真正的最大库长)。
               rho_sd    = { v <- cs[rho_var > 0][order(-LibSize)]
                             if (nrow(v)) sqrt(v$rho_var[1]) else NA_real_ },
               # Ye et al. 2015 风格：点估计和误差棒来自同一个固定库长(60% lib_max)
               # 的同一批 50 次抽样，跟上面 rho/rho_sd 那套拼凑法是两回事。
               lib_ye    = lib_ye,
               rho_ye    = mean(cs_ye$`driv:resp`, na.rm = TRUE),
               rho_ye_sd = sqrt(mean(cs_ye$`driv:resp_var`, na.rm = TRUE)),
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
res[, example_type := fifelse(meteo_stat %in% EX_CLEAN, "clean",
                       fifelse(meteo_stat %in% EX_AMBIG, "ambiguous", "null"))]

stamp <- format(Sys.time(), "%Y%m%d_%H%M")
safe_io(saveRDS(res, file.path(OUT_DIR, sprintf("dirscan_buf1000_%s.rds", stamp))), "save rds")
safe_io(fwrite(res,  file.path(OUT_DIR, sprintf("dirscan_buf1000_%s.csv", stamp))), "save csv")
log_msg("写出: ", OUT_DIR, "/dirscan_buf1000_", stamp, ".{rds,csv}  (", nrow(res), " 行)")

cat("\n=== 各类型/各方向 显著 tp 覆盖情况(宽松判据 p_surr<0.1 且 rho-rho_min>0) ===\n")
res[, sig := !is.na(p_surr) & p_surr < 0.1 & (rho - rho_min) > 0]
print(res[, .(n_station = uniqueN(meteo_stat),
              n_sig_combo = sum(sig, na.rm = TRUE),
              pct_sig = round(100 * mean(sig, na.rm = TRUE), 1)),
          by = .(example_type, direction)][order(example_type, direction)])
