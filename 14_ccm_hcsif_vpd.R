#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 14_ccm_hcsif_vpd.R
#
# 用 HCSIF (500 m, 8 天) 重做 VPD -> SIF 的 CCM 因果检验。
# 结果写入独立目录 data_proc/ccm_hcsif_500m/，不触碰任何既有结果。
#
# 与既有 01_ccm_vpd.R 的方法保持可比:
#   - 同样的 VPD 定义(Tetens 公式, 阈值 2.0 kPa)与热胁迫指标
#   - 同样的线性去趋势、EmbedDimension 选 E、CCM + 收敛趋势、SMap 取效应方向
#   - 同样的滞后扫描 tp = 0..8
#
# 三处必要的差异:
#   1. 时间步长由"周"改为 HCSIF 的 8 天合成窗口(DOY 60+8k)，VPD 按同窗口聚合
#   2. 因变量跑两套: SIF(站点像元点值) 与 SIF_buf750(750 m 邻域均值)
#      —— 点值空间精度高但城市站点缺失严重，邻域均值样本完整但混入周边植被
#   3. 滞后 lag 在年内进行(group_by(year))，避免跨越 10 月-次年 4 月的数据空档
#
# 用法:
#   Rscript 14_ccm_hcsif_vpd.R              # 全量
#   Rscript 14_ccm_hcsif_vpd.R --stations 60  # 抽样试跑
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station")
METEO    <- file.path(PROJ, "data_raw/meteo_data_1961-2023")
OUT_DIR  <- file.path(PROJ, "data_proc/ccm_hcsif_500m")
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")

VPD_THRESHOLD <- 2.0     # kPa，与既有分析一致
TP_SEQ        <- 0:8     # 滞后步数(每步 8 天)
MIN_PTS       <- 40      # CCM 最少样本量
MIN_DAYS_WIN  <- 5       # 一个 8 天窗口至少要有几天有效气象数据
N_CORES       <- max(1L, parallel::detectCores() - 2L)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

args     <- commandArgs(trailingOnly = TRUE)
n_sample <- if (length(args) >= 2 && args[1] == "--stations") as.integer(args[2]) else NA_integer_

# ===========================================================================
# 1. HCSIF 站点序列
# ===========================================================================
fs <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
if (!length(fs)) stop("没有 HCSIF 提取结果")
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[, .(meteo_stat, year, doy, date, SIF, SIF_buf750)]
log_msg("HCSIF: ", nrow(sif), " 条, ", uniqueN(sif$meteo_stat), " 站, ",
        uniqueN(sif$year), " 年, ", uniqueN(sif$doy), " 个 DOY 窗口")

win_doys <- sort(unique(sif$doy))

# ===========================================================================
# 2. 日尺度 VPD -> 8 天窗口聚合
# ===========================================================================
build_vpd <- function() {
  yrs  <- sort(unique(sif$year))
  ids  <- sort(unique(sif$meteo_stat))
  fl   <- file.path(METEO, paste0(ids, ".txt"))
  fl   <- fl[file.exists(fl)]
  log_msg("读取气象文件 ", length(fl), " 个 (缺失 ", length(ids) - length(fl), " 个站)")

  one <- function(p) {
    sid <- as.integer(sub("\\.txt$", "", basename(p)))
    d <- tryCatch(fread(p, skip = 1, showProgress = FALSE), error = function(e) NULL)
    if (is.null(d) || !all(c("date", "tavg", "RH") %in% names(d))) return(NULL)
    d <- d[, .(date, tavg, RH)]
    # 哨兵值 999999 -> NA
    d[tavg >= 999990, tavg := NA_real_]
    d[RH   >= 999990, RH   := NA_real_]
    d[, date := as.Date(date)]
    d[, year := as.integer(format(date, "%Y"))]
    d <- d[year %in% yrs]
    if (!nrow(d)) return(NULL)
    d[, doy := as.integer(format(date, "%j"))]
    # 与 HCSIF 相同的 8 天窗口: 起始 DOY = 60 + 8k
    d[, win := 60L + 8L * ((doy - 60L) %/% 8L)]
    d <- d[win %in% win_doys]
    if (!nrow(d)) return(NULL)
    # Tetens 公式，与 _targets.R 中的定义完全一致
    d[, svp := 0.6112 * exp((17.67 * tavg) / (tavg + 243.5))]
    d[, avp := (RH / 100) * svp]
    d[, vpd := svp - avp]
    d[vpd < 0 | vpd > 10 | is.na(tavg) | is.na(RH), vpd := NA_real_]
    d[, `:=`(is_heat = vpd > VPD_THRESHOLD,
             heat_int = pmax(vpd - VPD_THRESHOLD, 0))]
    d[!is.na(vpd), .(
      meteo_stat    = sid,
      n_days        = .N,
      vpd_mean      = mean(vpd),
      heat_over_sum = sum(heat_int),
      heat_freq     = mean(is_heat)
    ), by = .(year, doy = win)]
  }

  res <- parallel::mclapply(fl, one, mc.cores = N_CORES)
  rbindlist(res[!vapply(res, is.null, logical(1))])
}

if (file.exists(CACHE)) {
  vpd <- readRDS(CACHE); log_msg("VPD 8 天序列: 读取缓存 ", nrow(vpd), " 条")
} else {
  log_msg("构建 VPD 8 天序列(首次运行较慢)...")
  vpd <- build_vpd()
  saveRDS(vpd, CACHE)
  log_msg("VPD 8 天序列: ", nrow(vpd), " 条, 已缓存")
}
vpd <- vpd[n_days >= MIN_DAYS_WIN]

# ===========================================================================
# 3. 合并 + 去趋势
# ===========================================================================
d <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
log_msg("SIF-VPD 匹配后: ", nrow(d), " 条, ", uniqueN(d$meteo_stat), " 站")

safe_detrend <- function(x, t) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch(as.numeric(residuals(lm(x ~ t, na.action = na.exclude))),
           error = function(e) rep(NA_real_, length(x)))
}

# 真实时间坐标: 以 8 天为单位的绝对序号。
# 用它(而非行号)铺时间轴，才能让每年 10 月-次年 4 月的休眠期表现为真实的空档。
d[, idx8 := as.integer(round(as.numeric(
  as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]

# 去趋势用真实时间坐标，避免缺失把行号和实际间隔的对应关系拉偏
for (v in c("SIF", "SIF_buf750", "vpd_mean", "heat_over_sum")) {
  d[, paste0(v, "_dt") := safe_detrend(get(v), idx8), by = meteo_stat]
}

# ===========================================================================
# 4. CCM
# ===========================================================================
run_ccm <- function(sid, y_col, x_col, tp_x) {
  dd <- d[meteo_stat == sid]
  # 滞后在年内进行，不跨越 10 月-次年 4 月的空档
  dd[, x_lag := shift(get(x_col), n = tp_x, type = "lag"), by = year]

  # --- 缺失值处理 ---------------------------------------------------------
  # 关键: 不能"删掉缺失行再重新编号 1..n"。CCM 靠时间延迟嵌入重建吸引子，
  # 相邻性就是它的全部信息来源；删行重编号会把不相邻的观测当成相邻的。
  # 合成实验(已知因果, 20% 成片缺失, n=12): 删行重编号使 rho 系统性低估
  # 0.069(相对 14%)，而保留 NA 的偏差仅 +0.005。
  #
  # 做法: 铺一条完整的 8 天规则时间轴，观测放到各自的真实位置上，
  # 其余(含每年 10 月-次年 4 月的休眠期)留 NA。跨越 NA 的嵌入自然失效，
  # 既不会伪造相邻性，也不会跨年份拼接。
  dd <- dd[!is.na(idx8)]
  if (!nrow(dd)) return(NULL)
  tmp  <- dd[, .(idx8, y = get(y_col), x = x_lag)]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)

  n_valid_pairs <- full[!is.na(y) & !is.na(x), .N]
  if (n_valid_pairs < MIN_PTS) return(NULL)

  df <- data.frame(time = seq_len(nrow(full)), sif = full$y, heat = full$x)
  n <- nrow(df)

  tryCatch({
    # maxE / libSizes 都按"实际有效配对数"定，而不是补齐后的序列长度 n,
    # 否则 NA 占位会把 E 和文库规模撑得虚高。
    ed <- rEDM::EmbedDimension(dataFrame = df, columns = "sif", target = "sif",
                               lib = paste("1", n), pred = paste("1", n),
                               maxE = max(2, min(8, floor(n_valid_pairs / 10))),
                               showPlot = FALSE)
    best_E <- ed$E[which.max(ed$rho)]

    lib_max <- min(n - best_E, n_valid_pairs)
    cm <- rEDM::CCM(dataFrame = df, E = best_E, Tp = 0,
                    columns = "heat", target = "sif",
                    libSizes = paste(best_E + 2, lib_max,
                                     max(2, floor((lib_max - best_E - 2) / 10))),
                    sample = 50, random = TRUE, showPlot = FALSE)
    cs <- as.data.table(cm)[, .(rho_m = mean(`heat:sif`, na.rm = TRUE)), by = LibSize]
    setorder(cs, LibSize)

    sm <- rEDM::SMap(dataFrame = df, E = best_E, theta = 2,
                     lib = paste("1", n), pred = paste("1", n),
                     columns = "heat", target = "sif", embedded = FALSE)
    cc <- which(grepl("heat", colnames(sm$coefficients), ignore.case = TRUE))[1]
    if (is.na(cc)) cc <- 2L

    data.table(meteo_stat = sid, y_var = y_col, x_var = x_col, tp = tp_x,
               E = best_E, n_obs = n_valid_pairs, n_grid = n,
               rho       = cs$rho_m[nrow(cs)],
               rho_min   = cs$rho_m[1],
               trend     = cor(cs$LibSize, cs$rho_m),
               mean_coef = mean(sm$coefficients[, cc], na.rm = TRUE))
  }, error = function(e) NULL)
}

stations <- sort(unique(d$meteo_stat))
if (!is.na(n_sample) && n_sample < length(stations)) {
  set.seed(42); stations <- sort(sample(stations, n_sample))
  log_msg("抽样试跑: ", length(stations), " 站")
}

grid <- CJ(y_col = c("SIF_dt", "SIF_buf750_dt"),
           x_col = c("vpd_mean_dt"),          # 仅 VPD->SIF（热胁迫驱动本次不算）
           tp    = TP_SEQ, sorted = FALSE)
log_msg("CCM 规模: ", length(stations), " 站 x ", nrow(grid), " 组合 = ",
        length(stations) * nrow(grid), " 次")

# 防御式读取/写入：Dropbox 文件偶发 EINTR(中断的系统调用)，重试
safe_io <- function(expr, what, tries = 6) {
  for (k in seq_len(tries)) {
    r <- tryCatch(force(expr), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    log_msg("IO失败[", what, "]: ", conditionMessage(r), " —— 重试 ", k, "/", tries)
    Sys.sleep(3)
  }
  stop("多次重试仍失败: ", what)
}

# 检查点目录：每个组合算完立即落盘，崩溃不丢；重跑自动断点续算
PARTS_DIR <- file.path(OUT_DIR, "parts")
dir.create(PARTS_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  part_f <- file.path(PARTS_DIR, sprintf("part_%02d.rds", i))
  if (file.exists(part_f)) {                       # 断点续算：已完成组合跳过
    log_msg(sprintf("  [%2d/%2d] 已有检查点，跳过", i, nrow(grid)))
    next
  }
  g <- grid[i]
  r <- parallel::mclapply(stations, function(s) run_ccm(s, g$y_col, g$x_col, g$tp),
                          mc.cores = N_CORES)
  part <- rbindlist(r[!vapply(r, is.null, logical(1))])
  safe_io(saveRDS(part, part_f), paste0("checkpoint part_", i))  # ★ 立即落盘
  log_msg(sprintf("  [%2d/%2d] %s <- %s tp=%d : %d 站成功 → %s (累计 %.1f 分钟)",
                  i, nrow(grid), g$y_col, g$x_col, g$tp, nrow(part),
                  basename(part_f), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# 从检查点汇总（即使某次崩溃，重跑续算后仍能在此完整拼装）
parts <- sort(list.files(PARTS_DIR, pattern = "^part_[0-9]+\\.rds$", full.names = TRUE))
res <- rbindlist(lapply(parts, readRDS), fill = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M")
# ★ 先保存原始结果(不含meta)，确保 8 小时计算成果一定落盘
safe_io(saveRDS(res, file.path(OUT_DIR, paste0("ccm_hcsif_vpd_", stamp, "_raw.rds"))),
        "save raw rds")

# 附加坐标(防御式读取 + 容错：合并失败也不影响已保存的原始结果)
res <- tryCatch({
  meta <- unique(safe_io(fread(file.path(PROJ, "data_raw/hcsif/stations_924.csv")),
                         "read stations_924"),
                 by = "meteo_stat")
  merge(res, meta, by = "meteo_stat", all.x = TRUE)
}, error = function(e) { log_msg("meta合并失败: ", conditionMessage(e), " —— 用无坐标结果"); res })

safe_io(saveRDS(res, file.path(OUT_DIR, paste0("ccm_hcsif_vpd_", stamp, ".rds"))), "save rds")
safe_io(fwrite(res,  file.path(OUT_DIR, paste0("ccm_hcsif_vpd_", stamp, ".csv"))), "save csv")
log_msg("写出: ", OUT_DIR, "/ccm_hcsif_vpd_", stamp, ".{rds,csv}  (", nrow(res), " 行)")

cat("\n=== 结果速览 (rho > 0.1 且收敛 trend > 0 视为有因果证据) ===\n")
print(res[, .(n = .N,
              median_rho = round(median(rho, na.rm = TRUE), 3),
              causal_pct = round(100*mean(rho > 0.1 & trend > 0, na.rm = TRUE), 1)),
          by = .(y_var, x_var, tp)][order(y_var, x_var, tp)])
