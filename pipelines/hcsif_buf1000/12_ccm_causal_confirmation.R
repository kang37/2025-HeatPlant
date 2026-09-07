#!/usr/bin/env Rscript
# =============================================================================
# 12_ccm_causal_confirmation.R
#   项目里唯一的"CCM确认VPD->SIF因果耦合"判据，替代之前散落在各处、互不一致
#   的多套标准(主线宽松/严格显著、旧134站方向确认判据等)。三条标准依次筛选,
#   全部阈值/范围集中在下面, 改标准只改这几行, 不用碰下面的逻辑:
#     (1) 在 tp∈[TP_SIG_LO, TP_SIG_HI] 范围内, 存在某个 tp 同时满足
#         Δρ=rho-rho_min > DRHO_MIN  且  p_surr < PSURR_MAX
#     (2) 在 tp∈[TP_OPT_LO, TP_OPT_HI](更大范围, 含负tp) 里取全局
#         optimal_tp = argmax(rho)(不要求这个tp本身显著), 要求 optimal_tp >= TP_MIN
#   数据源必须是负tp扫描版本(tp=-8..8)，否则第(2)条会因为搜索范围本身不含
#   负值而恒真、变成空判据(教训见 docs/02_ccm.md 2026-09-04 日志)。
#   额外沿用旧 06_smap_bivar_134.R 的地表覆盖断点过滤(森林占比年际序列有结构性
#   跳变的站, CCM单吸引子假设不成立)——这条不在用户列的三条标准里, 是从旧脚本
#   继承的数据质量把关, 单独标注方便以后决定要不要保留。
#   产出: data_proc/ccm_hcsif_buf1000_causal_confirmed/
#     funnel.csv              每一层筛选剩多少站(方便以后改标准时对比)
#     stations_confirmed.csv  最终通过的站点清单(含optimal_tp/rho等)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000_causal_confirmed")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 阈值/范围：以后改标准只改这几行 ----------------------------------------
DRHO_MIN  <- 0        # 标准(1)之一: Δρ 严格大于此值
PSURR_MAX <- 0.05     # 标准(1)之一: p_surr 严格小于此值
TP_SIG_LO <- 0; TP_SIG_HI <- 8    # 标准(1)的 tp 搜索范围(正向滞后, 因果方向有意义的区间)
TP_OPT_LO <- -8; TP_OPT_HI <- 8   # 标准(2)算 optimal_tp 时的 tp 搜索范围(含负tp)
TP_MIN    <- 0        # 标准(2): optimal_tp >= 此值
EXCLUDE_LANDCOVER_BREAK <- TRUE   # 沿用旧脚本的地表覆盖断点过滤

SRC <- "data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv"
d <- fread(SRC)
stopifnot(min(d$tp) <= TP_OPT_LO, max(d$tp) >= TP_OPT_HI)
d[, drho := rho - rho_min]

n0 <- uniqueN(d$meteo_stat)

# ---- 可选: 地表覆盖断点过滤(继承自旧 06_smap_bivar_134.R) ------------------
LC_BREAK_F <- file.path(PROJ, "data_proc/output_hcsif_buf1000/landcover_break_stations.csv")
n0b <- n0
if (EXCLUDE_LANDCOVER_BREAK && file.exists(LC_BREAK_F)) {
  lcb  <- fread(LC_BREAK_F)
  excl <- lcb[has_break == TRUE]$stat_id
  d    <- d[!meteo_stat %in% excl]
  n0b  <- uniqueN(d$meteo_stat)
}

# ---- 标准(1)的 tp 搜索子集 --------------------------------------------------
sig_range <- d[tp >= TP_SIG_LO & tp <= TP_SIG_HI]
n_drho  <- sig_range[drho > DRHO_MIN, uniqueN(meteo_stat)]
n_psurr <- sig_range[p_surr < PSURR_MAX, uniqueN(meteo_stat)]
sig_range[, pass1 := (drho > DRHO_MIN) & (p_surr < PSURR_MAX)]
n1 <- sig_range[pass1 == TRUE, uniqueN(meteo_stat)]
stations_pass1 <- unique(sig_range[pass1 == TRUE]$meteo_stat)

# ---- 标准(2): 在更大范围(含负tp)里找全局 optimal_tp, 要求 >= TP_MIN --------
opt_range <- d[meteo_stat %in% stations_pass1 & tp >= TP_OPT_LO & tp <= TP_OPT_HI]
opt <- opt_range[order(meteo_stat, -rho)][, .SD[1], by = meteo_stat]
opt[, pass2 := tp >= TP_MIN]
n2 <- opt[pass2 == TRUE, .N]

final <- opt[pass2 == TRUE][, .(meteo_stat, optimal_tp = tp, optimal_rho = rho,
                                  optimal_drho = drho, optimal_p_surr = p_surr, E, n_obs)]

# 附加诊断列(不参与筛选, 供下游可视化用): tp∈[0,8]内有几个tp同时满足标准(1),
# 以及全部满足标准(1)的tp列表(含负tp, 供检查是否有可疑的负tp显著)
nsig_0to8 <- sig_range[pass1 == TRUE & tp %in% 0:8, .(nsig_0to8 = .N), by = meteo_stat]
nsig_all  <- sig_range[pass1 == TRUE, .(nsig_all = .N,
                                          tp_sig_list = paste(sort(tp), collapse = ";")), by = meteo_stat]
final <- merge(final, nsig_0to8, by = "meteo_stat", all.x = TRUE)
final <- merge(final, nsig_all,  by = "meteo_stat", all.x = TRUE)
final[is.na(nsig_0to8), nsig_0to8 := 0L]
setorder(final, meteo_stat)

# ---- 漏斗表 ------------------------------------------------------------------
funnel <- data.table(
  step = c("0_total_with_ccm_output", "0b_after_landcover_break_filter",
           "1a_drho_pass_in_0to8(info)", "1b_psurr_pass_in_0to8(info)",
           "1_drho_and_psurr_same_tp_in_0to8", "2_optimal_tp_in_neg8to8_pass"),
  n_stations = c(n0, n0b, n_drho, n_psurr, n1, n2),
  criterion  = c("有CCM输出(tp=-8..8全量)",
                 sprintf("排除地表覆盖断点站(开关=%s)", EXCLUDE_LANDCOVER_BREAK),
                 sprintf("[仅供参考,非gating] tp∈[%d,%d]内至少1个tp: drho>%.3g", TP_SIG_LO, TP_SIG_HI, DRHO_MIN),
                 sprintf("[仅供参考,非gating] tp∈[%d,%d]内至少1个tp: p_surr<%.3g", TP_SIG_LO, TP_SIG_HI, PSURR_MAX),
                 sprintf("tp∈[%d,%d]内至少1个tp同时满足 drho>%.3g 且 p_surr<%.3g", TP_SIG_LO, TP_SIG_HI, DRHO_MIN, PSURR_MAX),
                 sprintf("+ 在tp∈[%d,%d]范围内的全局optimal_tp(argmax rho) >= %.3g", TP_OPT_LO, TP_OPT_HI, TP_MIN))
)
cat("=== 因果确认漏斗 ===\n"); print(funnel)

fwrite(funnel, file.path(OUT, "funnel.csv"))
fwrite(final,  file.path(OUT, "stations_confirmed.csv"))
cat("\n最终确认因果耦合的站数:", n2, "/", n0, "\n")
cat("-> ", file.path(OUT, "funnel.csv"), "\n")
cat("-> ", file.path(OUT, "stations_confirmed.csv"), "\n")
