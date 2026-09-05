#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 54_reclassify.R — 三套改进的分类判据, 并排评估
#
# 诊断(53 号)已经定性: multi_flip 主要不是"真的来回摆动", 而是"系数太小,
# 符号被噪声推着变"。变号处较小侧 |coef| 中位数 0.0056, 只有全部系数中位数
# 的四分之一; 71% 的变号发生在 |coef| < 0.01 处。
# 所以改进方向不是放宽"允许几次变号", 而是不要在近零系数上强行读符号。
#
#   S0 严格(现状)   9 个符号至多变号一次, 否则 multi_flip
#   S1 零带         |coef| < tau 视为"未定", 只看非零部分的符号序列
#   S2 加权变点     在 18 个候选模型(2 个常数 + 2x8 个单变点)里选加权一致度
#                   最高者。权重 = |coef|, 于是近零滞后自动没有发言权。
#                   每个站都得到类型, 外加一个一致度 A 作可信度权重。
#
# S2 是主推方案: 不丢站、不设阈值、天然给出可信度, 且在严格站上应当能复现 S0。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
P <- cm[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)]
P <- P[sapply(co, length) == 9L & !sapply(co, anyNA)]
TYPE <- function(s1, flip) fifelse(!flip, fifelse(s1 < 0, "always_inhibit", "always_promote"),
                                   fifelse(s1 < 0, "inhibit_promote", "promote_inhibit"))

# ---- S0 严格 --------------------------------------------------------------
s0 <- function(x) { s <- sign(x); k <- which(diff(s) != 0)
  if (length(k) > 1) return("multi_flip")
  TYPE(s[1], length(k) == 1) }
P[, S0 := sapply(co, s0)]

# ---- S1 零带 --------------------------------------------------------------
# |coef| < tau 视为未定; 剩下的非零符号序列至多变号一次即可分类,
# 且要求非零滞后 >= 4, 翻转型两侧各 >= 2, 否则仍判为未定。
s1 <- function(x, tau) {
  keep <- abs(x) >= tau
  if (sum(keep) < 4) return("undetermined")
  s <- sign(x[keep]); k <- which(diff(s) != 0)
  if (length(k) > 1) return("multi_flip")
  if (length(k) == 1 && min(k, length(s) - k) < 2) return("undetermined")
  TYPE(s[1], length(k) == 1)
}
TAUS <- c(0.005, 0.01, 0.015, 0.02, 0.03)
for (tau in TAUS) P[, paste0("S1_", tau) := sapply(co, s1, tau = tau)]

# ---- S2 加权变点 ----------------------------------------------------------
# 候选: 常数(+/-) 与 单变点(变点位置 1..8, 起始符号 +/-)
# 目标: 最大化 A = sum_t |coef_t| * 1{sign_t == fit_t} / sum_t |coef_t|
s2 <- function(x) {
  w <- abs(x); s <- sign(x); n <- length(x); tot <- sum(w)
  best <- list(A = -1)
  for (s1v in c(-1, 1)) {
    fit <- rep(s1v, n)                                   # 常数
    A <- sum(w * (s == fit)) / tot
    if (A > best$A) best <- list(A = A, s1 = s1v, k = NA_integer_)
    for (k in 1:(n-1)) {                                 # 在第 k 步之后变号
      fit <- c(rep(s1v, k), rep(-s1v, n - k))
      A <- sum(w * (s == fit)) / tot
      if (A > best$A) best <- list(A = A, s1 = s1v, k = k)
    }
  }
  c(best$A, best$s1, if (is.na(best$k)) NA_real_ else best$k)
}
r2 <- t(sapply(P$co, s2))
P[, `:=`(A = r2[,1], s2_s1 = r2[,2], s2_k = r2[,3])]
P[, S2 := TYPE(s2_s1, !is.na(s2_k))]
P[, tp0_dir := sapply(co, function(x) fifelse(x[1] < 0, "inhibit", "promote"))]

# ===== 评估 ================================================================
log_msg("################ 1. 可分类站数 ################\n")
tab <- rbindlist(c(
  list(data.table(判据 = "S0 严格(现状)", 可分类 = P[S0 != "multi_flip", .N])),
  lapply(TAUS, function(tau) { v <- P[[paste0("S1_", tau)]]
    data.table(判据 = sprintf("S1 零带 tau=%.3f", tau),
               可分类 = sum(!v %in% c("multi_flip", "undetermined"))) }),
  list(data.table(判据 = "S2 加权变点", 可分类 = nrow(P)))))
tab[, `占 898 站` := sprintf("%.1f%%", 100*可分类/nrow(P))]
print(tab)

log_msg("\n################ 2. 与严格判据的一致性(仅看 S0 已定型的 461 站) ################\n")
chk <- P[S0 != "multi_flip"]
log_msg("S2 复现 S0 的比例 = ", sprintf("%.1f%%", 100*mean(chk$S2 == chk$S0)),
        "  (", sum(chk$S2 == chk$S0), "/", nrow(chk), ")")
bad <- chk[S2 != S0]
if (nrow(bad)) {
  log_msg("不一致 ", nrow(bad), " 站, 其一致度 A 中位 = ", round(median(bad$A), 3),
          " vs 一致站 A 中位 = ", round(median(chk[S2 == S0]$A), 3))
  print(bad[, .N, by = .(S0, S2)][order(-N)][1:min(5,.N)])
}
for (tau in TAUS) { v <- P[[paste0("S1_", tau)]]
  ok <- chk[[paste0("S1_", tau)]]
  m <- !ok %in% c("multi_flip", "undetermined")
  log_msg(sprintf("S1 tau=%.3f 在 461 站上: 仍可定型 %d 站, 其中与 S0 一致 %.1f%%",
                  tau, sum(m), 100*mean(chk$S0[m] == ok[m])))
}

log_msg("\n################ 3. 起始方向是否被改动 ################\n")
log_msg("S2 的起始方向与 tp=0 符号一致的站: ",
        sum((P$s2_s1 < 0) == (P$tp0_dir == "inhibit")), " / ", nrow(P),
        "  (", sprintf("%.1f%%", 100*mean((P$s2_s1 < 0) == (P$tp0_dir == "inhibit"))), ")")
dif <- P[(s2_s1 < 0) != (tp0_dir == "inhibit")]
log_msg("被改动的 ", nrow(dif), " 站: |coef(tp0)| 中位 = ",
        round(median(sapply(dif$co, function(x) abs(x[1]))), 4),
        " vs 全样本 ", round(median(sapply(P$co, function(x) abs(x[1]))), 4))

log_msg("\n################ 4. S2 的类型分布与可信度 ################\n")
print(P[, .(站数 = .N, 占比 = sprintf("%.1f%%", 100*.N/nrow(P)),
            `一致度A 中位` = round(median(A), 3)), by = .(S2类型 = S2)][order(-站数)])
log_msg("\n一致度 A 的分布: ", paste(sprintf("%s=%.2f", c("P10","P25","中位","P75","P90"),
        quantile(P$A, c(.1,.25,.5,.75,.9))), collapse = "  "))
log_msg("A >= 0.8 的站: ", P[A >= .8, .N], " (", sprintf("%.1f%%", 100*mean(P$A >= .8)), ")")
log_msg("原 multi_flip 的 437 站里 A >= 0.8 的: ", P[S0 == "multi_flip" & A >= .8, .N])

fwrite(P[, .(stat_id, S0, S2, A, s2_start = fifelse(s2_s1 < 0, "inhibit", "promote"),
             s2_k, tp0_dir,
             S1_0.01 = get("S1_0.01"), S1_0.02 = get("S1_0.02"))],
       file.path(OUT, "reclassify.csv"))
log_msg("\n已写出 reclassify.csv")
