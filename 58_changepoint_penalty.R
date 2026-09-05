#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 58_changepoint_penalty.R — S2 的一个缺陷与修法
#
# 57 号发现: S2 与"简单多数"的 364 处分歧 100% 落在翻转型站上。查具体案例:
#   站 50468  mean_coef = -0.005 +0.062 +0.110 ... +0.185 +0.126
# 简单多数判"促进"(8/9 为正), S2 判"起始抑制、tp=1 变点", 一致度 1.00。
# 但那个"抑制"证据只有一个 -0.005 —— 是噪声。S2 之所以照单全收, 是因为
# 它在 18 个候选里取 argmax(加权一致度) 而不给变点任何惩罚:
# 只要挪一个变点就能把这个噪声滞后"解释掉", A 就升到 1.00。
#
# 修法 S2b: 变点模型必须比最好的常数模型多解释 delta 的权重才被接受。
# delta = 0 即退回 S2; delta 越大越倾向判为"始终抑制/始终促进"。
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
rc <- fread(file.path(OUT, "reclassify.csv"))
P <- merge(P, rc[, .(stat_id, S0)], by = "stat_id")

fitboth <- function(x) {
  w <- abs(x); s <- sign(x); n <- length(x); tot <- sum(w)
  bc <- list(A = -1); bk <- list(A = -1)
  for (s1 in c(-1, 1)) {
    A <- sum(w * (s == rep(s1, n))) / tot
    if (A > bc$A) bc <- list(A = A, s1 = s1)
    for (k in 1:(n-1)) {
      A <- sum(w * (s == c(rep(s1, k), rep(-s1, n-k)))) / tot
      if (A > bk$A) bk <- list(A = A, s1 = s1, k = k)
    }
  }
  c(bc$A, bc$s1, bk$A, bk$s1, bk$k)
}
r <- t(sapply(P$co, fitboth))
P[, `:=`(A_const = r[,1], s_const = r[,2], A_flip = r[,3], s_flip = r[,4], k_flip = r[,5])]
P[, gain := A_flip - A_const]

log_msg("################ 变点带来的额外解释力有多大 ################\n")
log_msg("gain = 最好的变点模型 - 最好的常数模型 (加权一致度之差)")
log_msg("分位: ", paste(sprintf("%s=%.3f", c("P10","P25","中位","P75","P90"),
        quantile(P$gain, c(.1,.25,.5,.75,.9))), collapse = "  "))
log_msg("gain < 0.05 的站: ", P[gain < .05, .N], " (",
        sprintf("%.1f%%", 100*mean(P$gain < .05)), ") —— 变点几乎没多解释什么")
log_msg("gain < 0.02 的站: ", P[gain < .02, .N], " (",
        sprintf("%.1f%%", 100*mean(P$gain < .02)), ")")

TY <- function(s1, isflip) fifelse(!isflip, fifelse(s1 < 0, "always_inhibit", "always_promote"),
                                   fifelse(s1 < 0, "inhibit_promote", "promote_inhibit"))
log_msg("\n################ delta 扫描 ################\n")
DEL <- c(0, .02, .05, .10, .15, .20)
res <- rbindlist(lapply(DEL, function(dl) {
  isf <- P$gain > dl
  ty  <- TY(fifelse(isf, P$s_flip, P$s_const), isf)
  st  <- fifelse(isf, P$s_flip, P$s_const)
  chk <- P$S0 != "multi_flip"
  data.table(delta = dl, 翻转型站数 = sum(isf),
             `翻转型占比` = sprintf("%.1f%%", 100*mean(isf)),
             `复现 S0` = sprintf("%.1f%%", 100*mean(ty[chk] == P$S0[chk])),
             `最小类 n` = min(table(ty)),
             `始终抑制` = sum(ty == "always_inhibit"),
             `抑制转促进` = sum(ty == "inhibit_promote"),
             `始终促进` = sum(ty == "always_promote"),
             `促进转抑制` = sum(ty == "promote_inhibit"))
}))
print(res)

log_msg("\n################ 站 50468 在各 delta 下被判成什么 ################\n")
i <- which(P$stat_id == 50468)
if (length(i)) {
  log_msg("mean_coef = ", paste(sprintf("%+.3f", P$co[[i]]), collapse = " "))
  log_msg("A_const = ", round(P$A_const[i], 3), " (", ifelse(P$s_const[i] < 0, "抑制", "促进"), ")",
          " | A_flip = ", round(P$A_flip[i], 3), " | gain = ", round(P$gain[i], 3))
  for (dl in DEL) log_msg("  delta=", sprintf("%.2f", dl), " -> ",
    TY(fifelse(P$gain[i] > dl, P$s_flip[i], P$s_const[i]), P$gain[i] > dl))
}

# 推荐 delta: 取能保住 S0 复现率 >= 95% 的最大值
ok <- res[as.numeric(sub("%","",`复现 S0`)) >= 95]
log_msg("\n复现 S0 >= 95% 的最大 delta = ", max(ok$delta))
DL <- 0.05
isf <- P$gain > DL
P[, `:=`(S2b = TY(fifelse(isf, s_flip, s_const), isf),
         S2b_start = fifelse(fifelse(isf, s_flip, s_const) < 0, "inhibit", "promote"),
         S2b_k = fifelse(isf, k_flip, NA_real_),
         A_used = fifelse(isf, A_flip, A_const))]
fwrite(P[, .(stat_id, S0, S2b, S2b_start, S2b_k, A_const, A_flip, gain, A_used)],
       file.path(OUT, "reclassify_s2b.csv"))
log_msg("已按 delta = ", DL, " 写出 reclassify_s2b.csv")
