#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 59_changepoint_test.R — 给每个站一个"这次变号是不是真的"的置换检验
#
# 58 号暴露的问题: S2 取 argmax(加权一致度), 对变点不设门槛, 于是一个
# -0.005 的噪声滞后也能换来一个"变点"。而固定阈值 delta 是任意的。
# 更好的办法是让数据自己定门槛: 对每个站做置换检验。
#
#   统计量 gain = 最好的单变点模型 - 最好的常数模型 (加权一致度之差)
#   零假设 系数的时间顺序无关紧要(即不存在真实转变)
#   置换   打乱该站 9 个 mean_coef 的 tp 顺序, 幅度与符号的多重集不变,
#          只破坏时间结构; 重算 gain, 999 次
#   判定   p = (1 + #{gain_perm >= gain_obs}) / 1000; p < 0.05 才认变点
#
# 这样"要不要判成翻转型"由该站自身的信号强度决定, 强信号站容易过关,
# 弱信号站的变号自动被否掉, 不需要任何人为阈值。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
P <- cm[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)]
P <- P[sapply(co, length) == 9L & !sapply(co, anyNA)]
rc <- fread(file.path(OUT, "reclassify.csv"))
P <- merge(P, rc[, .(stat_id, S0)], by = "stat_id")

# 用累积和把 18 个候选一次算完: O(n) 而不是 O(n^2)
fit <- function(x) {
  w <- abs(x); pos <- x > 0; tot <- sum(w)
  Wp <- sum(w[pos])                       # 判"正"时的一致权重
  cp <- c(0, cumsum(w * pos))             # 前 k 步里"正"的权重
  cw <- c(0, cumsum(w))
  n <- length(x)
  # 常数
  Ac <- max(Wp, tot - Wp) / tot
  sc <- if (Wp >= tot - Wp) 1 else -1
  # 单变点: 前 k 步判 s1, 后 n-k 步判 -s1
  k <- 1:(n-1)
  agr_p <- cp[k+1] + ((tot - Wp) - (cw[k+1] - cp[k+1]))   # s1=+1
  agr_n <- (cw[k+1] - cp[k+1]) + (Wp - cp[k+1])           # s1=-1
  Af <- max(c(agr_p, agr_n)) / tot
  j <- which.max(c(agr_p, agr_n))
  s1 <- if (j <= length(k)) 1 else -1
  kk <- k[if (j <= length(k)) j else j - length(k)]
  c(Ac, sc, Af, s1, kk)
}
NPERM <- 999L
out <- t(sapply(P$co, function(x) {
  o <- fit(x); g <- o[3] - o[1]
  gp <- replicate(NPERM, { y <- x[sample.int(length(x))]; f <- fit(y); f[3] - f[1] })
  c(o, g, (1 + sum(gp >= g - 1e-12)) / (NPERM + 1))
}))
P[, `:=`(A_const = out[,1], s_const = out[,2], A_flip = out[,3],
         s_flip = out[,4], k_flip = out[,5], gain = out[,6], p_flip = out[,7])]

TY <- function(s1, isf) fifelse(!isf, fifelse(s1 < 0, "always_inhibit", "always_promote"),
                                fifelse(s1 < 0, "inhibit_promote", "promote_inhibit"))
P[, isf := p_flip < .05]
P[, `:=`(S3 = TY(fifelse(isf, s_flip, s_const), isf),
         S3_start = fifelse(fifelse(isf, s_flip, s_const) < 0, "inhibit", "promote"),
         S3_k = fifelse(isf, k_flip, NA_real_))]

log_msg("################ 置换检验结果(", NPERM, " 次/站) ################\n")
log_msg("判为存在真实变点(p<0.05)的站: ", P[isf == TRUE, .N], " / ", nrow(P),
        " (", sprintf("%.1f%%", 100*mean(P$isf)), ")")
log_msg("S2 无惩罚时判为翻转型的站: 660 (73.5%)  —— 置换检验否掉了 ",
        660 - P[isf == TRUE, .N], " 站的变点")
log_msg("\np 值分布: ", paste(sprintf("%s=%.3f", c("P10","P25","中位","P75"),
        quantile(P$p_flip, c(.1,.25,.5,.75))), collapse = "  "))

log_msg("\n-- 四类分布对比 --")
LB <- c(always_inhibit="始终抑制", inhibit_promote="抑制转促进",
        always_promote="始终促进", promote_inhibit="促进转抑制")
s2 <- fread(file.path(OUT, "reclassify.csv"))
cmpt <- merge(
  merge(P[, .(S3 = .N), by = .(类型 = LB[S3])],
        s2[, .(S2 = .N), by = .(类型 = LB[S2])], by = "类型", all = TRUE),
  s2[S0 != "multi_flip", .(S0严格 = .N), by = .(类型 = LB[S0])], by = "类型", all = TRUE)
for (j in c("S0严格","S2","S3")) set(cmpt, which(is.na(cmpt[[j]])), j, 0L)
print(cmpt[order(-S3), .(类型, S0严格, `S2无惩罚` = S2, `S3置换检验` = S3)])

log_msg("\n-- 站 50468 复核 --")
i <- which(P$stat_id == 50468)
log_msg("mean_coef = ", paste(sprintf("%+.3f", P$co[[i]]), collapse = " "))
log_msg("gain = ", round(P$gain[i], 4), " | 置换 p = ", round(P$p_flip[i], 3),
        " -> 判为 ", LB[P$S3[i]], "  (S2 无惩罚时判 抑制转促进)")

log_msg("\n-- 与 S0 的关系(仅 461 个严格站) --")
chk <- P[S0 != "multi_flip"]
log_msg("S3 与 S0 一致: ", sprintf("%.1f%%", 100*mean(chk$S3 == chk$S0)))
FL <- c("inhibit_promote","promote_inhibit")
log_msg("不一致的 ", sum(chk$S3 != chk$S0), " 站里:")
log_msg("  S0 判翻转型 而 S3 判恒定型: ", chk[S3 != S0 & S0 %in% FL & !S3 %in% FL, .N])
log_msg("  S0 判恒定型 而 S3 判翻转型: ", chk[S3 != S0 & !S0 %in% FL & S3 %in% FL, .N])
log_msg("  两者都是翻转型但方向不同  : ", chk[S3 != S0 & S0 %in% FL & S3 %in% FL, .N])
bad <- chk[S3 != S0]
log_msg("这些站的 gain 中位 = ", round(median(bad$gain), 4),
        " vs 一致站 ", round(median(chk[S3 == S0]$gain), 4))
log_msg("  -> S0 把它们判成翻转型, 是因为严格判据同样不给变点设门槛:")
log_msg("     只要恰好只变一次号, 哪怕变在 |coef|=0.005 处也算数。")

fwrite(P[, .(stat_id, S0, S3, S3_start, S3_k, A_const, A_flip, gain, p_flip)],
       file.path(OUT, "reclassify_s3.csv"))
log_msg("\n已写出 reclassify_s3.csv")
