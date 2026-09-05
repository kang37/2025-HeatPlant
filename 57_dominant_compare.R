#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 57_dominant_compare.R — "主导方向"到底怎么定? 三种候选口径的对照
#
#   M1 简单多数   9 个符号里哪边多算哪边(用户猜的口径)
#   M2 加权多数   sign(sum(mean_coef)), 大系数说话更响, 但仍是全局一票制
#   S2 加权变点   先拟合"常数 或 单变点"的阶梯, 起始方向 = 第一段的符号
#
# 关键差别有两处:
#   (1) M1 数个数, S2/M2 按 |coef| 加权 —— |coef|=0.10 的滞后应当比 0.005 的重 20 倍
#   (2) M1/M2 是全局汇总, S2 的"起始方向"只由变点之前那一段决定 ——
#       对翻转型的站, 后半段再长也不该改写"起始方向"
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
P <- merge(P, rc[, .(stat_id, S0, S2, A, s2_start, s2_k)], by = "stat_id")

P[, `:=`(
  n_neg = sapply(co, function(x) sum(x < 0)),
  M1    = sapply(co, function(x) fifelse(sum(x < 0) > 4.5, "inhibit", "promote")),
  M2    = sapply(co, function(x) fifelse(sum(x)   < 0,     "inhibit", "promote")),
  tp0   = sapply(co, function(x) fifelse(x[1]     < 0,     "inhibit", "promote")))]

log_msg("################ 四种口径两两一致率(898 站) ################\n")
K <- c("s2_start","M1","M2","tp0")
NM <- c(s2_start = "S2 加权变点", M1 = "M1 简单多数", M2 = "M2 加权多数", tp0 = "tp=0 符号")
m <- outer(K, K, Vectorize(function(a, b) mean(P[[a]] == P[[b]])))
dimnames(m) <- list(NM[K], NM[K])
print(round(100*m, 1))

log_msg("\n################ S2 与简单多数在哪儿分道扬镳 ################\n")
d <- P[s2_start != M1]
log_msg("不一致 ", nrow(d), " 站 (", sprintf("%.1f%%", 100*nrow(d)/nrow(P)), ")")
log_msg("其中属于翻转型(S2 判为有变点)的: ", d[!is.na(s2_k), .N],
        " 站 (", sprintf("%.0f%%", 100*mean(!is.na(d$s2_k))), ")")
log_msg("  -> 这些站上 M1 把「后半段更长」误读成了「主导方向」")
log_msg("其中属于常数型的: ", d[is.na(s2_k), .N], " 站")
log_msg("  -> 这些站上少数几个滞后的 |coef| 压倒了多数个小系数")

log_msg("\n-- 两类分歧各举 3 个例子 --")
show <- function(dd, tag) {
  if (!nrow(dd)) return(invisible())
  log_msg("\n【", tag, "】")
  for (i in seq_len(min(3, nrow(dd)))) {
    x <- dd$co[[i]]
    log_msg("站 ", dd$stat_id[i], "  mean_coef = ",
            paste(sprintf("%+.3f", x), collapse = " "))
    log_msg("        负号 ", sum(x < 0), "/9 -> 简单多数判「",
            ifelse(dd$M1[i] == "inhibit", "抑制", "促进"), "」",
            " | S2 判「", ifelse(dd$s2_start[i] == "inhibit", "抑制", "促进"), "」起始",
            ifelse(is.na(dd$s2_k[i]), "(无变点)", sprintf("(变点在 tp=%d)", dd$s2_k[i])),
            " | 一致度 ", sprintf("%.2f", dd$A[i]))
  }
}
show(d[!is.na(s2_k)][order(-A)], "翻转型：后半段更长，但起始方向不该被改写")
show(d[is.na(s2_k)][order(-A)],  "常数型：少数大系数压倒多数小系数")

log_msg("\n################ 哪个口径更可预测(用已算好的折外结果做参照) ################\n")
log_msg("说明: 起始方向层的空间分块 AUC —— tp=0 符号 0.789/0.827, S2 方向 0.775/0.837")
log_msg("      (logit/XGBoost, n=762)。三个口径差别在 0.01 量级, 预测力上分不出高下,")
log_msg("      选哪个取决于定义本身是否贴合研究问题, 不取决于 AUC。")
fwrite(P[, .(stat_id, S0, S2, A, s2_start, s2_k, M1, M2, tp0, n_neg)],
       file.path(OUT, "dominant_compare.csv"))
