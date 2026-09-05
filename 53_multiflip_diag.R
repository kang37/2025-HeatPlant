#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 53_multiflip_diag.R — 诊断: 为什么一半站被判成 multi_flip?
#
# 严格判据要求 9 个滞后的符号序列至多变号一次。纯噪声下 8 个相邻对各自
# 独立变号的概率 0.5, 期望变号 4 次, "至多一次"的概率只有 18/512 = 3.5%。
# 所以 51% 的站能通过, 说明信号确实存在; 但问题是: 剩下 49% 到底是
# "真的来回摆动"还是"系数太小被噪声推着变号"。
#
# 三个诊断:
#   1) 变号次数分布 vs 纯噪声基准
#   2) 变号次数与 |mean_coef| 量级的关系
#   3) multi_flip 站的变号发生在什么量级的系数上
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
pat <- cm[, .(co = list(mean_coef), rho = list(rho)), by = .(stat_id = meteo_stat)]
pat <- pat[sapply(co, length) == 9L & !sapply(co, anyNA)]
pat[, `:=`(nchg   = sapply(co, function(x) sum(diff(sign(x)) != 0)),
           amed   = sapply(co, function(x) median(abs(x))),
           amin   = sapply(co, function(x) min(abs(x))),
           a0     = sapply(co, function(x) abs(x[1])),
           rhomed = sapply(rho, median))]
log_msg("完整站 ", nrow(pat), "\n")

log_msg("-- 变号次数分布 vs 纯噪声基准 --")
nb <- sapply(0:8, function(k) choose(8, k) * 2 / 2^9)      # 纯噪声下的概率
d <- pat[, .(实测站数 = .N), by = .(变号次数 = nchg)][order(变号次数)]
d <- merge(data.table(变号次数 = 0:8, 噪声期望站数 = round(nb * nrow(pat), 1)),
           d, by = "变号次数", all.x = TRUE)
d[is.na(实测站数), 实测站数 := 0L]
d[, 实测占比 := sprintf("%.1f%%", 100*实测站数/nrow(pat))]
print(d)
log_msg("\n严格可分类(变号<=1): ", pat[nchg <= 1, .N], " 站 (",
        round(100*mean(pat$nchg <= 1), 1), "%) | 纯噪声下应为 3.5%")

log_msg("\n-- 变号次数 x |mean_coef| 中位数 --")
pat[, abin := cut(amed, c(0, .01, .02, .05, .1, 1),
                  labels = c("<0.01","0.01-0.02","0.02-0.05","0.05-0.10",">0.10"))]
print(dcast(pat[, .N, by = .(abin, multi = nchg > 1)], abin ~ multi, value.var = "N")[
  , .(`|coef|中位数` = abin, 严格可分类 = `FALSE`, 多次变号 = `TRUE`,
      多次变号占比 = sprintf("%.1f%%", 100*`TRUE`/(`TRUE` + `FALSE`)))])

log_msg("\n-- 多次变号站: 变号发生在什么量级的系数上 --")
# 对每个 multi_flip 站, 取所有变号处两侧系数绝对值的较小者
mf <- pat[nchg > 1]
sw <- unlist(lapply(mf$co, function(x) {
  k <- which(diff(sign(x)) != 0); pmin(abs(x[k]), abs(x[k+1])) }))
allc <- unlist(pat$co)
log_msg("变号处较小侧 |coef| 的中位数 = ", round(median(sw), 4),
        " | 全部系数 |coef| 的中位数 = ", round(median(abs(allc)), 4))
log_msg("变号处较小侧 |coef| < 0.01 的比例 = ", round(100*mean(sw < .01), 1), "%",
        " | 全部系数中该比例 = ", round(100*mean(abs(allc) < .01), 1), "%")
log_msg("变号处较小侧 |coef| < 0.02 的比例 = ", round(100*mean(sw < .02), 1), "%")

log_msg("\n-- 多次变号站与严格站的对比 --")
print(pat[, .(站数 = .N, `|coef|中位` = round(median(amed), 4),
              `rho 中位` = round(median(rhomed), 3),
              `|coef(tp0)| 中位` = round(median(a0), 4)),
          by = .(类型 = fifelse(nchg > 1, "多次变号", "严格可分类"))])
fwrite(pat[, .(stat_id, nchg, amed, amin, a0, rhomed)],
       file.path(OUT, "multiflip_diag.csv"))
