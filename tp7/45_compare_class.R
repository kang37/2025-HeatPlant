#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 45_compare_class.R — tp<=8 与 tp<=7 两种截断下的站点分类对照
#
# 只做分类, 不进回归。目的是先看清"截掉末步"改变了哪些站的归属:
#   末步(tp=8)才发生的变号, 在 tp7 口径下看不见 -> 该站变成"始终抑制/始终促进";
#   仅因末步多出一次变号而被判为多次变号的站, 在 tp7 口径下重新变成干净的单次翻转。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)

classify <- function(co) {
  s <- sign(co)
  if (anyNA(s) || any(s == 0)) return(list(ty = "undefined", ev = NA_integer_))
  k <- which(diff(s) != 0)
  if (!length(k)) return(list(ty = if (s[1] < 0) "1 始终抑制" else "3 始终促进",
                              ev = NA_integer_))
  if (length(k) > 1) return(list(ty = "X 多次变号", ev = NA_integer_))
  list(ty = if (s[1] < 0) "2 先抑制后促进" else "4 先促进后抑制", ev = as.integer(k))
}

grab <- function(M) {
  x <- cm[tp <= M]
  p <- x[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)]
  p <- p[sapply(co, length) == M + 1L]
  cl <- lapply(p$co, classify)
  p[, `:=`(ty = sapply(cl, `[[`, "ty"), ev = as.integer(sapply(cl, `[[`, "ev")))][, co := NULL]
  setnames(p, c("ty", "ev"), paste0(c("ty", "ev"), M))
  p[]
}
a8 <- grab(8L); a7 <- grab(7L)

cat("\n===== 站点分类构成 =====\n")
tb <- merge(a8[, .N, by = .(ty = ty8)], a7[, .N, by = .(ty = ty7)],
            by = "ty", all = TRUE, suffixes = c("_tp8", "_tp7"))
setnafill(tb, fill = 0L, cols = c("N_tp8", "N_tp7"))
tb[, 变化 := N_tp7 - N_tp8]
print(tb[order(ty)])
cat("\n可用于生存/回归的四类站合计: tp8 =", a8[ty8 != "X 多次变号" & ty8 != "undefined", .N],
    " tp7 =", a7[ty7 != "X 多次变号" & ty7 != "undefined", .N], "\n")

cat("\n===== 逐站归属的迁移矩阵 (行 = tp8 口径, 列 = tp7 口径) =====\n")
m <- merge(a8, a7, by = "stat_id")
print(dcast(m[, .N, by = .(ty8, ty7)], ty8 ~ ty7, value.var = "N", fill = 0L))
cat("\n归属改变的站:", m[ty8 != ty7, .N], "/", nrow(m), "\n")

cat("\n===== 翻转时间分布 (仅单次翻转站) =====\n")
print(dcast(rbind(a8[!is.na(ev8), .(口径 = "tp<=8", 起始 = substr(ty8, 1, 1), tp = ev8)],
                  a7[!is.na(ev7), .(口径 = "tp<=7", 起始 = substr(ty7, 1, 1), tp = ev7)])[
        , .N, by = .(口径, 起始, tp)], 口径 + 起始 ~ tp, value.var = "N", fill = 0L))

fwrite(m, file.path(OUT, "compare_class_tp8_vs_tp7.csv"))
cat("\n已写出 compare_class_tp8_vs_tp7.csv\n")
