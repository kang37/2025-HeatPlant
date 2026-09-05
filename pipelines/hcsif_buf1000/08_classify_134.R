#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 08_classify_134.R — 对134站(确认VPD->SIF方向)做"方向 x tp响应速度"分类
#
# 用户决定(2026-09-05): 方向就用 06_smap_bivar_134.R 的二变量S-map主线结果
# (Promote/Inhibit/Ambiguous)，不再因为07脚本发现的稳健性问题深挖；也不做
# k-means，只走方向x tp分箱的先验方案。
#
# 三个维度的处理方式:
#   方向(direction_2v)  -> 主轴，3档(Promote/Inhibit/Ambiguous)
#   optimal tp          -> 副轴，分3箱: Fast(0-2,8-16天) / Mid(3-5,24-40天) /
#                           Slow(6-8,48-64天)，对应短滞后气孔响应到长滞后累积
#                           水分胁迫的生理时间尺度
#   nsig(显著tp数量/强度) -> 不再开分类轴(134站体量小，再分会有格子<5站)，作为
#                           连续叠加变量放进可视化(点大小/颜色深浅)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
OUT_DIR <- file.path(PROJ, "data_proc/smap_bivar_134")

st <- readRDS(file.path(OUT_DIR, "smap_bivar_134.rds"))
stopifnot(nrow(st) > 0)
cat("站数 n =", nrow(st), "(过滤地表覆盖断点站后，脚本/文件名沿用原 134 命名)\n")

st[, tp_bin := fifelse(optimal_tp <= 2, "Fast (tp0-2)",
                 fifelse(optimal_tp <= 5, "Mid (tp3-5)", "Slow (tp6-8)"))]
st[, tp_bin := factor(tp_bin, levels = c("Fast (tp0-2)", "Mid (tp3-5)", "Slow (tp6-8)"))]
st[, direction_2v := factor(direction_2v, levels = c("Promote", "Inhibit", "Ambiguous"))]
st[, class_label := paste0(direction_2v, " / ", tp_bin)]

cat("=== 分类矩阵: 方向 x tp速度分箱 (N) ===\n")
tab_n <- dcast(st, direction_2v ~ tp_bin, value.var = "meteo_stat", fun.aggregate = length)
print(tab_n)

cat("\n=== 每格 |median_coef_2v| 中位数(强度) ===\n")
tab_strength <- dcast(st, direction_2v ~ tp_bin,
                       value.var = "median_coef_2v",
                       fun.aggregate = function(x) round(median(abs(x), na.rm = TRUE), 4))
print(tab_strength)

cat("\n=== 每格 nsig_0to8 中位数(显著tp广度) ===\n")
tab_nsig <- dcast(st, direction_2v ~ tp_bin,
                   value.var = "nsig_0to8",
                   fun.aggregate = function(x) round(median(x, na.rm = TRUE), 1))
print(tab_nsig)

cat("\n=== class_label 频数(排序) ===\n")
print(st[, .N, by = class_label][order(-N)])

saveRDS(st, file.path(OUT_DIR, "classify_134.rds"))
fwrite(st, file.path(OUT_DIR, "classify_134.csv"))
cat("\n已写出 ", file.path(OUT_DIR, "classify_134.{rds,csv}"), "\n")
