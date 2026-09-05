#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 44_class_data.R — 装配"符号分类"分析用的建模表并缓存
#
# 组别定义(已核对): 抑制组 <=> mean_coef 在 tp=0 为负。四类型以起始方向命名,
# 所以"起始符号"与"两组归属"完全等价, 在严格分类的站上不一致数为 0。
#   grp2  二分类:  0 = 促进组(tp0 系数为正), 1 = 抑制组(tp0 系数为负)
#   stype 四分类:  always_inhibit / inhibit_promote / always_promote / promote_inhibit
#                  (多次变号的站 stype = multi_flip, 只在样本 B 中出现)
#
# 样本 A = 严格四分类站(至多变号一次); 样本 B = 全部有完整 9 个滞后的站。
# 权重 w = |mean_coef(tp=0)|, 因为约一半站的 |系数| < 0.03, 符号近乎抛硬币,
# 加权比剔除低信号站更好(剔除会按结果筛样本, 引入选择偏误)。
#
# 输出: data_proc/output_hcsif_buf1000/class_model_data.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")

# ---- 复用 22 号脚本的协变量装配(执行到 CORE 定义为止) --------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(
  eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)
log_msg("协变量 ", length(PRED), " 个: ", paste(PRED, collapse = ", "))

# ---- 分类 ----------------------------------------------------------------
cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
N_TP <- length(unique(cm$tp))

classify_strict <- function(co) {
  s <- sign(co)
  if (anyNA(s) || any(s == 0)) return("undefined")
  k <- which(diff(s) != 0)
  if (length(k) == 0) return(if (s[1] < 0) "always_inhibit" else "always_promote")
  if (length(k) > 1)  return("multi_flip")
  if (s[1] < 0) "inhibit_promote" else "promote_inhibit"
}

pat <- cm[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)]
pat <- pat[sapply(co, length) == N_TP]
pat[, stype := sapply(co, classify_strict)]
pat[, coef0 := sapply(co, `[`, 1)]
pat <- pat[stype != "undefined"]
pat[, co := NULL]
pat[, `:=`(grp2 = as.integer(coef0 < 0),          # 1 = 抑制组
           w    = abs(coef0),
           strict = stype != "multi_flip")]

log_msg("\n-- 全部完整站 ", nrow(pat), " --")
print(pat[, .N, by = .(stype)][order(-N)])
log_msg("样本 A(严格): ", pat[strict == TRUE, .N],
        " 站 | 抑制 ", pat[strict == TRUE & grp2 == 1, .N],
        " / 促进 ", pat[strict == TRUE & grp2 == 0, .N])
log_msg("样本 B(全部): ", nrow(pat),
        " 站 | 抑制 ", pat[grp2 == 1, .N], " / 促进 ", pat[grp2 == 0, .N])
# 核对: 严格站上"起始方向命名"与"tp0 符号"是否一致
chk <- pat[strict == TRUE][, sum((stype %in% c("always_inhibit","inhibit_promote")) != (grp2 == 1))]
log_msg("严格站上两种组别定义不一致的站数 = ", chk, "  (应为 0)")

st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
D <- merge(merge(pat, dt[, c("stat_id", PRED, "crop", "forest", "koppen_group"), with = FALSE],
                 by = "stat_id"),
           st[, .(stat_id, longitude, latitude)], by = "stat_id")
D[, cc := complete.cases(D[, PRED, with = FALSE])]
log_msg("\n协变量完整: 样本 B ", D[cc == TRUE, .N], " 站 | 样本 A ",
        D[cc == TRUE & strict == TRUE, .N], " 站")
print(D[cc == TRUE, .N, by = .(样本 = ifelse(strict, "A严格", "仅B"),
                               组 = ifelse(grp2 == 1, "抑制", "促进"))])

saveRDS(list(D = D[cc == TRUE], PRED = PRED), file.path(OUT, "class_model_data.rds"))
log_msg("\n已写出 class_model_data.rds")
