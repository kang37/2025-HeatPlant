#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 55_reclass_eval.R — 新判据(S2 加权变点)到底有没有用?
#
# S2 已通过两道内部检验: 在 461 个严格站上 100% 复现原判据, 且把 898 站全部定型。
# 但"能定型"不等于"定得对"。真正的检验是外部的:
#   Q1 起始方向层: 用 S2 的方向做因变量, 空间分块 AUC 会不会掉?
#      (样本从 385 扩到 762, 若 AUC 保持, 说明新增的站不是噪声)
#   Q2 翻转层: 样本翻倍之后, "会不会翻转"是不是终于可预测了?
#      (这是四分类相对二分类唯一可能的增量)
#   Q3 类型层: 四分类的折外表现是否改善, 最小类是否还小到不可用
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(xgboost) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
PRED <- obj$PRED; FEAT <- c(PRED, "crop")
rc <- fread(file.path(OUT, "reclassify.csv"))
D <- merge(obj$D, rc[, .(stat_id, S0, S2, A, s2_start, s2_k)], by = "stat_id")
D[, `:=`(dir_new  = as.integer(s2_start == "inhibit"),
         flip_new = as.integer(!is.na(s2_k)),
         flip_old = as.integer(stype %in% c("inhibit_promote","promote_inhibit")))]
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }
auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0); if (!length(p) || !length(n)) return(NA_real_)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(w[i]*w[n]*(d>0)) +
                   .5*sum(w[i]*w[n]*(d==0)); den <- den + sum(w[i]*w[n]) }
  num/den }
PAR <- list(objective = "binary:logistic", eval_metric = "logloss", max_depth = 3,
            eta = .03, min_child_weight = 8, subsample = .8, colsample_bytree = .8,
            lambda = 3, nthread = 4)
z <- copy(D); for (v in FEAT) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
z[, fold := { set.seed(42); kmeans(scale(cbind(longitude, latitude)), 10,
                                   nstart = 25, iter.max = 100)$cluster }]
z[, wt := winz(w, .01)][, wt := wt/mean(wt)]

cvrun <- function(d, yv, tag, useW = FALSE) {
  d <- copy(d); d[, yy := d[[yv]]]; d[, ww := if (useW) wt else 1]
  pl <- px <- rep(NA_real_, nrow(d))
  f <- as.formula(paste("yy ~", paste(PRED, collapse = "+")))
  X <- as.matrix(d[, FEAT, with = FALSE])
  for (k in sort(unique(d$fold))) {
    ti <- which(d$fold != k); vi <- which(d$fold == k)
    if (!length(vi) || length(unique(d$yy[ti])) < 2) next
    pl[vi] <- suppressWarnings(predict(glm(f, binomial(), d[ti], weights = ww),
                                       d[vi], type = "response"))
    dtr <- xgb.DMatrix(X[ti, ], label = d$yy[ti], weight = d$ww[ti])
    inner <- lapply(sort(unique(d$fold[ti])), function(j) which(d$fold[ti] == j))
    nb <- max(10L, xgb.cv(PAR, dtr, nrounds = 600, folds = inner,
                          early_stopping_rounds = 40, verbose = 0)$early_stop$best_iteration)
    px[vi] <- predict(xgb.train(PAR, dtr, nrounds = nb, verbose = 0), X[vi, , drop = FALSE])
  }
  ok <- !is.na(pl)
  data.table(设定 = tag, n = sum(ok), 阳性率 = round(mean(d$yy[ok]), 3),
             AUC_logit = round(auc_w(d$yy[ok], pl[ok]), 3),
             AUC_xgb   = round(auc_w(d$yy[ok], px[ok]), 3))
}

log_msg("=== Q1 起始方向层: 旧口径 vs 新口径 ===\n")
q1 <- rbindlist(list(
  cvrun(z[S0 != "multi_flip"], "grp2",    "旧: 严格站 x tp=0 符号(385)", TRUE),
  cvrun(z,                      "grp2",    "旧: 全部站 x tp=0 符号(762)", TRUE),
  cvrun(z[S0 != "multi_flip"], "dir_new", "新: 严格站 x S2 方向(385)",  TRUE),
  cvrun(z,                      "dir_new", "新: 全部站 x S2 方向(762)",  TRUE)))
print(q1)

log_msg("\n=== Q2 翻转层: 样本翻倍后是否变得可预测 ===\n")
q2 <- rbindlist(list(
  cvrun(z[S0 != "multi_flip" & grp2 == 1], "flip_old", "旧: 抑制组内翻转(214)"),
  cvrun(z[S0 != "multi_flip" & grp2 == 0], "flip_old", "旧: 促进组内翻转(171)"),
  cvrun(z[dir_new == 1],                   "flip_new", "新: 抑制组内翻转"),
  cvrun(z[dir_new == 0],                   "flip_new", "新: 促进组内翻转"),
  cvrun(z,                                  "flip_new", "新: 不分组直接预测翻转")))
print(q2)

log_msg("\n=== Q3 四分类的类型分布 ===\n")
LB <- c(always_inhibit="始终抑制", inhibit_promote="抑制转促进",
        always_promote="始终促进", promote_inhibit="促进转抑制")
cmp <- merge(z[S0 != "multi_flip", .(旧 = .N), by = .(类型 = LB[S0])],
             z[, .(新 = .N), by = .(类型 = LB[S2])], by = "类型", all = TRUE)
cmp[is.na(旧), 旧 := 0L]
print(cmp[order(-新)])
log_msg("最小类: 旧 ", min(cmp$旧[cmp$旧 > 0]), " 站 -> 新 ", min(cmp$新), " 站")
fwrite(rbind(q1, q2, fill = TRUE), file.path(OUT, "reclass_eval.csv"))
