#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 51_flip_predict.R — 决定性检验: 给定起始方向, "会不会翻转"到底可不可预测?
#
# 四分类 = 起始方向 x 是否翻转。起始方向那一层已知可预测(空间分块 AUC 0.867)。
# 所以"四分类比二分类多出的信息"全部装在第二层里。这个脚本就量它:
# 对两个组分别做 flip 的空间分块 CV(logit 与 XGBoost 各一遍), 看 AUC 离 0.5 有多远。
# 若 AUC ~ 0.5, 则四分类多出的那一维不含可被协变量恢复的信息, 二分类严格更优。
#
# 另报一个信息分解: 用折外预测的对数似然, 把四分类的总信息拆成
#   方向层贡献 + 翻转层贡献, 直接看后者占多少。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(xgboost) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
PRED <- obj$PRED; FEAT <- c(PRED, "crop")
D <- obj$D[strict == TRUE]
D[, flip := as.integer(stype %in% c("inhibit_promote", "promote_inhibit"))]
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

# 全样本一套折, 子集沿用同一套(保证与二分类结果可比)
z <- copy(D)
for (v in FEAT) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
z[, fold := { set.seed(42); kmeans(scale(cbind(longitude, latitude)), 10,
                                   nstart = 25, iter.max = 100)$cluster }]

cvrun <- function(d, yv, tag) {
  d <- copy(d); d[, yy := d[[yv]]]
  pl <- px <- rep(NA_real_, nrow(d))
  f <- as.formula(paste("yy ~", paste(PRED, collapse = "+")))
  X <- as.matrix(d[, FEAT, with = FALSE])
  for (k in sort(unique(d$fold))) {
    ti <- which(d$fold != k); vi <- which(d$fold == k)
    if (!length(vi) || length(unique(d$yy[ti])) < 2) next
    pl[vi] <- suppressWarnings(predict(glm(f, binomial(), d[ti]), d[vi], type = "response"))
    dtr <- xgb.DMatrix(X[ti, ], label = d$yy[ti])
    inner <- lapply(sort(unique(d$fold[ti])), function(j) which(d$fold[ti] == j))
    nb <- max(10L, xgb.cv(PAR, dtr, nrounds = 600, folds = inner,
                          early_stopping_rounds = 40, verbose = 0)$early_stop$best_iteration)
    px[vi] <- predict(xgb.train(PAR, dtr, nrounds = nb, verbose = 0), X[vi, , drop = FALSE])
  }
  ok <- !is.na(pl)
  ll <- function(p) { p <- pmin(pmax(p, 1e-6), 1-1e-6)
    sum(d$yy[ok]*log(p[ok]) + (1-d$yy[ok])*log(1-p[ok])) }
  base <- mean(d$yy[ok])
  data.table(层 = tag, n = sum(ok), 阳性率 = round(base, 3),
             AUC_logit = round(auc_w(d$yy[ok], pl[ok]), 3),
             AUC_xgb   = round(auc_w(d$yy[ok], px[ok]), 3),
             `伪R2_logit` = round(1 - ll(pl)/ll(rep(base, nrow(d))), 3),
             `伪R2_xgb`   = round(1 - ll(px)/ll(rep(base, nrow(d))), 3))
}

log_msg("=== 折外预测: 每一层各自能不能预测(样本 A, 同一套空间折) ===\n")
r <- rbindlist(list(
  cvrun(z,             "grp2", "第1层 起始方向"),
  cvrun(z[grp2 == 1],  "flip", "第2层a 抑制组内翻转"),
  cvrun(z[grp2 == 0],  "flip", "第2层b 促进组内翻转"),
  cvrun(z,             "flip", "参考: 不分组, 直接预测是否翻转")))
print(r)
log_msg("\n伪 R2 = 1 - 折外对数似然/常数模型对数似然。<= 0 表示还不如直接用基础发生率。")
fwrite(r, file.path(OUT, "class_flip_predictability.csv"))

# ---- 样本量代价 ----
DB <- obj$D
log_msg("\n=== 样本量代价 ===")
log_msg("二分类可用站(样本 B): ", nrow(DB), "  |  四分类可用站(仅严格): ", nrow(D),
        "  |  丢弃 ", nrow(DB) - nrow(D), " 站 (", round(100*(1-nrow(D)/nrow(DB)), 1), "%)")
log_msg("四类最小类 n = ", min(table(D$stype)), " (", names(which.min(table(D$stype))), ")")
log_msg("按 10 折分, 该类平均每折仅 ", round(min(table(D$stype))/10, 1), " 站")
