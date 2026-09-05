#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 46_class_xgb.R — "抑制组 vs 促进组"的 XGBoost + TreeSHAP(探索层, 不是推断层)
#
# 定位: 用来发现非线性形状与交互, 发现之后回到 logistic 模型里加样条/交互项检验。
# 树模型不受成分数据共线之困, 故这里把耕地也放回来(9 个地类全入), 共 19 个特征。
#
# 三件必须做对的事:
#   1) 空间分块 CV: 折 = 经纬度 k-means 10 类, 与 45 号脚本同一套折, 可直接比。
#   2) 折外 SHAP: 每个站的 SHAP 由"没见过它的那个模型"算出, 避免过拟合归因。
#   3) nrounds 用嵌套空间 CV 定, 不用固定值; n 小特征多, 参数取强正则。
# TreeSHAP 用 xgboost 自带 predcontrib(精确解, 非抽样近似), 单位是对数几率。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(xgboost) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
D <- obj$D
FEAT <- c(obj$PRED, "crop")            # 树模型不怕成分共线, 把参照类耕地放回来

auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(w[i]*w[n]*(d > 0)) +
                   .5*sum(w[i]*w[n]*(d == 0)); den <- den + sum(w[i]*w[n]) }
  num/den
}
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }
make_folds <- function(lon, lat, k = 10) {
  set.seed(42); kmeans(scale(cbind(lon, lat)), centers = k, nstart = 25, iter.max = 100)$cluster }

PAR <- list(objective = "binary:logistic", eval_metric = "logloss",
            max_depth = 3, eta = 0.03, min_child_weight = 8,
            subsample = 0.8, colsample_bytree = 0.8,
            lambda = 3, alpha = 0, nthread = 4)

run_xgb <- function(tag, d, use_w = TRUE) {
  X <- as.matrix(d[, FEAT, with = FALSE])
  y <- d$grp2
  wt <- winz(d$w, .01); wt <- wt / mean(wt)
  W <- if (use_w) wt else rep(1, nrow(d))
  fold <- make_folds(d$longitude, d$latitude, 10)

  pred <- rep(NA_real_, nrow(d))
  shap <- matrix(NA_real_, nrow(d), length(FEAT) + 1,
                 dimnames = list(NULL, c(FEAT, "BIAS")))
  nr <- integer(0)
  for (k in sort(unique(fold))) {
    ti <- which(fold != k); vi <- which(fold == k)
    dtr <- xgb.DMatrix(X[ti, ], label = y[ti], weight = W[ti])
    # nrounds: 训练集内部再按空间块折一次(嵌套), 避免用测试折信息定轮数
    inner <- lapply(sort(unique(fold[ti])), function(j) which(fold[ti] == j))
    cvf <- xgb.cv(PAR, dtr, nrounds = 800, folds = inner, early_stopping_rounds = 40,
                  verbose = 0)
    nb <- max(10L, cvf$early_stop$best_iteration); nr <- c(nr, nb)
    bst <- xgb.train(PAR, dtr, nrounds = nb, verbose = 0)
    pred[vi] <- predict(bst, X[vi, , drop = FALSE])
    shap[vi, ] <- predict(bst, X[vi, , drop = FALSE], predcontrib = TRUE)
  }
  log_msg("\n#### 样本 ", tag, " | n = ", nrow(d), " | 权重: ",
          ifelse(use_w, "|coef0|", "无"), " | 各折最佳轮数 ",
          paste(range(nr), collapse = "-"), " (中位 ", median(nr), ") ####")
  log_msg("空间分块 CV AUC = ", round(auc_w(y, pred), 3),
          " | 加权 AUC = ", round(auc_w(y, pred, wt), 3))

  # 折外 SHAP 汇总
  sl <- as.data.table(shap[, FEAT, drop = FALSE])
  sl[, stat_id := d$stat_id]
  sl <- melt(sl, id.vars = "stat_id", variable.name = "var", value.name = "shap")
  fv <- melt(d[, c("stat_id", FEAT), with = FALSE], id.vars = "stat_id",
             variable.name = "var", value.name = "x")
  sl <- merge(sl, fv, by = c("stat_id", "var"))
  sl[, `:=`(sample = tag, grp2 = d$grp2[match(stat_id, d$stat_id)])]
  imp <- sl[, .(mean_abs = mean(abs(shap)),
                cor_xs = suppressWarnings(cor(x, shap, method = "spearman"))), by = var]
  imp[, rel := mean_abs / sum(mean_abs)][order(-mean_abs)]
  setorder(imp, -mean_abs)
  log_msg("-- SHAP 重要性(前 10, 单位=对数几率) --")
  print(imp[1:10, .(var, mean_abs = round(mean_abs, 3), rel = round(rel, 3),
                    方向 = ifelse(cor_xs > 0, "越大越偏抑制", "越大越偏促进"),
                    rho = round(cor_xs, 2))])

  # 全样本模型: 用于交互检测(逐折交互噪声太大)
  dall <- xgb.DMatrix(X, label = y, weight = W)
  inner <- lapply(sort(unique(fold)), function(j) which(fold == j))
  nbf <- max(10L, xgb.cv(PAR, dall, nrounds = 800, folds = inner,
                         early_stopping_rounds = 40, verbose = 0)$early_stop$best_iteration)
  bfull <- xgb.train(PAR, dall, nrounds = nbf, verbose = 0)
  ia <- predict(bfull, X, predinteraction = TRUE)          # n x (p+1) x (p+1)
  nm <- c(FEAT, "BIAS")
  M <- apply(abs(ia), c(2, 3), mean)
  dimnames(M) <- list(nm, nm)
  iu <- as.data.table(as.table(M))[V1 != "BIAS" & V2 != "BIAS" & V1 != V2]
  setnames(iu, c("v1", "v2", "val"))
  iu[, key := paste(pmin(as.character(v1), as.character(v2)),
                    pmax(as.character(v1), as.character(v2)))]
  iu <- iu[, .(strength = sum(val)), by = key][order(-strength)][1:10]
  log_msg("-- 最强的 10 个成对交互(|SHAP 交互值|均值之和) --")
  print(iu)

  list(shap = sl, imp = imp[, sample := tag], inter = iu[, sample := tag],
       cv = data.table(sample = tag, stat_id = d$stat_id, y = y, pred = pred,
                       w = wt, absco = d$w, fold = fold))
}

A <- run_xgb("A", D[strict == TRUE])
B <- run_xgb("B", D)

fwrite(rbind(A$shap, B$shap),  file.path(OUT, "class_xgb_shap.csv"))
fwrite(rbind(A$imp,  B$imp),   file.path(OUT, "class_xgb_imp.csv"))
fwrite(rbind(A$inter, B$inter),file.path(OUT, "class_xgb_inter.csv"))
fwrite(rbind(A$cv,   B$cv),    file.path(OUT, "class_xgb_cv.csv"))

cv <- rbind(A$cv, B$cv)
log_msg("\n-- 按 |coef0| 三分位分层的 AUC --")
cv[, tert := cut(absco, quantile(absco, 0:3/3), include.lowest = TRUE,
                 labels = c("小","中","大")), by = sample]
print(dcast(cv[, .(AUC = round(auc_w(y, pred), 3)), by = .(sample, tert)],
            sample ~ tert, value.var = "AUC"))
log_msg("\n完成")
