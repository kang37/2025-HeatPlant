#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 47_class_four.R — 四分类(始终抑制 / 抑制转促进 / 始终促进 / 促进转抑制)
#
# 四类只在严格分类站上有定义, 故本脚本只用样本 A(385 站)。
#
# 三条并行路线:
#   (1) 分解式: 四类 = 起始方向(2) x 后续是否翻转(2), 两层都是二元 logit,
#       都能上 Conley 空间 HAC 标准误。P(四类) = P(方向) x P(翻转|方向)。
#       这是唯一能给出可信 p 值的路线——始终抑制只有二十来站, 直接上多项
#       logit 会把 3x19 个参数压在这点样本上。
#   (2) 多项 logit(nnet::multinom): 直接对照, 报系数但标准误按站聚类无从修,
#       只作定性参考, 不作推断依据。
#   (3) XGBoost multi:softprob + 逐类 TreeSHAP: 看每一类各自被什么驱动。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(nnet); library(xgboost) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
PRED <- obj$PRED; FEAT <- c(PRED, "crop")
D <- obj$D[strict == TRUE]
LAB <- c(always_inhibit = "始终抑制", inhibit_promote = "抑制转促进",
         always_promote = "始终促进", promote_inhibit = "促进转抑制")
D[, stype := factor(stype, levels = names(LAB))]
D[, flip := as.integer(stype %in% c("inhibit_promote", "promote_inhibit"))]

winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }
z <- copy(D); for (v in FEAT) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
z[, wt := winz(w, .01)][, wt := wt/mean(wt)]
fold <- { set.seed(42); kmeans(scale(cbind(z$longitude, z$latitude)), 10,
                               nstart = 25, iter.max = 100)$cluster }

log_msg("-- 四类样本量(样本 A, 协变量完整) --")
print(z[, .(n = .N, 占比 = round(.N/nrow(z), 3)), by = .(类型 = LAB[as.character(stype)])])

conley_glm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit)
  u <- as.numeric(fit$prior.weights) * residuals(fit, type = "response")
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat)*pi/180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  K <- pmax(0, 1 - sqrt(dx^2 + dy^2)/cut); dim(K) <- c(length(lon), length(lon))
  br <- vcov(fit); V <- br %*% crossprod(S, K %*% S) %*% br
  V * (nrow(X)/(nrow(X) - ncol(X)))
}
auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0); if (!length(p) || !length(n)) return(NA_real_)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(w[i]*w[n]*(d>0)) +
                   .5*sum(w[i]*w[n]*(d==0)); den <- den + sum(w[i]*w[n]) }
  num/den }

# ===== (1) 分解式: 方向层 + 翻转层 ========================================
fit_layer <- function(dd, yv, tag, use_w) {
  dd <- copy(dd); dd[, yy := dd[[yv]]]
  dd[, ww := if (use_w) wt else 1]
  f <- suppressWarnings(glm(as.formula(paste("yy ~", paste(PRED, collapse = "+"))),
                            binomial(), dd, weights = ww))
  se <- sqrt(diag(conley_glm(f, dd$longitude, dd$latitude, 200)))
  o <- data.table(layer = tag, var = names(coef(f)), beta = coef(f), se_conley = se)
  o[, `:=`(p = 2*pnorm(-abs(beta/se_conley)), or = exp(beta))]
  ns <- o[var != "(Intercept)" & p < .05, .N]
  # 分离诊断: 阳性率极端 + 拟合概率贴边 + 系数爆炸 => 系数与 p 值不可信
  ph <- fitted(f); sep <- mean(ph < .001 | ph > .999)
  bal <- min(mean(dd$yy), 1 - mean(dd$yy))
  log_msg("\n### ", tag, " | n = ", nrow(dd), " | 阳性 ", sum(dd$yy),
          " (", round(100*mean(dd$yy), 1), "%) | Conley 显著 ", ns, " 个 ###")
  log_msg("诊断: 拟合概率贴边占比 ", round(sep, 3),
          " | 最大 |beta| = ", round(max(abs(o$beta[-1])), 2),
          " | 少数类占比 ", round(bal, 3),
          if (bal < .15) "   <<< 少数类过小" else "",
          if (sep > .10 && max(abs(o$beta[-1])) > 4) "   <<< 疑似准完全分离" else "")
  print(o[var != "(Intercept)"][order(p)][1:6,
        .(var, beta = round(beta,3), OR = round(or,2), p = signif(p,3))])
  o
}
# 稀疏度: 占比近乎恒为 0 的地类, 标准化后极端值支配拟合, 大 |beta| 多半是这个原因
spa <- data.table(var = PRED, 非零占比 = sapply(PRED, function(v) mean(D[[v]] != 0)))
log_msg("\n-- 自变量非零占比(仅列 < 0.5 的) --")
print(spa[非零占比 < .5][order(非零占比)])

L1 <- fit_layer(z, "grp2", "第1层 起始方向(1=抑制)", TRUE)
L2i <- fit_layer(z[grp2 == 1], "flip", "第2层a 抑制组内是否翻转", FALSE)
L2p <- fit_layer(z[grp2 == 0], "flip", "第2层b 促进组内是否翻转", FALSE)
fwrite(rbindlist(list(L1, L2i, L2p)), file.path(OUT, "class4_layer_coef.csv"))

# ===== (2) 多项 logit ======================================================
mn <- multinom(as.formula(paste("stype ~", paste(PRED, collapse = "+"))),
               z, weights = wt, trace = FALSE, maxit = 1000)
cf <- coef(mn); sem <- summary(mn)$standard.errors
mo <- rbindlist(lapply(rownames(cf), function(r)
  data.table(vs = r, var = colnames(cf), beta = cf[r, ], se = sem[r, ])))
mo[, p := 2*pnorm(-abs(beta/se))]
log_msg("\n### 多项 logit(参照类 = 始终抑制; 标准误未做空间修正, 仅作参考) ###")
print(mo[var != "(Intercept)"][order(p)][1:12,
      .(对比 = LAB[vs], var, beta = round(beta,2), p = signif(p,3))])
fwrite(mo, file.path(OUT, "class4_multinom_coef.csv"))

# ===== (3) XGBoost 四分类 + 逐类 SHAP ======================================
X <- as.matrix(z[, FEAT, with = FALSE]); ycls <- as.integer(z$stype) - 1L
PAR <- list(objective = "multi:softprob", num_class = 4, eval_metric = "mlogloss",
            max_depth = 3, eta = 0.05, min_child_weight = 8, subsample = .8,
            colsample_bytree = .8, lambda = 3, nthread = 4)
P <- matrix(NA_real_, nrow(z), 4)
SH <- array(NA_real_, c(nrow(z), length(FEAT)+1, 4))
for (k in sort(unique(fold))) {
  ti <- which(fold != k); vi <- which(fold == k)
  dtr <- xgb.DMatrix(X[ti, ], label = ycls[ti], weight = z$wt[ti])
  inner <- lapply(sort(unique(fold[ti])), function(j) which(fold[ti] == j))
  nb <- max(10L, xgb.cv(PAR, dtr, nrounds = 600, folds = inner,
                        early_stopping_rounds = 40, verbose = 0)$early_stop$best_iteration)
  bst <- xgb.train(PAR, dtr, nrounds = nb, verbose = 0)
  P[vi, ] <- predict(bst, X[vi, , drop = FALSE])
  # 多分类的 predcontrib 返回 n x num_class x (p+1), 不是 list
  ct <- predict(bst, X[vi, , drop = FALSE], predcontrib = TRUE)
  for (c4 in 1:4) SH[vi, , c4] <- ct[, c4, ]
}
acc <- mean(max.col(P) == (ycls + 1L))
base <- max(table(ycls))/length(ycls)
log_msg("\n### XGBoost 四分类(空间分块 CV) ###")
log_msg("总体准确率 ", round(acc,3), " | 全猜最大类的基准 ", round(base,3))
ova <- sapply(1:4, function(c4) auc_w(as.integer(ycls == c4-1), P[, c4]))
names(ova) <- LAB[levels(z$stype)]
log_msg("逐类 one-vs-rest AUC:")
print(round(ova, 3))
log_msg("宏平均 AUC = ", round(mean(ova), 3))
cm <- table(真实 = LAB[levels(z$stype)][ycls+1], 预测 = LAB[levels(z$stype)][max.col(P)])
log_msg("\n-- 混淆矩阵 --"); print(cm)

impc <- rbindlist(lapply(1:4, function(c4) {
  s <- SH[, seq_along(FEAT), c4, drop = TRUE]
  data.table(cls = LAB[levels(z$stype)][c4], var = FEAT,
             mean_abs = colMeans(abs(s)),
             cor_xs = sapply(seq_along(FEAT), function(j)
               suppressWarnings(cor(X[, j], s[, j], method = "spearman"))))
}))
impc[, rel := mean_abs/sum(mean_abs), by = cls]
log_msg("\n-- 每类 SHAP 重要性前 5 --")
for (cc in unique(impc$cls)) {
  log_msg("[", cc, "]")
  print(impc[cls == cc][order(-mean_abs)][1:5, .(var, mean_abs = round(mean_abs,3),
        方向 = ifelse(cor_xs > 0, "越大越像这类", "越大越不像这类"))])
}
shl <- rbindlist(lapply(1:4, function(c4) {
  s <- as.data.table(SH[, seq_along(FEAT), c4, drop = TRUE]); setnames(s, FEAT)
  s[, `:=`(stat_id = z$stat_id, cls = LAB[levels(z$stype)][c4])]
  melt(s, id.vars = c("stat_id","cls"), variable.name = "var", value.name = "shap") }))
fv <- melt(D[, c("stat_id", FEAT), with = FALSE], id.vars = "stat_id",
           variable.name = "var", value.name = "x")
shl <- merge(shl, fv, by = c("stat_id","var"))
fwrite(shl,  file.path(OUT, "class4_xgb_shap.csv"))
fwrite(impc, file.path(OUT, "class4_xgb_imp.csv"))
fwrite(data.table(stat_id = z$stat_id, stype = as.character(z$stype),
                  grp2 = z$grp2, flip = z$flip, fold = fold,
                  p1 = P[,1], p2 = P[,2], p3 = P[,3], p4 = P[,4]),
       file.path(OUT, "class4_xgb_cv.csv"))
log_msg("\n完成")
