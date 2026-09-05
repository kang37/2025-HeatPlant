#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 45_class_logit.R — "抑制组 vs 促进组"二分类的 logistic 回归(推断主力)
#
# 三个必须处理的问题:
#   1) 标签信噪比: |mean_coef(tp=0)| 中位仅 0.028, 约一半站符号近乎抛硬币。
#      用 w = |coef0| 作先验权重(99% 缩尾后归一到均值 1), 而不是剔除低信号站
#      ——剔除是按结果筛样本, 会引入选择偏误。
#   2) 空间相关 + 异方差: 标准误用 Conley 空间 HAC(Bartlett 核, 截断 200 km),
#      GLM 版把 OLS 的残差换成得分 X*w*(y-p), 面包换成 (X'WX)^-1。
#   3) 交叉验证必须空间分块: 对经纬度 k-means 聚成 10 折。随机分折会把
#      空间相邻站分到不同折, AUC 偏乐观约 0.05。
#
# 样本 A = 严格四分类站(385, 用户指定主口径); 样本 B = 全部完整站(762, 稳健性)。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
D <- obj$D; PRED <- obj$PRED
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

# ---- AUC(Mann-Whitney, 支持权重) -----------------------------------------
auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0)
  if (!length(p) || !length(n)) return(NA_real_)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(w[i] * w[n] * (d > 0)) +
                   .5 * sum(w[i] * w[n] * (d == 0)); den <- den + sum(w[i] * w[n]) }
  num / den
}

# ---- Conley 空间 HAC(GLM 版) ---------------------------------------------
conley_glm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit)
  u <- as.numeric(fit$prior.weights) * residuals(fit, type = "response")
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  # pmax(0, M) 会丢掉 dim 属性(属性只从第一个参数复制), 必须补回来
  K  <- pmax(0, 1 - sqrt(dx^2 + dy^2) / cut)      # Bartlett, 截断外为 0
  dim(K) <- c(length(lon), length(lon))
  br <- vcov(fit)
  V  <- br %*% (crossprod(S, K %*% S)) %*% br
  n <- nrow(X); k <- ncol(X)
  V * (n / (n - k))
}

# ---- 空间分块折 -----------------------------------------------------------
make_folds <- function(lon, lat, k = 10) {
  set.seed(42)
  kmeans(scale(cbind(lon, lat)), centers = k, nstart = 25, iter.max = 100)$cluster
}

prep <- function(d) {
  z <- copy(d)
  for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  z[, wt := winz(w, .01)][, wt := wt / mean(wt)]
  z
}

fml <- as.formula(paste("grp2 ~", paste(PRED, collapse = "+")))

run_sample <- function(tag, d) {
  z <- prep(d)
  fold <- make_folds(z$longitude, z$latitude, 10)
  ess <- sum(z$wt)^2 / sum(z$wt^2)
  log_msg("\n################ 样本 ", tag, " | n = ", nrow(z),
          " (抑制 ", sum(z$grp2), " / 促进 ", sum(1 - z$grp2), ") ################")
  log_msg("权重有效样本量 ESS = ", round(ess, 1), " (", round(100*ess/nrow(z), 1), "% of n)")

  # ---- VIF ----
  vif <- sapply(PRED, function(v) {
    r2 <- summary(lm(as.formula(paste(v, "~", paste(setdiff(PRED, v), collapse = "+"))),
                     z))$r.squared; 1/(1-r2) })
  log_msg("最大 VIF = ", round(max(vif), 2), " (", names(which.max(vif)), ")")

  res <- list()
  for (wmode in c("unw", "wgt")) {
    z[, ww := if (wmode == "wgt") wt else 1]
    fit <- suppressWarnings(glm(fml, binomial(), z, weights = ww))
    b  <- coef(fit)
    se_n <- sqrt(diag(vcov(fit)))
    se_c <- sqrt(diag(conley_glm(fit, z$longitude, z$latitude, 200)))
    o <- data.table(sample = tag, wmode = wmode, var = names(b), beta = b,
                    se_naive = se_n, se_conley = se_c)
    o[, `:=`(z_c = beta/se_conley, p_naive = 2*pnorm(-abs(beta/se_naive)),
             p_conley = 2*pnorm(-abs(beta/se_conley)),
             or = exp(beta), or_lo = exp(beta - 1.96*se_conley),
             or_hi = exp(beta + 1.96*se_conley))]
    res[[wmode]] <- o
    log_msg("\n--- ", ifelse(wmode == "wgt", "加权(w=|coef0|)", "不加权"),
            " logit | 显著数: 朴素 ", o[var != "(Intercept)" & p_naive < .05, .N],
            " | Conley200 ", o[var != "(Intercept)" & p_conley < .05, .N], " ---")
    pr <- o[var != "(Intercept)"][order(p_conley)][1:8]
    print(pr[, .(var, beta = round(beta, 3), OR = round(or, 2),
                 se_n = round(se_naive, 3), se_C = round(se_conley, 3),
                 p_C = signif(p_conley, 3))])
  }

  # ---- 空间分块 CV ----
  cvres <- rbindlist(lapply(c("unw", "wgt"), function(wmode) {
    pr <- rep(NA_real_, nrow(z))
    for (k in unique(fold)) {
      tr <- z[fold != k]; te <- which(fold == k)
      tr[, ww := if (wmode == "wgt") wt else 1]
      f <- suppressWarnings(glm(fml, binomial(), tr, weights = ww))
      pr[te] <- predict(f, z[te], type = "response")
    }
    data.table(sample = tag, model = paste0("logit_", wmode), pred = pr,
               y = z$grp2, w = z$wt, absco = z$w, fold = fold)
  }))
  list(coef = rbindlist(res), cv = cvres, z = z, fold = fold)
}

A <- run_sample("A", D[strict == TRUE])
B <- run_sample("B", D)

cv <- rbind(A$cv, B$cv)
fwrite(rbind(A$coef, B$coef), file.path(OUT, "class_logit_coef.csv"))
fwrite(cv, file.path(OUT, "class_logit_cv.csv"))

log_msg("\n================ 空间分块 CV 的 AUC ================")
s <- cv[, .(AUC = round(auc_w(y, pred), 3),
            AUC加权 = round(auc_w(y, pred, w), 3)), by = .(sample, model)]
print(s)
log_msg("\n-- 按 |coef0| 三分位分层的 AUC(不加权) --")
cv[, tert := cut(absco, quantile(absco, 0:3/3), include.lowest = TRUE,
                 labels = c("小","中","大")), by = .(sample, model)]
print(dcast(cv[, .(AUC = round(auc_w(y, pred), 3)), by = .(sample, model, tert)],
            sample + model ~ tert, value.var = "AUC"))
log_msg("\n完成")
