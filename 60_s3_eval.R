#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 60_s3_eval.R — S3(置换检验判据)的外部验证
# 方向层还准不准? 翻转层在"只保留真变点"之后是不是终于可预测?
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(xgboost) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")
obj <- readRDS(file.path(OUT, "class_model_data.rds"))
PRED <- obj$PRED; FEAT <- c(PRED, "crop")
s3 <- fread(file.path(OUT, "reclassify_s3.csv"))
D <- merge(obj$D, s3, by = "stat_id")
D[, `:=`(dir3 = as.integer(S3_start == "inhibit"), flip3 = as.integer(!is.na(S3_k)))]
winz <- function(x, p = .01) { q <- quantile(x, c(p,1-p), na.rm=TRUE); pmin(pmax(x,q[1]),q[2]) }
auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y==1); n <- which(y==0); if(!length(p)||!length(n)) return(NA_real_)
  num<-0; den<-0
  for (i in p){d<-sign(s[i]-s[n]); num<-num+sum(w[i]*w[n]*(d>0))+.5*sum(w[i]*w[n]*(d==0))
               den<-den+sum(w[i]*w[n])}; num/den }
PAR <- list(objective="binary:logistic", eval_metric="logloss", max_depth=3, eta=.03,
            min_child_weight=8, subsample=.8, colsample_bytree=.8, lambda=3, nthread=4)
z <- copy(D); for (v in FEAT) set(z, j=v, value=as.numeric(scale(winz(as.numeric(z[[v]])))))
z[, fold := { set.seed(42); kmeans(scale(cbind(longitude, latitude)), 10,
                                   nstart=25, iter.max=100)$cluster }]
z[, wt := winz(w,.01)][, wt := wt/mean(wt)]
cvrun <- function(d, yv, tag, useW=FALSE) {
  d <- copy(d); d[, yy := d[[yv]]]; d[, ww := if (useW) wt else 1]
  pl <- px <- rep(NA_real_, nrow(d))
  f <- as.formula(paste("yy ~", paste(PRED, collapse="+")))
  X <- as.matrix(d[, FEAT, with=FALSE])
  for (k in sort(unique(d$fold))) {
    ti <- which(d$fold!=k); vi <- which(d$fold==k)
    if (!length(vi) || length(unique(d$yy[ti]))<2) next
    pl[vi] <- suppressWarnings(predict(glm(f, binomial(), d[ti], weights=ww), d[vi], type="response"))
    dtr <- xgb.DMatrix(X[ti,], label=d$yy[ti], weight=d$ww[ti])
    inner <- lapply(sort(unique(d$fold[ti])), function(j) which(d$fold[ti]==j))
    nb <- max(10L, xgb.cv(PAR, dtr, nrounds=600, folds=inner,
                          early_stopping_rounds=40, verbose=0)$early_stop$best_iteration)
    px[vi] <- predict(xgb.train(PAR, dtr, nrounds=nb, verbose=0), X[vi,,drop=FALSE])
  }
  ok <- !is.na(pl)
  data.table(设定=tag, n=sum(ok), 阳性率=round(mean(d$yy[ok]),3),
             AUC_logit=round(auc_w(d$yy[ok], pl[ok]),3),
             AUC_xgb=round(auc_w(d$yy[ok], px[ok]),3)) }
log_msg("=== 方向层 ===\n")
print(rbindlist(list(
  cvrun(z, "grp2", "tp=0 符号(762)", TRUE),
  cvrun(z, "dir3", "S3 方向(762)",  TRUE))))
log_msg("\n=== 翻转层: 只保留置换检验通过的变点 ===\n")
print(rbindlist(list(
  cvrun(z[dir3==1], "flip3", "抑制支内是否有真变点"),
  cvrun(z[dir3==0], "flip3", "促进支内是否有真变点"),
  cvrun(z,          "flip3", "不分组: 该站有没有真变点"))))
log_msg("\n方向一致性: S3 方向 与 tp=0 符号一致 ",
        sprintf("%.1f%%", 100*mean(D$dir3 == D$grp2)), " (", sum(D$dir3==D$grp2), "/", nrow(D), ")")
