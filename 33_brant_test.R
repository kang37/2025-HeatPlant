#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 33_brant_test.R — 有序 logit 的比例优势假设检验(Brant)
#
# 有序 logit 相对 OLS 的优势是不假设各等级等距, 但它有自己的关键假设:
# 比例优势(proportional odds)——每个自变量的效应在所有切点上相同。
# 该假设若不成立, 有序 logit 的单一系数就无法解释。
#
# brant / VGAM 未安装, 故手写: 在各切点拟合二元 logit, 用嵌套事件的
# 联合协方差构造 Wald 统计量检验系数跨切点相等。
#   嵌套性 {Y>c_1} ⊇ {Y>c_2} (c_1<c_2) 使 P(both) = P(较严事件),
#   故 Cov 的中间项对角元为 min(pi_j, pi_k) - pi_j*pi_k。
# ---------------------------------------------------------------------------

Sys.setlocale("LC_ALL","en_US.UTF-8")
suppressPackageStartupMessages(library(data.table))
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
M <- 8L
env <- new.env(); L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1]-1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text=paste(src,collapse="\n")), envir=env)))
pat <- env$pat; dt <- env$dt; PRED <- env$PRED; winz <- env$winz
pat[, score := fifelse(start_dir=="inhibit_first",
       fifelse(event==1,(M+2)-event_tp,1), fifelse(event==1,as.numeric(event_tp),M+1))]
d0 <- merge(pat[,.(stat_id,start_dir,score)], dt[,c("stat_id",PRED),with=FALSE], by="stat_id")

# Brant 检验: 在各切点拟合二元 logit, 检验系数是否跨切点相等
# 嵌套事件 {Y>c} 使 Cov(b_j,b_k) = (X'W_j X)^-1 X'W_jk X (X'W_k X)^-1,
# 其中 W_jk 的对角元为 pi_max(j,k) - pi_j*pi_k
brant <- function(d, vars, cuts) {
  X <- cbind(1, as.matrix(d[, ..vars])); k <- ncol(X); J <- length(cuts)
  fits <- lapply(cuts, function(c) glm(as.numeric(d$score > c) ~ X - 1, family=binomial))
  b <- lapply(fits, coef); pi <- lapply(fits, fitted)
  Wm <- function(j,k2) {
    pj <- pi[[j]]; pk <- pi[[k2]]; pmax_ <- pmin(pj, pk)   # 嵌套: P(both)=P(较严的)
    diag(pmax_ - pj*pk) }
  br <- lapply(fits, function(f) solve(crossprod(X, diag(fitted(f)*(1-fitted(f))) %*% X)))
  V <- matrix(0, J*k, J*k)
  for (i in 1:J) for (j2 in 1:J) {
    blk <- br[[i]] %*% crossprod(X, Wm(i,j2) %*% X) %*% br[[j2]]
    V[((i-1)*k+1):(i*k), ((j2-1)*k+1):(j2*k)] <- blk }
  bb <- unlist(b)
  # 对每个自变量: 检验其在 J 个切点上的系数全部相等
  res <- rbindlist(lapply(2:k, function(p) {
    D <- matrix(0, J-1, J*k)
    for (i in 2:J) { D[i-1, p] <- 1; D[i-1, (i-1)*k+p] <- -1 }
    s <- as.numeric(D %*% bb); Vd <- D %*% V %*% t(D)
    chi <- tryCatch(as.numeric(t(s) %*% solve(Vd) %*% s), error=function(e) NA)
    data.table(var=vars[p-1], chi2=chi, df=J-1, p=pchisq(chi, J-1, lower.tail=FALSE),
               系数范围=sprintf("%.2f ~ %.2f", min(sapply(b,`[`,p)), max(sapply(b,`[`,p)))) }))
  # 总体检验
  D <- matrix(0, (J-1)*(k-1), J*k); r <- 0
  for (p in 2:k) for (i in 2:J) { r <- r+1; D[r,p] <- 1; D[r,(i-1)*k+p] <- -1 }
  s <- as.numeric(D %*% bb)
  omni <- tryCatch({ chi <- as.numeric(t(s) %*% solve(D %*% V %*% t(D)) %*% s)
    c(chi, (J-1)*(k-1), pchisq(chi,(J-1)*(k-1),lower.tail=FALSE)) }, error=function(e) rep(NA,3))
  list(byvar=res, omni=omni)
}

KEY <- c("urban_rate","mixedleaf","tavg","invest","water","imperv","elev","rh")
for (g in c("inhibit_first","promote_first")) {
  d <- d0[start_dir==g]; d <- d[complete.cases(d[,c("score",PRED),with=FALSE])]
  for (v in PRED) set(d, j=v, value=as.numeric(scale(winz(as.numeric(d[[v]])))))
  cuts <- c(3,5,7)     # 三个切点, 兼顾稳定与覆盖
  cat("\n##########", ifelse(g=="inhibit_first","抑制组","促进组"),
      "| n =", nrow(d), "| 切点 score>", paste(cuts,collapse=","), "##########\n")
  cat("各切点的事件比例:", paste(sprintf("%.2f", sapply(cuts, function(c) mean(d$score>c))), collapse=" "), "\n\n")
  r <- brant(d, KEY, cuts)
  cat("--- 逐变量的比例优势假设检验 ---\n")
  print(r$byvar[order(p), .(var, chi2=round(chi2,2), df, p=signif(p,3),
                            违反=fifelse(p<.05,"*",""), 系数范围)])
  cat(sprintf("\n总体检验: chi2 = %.1f, df = %d, p = %.4g  %s\n", r$omni[1], r$omni[2], r$omni[3],
      ifelse(!is.na(r$omni[3]) && r$omni[3]<.05, "<<< 拒绝比例优势假设", "未拒绝")))
}
