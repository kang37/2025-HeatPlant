#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 29_robust_se.R — 三种标准误的对照: 普通 / HC3 / Conley 空间 HAC
#
# 28 号诊断查出两组都有异方差, 促进组还有残差空间自相关。这两种违反需要
# 不同的修正:
#
#   HC3    只修异方差。把误差方差随观测变化考虑进去, 并按杠杆值 (1-h_i)^2
#          放大, 小样本下比 HC0/HC1 保守。仍假设观测之间互相独立。
#   Conley 同时修异方差与空间相关。用距离核给邻近站点的残差乘积加权,
#          相当于二维版的 Newey-West。截断距离外权重为 0。
#          截断距离 -> 0 时 Conley 退化为 HC0, 故它是 HC3 的推广。
#
# 系数在三种方法下完全相同, 变的只有标准误与 p 值。
#
# 用法: Rscript 29_robust_se.R [--cut 200]   (Conley 截断距离, km)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(sandwich) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
aa <- commandArgs(trailingOnly = TRUE)
CUTS <- if (any(aa == "--cut")) as.numeric(aa[which(aa == "--cut") + 1]) else c(100, 200, 400)
M <- 8L
log_msg <- function(...) cat(..., "\n", sep = "")

env <- new.env(); L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
pat <- env$pat; dt <- env$dt; PRED <- env$PRED; winz <- env$winz
pat[, score := fifelse(start_dir == "inhibit_first",
                       fifelse(event == 1, (M + 2) - event_tp, 1),
                       fifelse(event == 1, as.numeric(event_tp), M + 1))]
st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
d0 <- merge(merge(pat[, .(stat_id, start_dir, score)],
                  dt[, c("stat_id", PRED), with = FALSE], by = "stat_id"),
            st[, .(stat_id, longitude, latitude)], by = "stat_id")

# Conley 空间 HAC: Bartlett 核, 截断距离 cut(km)
conley_vcov <- function(m, lon, lat, cut) {
  X <- model.matrix(m); u <- residuals(m); n <- nrow(X)
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  D  <- sqrt(dx^2 + dy^2)
  K  <- pmax(0, 1 - D / cut)                    # Bartlett 权重, 截断外为 0
  meat <- crossprod(X, (K * outer(u, u)) %*% X)
  br <- solve(crossprod(X))
  V <- br %*% meat %*% br
  V * (n / (n - ncol(X)))                       # 小样本校正
}

res <- list()
for (g in c("inhibit_first", "promote_first")) {
  d <- d0[start_dir == g]
  d <- d[complete.cases(d[, c("score", PRED), with = FALSE])]
  z <- copy(d)[, c("score", PRED), with = FALSE]
  for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  m <- lm(as.formula(paste("score ~", paste(PRED, collapse = "+"))), z)
  b <- coef(m)[PRED]
  se0 <- sqrt(diag(vcov(m)))[PRED]
  se3 <- sqrt(diag(vcovHC(m, type = "HC3")))[PRED]
  out <- data.table(grp = g, var = PRED, beta = b, se_ols = se0, se_hc3 = se3)
  for (cut in CUTS) {
    sc <- sqrt(diag(conley_vcov(m, d$longitude, d$latitude, cut)))[PRED]
    out[, paste0("se_conley", cut) := sc]
  }
  pv <- function(s) 2 * pnorm(-abs(b / s))
  out[, `:=`(p_ols = pv(se0), p_hc3 = pv(se3))]
  for (cut in CUTS) out[, paste0("p_conley", cut) := pv(out[[paste0("se_conley", cut)]])]
  res[[g]] <- out

  lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")
  log_msg("\n############ ", lab, " | n = ", nrow(z), " ############")
  cc <- paste0("p_conley", CUTS)
  pr <- out[order(p_hc3), c("var", "beta", "p_ols", "p_hc3", cc), with = FALSE]
  pr[, (c("beta")) := round(beta, 3)]
  for (v in c("p_ols", "p_hc3", cc)) pr[, (v) := signif(get(v), 3)]
  print(pr[1:10])
  log_msg("显著数: 普通 ", out[p_ols < .05, .N], " | HC3 ", out[p_hc3 < .05, .N],
          " | Conley ", paste(sapply(CUTS, function(k)
            sprintf("%dkm=%d", k, out[get(paste0("p_conley", k)) < .05, .N])), collapse = " "))
}
fwrite(rbindlist(res), file.path(OUT, sprintf("robust_se_tp%d.csv", M)))
log_msg("\n已写出 robust_se_tp", M, ".csv")
