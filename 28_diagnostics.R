#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 28_diagnostics.R — 得分回归的共线性与模型假设检验
#
# 检查项:
#   共线性   VIF、设计矩阵条件数、最大成对相关
#   线性性   RESET 检验(遗漏非线性项)
#   同方差   Breusch-Pagan
#   正态性   Shapiro-Wilk + 偏度峰度
#   独立性   残差的 Moran's I(站点数据必查, 手写 k 近邻权重 + 置换检验)
#   强影响点 Cook 距离
#
# Moran's I 自己实现(未装 ape/spdep):
#   I = (n/S0) * sum_ij w_ij (x_i-xbar)(x_j-xbar) / sum_i (x_i-xbar)^2
#   权重用 k 近邻(行标准化), 显著性用置换检验。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(lmtest) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")
M <- 8L

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

morans_I <- function(resid, lon, lat, k = 8, nperm = 999) {
  n <- length(resid)
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  D <- sqrt(dx^2 + dy^2); diag(D) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, order(D[i, ])[1:k]] <- 1
  W <- W / rowSums(W)                       # 行标准化
  z <- resid - mean(resid); S0 <- sum(W)
  I <- (n / S0) * sum(W * outer(z, z)) / sum(z^2)
  perm <- replicate(nperm, { zp <- sample(z)
    (n / S0) * sum(W * outer(zp, zp)) / sum(zp^2) })
  list(I = I, E = -1/(n-1), p = (sum(abs(perm) >= abs(I)) + 1) / (nperm + 1))
}

for (g in c("inhibit_first", "promote_first")) {
  d <- d0[start_dir == g]
  d <- d[complete.cases(d[, c("score", PRED), with = FALSE])]
  lon <- d$longitude; lat <- d$latitude
  z <- copy(d)[, c("score", PRED), with = FALSE]
  for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  f <- as.formula(paste("score ~", paste(PRED, collapse = " + ")))
  m <- lm(f, z)
  lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")
  log_msg("\n#################### ", lab, " | n = ", nrow(z), " ####################")

  log_msg("\n--- 1. 共线性 ---")
  vif <- sapply(PRED, function(v) 1 / (1 - summary(lm(
    as.formula(paste(v, "~", paste(setdiff(PRED, v), collapse = "+"))), z))$r.squared))
  vs <- sort(vif, decreasing = TRUE)
  print(round(vs[1:6], 2))
  log_msg("最大 VIF ", round(max(vif), 2), " (", names(which.max(vif)), ")",
          if (max(vif) > 10) "  <<< 严重" else if (max(vif) > 5) "  (>5 需留意)" else "  正常")
  X <- model.matrix(m)[, -1]
  cn <- sqrt(max(eigen(cor(X))$values) / min(eigen(cor(X))$values))
  log_msg("设计矩阵条件数 ", round(cn, 1), if (cn > 30) "  <<< >30 提示共线" else "  (<30 可接受)")
  cr <- cor(z[, ..PRED]); diag(cr) <- 0
  w <- which(abs(cr) == max(abs(cr)), arr.ind = TRUE)[1, ]
  log_msg("最大成对相关 ", round(cr[w[1], w[2]], 3), " (",
          PRED[w[1]], " ~ ", PRED[w[2]], ")")

  log_msg("\n--- 2. 线性性 (RESET) ---")
  rs <- resettest(m, power = 2:3, type = "fitted")
  log_msg(sprintf("F = %.2f, p = %.4f  %s", rs$statistic, rs$p.value,
                  ifelse(rs$p.value < .05, "<<< 拒绝线性, 存在遗漏的非线性", "未拒绝")))

  log_msg("\n--- 3. 同方差 (Breusch-Pagan) ---")
  bp <- bptest(m)
  log_msg(sprintf("BP = %.2f, df = %d, p = %.4f  %s", bp$statistic, bp$parameter,
                  bp$p.value, ifelse(bp$p.value < .05, "<<< 异方差", "未拒绝同方差")))

  log_msg("\n--- 4. 残差正态性 ---")
  r <- residuals(m); sw <- shapiro.test(r)
  sk <- mean((r - mean(r))^3) / sd(r)^3; ku <- mean((r - mean(r))^4) / sd(r)^4 - 3
  log_msg(sprintf("Shapiro-Wilk W = %.4f, p = %.5f  %s", sw$statistic, sw$p.value,
                  ifelse(sw$p.value < .05, "<<< 拒绝正态", "未拒绝")))
  log_msg(sprintf("偏度 %.3f | 超额峰度 %.3f", sk, ku))

  log_msg("\n--- 5. 残差空间自相关 (Moran's I, k=8 近邻, 999 次置换) ---")
  mi <- morans_I(r, lon, lat)
  log_msg(sprintf("I = %.4f (期望 %.4f), p = %.4f  %s", mi$I, mi$E, mi$p,
                  ifelse(mi$p < .05, "<<< 残差存在空间自相关, 违反独立性", "未检出")))
  yi <- morans_I(z$score, lon, lat)
  log_msg(sprintf("(对照: 因变量本身 I = %.4f, p = %.4f)", yi$I, yi$p))

  log_msg("\n--- 6. 强影响点 (Cook 距离) ---")
  ck <- cooks.distance(m); thr <- 4 / nrow(z)
  log_msg(sprintf("超过 4/n=%.4f 的点: %d (%.1f%%) | 最大 %.3f | 超过 1 的: %d",
                  thr, sum(ck > thr), 100 * mean(ck > thr), max(ck), sum(ck > 1)))
}
