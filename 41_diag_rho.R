#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 41_diag_rho.R — 因果强度回归(rho_tp0 / rho_max / coef_tp0)的模型假设诊断
#
# 背景: 转变时间得分那一支已做过完整诊断与修正(28/29/31/32), 但因果强度这一支
# (20_drivers_all.R / 22_varpart_3dim.R)至今只用普通 OLS 标准误, 从未检验过
# 同方差、残差独立性与函数形式。本脚本把 28 号那套诊断原样搬过来。
#
# 检查项与 28_diagnostics.R 一致, 另加两项:
#   White 检验   BP 只对拟合值做辅助回归容易漏掉与 x 的交互型异方差
#   影响点重拟合 只报 Cook 距离个数说明不了问题, 直接剔除后重跑, 看系数与
#                显著性是否变化——"强影响点检验"的落点应该是结论稳不稳
#   非线性定位   RESET 拒绝之后逐变量加二次项, 找出非线性到底在哪个变量上,
#                并检验加入后线性项的结论是否改变
#
# 设定与 22_varpart_3dim.R 的 run(CORE, ...) 完全一致:
#   同一 18 变量 CORE 集、同样 1% 双侧缩尾、自变量标准化而因变量保持原尺度,
#   所以点估计与 22 的输出逐位相同, 本脚本只是在同一个模型上做诊断。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(lmtest) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

# ===== 复用 22 号的数据装配(执行到 winz 定义为止, 不跑后面的方差分解) =====
env <- new.env(); L <- readLines("22_varpart_3dim.R")
src <- L[1:grep("^winz <- function", L)[1]]
invisible(capture.output(suppressWarnings(
  eval(parse(text = paste(src, collapse = "\n")), envir = env))))
dt <- env$dt; CORE <- env$CORE; winz <- env$winz; st <- env$st
VARS <- unlist(CORE, use.names = FALSE)
stopifnot(length(VARS) == 18L)
dt <- merge(dt, st[, .(stat_id, longitude, latitude)], by = "stat_id", all.x = TRUE)

RESP <- c(rho_tp0 = "rho_tp0（tp=0 因果强度）",
          rho_max = "rho_max（跨滞后最强因果强度）",
          coef_tp0 = "coef_tp0（tp=0 的 S-map 偏导数）")

distmat <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2)
}
# Moran's I: k 近邻行标准化权重 + 置换检验。与 28 号同一估计量,
# 只是把 sum(W * outer(z,z)) 换成等价的二次型, n=768 下快两个量级。
morans_I <- function(r, lon, lat, k = 8, nperm = 999) {
  n <- length(r)
  D <- distmat(lon, lat); diag(D) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, order(D[i, ])[1:k]] <- 1
  W <- W / rowSums(W); S0 <- sum(W)
  z <- r - mean(r)
  stat <- function(v) (n / S0) * as.numeric(crossprod(v, W %*% v)) / sum(v^2)
  I <- stat(z)
  perm <- replicate(nperm, stat(sample(z)))
  list(I = I, E = -1 / (n - 1), p = (sum(abs(perm) >= abs(I)) + 1) / (nperm + 1))
}

summ <- list()
for (yv in names(RESP)) {
  cols <- c(yv, VARS)
  d <- dt[complete.cases(dt[, ..cols])]
  z <- d[, ..cols]
  for (v in VARS) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  f <- as.formula(paste(yv, "~", paste(VARS, collapse = " + ")))
  m <- lm(f, z)
  n <- nrow(z)
  log_msg("\n#################### ", RESP[yv], " | n = ", n, " ####################")
  log_msg("R2 = ", round(summary(m)$r.squared, 4),
          " | 调整 R2 = ", round(summary(m)$adj.r.squared, 4))

  log_msg("\n--- 1. 共线性 ---")
  vif <- sapply(VARS, function(v) 1 / (1 - summary(lm(
    as.formula(paste(v, "~", paste(setdiff(VARS, v), collapse = "+"))), z))$r.squared))
  print(round(sort(vif, decreasing = TRUE)[1:6], 2))
  log_msg("最大 VIF ", round(max(vif), 2), " (", names(which.max(vif)), ")",
          if (max(vif) > 10) "  <<< 严重" else if (max(vif) > 5) "  (>5 需留意)" else "  正常")
  X <- model.matrix(m)[, -1]
  ev <- eigen(cor(X))$values
  cn <- sqrt(max(ev) / min(ev))
  log_msg("设计矩阵条件数 ", round(cn, 1), if (cn > 30) "  <<< >30 提示共线" else "  (<30 可接受)")
  cr <- cor(z[, ..VARS]); diag(cr) <- 0
  w <- which(abs(cr) == max(abs(cr)), arr.ind = TRUE)[1, ]
  log_msg("最大成对相关 ", round(cr[w[1], w[2]], 3), " (", VARS[w[1]], " ~ ", VARS[w[2]], ")")

  log_msg("\n--- 2. 线性性 (RESET) ---")
  rs <- resettest(m, power = 2:3, type = "fitted")
  log_msg(sprintf("F = %.2f, df = (%d, %d), p = %.4g  %s", rs$statistic,
                  rs$parameter[1], rs$parameter[2], rs$p.value,
                  ifelse(rs$p.value < .05, "<<< 拒绝线性, 存在遗漏的非线性", "未拒绝")))

  log_msg("\n--- 3. 同方差 ---")
  bp <- bptest(m)                                        # 辅助回归用全部自变量
  wh <- bptest(m, ~ fitted(m) + I(fitted(m)^2))          # White 型(拟合值二次)
  log_msg(sprintf("Breusch-Pagan  BP = %.2f, df = %d, p = %.4g  %s",
                  bp$statistic, bp$parameter, bp$p.value,
                  ifelse(bp$p.value < .05, "<<< 异方差", "未拒绝同方差")))
  log_msg(sprintf("White(拟合值型) BP = %.2f, df = %d, p = %.4g  %s",
                  wh$statistic, wh$parameter, wh$p.value,
                  ifelse(wh$p.value < .05, "<<< 异方差", "未拒绝同方差")))

  log_msg("\n--- 4. 残差正态性 ---")
  r <- residuals(m); sw <- shapiro.test(r)
  sk <- mean((r - mean(r))^3) / sd(r)^3; ku <- mean((r - mean(r))^4) / sd(r)^4 - 3
  log_msg(sprintf("Shapiro-Wilk W = %.4f, p = %.4g  %s", sw$statistic, sw$p.value,
                  ifelse(sw$p.value < .05, "<<< 拒绝正态", "未拒绝")))
  log_msg(sprintf("偏度 %.3f | 超额峰度 %.3f", sk, ku))

  log_msg("\n--- 5. 残差空间自相关 (Moran's I, k=8 近邻, 999 次置换) ---")
  mi <- morans_I(r, d$longitude, d$latitude)
  log_msg(sprintf("I = %.4f (期望 %.4f), p = %.4f  %s", mi$I, mi$E, mi$p,
                  ifelse(mi$p < .05, "<<< 残差存在空间自相关, 违反独立性", "未检出")))
  yi <- morans_I(z[[yv]], d$longitude, d$latitude)
  log_msg(sprintf("(对照: 因变量本身 I = %.4f, p = %.4f)", yi$I, yi$p))

  log_msg("\n--- 6. 强影响点 (Cook 距离) 与剔除后重拟合 ---")
  ck <- cooks.distance(m); thr <- 4 / n
  log_msg(sprintf("超过 4/n=%.4f 的点: %d (%.1f%%) | 最大 %.3f | 超过 1 的: %d",
                  thr, sum(ck > thr), 100 * mean(ck > thr), max(ck), sum(ck > 1)))
  m2 <- lm(f, z[ck <= thr])
  b1 <- coef(m)[VARS]; b2 <- coef(m2)[VARS]
  p1 <- coef(summary(m))[VARS, 4]; p2 <- coef(summary(m2))[VARS, 4]
  sd1 <- coef(summary(m))[VARS, 2]
  flip <- VARS[(p1 < .05) != (p2 < .05)]
  log_msg(sprintf("剔除后 n = %d | 系数最大变动 %.4f (=%.2f 个原标准误, %s)",
                  nrow(m2$model), max(abs(b2 - b1)),
                  max(abs(b2 - b1) / sd1), VARS[which.max(abs(b2 - b1))]))
  log_msg("显著性翻转的变量: ", if (length(flip)) paste(flip, collapse = ", ") else "无",
          "  (原显著 ", sum(p1 < .05), " -> 剔除后 ", sum(p2 < .05), ")")


  log_msg("\n--- 7. 函数形式: 非线性藏在哪个变量上 ---")
  # RESET 只说"存在遗漏的非线性", 不说在哪。逐个变量加二次项做 t 检验定位,
  # 再把显著的二次项一起放进去重跑 RESET, 看非线性是否被吸收掉。
  q <- rbindlist(lapply(VARS, function(v) {
    mq <- lm(update(f, sprintf("~ . + I(%s^2)", v)), z)
    cf <- coef(summary(mq))
    data.table(var = v, b2 = cf[nrow(cf), 1], p2 = cf[nrow(cf), 4],
               dR2 = summary(mq)$r.squared - summary(m)$r.squared) }))
  setorder(q, p2)
  print(q[1:6, .(var, b2 = signif(b2, 3), p2 = signif(p2, 3), dR2 = round(dR2, 4))])
  qv <- q[p2 < .05 / length(VARS), var]              # Bonferroni, 避免 18 次检验的假阳性
  log_msg("Bonferroni(p<", signif(.05 / length(VARS), 2), ") 下显著的二次项: ",
          if (length(qv)) paste(qv, collapse = ", ") else "无")
  if (length(qv)) {
    mq <- lm(update(f, paste("~ . +", paste(sprintf("I(%s^2)", qv), collapse = "+"))), z)
    rq <- resettest(mq, power = 2:3, type = "fitted")
    log_msg(sprintf("加入这些二次项后: R2 %.4f -> %.4f | RESET p %.3g -> %.3g  %s",
                    summary(m)$r.squared, summary(mq)$r.squared, rs$p.value, rq$p.value,
                    ifelse(rq$p.value < .05, "<<< 仍拒绝, 非线性未被吸收完",
                           "非线性已被吸收")))
    bq <- coef(mq)[VARS]; pq <- coef(summary(mq))[VARS, 4]
    fl <- VARS[(p1 < .05) != (pq < .05)]
    log_msg("线性项显著性翻转的变量: ",
            if (length(fl)) paste(fl, collapse = ", ") else "无",
            "  (原 ", sum(p1 < .05), " -> 加二次项后 ", sum(pq < .05), ")")
    log_msg("线性项最大变动 ", round(max(abs(bq - b1)), 4),
            " (=", round(max(abs(bq - b1) / sd1), 2), " 个原标准误)")
    n_qflip <- length(fl); reset_q <- rq$p.value
  } else { n_qflip <- NA_integer_; reset_q <- NA_real_ }

  summ[[yv]] <- data.table(
    resp = yv, n = n, R2 = summary(m)$r.squared, maxVIF = max(vif), cond = cn,
    reset_F = as.numeric(rs$statistic), reset_p = rs$p.value,
    bp = as.numeric(bp$statistic), bp_p = bp$p.value,
    white_p = wh$p.value, shapiro_p = sw$p.value, skew = sk, kurt = ku,
    moran_I = mi$I, moran_p = mi$p, moran_I_y = yi$I,
    cook_gt4n = sum(ck > thr), cook_max = max(ck),
    drop_maxdb_se = max(abs(b2 - b1) / sd1), n_flip = length(flip),
    n_quad = length(qv), reset_p_quad = reset_q, n_quad_flip = n_qflip)
}

S <- rbindlist(summ)
fwrite(S, file.path(OUT, "diag_rho_summary.csv"))
log_msg("\n==================== 汇总 ====================")
print(S[, .(响应 = resp, n, R2 = round(R2, 3), 最大VIF = round(maxVIF, 2),
            RESET_p = signif(reset_p, 3), BP_p = signif(bp_p, 3),
            White_p = signif(white_p, 3), Moran_I = round(moran_I, 4),
            Moran_p = signif(moran_p, 3), Cook超标 = cook_gt4n, 显著翻转 = n_flip,
            二次项数 = n_quad, RESET_p_加二次 = signif(reset_p_quad, 3))])
log_msg("\n已写出 diag_rho_summary.csv")
