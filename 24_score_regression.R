#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 24_score_regression.R — 用"转变时间得分"作因变量的回归
#
# 与 23 号脚本的生存框架不同: 这里把每站压成一个有序得分, 回到普通回归。
#
# 编码(设最大滞后为 M):
#   抑制组(始终抑制 + 抑制后促进)
#       始终抑制            -> 1
#       在 tp = k 处翻转    -> (M + 1) - k     k=M 得 1, k=1 得 M
#       数值越大 = 越早脱离抑制
#   促进组(始终促进 + 先促进后抑制)
#       在 tp = k 处翻转    -> k               k=1 得 1, k=M 得 M
#       始终促进            -> M + 1
#       数值越大 = 促进状态维持越久
#
#   两组的共同含义: 得分越高, 系统停留在"促进"状态的时间越长。
#
# 注意这个编码把"始终抑制"与"在最末步翻转"合并为同一分值(均为 1),
# 即右删失信息被抹去——这是它相对生存模型的代价, 换来的是单一有序因变量。
#
# 模型: OLS(标准化系数) + 有序 logit(得分本质是有序而非等距) 互为稳健性检验。
#
# 用法: Rscript 24_score_regression.R [--maxtp 7]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(MASS) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")
aa <- commandArgs(trailingOnly = TRUE)
M <- if (any(aa == "--maxtp")) as.integer(aa[which(aa == "--maxtp") + 1]) else 7L
log_msg("最大滞后 M = ", M)

# ---- 复用 23 号脚本的分类与协变量(执行到 pat/pp 建好) --------------------
env <- new.env()
L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
pat <- env$pat; dt <- env$dt; PRED <- env$PRED; winz <- env$winz

# ---- 构造得分 -------------------------------------------------------------
pat[, score := fifelse(
  start_dir == "inhibit_first",
  fifelse(event == 1, (M + 1) - event_tp, 1),        # 抑制组
  fifelse(event == 1, as.numeric(event_tp), M + 1))] # 促进组

log_msg("\n-- 得分分布 --")
for (g in c("inhibit_first", "promote_first")) {
  s <- pat[start_dir == g]
  log_msg("\n[", g, "]  n = ", nrow(s))
  print(s[, .(站数 = .N), by = .(得分 = score, 类型 = stype,
                                 翻转tp = event_tp)][order(得分)])
}

d0 <- merge(pat[, .(stat_id, start_dir, score, stype, event_tp)],
            dt[, c("stat_id", PRED, "koppen_group"), with = FALSE], by = "stat_id")

# ---- 回归 -----------------------------------------------------------------
fit_score <- function(g, label) {
  d <- d0[start_dir == g, c("score", PRED), with = FALSE][complete.cases(
       d0[start_dir == g, c("score", PRED), with = FALSE])]
  log_msg("\n================ ", label, " | n = ", nrow(d), " ================")
  log_msg("得分 均值 ", round(mean(d$score), 2), " | 中位 ", median(d$score),
          " | 标准差 ", round(sd(d$score), 2))
  z <- copy(d); for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  f <- as.formula(paste("score ~", paste(PRED, collapse = " + ")))

  m <- lm(f, z); co <- coef(summary(m))
  tab <- data.table(var = rownames(co), beta = co[, 1], se = co[, 2], p = co[, 4])[var %in% PRED]

  # 有序 logit: 得分是有序类别, 不必等距
  z2 <- copy(z); z2[, score := factor(score, ordered = TRUE)]
  om <- tryCatch(polr(f, data = z2, Hess = TRUE), error = function(e) NULL)
  if (!is.null(om)) {
    oc <- coef(summary(om))
    ot <- data.table(var = rownames(oc), or_coef = oc[, 1],
                     p_ord = 2 * pnorm(-abs(oc[, 3])))[var %in% PRED]
    tab <- merge(tab, ot, by = "var", all.x = TRUE)
  }
  tab <- tab[order(p)]
  log_msg("OLS R2 = ", round(summary(m)$r.squared, 4),
          " | 调整 R2 = ", round(summary(m)$adj.r.squared, 4))
  print(tab[, .(var, beta = round(beta, 3), p = signif(p, 3),
                sig = fifelse(p < .05, "*", ""),
                ord_coef = round(or_coef, 3), p_ord = signif(p_ord, 3),
                ord_sig = fifelse(p_ord < .05, "*", ""))])
  fwrite(tab, file.path(OUT, sprintf("score_tp%d_%s.csv", M, g)))
  invisible(tab)
}
fit_score("inhibit_first", "抑制组: 得分越高 = 越早脱离抑制")
fit_score("promote_first", "促进组: 得分越高 = 促进维持越久")

# ---- 分气候区 -------------------------------------------------------------
PRED_S <- c("imperv", "grass", "elev", "rsds_sd", "ntl", "urban_rate")
log_msg("\n\n############ 分气候区(精简 6 变量) ############")
for (kg in sort(unique(na.omit(d0$koppen_group)))) {
  for (g in c("inhibit_first", "promote_first")) {
    d <- d0[koppen_group == kg & start_dir == g, c("score", PRED_S), with = FALSE]
    d <- d[complete.cases(d)]
    lab <- sprintf("%s 区 | %s", kg, ifelse(g == "inhibit_first", "抑制组", "促进组"))
    if (nrow(d) < 10 * length(PRED_S)) {
      log_msg(sprintf("\n[跳过] %s: n = %d (低于 %d)", lab, nrow(d), 10 * length(PRED_S))); next }
    if (uniqueN(d$score) < 3) { log_msg("\n[跳过] ", lab, ": 得分取值不足"); next }
    for (v in PRED_S) set(d, j = v, value = as.numeric(scale(winz(as.numeric(d[[v]])))))
    m <- lm(as.formula(paste("score ~", paste(PRED_S, collapse = "+"))), d)
    co <- coef(summary(m))
    t <- data.table(var = rownames(co), beta = co[, 1], p = co[, 4])[var %in% PRED_S]
    log_msg(sprintf("\n===== %s | n = %d | R2 = %.3f =====", lab, nrow(d), summary(m)$r.squared))
    print(t[order(p), .(var, beta = round(beta, 3), p = signif(p, 3),
                        sig = fifelse(p < .05, "*", ""))])
  }
}
