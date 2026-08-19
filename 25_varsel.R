#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 25_varsel.R — 转变时间得分回归的变量精简与滞后截断敏感性
#
# 两个问题合并在一个脚本里回答:
#   (1) 18 个自变量对 200 余个观测偏多(调整 R2 明显低于 R2), 精简后结论稳不稳;
#   (2) 观测窗口截到 tp<=8 / 7 / 6, 样本量与结论如何变化。
#
# 三套自变量方案:
#   full   18 个全上
#   lasso  glmnet 交叉验证选变量(lambda.1se), 再对选中者重拟合 OLS
#   dim3   每个维度按单变量相关强度各留 3 个
#
# 重要提醒: lasso 与 dim3 都在同一份数据上选变量再检验, p 值偏乐观,
#   不能当作正式推断。这里的用途是看"哪些变量在不同方案下反复出现",
#   即稳定性筛查, 而非显著性检验。
#
# 用法: Rscript 25_varsel.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(MASS); library(glmnet) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

DIMS <- list(
  landuse = c("imperv","grass","water","ever_needle","deci_needle",
              "ever_broad","deci_broad","mixedleaf"),
  climate = c("tavg","rh","cloud","precip","rsds_mean","rsds_sd","elev"),
  socioec = c("ntl","invest","urban_rate"))
PRED <- unlist(DIMS, use.names = FALSE)
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

build <- function(M) {
  env <- new.env()
  L <- readLines("23_hazard_all.R")
  src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
  src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
  src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
  invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
  pat <- env$pat
  pat[, score := fifelse(start_dir == "inhibit_first",
                         fifelse(event == 1, (M + 1) - event_tp, 1),
                         fifelse(event == 1, as.numeric(event_tp), M + 1))]
  merge(pat[, .(stat_id, start_dir, score)],
        env$dt[, c("stat_id", PRED), with = FALSE], by = "stat_id")
}

res <- list()
for (M in c(8L, 7L, 6L)) {
  d0 <- build(M)
  log_msg("\n\n##################################################")
  log_msg("############  最大滞后 M = ", M, "  ############")
  log_msg("##################################################")
  log_msg("抑制组 ", d0[start_dir == "inhibit_first", .N], " 站 | 促进组 ",
          d0[start_dir == "promote_first", .N], " 站 | 合计 ", nrow(d0))

  for (g in c("inhibit_first", "promote_first")) {
    d <- d0[start_dir == g, c("score", PRED), with = FALSE]
    d <- d[complete.cases(d)]
    for (v in PRED) set(d, j = v, value = as.numeric(scale(winz(as.numeric(d[[v]])))))
    y <- d$score; X <- as.matrix(d[, ..PRED])
    lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")
    log_msg("\n======== M=", M, " | ", lab, " | n = ", nrow(d), " ========")

    # --- full ---
    m_full <- lm(score ~ ., d)
    s_full <- summary(m_full)

    # --- lasso ---
    cvf <- cv.glmnet(X, y, alpha = 1, nfolds = 10)
    cf <- coef(cvf, s = "lambda.1se")
    sel <- rownames(cf)[which(cf[, 1] != 0)]; sel <- setdiff(sel, "(Intercept)")
    if (!length(sel)) {                       # 1se 过严时退回 min
      cf <- coef(cvf, s = "lambda.min")
      sel <- setdiff(rownames(cf)[which(cf[, 1] != 0)], "(Intercept)")
    }
    log_msg("LASSO 选中 ", length(sel), " 个: ",
            if (length(sel)) paste(sel, collapse = ", ") else "(无)")
    m_las <- if (length(sel)) lm(as.formula(paste("score ~", paste(sel, collapse = "+"))), d) else NULL

    # --- dim3: 各维度按单变量相关强度取前 3 ---
    d3 <- unlist(lapply(DIMS, function(vs) {
      r <- sapply(vs, function(v) abs(cor(d[[v]], y)))
      names(sort(r, decreasing = TRUE))[1:min(3, length(vs))] }), use.names = FALSE)
    m_d3 <- lm(as.formula(paste("score ~", paste(d3, collapse = "+"))), d)
    log_msg("dim3 选中: ", paste(d3, collapse = ", "))

    log_msg(sprintf("R2/调整R2  full %.3f/%.3f | lasso %s | dim3 %.3f/%.3f",
      s_full$r.squared, s_full$adj.r.squared,
      if (is.null(m_las)) "-" else sprintf("%.3f/%.3f", summary(m_las)$r.squared,
                                           summary(m_las)$adj.r.squared),
      summary(m_d3)$r.squared, summary(m_d3)$adj.r.squared))

    grab <- function(m, tag) {
      if (is.null(m)) return(NULL)
      co <- coef(summary(m))
      data.table(M = M, grp = lab, method = tag, var = rownames(co),
                 beta = co[, 1], p = co[, 4])[var %in% PRED]
    }
    tt <- rbindlist(list(grab(m_full, "full"), grab(m_las, "lasso"), grab(m_d3, "dim3")))
    res[[paste(M, g)]] <- tt
    log_msg("-- 各方案下显著(p<.05)的变量 --")
    w <- tt[p < .05]
    if (nrow(w)) print(dcast(w, var ~ method, value.var = "beta",
                             fun.aggregate = function(x) round(x[1], 3), fill = NA))
    else log_msg("(无)")
  }
}

all <- rbindlist(res)
fwrite(all, file.path(OUT, "varsel_score_all.csv"))
log_msg("\n\n############ 稳定性汇总: 在多少个(M x 方案)组合里显著 ############")
for (gg in c("抑制组", "促进组")) {
  n_comb <- all[grp == gg, uniqueN(paste(M, method))]
  s <- all[grp == gg & p < .05, .(显著次数 = .N,
                                  符号 = paste(unique(sign(beta)), collapse = "/"),
                                  平均beta = round(mean(beta), 3)), by = var]
  setorder(s, -显著次数)
  log_msg("\n[", gg, "]  共 ", n_comb, " 个组合")
  print(s)
}
