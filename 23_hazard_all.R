#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 23_hazard_all.R — 按"抑制组 / 促进组"分别做转变时间的离散时间风险模型
#
# 沿用 15_hcsif_new_strategy.R 的分类口径(逐字复现 find_transition):
#   起始方向由 tp=0 的因果系数符号决定;
#   翻转须"持续"——某步变号后, 其后各步中同号占比 >= 50% 才算真翻转;
#   抑制组 = always_inhibit + inhibit_promote
#   促进组 = always_promote + promote_inhibit
#   生存框架: time = 翻转所在 tp, 未翻转者右删失于 tp=8; event = 是否翻转。
#
# 与原脚本的两点不同:
#   1) 协变量换成本轮装配的三维度全集(以耕地为参照类, maxVIF 5.7);
#   2) 补按站点聚类的稳健标准误——person-period 展开后同一站有多行,
#      普通 GLM 把它们当独立观测, 标准误偏小、p 值偏乐观。
#
# 输出: data_proc/output_hcsif_buf1000/hazard_all_<组>.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")

# ---- 复用 22 号脚本的协变量装配(执行到 dt 建好为止) ----------------------
env <- new.env()
lines <- readLines("22_varpart_3dim.R")
cut <- grep("^CORE <- lapply", lines)[1]   # 截到 CORE 定义完为止
invisible(capture.output(
  eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)
log_msg("协变量 ", length(PRED), " 个: ", paste(PRED, collapse = ", "))

# ---- 分类: 逐字复现原脚本的持续性翻转判据 --------------------------------
cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
N_TP <- length(unique(cm$tp))

find_transition <- function(co, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(co) < 2) return(NA_integer_)
  for (i in 2:length(co)) {
    if (!is.na(co[i]) && cond(co[i])) {
      tail_ok <- mean(sapply(co[i:length(co)], cond), na.rm = TRUE)
      if (!is.na(tail_ok) && tail_ok >= 0.5) return(as.integer(i - 1L))
    }
  }
  NA_integer_
}

pat <- cm[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)]
pat <- pat[sapply(co, length) == N_TP]
pat[, `:=`(
  tp0_pos = sapply(co, function(x) !is.na(x[1]) && x[1] > 0),
  HTW     = sapply(co, find_transition, to_neg = TRUE),
  ITW     = sapply(co, find_transition, to_neg = FALSE))]
pat[, stype := fcase(
  tp0_pos  & is.na(HTW),  "always_promote",
  !tp0_pos & is.na(ITW),  "always_inhibit",
  tp0_pos  & !is.na(HTW), "promote_inhibit",
  !tp0_pos & !is.na(ITW), "inhibit_promote",
  default = "always_inhibit")]
pat[, `:=`(
  start_dir = fifelse(tp0_pos, "promote_first", "inhibit_first"),
  flipped   = stype %in% c("inhibit_promote", "promote_inhibit"))]
pat[, event_tp := fcase(stype == "inhibit_promote", ITW,
                        stype == "promote_inhibit", HTW,
                        default = NA_integer_)]
pat[, `:=`(time = fifelse(flipped, as.numeric(event_tp), max(cm$tp)),
           event = as.integer(flipped))]
pat[, co := NULL]

log_msg("\n-- 站点类型分布 --")
print(pat[, .(n = .N, 占比 = round(.N / nrow(pat), 3)), by = stype][order(-n)])
log_msg("\n-- 起始方向 x 是否翻转 --")
print(dcast(pat[, .N, by = .(start_dir, flipped)], start_dir ~ flipped, value.var = "N"))
log_msg("\n-- 翻转时间(tp, 每步8天)分布, 仅已翻转者 --")
print(pat[event == 1, .(n = .N, 中位 = median(time), 均值 = round(mean(time), 2)), by = start_dir])

# ---- person-period 展开 ---------------------------------------------------
pp <- pat[, .(tp = seq_len(max(1, time))), by = .(stat_id, start_dir, time, event)]
pp[, ev := as.integer(tp == time & event == 1)]
log_msg("\nperson-period: ", nrow(pp), " 行 = 站点 x 滞后步")

pp <- merge(pp, dt[, c("stat_id", PRED), with = FALSE], by = "stat_id")
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

# ---- 离散时间风险模型 + 按站聚类稳健标准误 -------------------------------
cluster_se <- function(fit, cl) {
  X <- model.matrix(fit)
  u <- residuals(fit, type = "response")            # y - p, 二项 GLM 的得分权重
  br <- vcov(fit)
  sc <- X * u
  meat <- crossprod(rowsum(sc, cl))
  G <- length(unique(cl)); n <- nrow(X); k <- ncol(X)
  adj <- (G / (G - 1)) * ((n - 1) / (n - k))        # 小样本校正
  sqrt(diag(br %*% meat %*% br) * adj)
}

fit_group <- function(g, label) {
  d <- pp[start_dir == g, c("ev", "tp", "stat_id", PRED), with = FALSE]
  d <- d[complete.cases(d)]
  for (v in PRED) set(d, j = v, value = as.numeric(scale(winz(as.numeric(d[[v]])))))
  n_ev <- sum(d$ev)
  log_msg("\n================ ", label, " ================")
  log_msg("站点 ", uniqueN(d$stat_id), " | person-period 行 ", nrow(d), " | 翻转事件 ", n_ev)
  if (n_ev < 10) { log_msg("事件过少, 跳过"); return(invisible(NULL)) }

  # tp 作为基线风险(因子), 与原脚本一致
  f <- as.formula(paste("ev ~ factor(tp) +", paste(PRED, collapse = " + ")))
  fit <- glm(f, data = d, family = binomial)
  co <- summary(fit)$coefficients
  cse <- cluster_se(fit, d$stat_id)
  tab <- data.table(var = rownames(co), coef = co[, 1],
                    se_naive = co[, 2], se_clust = cse[rownames(co)])
  tab <- tab[var %in% PRED]
  tab[, `:=`(z = coef / se_clust)]
  tab[, p_clust := 2 * pnorm(-abs(z))]
  tab[, p_naive := 2 * pnorm(-abs(coef / se_naive))]
  tab <- tab[order(-abs(z))]

  log_msg("系数>0 → 每步翻转风险更高 → 翻转更早")
  print(tab[, .(var, coef = round(coef, 3), se_clust = round(se_clust, 3),
                p_clust = signif(p_clust, 3), p_naive = signif(p_naive, 3),
                sig = fifelse(p_clust < .05, "*", ""))])
  log_msg("聚类后仍显著: ", tab[p_clust < .05, .N], " 个 | 若不聚类会有 ",
          tab[p_naive < .05, .N], " 个")
  fwrite(tab, file.path(OUT, sprintf("hazard_all_%s.csv", g)))
  invisible(tab)
}

fit_group("inhibit_first", "抑制组 (全抑制 + 先抑制后促进) → 翻转为促进")
fit_group("promote_first", "促进组 (全促进 + 先促进后抑制) → 翻转为抑制")
log_msg("\n完成")
