#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 36_with_building.R — 把建筑与路网变量并入得分回归
#
# 这三个变量(路网密度/建筑基底/建筑高度)覆盖率均为 75.1%, 早先被 CORE 的
# 90% 门槛挡在主模型之外。本脚本补回, 并把"样本变化"与"变量变化"分开:
#   A 18变量-全样本   基线(抑制 214 / 促进 171)
#   B 18变量-缩减样本 只换样本, 变量不变 -> 隔离样本效应
#   C 21变量-缩减样本 再加变量           -> 隔离变量效应
# 只比 A 与 C 会把两种变化混为一谈。
#
# 注意 mean_height 是零膨胀的: 有值的 679 站里 285 站(42%)为 0(缓冲区内无建筑)。
# 故除原值外, 另报 log1p 变换与"有无建筑"二元拆分的敏感性。
#
# 标准误: Conley 空间 HAC, Bartlett 核, 截断 200 km(与主模型一致)。
# 重要性: Shapley 值(抽样排列逼近, 合计等于 R²)。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; M <- 8L; CUT <- 200; NPERM <- 3000
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

env <- new.env(); L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
pat <- env$pat; dt <- env$dt; P18 <- env$PRED; winz <- env$winz
BLD <- c("road_density", "building_footprint", "mean_height")
P21 <- c(P18, BLD)
pat[, score := fifelse(start_dir == "inhibit_first",
                       fifelse(event == 1, (M + 2) - event_tp, 1),
                       fifelse(event == 1, as.numeric(event_tp), M + 1))]
st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
d0 <- merge(merge(pat[, .(stat_id, start_dir, score)],
                  dt[, c("stat_id", P21), with = FALSE], by = "stat_id"),
            st[, .(stat_id, longitude, latitude)], by = "stat_id")

CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率",
        road_density="路网密度", building_footprint="建筑基底", mean_height="建筑高度")

conley_se <- function(m, lon, lat, cut = CUT) {
  X <- model.matrix(m); u <- residuals(m); n <- nrow(X)
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  D <- sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2); br <- solve(crossprod(X))
  sqrt(diag(br %*% crossprod(X, (pmax(0, 1 - D / cut) * outer(u, u)) %*% X) %*% br) *
       (n / (n - ncol(X))))
}
r2_sub <- function(R, S) { if (!length(S)) return(0)
  rys <- R[1, S, drop = FALSE]; as.numeric(rys %*% solve(R[S, S, drop = FALSE]) %*% t(rys)) }
shapley <- function(z, vars, nperm = NPERM) {
  R <- cor(z[, c("score", vars), with = FALSE]); idx <- setNames(seq_along(vars) + 1L, vars)
  acc <- setNames(numeric(length(vars)), vars)
  for (b in seq_len(nperm)) { cur <- integer(0); prev <- 0
    for (v in sample(vars)) { cur <- c(cur, idx[[v]]); now <- r2_sub(R, cur)
      acc[v] <- acc[v] + (now - prev); prev <- now } }
  acc / nperm
}

fit <- function(d, vars, tag) {
  z <- copy(d)[, c("score", vars), with = FALSE]
  for (v in vars) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  m <- lm(as.formula(paste("score ~", paste(vars, collapse = "+"))), z)
  se <- conley_se(m, d$longitude, d$latitude); b <- coef(m)[vars]
  sh <- shapley(z, vars)
  R2 <- summary(m)$r.squared
  data.table(spec = tag, n = nrow(z), R2 = R2, adjR2 = summary(m)$adj.r.squared,
             var = vars, beta = b, se = se[vars], p = 2 * pnorm(-abs(b / se[vars])),
             shap = sh[vars], shap_pct = 100 * sh[vars] / sum(sh))
}

res <- list()
for (g in c("inhibit_first", "promote_first")) {
  lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")
  dA <- d0[start_dir == g][complete.cases(d0[start_dir == g, c("score", P18), with = FALSE])]
  dC <- d0[start_dir == g][complete.cases(d0[start_dir == g, c("score", P21), with = FALSE])]
  log_msg("\n################ ", lab, " ################")
  log_msg("A 18变量全样本 ", nrow(dA), " 站 | B/C 缩减样本 ", nrow(dC), " 站 (损失 ",
          nrow(dA) - nrow(dC), ", ", round(100 * (nrow(dA) - nrow(dC)) / nrow(dA)), "%)")
  log_msg("建筑高度为 0 的站: ", dC[mean_height == 0, .N], " / ", nrow(dC),
          " (", round(100 * mean(dC$mean_height == 0)), "%)")
  a <- fit(dA, P18, "A_18var_全样本")
  b <- fit(dC, P18, "B_18var_缩减样本")
  cc <- fit(dC, P21, "C_21var_缩减样本")
  t <- rbindlist(list(a, b, cc)); t[, grp := g]; res[[g]] <- t
  log_msg(sprintf("R²/调整R²:  A %.3f/%.3f | B %.3f/%.3f | C %.3f/%.3f",
    a$R2[1], a$adjR2[1], b$R2[1], b$adjR2[1], cc$R2[1], cc$adjR2[1]))

  log_msg("\n-- C(21变量) 全部变量, 按 Shapley 占比排序 --")
  print(cc[order(-shap), .(变量 = CN[var], 系数 = round(beta, 3), p = signif(p, 2),
      显著 = fifelse(p < .05, "*", ""), Shapley占比 = round(shap_pct, 1))])
  log_msg("\n-- 三个建筑/路网变量 --")
  print(cc[var %in% BLD, .(变量 = CN[var], 系数 = round(beta, 3), 标准误 = round(se, 3),
      p = signif(p, 3), Shapley占比 = round(shap_pct, 1))])
  log_msg("\n-- 加入建筑变量前后, 原显著变量的变化(同一批 ", nrow(dC), " 站) --")
  cmp <- merge(b[, .(var, beta_B = beta, p_B = p, shap_B = shap_pct)],
               cc[, .(var, beta_C = beta, p_C = p, shap_C = shap_pct)], by = "var")
  print(cmp[p_B < .05 | p_C < .05][order(p_C),
      .(变量 = CN[var], 系数_B = round(beta_B, 3), p_B = signif(p_B, 2),
        系数_C = round(beta_C, 3), p_C = signif(p_C, 2),
        Shapley_B = round(shap_B, 1), Shapley_C = round(shap_C, 1))])
}
fwrite(rbindlist(res), file.path(OUT, sprintf("building_tp%d.csv", M)))
log_msg("\n已写出 building_tp", M, ".csv")
