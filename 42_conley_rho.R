#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 42_conley_rho.R — 因果强度回归的稳健标准误: 普通 / HC3 / Conley 空间 HAC
#
# 41 号诊断查出三个响应变量全部拒绝同方差, 且残差 Moran's I 达 0.16-0.29
# (远高于得分回归促进组的 0.102), 所以这一支比得分那一支更需要空间修正。
#
# 截断距离不直接照搬得分回归的 200 km: 那个数是在 n=214/171 的子样本上定的,
# 本支样本 n=768、空间覆盖不同, 故按 31_conley_cutoff.R 的三条依据重新验证:
#   (1) 残差相关图: 空间相关在哪个距离带落回零
#   (2) 标准误对截断距离的敏感曲线: 应落在平台段
#   (3) 有效独立样本量 n / 平均邻居数: 太小则 meat 矩阵不稳
#
# 点估计三种方法完全相同, 变的只有标准误与 p 值。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(sandwich); library(ggplot2)
  library(patchwork); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"
set.seed(42)
CUTS <- c(50, 75, 100, 150, 200, 300, 400, 600, 800)
CUT  <- 200
log_msg <- function(...) cat(..., "\n", sep = "")

env <- new.env(); L <- readLines("22_varpart_3dim.R")
src <- L[1:grep("^winz <- function", L)[1]]
invisible(capture.output(suppressWarnings(
  eval(parse(text = paste(src, collapse = "\n")), envir = env))))
dt <- env$dt; CORE <- env$CORE; winz <- env$winz; st <- env$st
VARS <- unlist(CORE, use.names = FALSE)
dt <- merge(dt, st[, .(stat_id, longitude, latitude)], by = "stat_id", all.x = TRUE)

RESP <- c(rho_tp0 = "rho_tp0", rho_max = "rho_max", coef_tp0 = "coef_tp0")
RCN  <- c(rho_tp0 = "同期因果强度 rho(tp=0)", rho_max = "最强因果强度 rho_max",
          coef_tp0 = "因果方向 S-map 偏导数(tp=0)")
CN <- c(imperv="不透水面", grass="草地", water="水体",
        ever_needle="常绿针叶林", deci_needle="落叶针叶林",
        ever_broad="常绿阔叶林", deci_broad="落叶阔叶林", mixedleaf="混交林",
        tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
DIMS <- list(`土地利用与空间结构` = CORE$landuse, `气候与自然环境` = CORE$climate,
             `社会经济与人类活动` = CORE$socioec)
dim_of <- setNames(rep(names(DIMS), lengths(DIMS)), unlist(DIMS))

distmat <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2)
}
# Conley 空间 HAC: Bartlett 核, 截断距离外权重为 0
conley_vcov <- function(m, D, cut) {
  X <- model.matrix(m); u <- residuals(m); n <- nrow(X)
  K <- pmax(0, 1 - D / cut)
  br <- solve(crossprod(X))
  br %*% crossprod(X, (K * outer(u, u)) %*% X) %*% br * (n / (n - ncol(X)))
}
band_I <- function(r, D, lo, hi, nperm = 499) {
  W <- (D > lo & D <= hi) * 1; diag(W) <- 0
  if (sum(W) < 20) return(c(NA, NA, sum(W) / 2))
  z <- r - mean(r); n <- length(r)
  s <- function(v) (n / sum(W)) * as.numeric(crossprod(v, W %*% v)) / sum(v^2)
  I <- s(z); pm <- replicate(nperm, s(sample(z)))
  c(I, (sum(abs(pm) >= abs(I)) + 1) / (nperm + 1), sum(W) / 2)
}
BANDS <- c(0, 100, 200, 300, 400, 600, 800, 1200, 2000)

fit_one <- function(yv) {
  cols <- c(yv, VARS)
  d <- dt[complete.cases(dt[, ..cols])]
  z <- d[, ..cols]
  for (v in VARS) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  list(m = lm(as.formula(paste(yv, "~", paste(VARS, collapse = "+"))), z),
       d = d, z = z)
}

cor_all <- list(); se_all <- list(); res <- list()
for (yv in names(RESP)) {
  fo <- fit_one(yv); m <- fo$m; d <- fo$d
  D <- distmat(d$longitude, d$latitude); r <- residuals(m)
  log_msg("\n#################### ", RCN[yv], " | n = ", nrow(d), " ####################")
  log_msg("站点两两距离: 中位 ", round(median(D[upper.tri(D)])), " km | 最大 ", round(max(D)), " km")

  log_msg("\n--- (1) 残差相关图: 各距离带的 Moran's I ---")
  cg <- rbindlist(lapply(seq_len(length(BANDS) - 1), function(i) {
    v <- band_I(r, D, BANDS[i], BANDS[i + 1])
    data.table(resp = yv, lo = BANDS[i], hi = BANDS[i + 1],
               I = v[1], p = v[2], npair = v[3]) }))
  print(cg[, .(距离带 = paste0(lo, "-", hi, "km"), I = round(I, 4), p = signif(p, 3),
               站对数 = npair, 显著 = fifelse(!is.na(p) & p < .05, "*", ""))])
  cor_all[[yv]] <- cg

  log_msg("\n--- (2)(3) 标准误敏感曲线与有效样本量 ---")
  se0 <- sqrt(diag(vcov(m)))[VARS]
  st_tab <- rbindlist(lapply(CUTS, function(k) {
    s <- sqrt(diag(conley_vcov(m, D, k)))[VARS]
    nb <- mean(rowSums(D <= k)) - 1
    data.table(resp = yv, cut = k, var = VARS, se = s, nb = nb,
               eff_n = nrow(d) / max(1, nb)) }))
  st_tab[, ratio := se / se0[var]]
  print(st_tab[, .(平均邻居 = round(nb[1], 1), 有效样本 = round(eff_n[1], 1),
                   SE比中位 = round(median(ratio), 3),
                   SE比范围 = sprintf("%.2f-%.2f", min(ratio), max(ratio))), by = cut])
  se_all[[yv]] <- st_tab

  # --- 三种标准误对照 ---
  b <- coef(m)[VARS]
  out <- data.table(resp = yv, var = VARS, beta = b, se_ols = se0,
                    se_hc3 = sqrt(diag(vcovHC(m, type = "HC3")))[VARS])
  for (k in c(100, 200, 400))
    out[, paste0("se_conley", k) := sqrt(diag(conley_vcov(m, D, k)))[VARS]]
  pv <- function(s) 2 * pnorm(-abs(b / s))
  out[, `:=`(p_ols = pv(se_ols), p_hc3 = pv(se_hc3))]
  for (k in c(100, 200, 400)) out[, paste0("p_conley", k) := pv(out[[paste0("se_conley", k)]])]
  res[[yv]] <- out
}

R <- rbindlist(res)
fwrite(R, file.path(OUT, "robust_se_rho.csv"))
fwrite(rbindlist(cor_all), file.path(OUT, "conley_correlogram_rho.csv"))
SE <- rbindlist(se_all); fwrite(SE, file.path(OUT, "conley_se_by_cutoff_rho.csv"))

# ===== 结论是否变化 =====
log_msg("\n==================== 普通 OLS vs Conley ", CUT, "km: 结论是否变化 ====================")
R[, `:=`(sig_ols = p_ols < .05, sig_con = get(paste0("p_conley", CUT)) < .05)]
print(R[, .(显著_OLS = sum(sig_ols), 显著_HC3 = sum(p_hc3 < .05),
            显著_Conley100 = sum(p_conley100 < .05),
            显著_Conley200 = sum(p_conley200 < .05),
            显著_Conley400 = sum(p_conley400 < .05),
            SE倍数_HC3 = round(median(se_hc3 / se_ols), 3),
            SE倍数_Conley200 = round(median(se_conley200 / se_ols), 3)), by = resp])
log_msg("\n-- 由显著转为不显著的变量(OLS 显著, Conley ", CUT, "km 不显著) --")
lost <- R[sig_ols & !sig_con, .(响应 = resp, 变量 = CN[var], 维度 = dim_of[var],
                                系数 = round(beta, 4), p_OLS = signif(p_ols, 3),
                                p_HC3 = signif(p_hc3, 3),
                                p_Conley = signif(get(paste0("p_conley", CUT)), 3))]
if (nrow(lost)) print(lost) else log_msg("无")
log_msg("\n-- 反向: OLS 不显著而 Conley 显著 --")
gain <- R[!sig_ols & sig_con, .(响应 = resp, 变量 = CN[var], 系数 = round(beta, 4),
                                p_OLS = signif(p_ols, 3),
                                p_Conley = signif(get(paste0("p_conley", CUT)), 3))]
if (nrow(gain)) print(gain) else log_msg("无")

log_msg("\n-- Conley ", CUT, " km 下仍显著的变量(按响应) --")
for (yv in names(RESP)) {
  x <- R[resp == yv & sig_con][order(-abs(beta))]
  log_msg("\n[", RCN[yv], "]")
  if (nrow(x)) print(x[, .(变量 = CN[var], 维度 = dim_of[var], 系数 = round(beta, 4),
                           Conley_SE = round(get(paste0("se_conley", CUT)), 4),
                           p = signif(get(paste0("p_conley", CUT)), 3))])
  else log_msg("无")
}

# ===== 图 =====
cg <- rbindlist(cor_all); cg[, mid := (lo + hi) / 2]
cg[, rcn := factor(RCN[resp], levels = RCN)]
PALR <- setNames(c("#1C6B66", "#C2703A", "#4B6BA8"), RCN)
p1 <- ggplot(cg, aes(mid, I, colour = rcn)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_vline(xintercept = CUT, linetype = "22", colour = "grey45") +
  geom_line(linewidth = .8) +
  geom_point(aes(shape = !is.na(p) & p < .05), size = 2.3) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PALR, name = NULL) +
  labs(title = "A  残差相关图：因果强度回归的空间相关有多远",
       subtitle = "实心点 = 该距离带 Moran's I 显著（置换检验 p<0.05）；虚线为采用的 200 km",
       x = "距离带中点（km）", y = "Moran's I") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

MTH <- c(se_ols = "普通 OLS", se_hc3 = "HC3（修异方差）",
         se_conley200 = "Conley 200 km（异方差＋空间相关）")
dd <- melt(R[, .(resp, var, beta, se_ols, se_hc3, se_conley200)],
           id.vars = c("resp", "var", "beta"), variable.name = "method", value.name = "se")
dd[, `:=`(method = factor(MTH[as.character(method)], levels = MTH),
          lo = beta - 1.96 * se, hi = beta + 1.96 * se,
          cn = CN[var], rcn = factor(RCN[resp], levels = RCN))]
dd[, sig := (lo > 0 | hi < 0)]
ord <- dd[resp == "rho_tp0" & method == MTH[1]][order(beta)]
dd[, cn := factor(cn, levels = ord$cn)]
PAL <- c(`普通 OLS` = "#9BA8A5", `HC3（修异方差）` = "#C2703A",
         `Conley 200 km（异方差＋空间相关）` = "#1C6B66")
p2 <- ggplot(dd, aes(beta, cn, colour = method)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = .65,
                 position = position_dodge(width = .68)) +
  geom_point(size = 1.4, position = position_dodge(width = .68)) +
  geom_point(data = dd[sig == TRUE], size = 2.4, shape = 21, fill = "white",
             stroke = .75, position = position_dodge(width = .68)) +
  facet_wrap(~ rcn, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = PAL, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) +
  labs(title = "B  标准误方法对推断的影响：点估计不变，区间在变",
       subtitle = "空心圈标记该方法下 95% 区间不跨零；三个响应变量量纲不同，横轴各自缩放",
       x = "回归系数（自变量已标准化，95% 置信区间）", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = .35),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10.5, hjust = 0),
        legend.position = "top", legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2))

ggsave(file.path(OUT, "fig_robust_se_rho.png"),
       p1 / p2 + plot_layout(heights = c(1, 1.9)), width = 11.5, height = 10.5, dpi = 300)
log_msg("\n已输出 robust_se_rho.csv / conley_correlogram_rho.csv / ",
        "conley_se_by_cutoff_rho.csv / fig_robust_se_rho.png")


# ===== 与非线性一起看: 二次项 + Conley 标准误 =====
# 41 号的 RESET 在三个响应上全部拒绝, 说明除了标准误, 函数形式本身也有问题。
# 这里把两处修正叠加: 按 Bonferroni 选入显著的二次项后再算 Conley 标准误。
# 自变量已标准化(均值 0), 故二次模型里线性项恰好等于平均边际效应,
# 与线性模型的系数可直接比较——这一点是下面"是否翻转"能成立的前提。
log_msg("\n==================== 叠加非线性修正后 结论是否还站得住 ====================")
rob <- list()
for (yv in names(RESP)) {
  fo <- fit_one(yv); m <- fo$m; d <- fo$d; z <- fo$z
  D <- distmat(d$longitude, d$latitude)
  f <- formula(m)
  q <- sapply(VARS, function(v) {
    cf <- coef(summary(lm(update(f, sprintf("~ . + I(%s^2)", v)), z))); cf[nrow(cf), 4] })
  qv <- names(q)[q < .05 / length(VARS)]
  mq <- lm(update(f, paste("~ . +", paste(sprintf("I(%s^2)", qv), collapse = "+"))), z)
  pq <- 2 * pnorm(-abs(coef(mq)[VARS] / sqrt(diag(conley_vcov(mq, D, CUT)))[VARS]))
  base <- R[resp == yv]
  rob[[yv]] <- data.table(resp = yv, var = VARS, beta_lin = base$beta[match(VARS, base$var)],
                          p_con_lin = base[[paste0("p_conley", CUT)]][match(VARS, base$var)],
                          beta_quad = coef(mq)[VARS], p_con_quad = pq)
  # 二次项与原变量高度相关(右偏变量尤甚), 必须报共线, 否则"翻转"可能只是方差膨胀
  Xq <- model.matrix(mq)[, -1]
  vq <- sapply(seq_len(ncol(Xq)), function(j)
    1 / (1 - summary(lm(Xq[, j] ~ Xq[, -j]))$r.squared))
  names(vq) <- colnames(Xq)
  log_msg("\n[", RCN[yv], "] 选入二次项: ", paste(qv, collapse = ", "),
          " | R2 ", round(summary(m)$r.squared, 3), " -> ", round(summary(mq)$r.squared, 3))
  log_msg("  二次模型最大 VIF ", round(max(vq), 1), " (", names(which.max(vq)), ")",
          if (max(vq) > 10) "  <<< 超过线性模型的 5.7, 该式只作稳健性参考, 不宜当最终设定"
          else "  [<10]")
}
RB <- rbindlist(rob)
RB[, `:=`(s_lin = p_con_lin < .05, s_quad = p_con_quad < .05)]
fwrite(RB, file.path(OUT, "robust_se_rho_quad.csv"))
log_msg("\n-- 线性 vs 二次(均用 Conley ", CUT, " km)显著变量数 --")
print(RB[, .(线性模型 = sum(s_lin), 二次模型 = sum(s_quad),
             两者皆显著 = sum(s_lin & s_quad)), by = resp])
log_msg("\n-- 加入二次项后显著性翻转的变量 --")
fl <- RB[s_lin != s_quad, .(响应 = resp, 变量 = CN[var],
                            线性系数 = round(beta_lin, 4), p_线性 = signif(p_con_lin, 3),
                            二次模型线性项 = round(beta_quad, 4), p_二次 = signif(p_con_quad, 3),
                            方向 = fifelse(s_lin, "显著->不显著", "不显著->显著"))]
if (nrow(fl)) print(fl) else log_msg("无")
log_msg("\n-- 同时经受住空间 HAC 与非线性两重修正的变量 --")
print(RB[s_lin & s_quad, .(响应 = resp, 变量 = CN[var], 维度 = dim_of[var],
                           系数 = round(beta_lin, 4))][order(响应)])
log_msg("\n已写出 robust_se_rho_quad.csv")

# 与 22 号原有 OLS 结果核对: 点估计应逐位相同
chk <- merge(R[resp == "rho_tp0", .(var, beta)],
             fread(file.path(OUT, "varpart3_CORE_rho_tp0_coef.csv"))[, .(var, beta0 = beta)],
             by = "var")
log_msg("与 22_varpart_3dim.R 的 OLS 点估计最大绝对差: ",
        signif(max(abs(chk$beta - chk$beta0)), 3))
