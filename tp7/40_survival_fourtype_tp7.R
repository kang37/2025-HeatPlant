#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 40_survival_fourtype.R — "四类站"口径的生存分析(全流程统一 n = 385)
#
# 样本 = 严格判据下的四类站, 且 18 个协变量齐全:
#   1 始终抑制      9 个滞后全为负                  -> 右删失于 tp8
#   2 先抑制后促进  起始为负, 变号一次后其后全为正  -> 事件
#   3 始终促进      9 个滞后全为正                  -> 右删失于 tp8
#   4 先促进后抑制  起始为正, 变号一次后其后全为负  -> 事件
#   多次变号的 437 站不入样。
#
# 四类站共 461, 协变量齐全者 385(缺: invest 49, urban_rate 42, cloud 8)。
# KM/log-rank 与回归全部用同一批 385 站, 避免不同结论建立在不同样本上。
#
# 注意: 该口径下不存在竞争事件(每站要么删失、要么恰好一次持续翻转),
#       故 Fine-Gray 子分布风险退化为普通 Cox, 竞争风险分析在此不适用。
#       同一样本上可做的"两模型比对"是抑制组与促进组的系数对照(图 c)。
#
# 输出: surv4_*.csv, fig_survival_fourtype.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
  library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; CUT <- 200; MAXTP <- 7L
log_msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) log_msg("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))

CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率",
        elev="海拔", ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
GL <- c(inhibit_first = "抑制组", promote_first = "促进组")
CL <- setNames(c("#2166AC", "#B2182B"), GL)

dist_km <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2) }
conley_cox <- function(fit, lon, lat, cut) {
  Dfb <- as.matrix(residuals(fit, type = "dfbeta"))
  DM <- dist_km(lon, lat); K <- DM; K[] <- pmax(0, 1 - DM / cut)
  list(se = sqrt(pmax(diag(t(Dfb) %*% K %*% Dfb), 0)),
       eff_n = nrow(DM) / max(1, mean(rowSums(DM <= cut)))) }
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

# ===========================================================================
hdr("0. 样本")
env <- new.env(); L <- readLines("22_varpart_3dim.R"); k <- grep("^CORE <- lapply", L)[1]
invisible(capture.output(eval(parse(text = paste(L[1:k], collapse = "\n")), envir = env)))
dtc <- env$dt; PRED <- unlist(env$CORE, use.names = FALSE)

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
cm <- cm[tp <= MAXTP]; setorder(cm, meteo_stat, tp)
pat <- cm[, .(co = list(mean_coef)), by = .(stat_id = meteo_stat)][sapply(co, length) == MAXTP + 1L]

cls <- function(x) {
  s <- sign(x)
  if (anyNA(s) || any(s == 0)) return(list(ty = "undef", t = NA_real_, e = NA_integer_, d = NA_character_))
  kk <- which(diff(s) != 0); d <- if (s[1] < 0) "inhibit_first" else "promote_first"
  if (!length(kk))  return(list(ty = if (s[1] < 0) "1 始终抑制" else "3 始终促进",
                                t = as.numeric(MAXTP), e = 0L, d = d))
  if (length(kk) > 1) return(list(ty = "X 多次变号", t = NA_real_, e = NA_integer_, d = d))
  list(ty = if (s[1] < 0) "2 先抑制后促进" else "4 先促进后抑制",
       t = as.numeric(kk), e = 1L, d = d)
}
b <- lapply(pat$co, cls)
pat[, `:=`(ty = sapply(b, `[[`, "ty"), time = sapply(b, `[[`, "t"),
           event = as.integer(sapply(b, `[[`, "e")),
           start_dir = sapply(b, `[[`, "d"))][, co := NULL]
log_msg("完整 ", MAXTP + 1L, " 个滞后的站 ", nrow(pat), " -> 四类站 ", pat[ty != "X 多次变号", .N],
        " (多次变号 ", pat[ty == "X 多次变号", .N], " 站不入样)")

st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
S <- merge(merge(pat[ty != "X 多次变号"],
                 dtc[, c("stat_id", "koppen_group", PRED), with = FALSE], by = "stat_id"),
           st[, .(stat_id, longitude, latitude)], by = "stat_id")
n461 <- nrow(S)
S <- S[complete.cases(S[, PRED, with = FALSE])]
log_msg("协变量齐全 -> 分析样本 n = ", nrow(S), " (缺失掉 ", n461 - nrow(S), " 站)")
log_msg("\n-- 四类构成(n = ", nrow(S), ") --")
print(S[, .N, by = ty][order(ty)])
print(dcast(S[, .N, by = .(start_dir, event)], start_dir ~ event, value.var = "N"))
for (v in PRED) set(S, j = v, value = as.numeric(scale(winz(as.numeric(S[[v]])))))
fwrite(S[, .(stat_id, ty, start_dir, time, event, koppen_group)],
       file.path(OUT, "surv4_sample_tp7.csv"))

# ===========================================================================
hdr("1. 前提检验")
log_msg("-- 每自变量事件数 EPV --")
print(S[, .(n = .N, 事件 = sum(event), 删失 = sum(event == 0),
            EPV = round(sum(event) / length(PRED), 1)), by = start_dir])
vif <- sapply(PRED, function(v) 1 / (1 - summary(lm(as.formula(paste(v, "~",
        paste(setdiff(PRED, v), collapse = "+"))), S))$r.squared))
log_msg("\n最大 VIF = ", round(max(vif), 2), " (", CN[names(which.max(vif))], ")",
        if (max(vif) < 10) "  [<10, 可解释]" else "  <<< 共线警告")
print(round(sort(vif, decreasing = TRUE)[1:5], 2))

log_msg("\n-- 事件时间的并结 --")
print(S[event == 1, .N, by = time][order(time)])
log_msg("时间仅 ", S[event == 1, uniqueN(time)],
        " 个离散取值, 并结重 -> Cox 用 Efron 近似, 另以离散时间 cloglog 交叉核对")

m0 <- coxph(as.formula(paste("Surv(time, event) ~ strata(start_dir) +",
            paste(PRED, collapse = "+"))), data = S, ties = "efron")
zp <- cox.zph(m0)
log_msg("\n-- 比例风险假设(Schoenfeld 残差) --")
zt <- data.table(var = rownames(zp$table), chisq = zp$table[, "chisq"], p = zp$table[, "p"])
print(zt[order(p)][1:6, .(变量 = fifelse(var == "GLOBAL", "全局", CN[var]),
                          chisq = round(chisq, 2), p = signif(p, 3),
                          sig = fifelse(p < .05, "违反", ""))])
log_msg("全局 p = ", signif(zt[var == "GLOBAL", p], 3),
        if (zt[var == "GLOBAL", p] < .05) "  <<< 被拒" else "  [不拒绝]")
fwrite(zt, file.path(OUT, "surv4_phtest_tp7.csv"))

mres <- residuals(m0, type = "martingale")
DM <- dist_km(S$longitude, S$latitude); W <- (DM > 0 & DM <= CUT) * 1
z <- mres - mean(mres); I <- (length(z) / sum(W)) * sum(W * outer(z, z)) / sum(z^2)
set.seed(42)
Ip <- replicate(999, { zz <- sample(z); (length(z)/sum(W)) * sum(W * outer(zz, zz)) / sum(zz^2) })
log_msg("\n-- 残差空间自相关 --")
log_msg("Moran's I = ", round(I, 4), " (", CUT, " km 邻域), 置换 p = ",
        signif((1 + sum(abs(Ip) >= abs(I))) / 1000, 3),
        "  -> 标准误用 Conley 空间 HAC")
log_msg("\n-- 不适用的检验 --")
log_msg("正态性/等方差: KM 与 Cox 不对结局或残差作分布假定, 不检验。")
log_msg("非信息删失: 删失 = 观测窗口到 tp", MAXTP, ", 由设计决定, 假定成立。")

# ===========================================================================
hdr("2. Kaplan-Meier 与 log-rank")
km <- survfit(Surv(time, event) ~ start_dir, data = S)
kt <- data.table(time = km$time, n_risk = km$n.risk, n_event = km$n.event,
                 surv = round(km$surv, 4), lo = round(km$lower, 4), hi = round(km$upper, 4),
                 grp = rep(names(km$strata), km$strata))
print(kt); fwrite(kt, file.path(OUT, "surv4_km_tp7.csv"))
log_msg("\n翻转时间四分位(步, 每步 8 天):")
print(round(quantile(km, probs = c(.25, .5, .75))$quantile, 2))
log_msg("到 tp", MAXTP, " 仍未翻转: ",
        paste(sprintf("%s %.1f%%", GL[c("inhibit_first","promote_first")],
              100 * sapply(c("inhibit_first","promote_first"),
                           function(g) S[start_dir == g, 1 - mean(event)])), collapse = " | "))
s0 <- survdiff(Surv(time, event) ~ start_dir, data = S)
s1 <- survdiff(Surv(time, event) ~ start_dir, data = S, rho = 1)
p_lr <- pchisq(s0$chisq, 1, lower.tail = FALSE)
log_msg("\nlog-rank       chisq = ", round(s0$chisq, 2), ", p = ", signif(p_lr, 3))
log_msg("Peto-Wilcoxon  chisq = ", round(s1$chisq, 2), ", p = ",
        signif(pchisq(s1$chisq, 1, lower.tail = FALSE), 3))
mg <- coxph(Surv(time, event) ~ start_dir, data = S, ties = "efron")
cg <- conley_cox(mg, S$longitude, S$latitude, CUT)
p_sp <- 2 * pnorm(-abs(coef(mg) / cg$se))
log_msg("空间稳健(Conley ", CUT, "km): beta = ", round(coef(mg), 4), ", p = ", signif(p_sp, 3),
        " | 有效独立样本 ~", round(cg$eff_n, 1), " (名义 ", nrow(S), ")")

log_msg("\n-- 分气候区 log-rank(组内, 仅站数>=20 的区) --")
S[koppen_group == "", koppen_group := NA_character_]   # 空串也是缺失
for (g in c("inhibit_first", "promote_first")) {
  d <- S[start_dir == g & !is.na(koppen_group)]
  kp <- d[, .N, by = koppen_group][N >= 20, koppen_group]; d <- d[koppen_group %in% kp]
  if (length(kp) < 2) {
    log_msg(sprintf("%s: 仅 %s 区站数达到 20, 无法做组间比较(其余区: %s)", GL[g],
      paste(sort(kp), collapse = "/"),
      paste(S[start_dir == g & !is.na(koppen_group), .N, by = koppen_group][
        N < 20, sprintf("%s=%d", koppen_group, N)], collapse = ", ")))
    next }
  sk <- survdiff(Surv(time, event) ~ koppen_group, data = d)
  log_msg(sprintf("%s: n=%d, 区=%s, chisq=%.2f, p=%.3g", GL[g], nrow(d),
                  paste(sort(kp), collapse = "/"), sk$chisq,
                  pchisq(sk$chisq, length(sk$n) - 1, lower.tail = FALSE)))
}

# ===========================================================================
hdr("3. Cox 比例风险模型 (18 协变量, Conley 200km)")
cox_tab <- function(d, vars, tag, strat = TRUE) {
  f <- as.formula(paste("Surv(time, event) ~", if (strat) "strata(start_dir) +" else "",
                        paste(vars, collapse = "+")))
  m <- coxph(f, data = d, ties = "efron")
  cs <- conley_cox(m, d$longitude, d$latitude, CUT); bb <- coef(m)
  zz <- cox.zph(m)$table
  o <- data.table(tag = tag, var = vars, cn = CN[vars], beta = bb[vars],
                  se_naive = sqrt(diag(vcov(m)))[vars], se = cs$se[seq_along(vars)],
                  n = m$n, nev = m$nevent, eff_n = cs$eff_n, ph_global = zz["GLOBAL", "p"])
  o[, `:=`(HR = exp(beta), p = 2 * pnorm(-abs(beta / se)),
           p_naive = 2 * pnorm(-abs(beta / se_naive)),
           lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se))]
  o[, ph_p := zz[match(var, rownames(zz)), "p"]][]
}
res <- list(); res$all <- cox_tab(S, PRED, "全体")
log_msg("n = ", res$all$n[1], ", 事件 ", res$all$nev[1], ", 有效独立样本 ~",
        round(res$all$eff_n[1], 1), ", PH 全局 p = ", signif(res$all$ph_global[1], 3))
log_msg("HR > 1 = 该变量越大, 每步翻转风险越高 -> 翻转越早")
print(res$all[order(p)][1:10, .(cn, HR = round(HR, 3), CI = sprintf("[%.2f,%.2f]", lo, hi),
      p = signif(p, 3), p未修正 = signif(p_naive, 3), sig = fifelse(p < .05, "*", ""))])
log_msg("Conley 下显著 ", res$all[p < .05, .N], " 个 | 不修正会有 ",
        res$all[p_naive < .05, .N], " 个")

for (g in c("inhibit_first", "promote_first")) {
  o <- cox_tab(S[start_dir == g], PRED, GL[g], strat = FALSE); res[[g]] <- o
  log_msg("\n#### ", GL[g], " | n = ", o$n[1], ", 事件 ", o$nev[1],
          ", EPV ", round(o$nev[1] / length(PRED), 1),
          ", 有效独立样本 ~", round(o$eff_n[1], 1),
          ", PH 全局 p = ", signif(o$ph_global[1], 3))
  print(o[order(p)][1:6, .(cn, HR = round(HR, 3), CI = sprintf("[%.2f,%.2f]", lo, hi),
        p = signif(p, 3), sig = fifelse(p < .05, "*", ""))])
}
fwrite(rbindlist(res), file.path(OUT, "surv4_cox_coef_tp7.csv"))

log_msg("\n-- 3.1 与离散时间模型交叉核对(同一批站) --")
pp <- S[, .(tp = seq_len(time)), by = .(stat_id, start_dir, time, event)]
pp[, ev := as.integer(tp == time & event == 1)]
pp <- merge(pp, S[, c("stat_id", PRED), with = FALSE], by = "stat_id")
f_d <- as.formula(paste("ev ~ factor(tp) + start_dir +", paste(PRED, collapse = "+")))
g_cll <- glm(f_d, data = pp, family = binomial(link = "cloglog"))
g_lgt <- glm(f_d, data = pp, family = binomial(link = "logit"))
log_msg("Cox(Efron) 与离散 cloglog 系数相关 r = ",
        round(cor(coef(m0)[PRED], coef(g_cll)[PRED]), 4),
        " | 与 logit r = ", round(cor(coef(m0)[PRED], coef(g_lgt)[PRED]), 4))

log_msg("\n-- 3.2 变量集敏感性(固定同一批站, 只改变量集) --")
SETS <- list("18变量(全)" = PRED, "去绿地投资" = setdiff(PRED, "invest"),
             "去绿地投资+城镇化率" = setdiff(PRED, c("invest", "urban_rate")),
             "去五类森林" = setdiff(PRED, c("ever_needle","deci_needle","ever_broad",
                                            "deci_broad","mixedleaf")))
sens <- rbindlist(lapply(names(SETS), function(nm) cox_tab(S, SETS[[nm]], nm)))
for (nm in names(SETS))
  log_msg(sprintf("[%s] 变量 %d 个, 显著 %d 个: %s", nm, length(SETS[[nm]]),
    sens[tag == nm & p < .05, .N],
    paste(sens[tag == nm & p < .05][order(p), sprintf("%s(HR=%.2f)", cn, HR)], collapse = ", ")))
shared <- Reduce(intersect, SETS)
cmw <- dcast(sens[var %in% shared], cn ~ tag, value.var = "beta")
log_msg("(共有变量 ", length(shared), " 个)")
log_msg("\n各变量集间共有变量的系数相关:")
print(round(cor(as.matrix(cmw[, setdiff(names(cmw), "cn"), with = FALSE])), 3))
fwrite(sens, file.path(OUT, "surv4_varset_sens_tp7.csv"))

# ===========================================================================
hdr("4. 两组的驱动因素是否相同")
w <- merge(res$inhibit_first[, .(var, cn, b_inh = beta, p_inh = p, lo_i = lo, hi_i = hi)],
           res$promote_first[, .(var, b_pro = beta, p_pro = p)], by = "var")
w[, `:=`(同向 = sign(b_inh) == sign(b_pro), 差 = b_inh - b_pro)]
log_msg("18 个变量中方向一致的有 ", w[同向 == TRUE, .N], " 个; 两组系数相关 r = ",
        round(cor(w$b_inh, w$b_pro), 3))
log_msg("\n-- 任一组显著的变量 --")
print(w[p_inh < .05 | p_pro < .05][order(p_inh), .(cn,
      抑制组 = round(b_inh, 3), p1 = signif(p_inh, 2),
      促进组 = round(b_pro, 3), p2 = signif(p_pro, 2),
      方向 = fifelse(同向, "一致", "相反"))])
# 交互检验: 效应是否真的两组不同
log_msg("\n-- 组 x 变量交互的 LRT(检验'两组效应相同'的原假设) --")
base <- coxph(as.formula(paste("Surv(time, event) ~ strata(start_dir) +",
              paste(PRED, collapse = "+"))), data = S, ties = "efron")
het <- rbindlist(lapply(PRED, function(v) {
  m2 <- coxph(as.formula(paste("Surv(time, event) ~ strata(start_dir) +",
        paste(PRED, collapse = "+"), "+", v, ":start_dir")), data = S, ties = "efron")
  data.table(var = v, cn = CN[v], lrt = 2 * (m2$loglik[2] - base$loglik[2]), df = 1) }))
het[, p := pchisq(lrt, df, lower.tail = FALSE)][, p_bonf := pmin(1, p * .N)]
print(het[order(p)][1:6, .(cn, lrt = round(lrt, 2), p = signif(p, 3), p_bonf = signif(p_bonf, 3))])
log_msg("Bonferroni 后两组效应确有差异的变量: ",
        if (het[p_bonf < .05, .N]) paste(het[p_bonf < .05, cn], collapse = ", ") else "无")
fwrite(het, file.path(OUT, "surv4_group_interaction_tp7.csv"))

# ===========================================================================
# 图
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), legend.key.size = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))
kd <- data.table(time = km$time, surv = km$surv, lo = km$lower, hi = km$upper,
                 grp = GL[sub("start_dir=", "", rep(names(km$strata), km$strata))])
kd <- rbind(data.table(time = 0, surv = 1, lo = 1, hi = 1, grp = unique(kd$grp)), kd)
setorder(kd, grp, time)
pa <- ggplot(kd, aes(time, surv, colour = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_step(linewidth = .8) + geom_hline(yintercept = .5, linetype = 3, colour = "grey50") +
  scale_colour_manual(values = CL) + scale_fill_manual(values = CL) +
  scale_x_continuous(breaks = 0:MAXTP) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "(a) Kaplan-Meier: 尚未翻转的站点比例",
       subtitle = sprintf("n = %d | log-rank χ²=%.1f, p=%.2g | 空间稳健 p=%.2g",
                          nrow(S), s0$chisq, p_lr, p_sp),
       x = sprintf("滞后 tp(每步 8 天), 截断至 tp%d", MAXTP), y = "S(t)") + th

fp <- rbindlist(list(res$inhibit_first, res$promote_first))
ordv <- res$inhibit_first[order(beta), var]
fp[, lab := factor(CN[var], levels = CN[ordv])]
fp[, sig := fifelse(p < .05, "p<0.05", "n.s.")]
pb <- ggplot(fp, aes(HR, lab, colour = tag, shape = sig)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 position = position_dodge(width = .6), linewidth = .45) +
  geom_point(position = position_dodge(width = .6), size = 1.9) +
  scale_colour_manual(values = CL) +
  scale_shape_manual(values = c("p<0.05" = 16, "n.s." = 1)) + scale_x_log10() +
  labs(title = "(b) Cox 风险比", subtitle = sprintf("每 1 个标准差; 区间为 Conley 200km; HR>1 = 翻转更早; tp = 0-%d", MAXTP),
       x = "风险比 HR(对数轴)", y = NULL) + th

# (c) 两组系数散点: 同一样本、同一变量集, 两个模型的结果放在一张图上
w[, 显著 := fcase(p_inh < .05 & p_pro < .05, "两组都显著",
                  p_inh < .05, "仅抑制组显著", p_pro < .05, "仅促进组显著",
                  default = "均不显著")]
rr <- range(c(w$b_inh, w$b_pro)) * 1.08
pc <- ggplot(w, aes(b_inh, b_pro)) +
  annotate("rect", xmin = 0, xmax = rr[2], ymin = rr[1], ymax = 0, alpha = .05, fill = "red") +
  annotate("rect", xmin = rr[1], xmax = 0, ymin = 0, ymax = rr[2], alpha = .05, fill = "red") +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") +
  geom_point(aes(fill = 显著), shape = 21, size = 2.8, colour = "grey20", stroke = .3) +
  ggrepel::geom_text_repel(aes(label = cn), size = 2.6, max.overlaps = 30,
                           min.segment.length = 0, segment.colour = "grey70") +
  annotate("text", x = rr[2], y = rr[1], hjust = 1, vjust = -0.6, size = 2.7,
           colour = "grey35", label = "方向相反") +
  annotate("text", x = rr[1], y = rr[2], hjust = 0, vjust = 1.4, size = 2.7,
           colour = "grey35", label = "方向相反") +
  scale_fill_manual(values = c("两组都显著" = "#6A3D9A", "仅抑制组显著" = "#2166AC",
                               "仅促进组显著" = "#B2182B", "均不显著" = "white")) +
  coord_equal(xlim = rr, ylim = rr) +
  labs(title = "(c) 两组的驱动因素是否相同",
       subtitle = sprintf("虚线 = 效应相同; 两组系数相关 r = %.2f; 交互 LRT 经 Bonferroni 后 %d 个变量两组效应确有差异(%s)",
                          cor(w$b_inh, w$b_pro), het[p_bonf < .05, .N],
                          paste(het[p_bonf < .05][order(p), cn], collapse = "、")),
       x = "抑制组 β(始终抑制 + 先抑制后促进)",
       y = "促进组 β(始终促进 + 先促进后抑制)") + th

ggsave(file.path(OUT, "fig_survival_fourtype_tp7.png"),
       (pa | pb) / pc + plot_layout(heights = c(1, 1.15)), width = 12, height = 11, dpi = 300)
log_msg("\n已写出 fig_survival_fourtype_tp7.png")
