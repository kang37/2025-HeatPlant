#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 39_survival_exclude_transient.R — 剔除"瞬时翻转"站的敏感性分析
#
# 同一个结局(持续翻转), 三种对待竞争事件(status==2)的方式, 并排对比:
#
#   A 剔除   把这些站整个从样本里拿掉 (= 23_hazard_all.R 的旧口径)
#            它们变号前那几步"还没翻转"的信息也一并丢掉。
#   B 删失   保留在样本里, 在其首次变号的那一步删失(cause-specific 风险)
#            变号前的人时(person-time)被用上, 变号后不再跟踪。
#   C 保留   Fine-Gray 子分布风险, 变号后仍留在风险集里。
#
# A 与 B 的差别正是问题所在: 一个站要不要被删, 取决于它**首次变号之后**
# 的走势。用未来决定现在是否入样, 是典型的信息性选择。
#
# 输出: surv_exclude_compare.csv / fig_exclude_compare.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
  library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; CUT <- 200; MAXTP <- 8L
log_msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) log_msg("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))

D <- readRDS(file.path(OUT, "surv_base.rds"))
PRED <- c("imperv","grass","water","ever_needle","deci_needle","ever_broad",
          "deci_broad","mixedleaf","tavg","rh","cloud","precip","rsds_mean",
          "rsds_sd","elev","ntl","invest","urban_rate")
CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率",
        elev="海拔", ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
GL <- c(inhibit_first = "抑制组", promote_first = "促进组")

dist_km <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2) }
conley_cox <- function(fit, lon, lat, cut, id = NULL, key = NULL) {
  if (is.null(id)) Dfb <- as.matrix(residuals(fit, type = "dfbeta")) else {
    Dfb <- as.matrix(residuals(fit, type = "dfbeta", collapse = id))
    ord <- match(rownames(Dfb), as.character(key)); stopifnot(!anyNA(ord))
    lon <- lon[ord]; lat <- lat[ord] }
  DM <- dist_km(lon, lat); K <- DM; K[] <- pmax(0, 1 - DM / cut)
  list(se = sqrt(pmax(diag(t(Dfb) %*% K %*% Dfb), 0)),
       eff_n = nrow(DM) / max(1, mean(rowSums(DM <= cut)))) }

hdr("0. 三种设定下的样本")
cmpn <- D[, .(全部 = .N, 删失0 = sum(status == 0), 事件1 = sum(status == 1),
              瞬时2 = sum(status == 2)), by = start_dir]
cmpn[, `A剔除后样本` := 全部 - 瞬时2]
print(cmpn)
log_msg("\nA 剔除: n = ", D[status != 2, .N], " 站 (旧口径, 与 23_hazard_all.R 同)")
log_msg("B 删失 / C 保留: n = ", nrow(D), " 站")
log_msg("三者的**事件数完全相同**(持续翻转 ", D[status == 1, .N],
        " 例) —— 差别只在'谁留在风险集里、留多久'")

# 剔除掉的人时有多少?
pt_all <- D[, sum(time)]; pt_ex <- D[status != 2, sum(time)]
log_msg("\n人时(person-period 总步数): 全样本 ", pt_all, " 步 | 剔除后 ", pt_ex,
        " 步, 丢掉 ", pt_all - pt_ex, " 步 (",
        round(100 * (pt_all - pt_ex) / pt_all, 1), "%)")
log_msg("这些被丢掉的步里, 站点确实处于'尚未持续翻转'的状态, 本应进入风险集分母。")

# ===========================================================================
hdr("1. 描述: 三种设定给出的累积发生率")
res_cif <- list()
for (g in c("inhibit_first", "promote_first")) {
  d <- D[start_dir == g]
  kA <- survfit(Surv(time, status == 1L) ~ 1, data = d[status != 2])   # 剔除
  kB <- survfit(Surv(time, status == 1L) ~ 1, data = d)                # 删失
  aj <- survfit(Surv(time, status_f) ~ 1, data = d)                    # 竞争风险
  fA <- 1 - summary(kA, times = MAXTP)$surv
  fB <- 1 - summary(kB, times = MAXTP)$surv
  fC <- aj$pstate[nrow(aj$pstate), "sustained"]
  log_msg(sprintf("\n%s: tp8 处'持续翻转'的累积发生率", GL[g]))
  log_msg(sprintf("  A 剔除 (1-KM, n=%d)      %.1f%%", d[status != 2, .N], 100 * fA))
  log_msg(sprintf("  B 删失 (1-KM, n=%d)      %.1f%%", nrow(d), 100 * fB))
  log_msg(sprintf("  C 竞争风险 (AJ, n=%d)    %.1f%%   <-- 唯一无偏的那个", nrow(d), 100 * fC))
  log_msg(sprintf("  经验占比 事件1/全部 = %d/%d = %.1f%%",
                  d[status == 1, .N], nrow(d), 100 * d[status == 1, .N] / nrow(d)))
  res_cif[[g]] <- data.table(grp = GL[g], A_剔除 = fA, B_删失 = fB, C_竞争风险 = fC)
}
log_msg("\nA 和 B 都是 1-KM, 都把'不是持续翻转'当成了'还有机会持续翻转', 故都高估;")
log_msg("A 还额外丢掉了人时, 所以两者高估的量还不一样。C 与经验占比一致。")

# ===========================================================================
hdr("2. 回归: 三种设定的系数并排")
fit_one <- function(d, mode, g) {
  if (mode == "A") d <- d[status != 2]
  if (mode %in% c("A", "B")) {
    m <- coxph(as.formula(paste("Surv(time, status == 1L) ~", paste(PRED, collapse = "+"))),
               data = d, ties = "efron")
    cs <- conley_cox(m, d$longitude, d$latitude, CUT)
    b <- coef(m); se <- cs$se; ph <- cox.zph(m)$table["GLOBAL", "p"]
    n <- m$n; nev <- m$nevent; en <- cs$eff_n
  } else {
    fg <- finegray(Surv(time, status_f) ~ ., data = d, etype = "sustained")
    m <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
               paste(PRED, collapse = "+"))), weights = fgwt, data = fg,
               robust = TRUE, id = stat_id)
    cs <- conley_cox(m, d$longitude, d$latitude, CUT, id = fg$stat_id, key = d$stat_id)
    b <- coef(m); se <- cs$se; ph <- NA_real_
    n <- nrow(d); nev <- sum(d$status == 1L); en <- cs$eff_n
  }
  data.table(grp = GL[g], mode = mode, var = PRED, beta = b[PRED], se = se[seq_along(PRED)],
             n = n, nev = nev, eff_n = en, ph_global = ph)
}
MODE <- c(A = "A 剔除(旧口径)", B = "B 删失(病因别)", C = "C 保留(Fine-Gray)")
all_res <- rbindlist(lapply(c("inhibit_first", "promote_first"), function(g)
  rbindlist(lapply(c("A", "B", "C"), function(m) fit_one(D[start_dir == g], m, g)))))
all_res[, `:=`(HR = exp(beta), p = 2 * pnorm(-abs(beta / se)),
               lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se))]
all_res[, mode_lab := MODE[mode]]
fwrite(all_res, file.path(OUT, "surv_exclude_compare.csv"))

for (g in unique(all_res$grp)) {
  s <- all_res[grp == g]
  log_msg("\n########## ", g, " ##########")
  print(unique(s[, .(mode = mode_lab, n, 事件 = nev, 有效独立样本 = round(eff_n, 1),
                     PH全局p = signif(ph_global, 3))]))
  w <- dcast(s, var ~ mode, value.var = c("beta", "p"))
  w[, cn := CN[var]]
  setorder(w, p_A)
  print(w[1:8, .(cn, A_beta = round(beta_A, 3), A_p = signif(p_A, 2),
                 B_beta = round(beta_B, 3), B_p = signif(p_B, 2),
                 C_beta = round(beta_C, 3), C_p = signif(p_C, 2))])
  log_msg("显著变量数: A ", w[p_A < .05, .N], " | B ", w[p_B < .05, .N],
          " | C ", w[p_C < .05, .N])
  chg <- w[(p_A < .05) != (p_B < .05)]
  log_msg("A 与 B 显著性结论不同的变量: ",
          if (nrow(chg)) paste(sprintf("%s(A p=%.3g, B p=%.3g)", chg$cn, chg$p_A, chg$p_B),
                               collapse = "; ") else "无")
  sg <- w[sign(beta_A) != sign(beta_B)]
  log_msg("A 与 B 符号相反的变量: ",
          if (nrow(sg)) paste(sg$cn, collapse = ", ") else "无")
  log_msg("A 与 B 系数相关 r = ", round(cor(w$beta_A, w$beta_B), 3),
          " | 最大绝对变动 ", round(max(abs(w$beta_A - w$beta_B)), 3))
}

# ===========================================================================
hdr("3. 剔除是否有偏: 被剔者与保留者的协变量差异")
D[, dropped := status == 2L]
bal <- rbindlist(lapply(PRED, function(v) {
  t <- t.test(D[[v]] ~ D$dropped)
  data.table(var = v, cn = CN[v], 保留 = t$estimate[1], 被剔 = t$estimate[2],
             标准化差 = (t$estimate[2] - t$estimate[1]) / sd(D[[v]]), p = t$p.value) }))
bal[, absd := abs(标准化差)]; setorder(bal, -absd); bal[, absd := NULL]
print(bal[1:8, .(cn, 保留 = round(保留, 3), 被剔 = round(被剔, 3),
                 标准化差 = round(标准化差, 3), p = signif(p, 3))])
log_msg("\n标准化差 >0.1 的变量数: ", bal[abs(标准化差) > .1, .N], " / ", nrow(bal))
log_msg("(倾向得分文献常用 0.1 作'两组不可比'的阈值)")
fwrite(bal, file.path(OUT, "surv_exclude_balance.csv"))

# ===========================================================================
# 图: 三种设定的系数对照
CLM <- c("A 剔除(旧口径)" = "#B2182B", "B 删失(病因别)" = "#2166AC",
         "C 保留(Fine-Gray)" = "#1A7F37")
pf <- copy(all_res)
ordv <- pf[grp == "抑制组" & mode == "B"][order(beta), var]
pf[, lab := factor(CN[var], levels = CN[ordv])]
pf[, sig := fifelse(p < .05, "p<0.05", "n.s.")]
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), legend.key.size = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

pp <- ggplot(pf, aes(HR, lab, colour = mode_lab, shape = sig)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 position = position_dodge(width = .7), linewidth = .45) +
  geom_point(position = position_dodge(width = .7), size = 1.8) +
  scale_colour_manual(values = CLM) +
  scale_shape_manual(values = c("p<0.05" = 16, "n.s." = 1)) +
  scale_x_log10() + facet_wrap(~ grp) +
  labs(title = "剔除 / 删失 / 保留: 同一结局(持续翻转)的三种处理",
       subtitle = sprintf("每 1 个标准差, 区间为 Conley %dkm。三者事件数相同, 差别只在风险集", CUT),
       x = "风险比 HR(对数轴)", y = NULL) + th
ggsave(file.path(OUT, "fig_exclude_compare.png"), pp, width = 11, height = 7, dpi = 300)
log_msg("\n已写出 fig_exclude_compare.png")
