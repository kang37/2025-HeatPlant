#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 56_fig_reclass.R — multi_flip 的诊断与新判据的效果图
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(sysfonts)
  library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"
INH <- "#1C6B66"; PRO <- "#C2703A"
th <- function(b = 10) theme_minimal(base_size = b, base_family = "heiti") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = .3),
        strip.text = element_text(face = "bold", hjust = 0, size = b - .5),
        plot.title = element_text(face = "bold", size = b + 2),
        plot.subtitle = element_text(colour = "grey35", size = b - 1.5, lineheight = 1.2))

dg <- fread(file.path(OUT, "multiflip_diag.csv"))
rc <- fread(file.path(OUT, "reclassify.csv"))
cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)

# --- A. 变号次数分布 vs 纯噪声 ---------------------------------------------
n <- nrow(dg)
obs <- dg[, .N, by = .(k = nchg)]
nb <- data.table(k = 0:8, 噪声 = choose(8, 0:8) * 2 / 2^9 * n)
a <- merge(nb, obs, by = "k", all.x = TRUE); a[is.na(N), N := 0]
am <- melt(a, id.vars = "k", variable.name = "src", value.name = "v")
am[, src := factor(ifelse(src == "N", "实测", "纯噪声基准"), c("实测","纯噪声基准"))]
pA <- ggplot(am, aes(factor(k), v, fill = src)) +
  geom_col(position = position_dodge(.75), width = .68) +
  annotate("rect", xmin = .4, xmax = 2.6, ymin = -Inf, ymax = Inf,
           fill = INH, alpha = .06) +
  annotate("text", x = 1.5, y = max(am$v)*1.12, label = "严格判据接受区",
           family = "heiti", size = 2.9, colour = INH) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  scale_fill_manual(values = c(实测 = INH, 纯噪声基准 = "grey72"), name = NULL) +
  labs(subtitle = "A　9 个滞后里的变号次数：实测 vs 纯噪声",
       x = "变号次数", y = "站数") + th() + theme(legend.position = "top")

# --- B. multi_flip 占比 vs |coef| ------------------------------------------
dg[, abin := cut(amed, c(0,.01,.02,.05,.1,1),
                 labels = c("<0.01","0.01–0.02","0.02–0.05","0.05–0.10",">0.10"))]
b <- dg[, .(占比 = mean(nchg > 1), n = .N), by = abin][order(abin)]
pB <- ggplot(b, aes(abin, 占比)) +
  geom_col(fill = PRO, width = .62) +
  geom_text(aes(label = sprintf("%.0f%%\nn=%d", 100*占比, n)), vjust = -.25,
            size = 2.7, family = "heiti", lineheight = .95) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.02),
                     expand = expansion(mult = c(0, .1))) +
  labs(subtitle = "B　多次变号的比例随 |mean_coef| 单调下降",
       x = "该站 9 个滞后的 |mean_coef| 中位数", y = "被判为多次变号的比例") + th()

# --- C. 示例: 变号发生在近零系数上 -----------------------------------------
ex <- rc[S0 == "multi_flip" & A >= .92][order(-A)][seq(1, .N, length.out = 6)]$stat_id
E <- cm[meteo_stat %in% ex, .(stat_id = meteo_stat, tp, co = mean_coef)]
E <- merge(E, rc[, .(stat_id, A, s2_start, s2_k, S2)], by = "stat_id")
E[, fit := { s <- ifelse(s2_start == "inhibit", -1, 1)
             ifelse(is.na(s2_k) | tp < s2_k, s, -s) }, by = stat_id]
E[, lab := sprintf("站 %d ｜ 一致度 %.2f", stat_id, A)]
seg <- E[, .(tp = c(tp, max(tp) + .5), fit = c(fit, fit[.N])), by = lab]
pC <- ggplot(E, aes(tp, co)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = .35) +
  geom_col(aes(fill = co < 0), width = .68) +
  geom_step(data = seg, aes(tp - .5, fit * max(abs(E$co)) * .92), colour = "grey25",
            linewidth = .55, linetype = 2) +
  facet_wrap(~ lab, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = INH, `FALSE` = PRO),
                    labels = c("正（促进）","负（抑制）"), name = NULL) +
  labs(subtitle = paste0("C　六个原本被丢弃的「多次变号」站：柱是 mean_coef，",
                         "虚线是加权变点模型拟合出的符号"),
       x = "滞后 tp（每步 8 天）", y = "mean_coef") +
  th() + theme(legend.position = "top")

# --- D. 判据对比 -----------------------------------------------------------
d <- data.table(
  判据 = c("S0 严格（现状）","S1 零带 τ=0.01","S2 加权变点"),
  可分类 = c(461, 543, 898),
  复现率 = c(1, .901, 1))
d[, 判据 := factor(判据, levels = 判据)]
pD <- ggplot(d, aes(判据, 可分类, fill = 判据)) +
  geom_col(width = .6) +
  geom_hline(yintercept = 898, linetype = 2, colour = "grey45", linewidth = .35) +
  geom_text(aes(label = sprintf("%d 站\n(%.0f%%)\n复现旧标签 %.0f%%",
                                可分类, 100*可分类/898, 100*复现率)),
            vjust = -.15, size = 2.7, family = "heiti", lineheight = 1) +
  scale_fill_manual(values = c("grey72", PRO, INH), guide = "none") +
  scale_y_continuous(limits = c(0, 1150), expand = expansion(mult = c(0, .02))) +
  labs(subtitle = "D　可分类站数：虚线为全部 898 站", x = NULL, y = "可分类站数") + th()

# --- E. 折外 AUC 对比 -------------------------------------------------------
ev <- fread(file.path(OUT, "reclass_eval.csv"))
ev <- ev[设定 %like% "方向|符号|翻转"]
ev[, 层 := ifelse(设定 %like% "翻转", "翻转层", "起始方向层")]
ev[, 口径 := ifelse(设定 %like% "^新", "新判据 S2", "旧判据")]
em <- melt(ev, id.vars = c("设定","层","口径","n"),
           measure.vars = c("AUC_logit","AUC_xgb"),
           variable.name = "模型", value.name = "AUC")
em[, 模型 := ifelse(模型 == "AUC_logit", "logit", "XGBoost")]
em[, 设定 := factor(设定, levels = ev$设定)]
pE <- ggplot(em, aes(AUC, 设定, colour = 层, shape = 模型)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey45", linewidth = .35) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(起始方向层 = INH, 翻转层 = PRO), name = NULL) +
  scale_shape_manual(values = c(logit = 19, XGBoost = 17), name = NULL) +
  scale_y_discrete(limits = rev) +
  labs(subtitle = "E　空间分块折外 AUC：换判据前后，方向层没掉、翻转层仍然弱",
       x = "AUC（0.5 = 抛硬币）", y = NULL) +
  th() + theme(legend.position = "top", axis.text.y = element_text(size = 8))

p <- (pA | pB) / pC / (pD | pE) + plot_layout(heights = c(1, 1.15, 1)) +
  plot_annotation(
    title = "一半站被判成「多次变号」不是站的问题，是判据的问题",
    subtitle = paste0("变号处的 |mean_coef| 中位数只有 0.0056，是全部系数中位数（0.0236）的四分之一；",
                      "71% 的变号发生在 |coef| < 0.01 处。\n",
                      "加权变点判据（S2）在近零滞后上不强行读符号：",
                      "它把 898 站全部定型，且在 461 个原本可定型的站上 100% 复现旧标签。"),
    theme = th(11))
ggsave(file.path(OUT, "fig_reclass.png"), p, width = 12.5, height = 13.5, dpi = 300)
cat("已输出 fig_reclass.png\n")
