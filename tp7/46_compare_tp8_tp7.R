#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 46_compare_tp8_tp7.R — 末步截断的敏感性: tp<=8 与 tp<=7 两个口径并排对照
#
# 两张图:
#   fig_compare_tp8_tp7_score.png  得分回归(OLS + Conley 200km)的 18 个系数
#   fig_compare_tp8_tp7_cox.png    Cox 风险比(同样 18 个变量, 分抑制/促进组)
# 一张表:
#   compare_tp8_tp7_summary.csv    两口径的系数、CI、p 与显著性变化
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"

CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率",
        elev="海拔", ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
PAL <- c(`tp ≤ 8（原口径）` = "#B2182B", `tp ≤ 7（末步截断）` = "#2166AC")
th <- theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = .3),
        panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 10, hjust = 0),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9, lineheight = 1.2))

# ---- 1. 得分回归 ----------------------------------------------------------
rd <- function(M) {
  x <- fread(file.path(OUT, sprintf("final_conley200_tp%d.csv", M)))
  x[, .(spec = sprintf("tp ≤ %d%s", M, ifelse(M == 8, "（原口径）", "（末步截断）")),
        grp = fifelse(grepl("^抑制", 组), "inhibit_first", "promote_first"),
        var = 英文名, beta = 标准化系数, lo = CI下限, hi = CI上限, p = p值,
        n = as.integer(sub(".*n = (\\d+).*", "\\1", 组)))]
}
S <- rbind(rd(8L), rd(7L))
NS <- dcast(unique(S[, .(grp, spec, n)]), grp ~ spec, value.var = "n")
gn <- function(g, sp) NS[grp == g][[sp]]        # 按列名取, 避免 dcast 的列序不定
S[, `:=`(spec = factor(spec, levels = names(PAL)), cn = CN[var], sig = p < .05)]
lv <- S[spec == levels(spec)[1] & grp == "inhibit_first"][order(beta), cn]
S[, cn := factor(cn, levels = lv)]
S[, grp_cn := factor(fifelse(grp == "inhibit_first",
      sprintf("抑制组：得分高 = 越早脱离抑制（n: %d → %d）",
              gn("inhibit_first", names(PAL)[1]), gn("inhibit_first", names(PAL)[2])),
      sprintf("促进组：得分高 = 促进维持越久（n: %d → %d）",
              gn("promote_first", names(PAL)[1]), gn("promote_first", names(PAL)[2]))))]
S[, grp_cn := factor(grp_cn, levels = c(grep("^抑制", levels(grp_cn), value = TRUE),
                                        grep("^促进", levels(grp_cn), value = TRUE)))]

p1 <- ggplot(S, aes(beta, cn, colour = spec, shape = sig)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = .55, position = position_dodge(width = .65)) +
  geom_point(aes(fill = spec), size = 2.2, stroke = .7,
             position = position_dodge(width = .65)) +
  facet_wrap(~ grp_cn, nrow = 1) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL, guide = "none") +
  labs(title = "转变时间得分回归：末步截断前后的对照",
       subtitle = "OLS 标准化系数 + Conley 空间 HAC 95% 置信区间（截断 200 km）；实心 = p<0.05",
       x = "标准化回归系数", y = NULL) + th
ggsave(file.path(OUT, "fig_compare_tp8_tp7_score.png"), p1,
       width = 10.5, height = 6.4, dpi = 300)

# ---- 2. Cox 风险比 --------------------------------------------------------
rc <- function(M) {
  f <- if (M == 8L) "surv4_cox_coef.csv" else sprintf("surv4_cox_coef_tp%d.csv", M)
  x <- fread(file.path(OUT, f))[tag %in% c("抑制组", "促进组")]
  x[, .(spec = sprintf("tp ≤ %d%s", M, ifelse(M == 8, "（原口径）", "（末步截断）")),
        grp = tag, var, cn, HR, lo, hi, p, n, nev)]
}
K <- rbind(rc(8L), rc(7L))
NK <- dcast(unique(K[, .(grp, spec, n, nev)])[, .(grp, spec, lab = sprintf("%d站/%d事件", n, nev))],
            grp ~ spec, value.var = "lab")
gk <- function(g, sp) NK[grp == g][[sp]]        # 同上, 按列名取
K[, `:=`(spec = factor(spec, levels = names(PAL)), sig = p < .05)]
lvk <- K[spec == levels(spec)[1] & grp == "抑制组"][order(HR), cn]
K[, cn := factor(cn, levels = lvk)]
K[, grp_cn := factor(fifelse(grp == "抑制组",
      sprintf("抑制组 → 翻转为促进（%s → %s）",
              gk("抑制组", names(PAL)[1]), gk("抑制组", names(PAL)[2])),
      sprintf("促进组 → 翻转为抑制（%s → %s）",
              gk("促进组", names(PAL)[1]), gk("促进组", names(PAL)[2]))))]
K[, grp_cn := factor(grp_cn, levels = c(grep("^抑制", levels(grp_cn), value = TRUE),
                                        grep("^促进", levels(grp_cn), value = TRUE)))]

p2 <- ggplot(K, aes(HR, cn, colour = spec, shape = sig)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = .35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = .55, position = position_dodge(width = .65)) +
  geom_point(aes(fill = spec), size = 2.2, stroke = .7,
             position = position_dodge(width = .65)) +
  facet_wrap(~ grp_cn, nrow = 1, scales = "free_x") +
  scale_x_log10() +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL, guide = "none") +
  labs(title = "Cox 风险比：末步截断前后的对照",
       subtitle = "每 1 个标准差；区间为 Conley 200 km；HR>1 = 翻转更早；实心 = p<0.05",
       x = "风险比 HR（对数轴）", y = NULL) + th
ggsave(file.path(OUT, "fig_compare_tp8_tp7_cox.png"), p2,
       width = 10.5, height = 6.4, dpi = 300)

# ---- 3. 汇总表 ------------------------------------------------------------
w <- merge(S[spec == levels(spec)[1], .(grp, var, cn, b8 = beta, p8 = p)],
           S[spec == levels(spec)[2], .(grp, var, b7 = beta, p7 = p)],
           by = c("grp", "var"))
w[, 显著性变化 := fcase(p8 < .05 & p7 < .05, "两口径均显著",
                        p8 < .05 & p7 >= .05, "仅原口径显著（截断后消失）",
                        p8 >= .05 & p7 < .05, "仅截断后显著（新出现）",
                        default = "两口径均不显著")]
w[, `:=`(组 = fifelse(grp == "inhibit_first", "抑制组", "促进组"), 变量 = cn,
         系数_tp8 = round(b8, 3), p_tp8 = signif(p8, 3),
         系数_tp7 = round(b7, 3), p_tp7 = signif(p7, 3),
         系数变化 = round(b7 - b8, 3), 符号一致 = sign(b7) == sign(b8))]
tab <- w[order(组, -abs(b8)), .(组, 变量, 系数_tp8, p_tp8, 系数_tp7, p_tp7,
                                系数变化, 符号一致, 显著性变化)]
fwrite(tab, file.path(OUT, "compare_tp8_tp7_summary.csv"))
cat("\n===== 得分回归: 两口径系数对照 =====\n"); print(tab)
cat("\n两口径系数相关 r =",
    round(cor(w$b8, w$b7), 3), "| 符号一致", w[符号一致 == TRUE, .N], "/", nrow(w), "\n")
cat("\n显著性变化计数:\n"); print(w[, .N, by = .(组, 显著性变化)][order(组, -N)])
cat("\n已输出 fig_compare_tp8_tp7_score.png / fig_compare_tp8_tp7_cox.png / compare_tp8_tp7_summary.csv\n")
