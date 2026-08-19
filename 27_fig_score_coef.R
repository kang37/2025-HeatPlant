#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 27_fig_score_coef.R — 转变时间得分回归的系数表与系数图
#
# 输入: score_tp<M>_{inhibit,promote}_first.csv (由 24_score_regression.R 产出)
# 输出: score_tp<M>_coef_table.csv   完整系数表(两组 x 18 变量)
#       fig_score_coef_tp<M>.png     系数图(点=标准化系数, 线=95% CI)
#
# 系数为标准化回归系数, 故各变量可直接横向比较。
# 实心点 = OLS p<0.05; 方框 = 有序 logit 也显著(两法一致)。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"
aa <- commandArgs(trailingOnly = TRUE)
M <- if (any(aa == "--maxtp")) as.integer(aa[which(aa == "--maxtp") + 1]) else 8L

DIMS <- list(
  `土地利用与空间结构` = c("imperv","grass","water","ever_needle","deci_needle",
                           "ever_broad","deci_broad","mixedleaf"),
  `气候与自然环境`     = c("tavg","rh","cloud","precip","rsds_mean","rsds_sd","elev"),
  `社会经济与人类活动` = c("ntl","invest","urban_rate"))
CN <- c(imperv="不透水面", grass="草地", water="水体",
        ever_needle="常绿针叶林", deci_needle="落叶针叶林",
        ever_broad="常绿阔叶林", deci_broad="落叶阔叶林", mixedleaf="混交林",
        tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
dim_of <- setNames(rep(names(DIMS), lengths(DIMS)), unlist(DIMS))

grp_lab <- c(inhibit_first = "抑制组：得分高 = 越早脱离抑制",
             promote_first = "促进组：得分高 = 促进维持越久")
d <- rbindlist(lapply(names(grp_lab), function(g) {
  x <- fread(file.path(OUT, sprintf("score_tp%d_%s.csv", M, g)))
  x[, grp := grp_lab[g]][] }))
d[, grp := factor(grp, levels = grp_lab)]   # 抑制组在左
d[, `:=`(dimen = factor(dim_of[var], levels = names(DIMS)),
         cn = CN[var],
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         sig = p < .05, sig2 = p < .05 & p_ord < .05)]

# 按维度分组、组内按抑制组的系数排序，两个面板共用同一纵轴顺序
ord <- d[grp == grp_lab[1]][order(dimen, beta)]
d[, cn := factor(cn, levels = ord$cn)]

fwrite(d[order(grp, dimen, -abs(beta)),
         .(组 = grp, 维度 = dimen, 变量 = cn, 英文名 = var,
           标准化系数 = round(beta, 4), 标准误 = round(se, 4),
           CI下限 = round(lo, 4), CI上限 = round(hi, 4), p值 = signif(p, 3),
           有序logit系数 = round(or_coef, 4), 有序logit_p = signif(p_ord, 3))],
       file.path(OUT, sprintf("score_tp%d_coef_table.csv", M)))

PAL <- c(`土地利用与空间结构` = "#8C6D3F",
         `气候与自然环境`     = "#2E7D74",
         `社会经济与人类活动` = "#7A4A78")
p <- ggplot(d, aes(beta, cn, colour = dimen)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .55, alpha = .85) +
  geom_point(aes(shape = sig, fill = dimen), size = 2.4, stroke = .7) +
  geom_point(data = d[sig2 == TRUE], shape = 22, size = 4.4,
             fill = NA, colour = "grey25", stroke = .45) +
  facet_wrap(~ grp, nrow = 1) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(title = sprintf("转变时间得分的影响因素（严格分类，tp = 0–%d）", M),
       subtitle = paste0("点为标准化回归系数，横线为 95% 置信区间；",
                         "实心 = OLS 显著，外框 = 有序 logit 亦显著\n",
                         "抑制组 n = 214，促进组 n = 171"),
       x = "标准化回归系数", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = .35),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11, hjust = 0),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2))
ggsave(file.path(OUT, sprintf("fig_score_coef_tp%d.png", M)), p,
       width = 10, height = 6.4, dpi = 300)
cat("已输出 score_tp", M, "_coef_table.csv 与 fig_score_coef_tp", M, ".png\n", sep = "")
print(d[order(grp, -abs(beta)), .(grp, cn, beta = round(beta,3),
        CI = sprintf("[%.2f, %.2f]", lo, hi), p = signif(p,3),
        p_ord = signif(p_ord,3))][1:12])
