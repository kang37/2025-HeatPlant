#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 32_final_conley.R — 最终回归结果: 转变时间得分 ~ 三维度协变量
#                     标准误统一用 Conley 空间 HAC, 截断 200 km
#
# 截断距离的依据(见 31_conley_cutoff.R):
#   残差相关图显示空间相关集中在 200 km 内; 标准误敏感曲线在 150-200 km
#   出现平台; 有效样本量(n/平均邻居)在 200 km 时为 12, 再大即不稳。
#   三条依据的交点即 200 km。
#
# 点估计与普通 OLS 完全相同(已验证最大绝对差=0), 变的只有标准误与 p 值。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; M <- 8L; CUT <- 200

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

r <- fread(file.path(OUT, sprintf("robust_se_tp%d.csv", M)))
d <- r[, .(grp, var, beta, se = se_conley200)]
d[, `:=`(z = beta / se)]
d[, `:=`(p = 2 * pnorm(-abs(z)),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         dimen = factor(dim_of[var], levels = names(DIMS)),
         cn = CN[var],
         grp_cn = factor(fifelse(grp == "inhibit_first",
                                 "抑制组：得分高 = 越早脱离抑制（n = 214）",
                                 "促进组：得分高 = 促进维持越久（n = 171）"),
                levels = c("抑制组：得分高 = 越早脱离抑制（n = 214）",
                           "促进组：得分高 = 促进维持越久（n = 171）")))]
d[, sig := p < .05]
d[, star := fcase(p < .001, "***", p < .01, "**", p < .05, "*", default = "")]

fwrite(d[order(grp_cn, dimen, -abs(beta)),
         .(组 = grp_cn, 维度 = dimen, 变量 = cn, 英文名 = var,
           标准化系数 = round(beta, 4), Conley标准误 = round(se, 4),
           z值 = round(z, 3), p值 = signif(p, 4), 显著性 = star,
           CI下限 = round(lo, 4), CI上限 = round(hi, 4))],
       file.path(OUT, sprintf("final_conley%d_tp%d.csv", CUT, M)))

ord <- d[grp == "inhibit_first"][order(dimen, beta)]
d[, cn := factor(cn, levels = ord$cn)]
PAL <- c(`土地利用与空间结构` = "#8C6D3F", `气候与自然环境` = "#2E7D74",
         `社会经济与人类活动` = "#7A4A78")
p <- ggplot(d, aes(beta, cn, colour = dimen)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = .6) +
  geom_point(aes(shape = sig, fill = dimen), size = 2.5, stroke = .75) +
  geom_text(data = d[sig == TRUE], aes(x = hi, label = star), hjust = -0.35,
            size = 3.4, show.legend = FALSE) +
  facet_wrap(~ grp_cn, nrow = 1) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.05, .12))) +
  labs(title = "转变时间得分的影响因素",
       subtitle = paste0("严格分类，tp = 0–8；标准化回归系数与 95% 置信区间；",
                         "标准误为 Conley 空间 HAC（截断 200 km）\n",
                         "* p<0.05  ** p<0.01  *** p<0.001"),
       x = "标准化回归系数", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = .35),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10.5, hjust = 0),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2))
ggsave(file.path(OUT, sprintf("fig_final_conley%d.png", CUT)), p,
       width = 10.2, height = 6.4, dpi = 300)

cat("=== 抑制组（按 |系数| 排序）===\n")
print(d[grp == "inhibit_first"][order(-abs(beta)),
      .(变量 = cn, 维度 = dimen, 系数 = round(beta,3),
        CI = sprintf("[%.2f, %.2f]", lo, hi), p = signif(p,3), 显著性 = star)])
cat("\n=== 促进组（按 |系数| 排序）===\n")
print(d[grp == "promote_first"][order(-abs(beta)),
      .(变量 = cn, 维度 = dimen, 系数 = round(beta,3),
        CI = sprintf("[%.2f, %.2f]", lo, hi), p = signif(p,3), 显著性 = star)])
cat("\n显著变量数: 抑制组", d[grp=="inhibit_first" & sig, .N],
    "| 促进组", d[grp=="promote_first" & sig, .N], "\n")
