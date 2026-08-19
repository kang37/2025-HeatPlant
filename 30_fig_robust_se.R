#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 30_fig_robust_se.R — 三种标准误的置信区间对照图
#
# 三种方法的点估计完全相同(已验证最大绝对差 = 0), 差别只在区间宽度。
# 因此图的设计是: 每个变量一个点, 叠加三条不同颜色的置信区间(纵向错开)。
# 这样"点不动、线在变"一眼可见, 比并排三张系数图更能说明问题。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"

CN <- c(imperv="不透水面", grass="草地", water="水体",
        ever_needle="常绿针叶林", deci_needle="落叶针叶林",
        ever_broad="常绿阔叶林", deci_broad="落叶阔叶林", mixedleaf="混交林",
        tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
MTH <- c(se_ols = "普通 OLS", se_hc3 = "HC3（修异方差）",
         se_conley200 = "Conley 200 km（异方差＋空间相关）")

r <- fread(file.path(OUT, "robust_se_tp8.csv"))
d <- melt(r[, .(grp, var, beta, se_ols, se_hc3, se_conley200)],
          id.vars = c("grp", "var", "beta"),
          variable.name = "method", value.name = "se")
d[, `:=`(method = factor(MTH[as.character(method)], levels = MTH),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         cn = CN[var],
         grp = factor(ifelse(grp == "inhibit_first",
                             "抑制组（n = 214）", "促进组（n = 171）"),
                      levels = c("抑制组（n = 214）", "促进组（n = 171）")))]
d[, sig := (lo > 0 | hi < 0)]
ord <- d[grp == levels(grp)[1] & method == MTH[1]][order(beta)]
d[, cn := factor(cn, levels = ord$cn)]

PAL <- c("普通 OLS" = "#9BA8A5",
         "HC3（修异方差）" = "#C2703A",
         "Conley 200 km（异方差＋空间相关）" = "#1C6B66")

p <- ggplot(d, aes(beta, cn, colour = method)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .7,
                 position = position_dodge(width = .68)) +
  geom_point(size = 1.5, position = position_dodge(width = .68)) +
  geom_point(data = d[sig == TRUE], size = 2.6, shape = 21, fill = "white",
             stroke = .8, position = position_dodge(width = .68)) +
  facet_wrap(~ grp, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) +
  labs(title = "标准误方法对推断的影响：点估计不变，区间在变",
       subtitle = paste0("三种方法的标准化回归系数完全相同（最大绝对差 = 0），",
                         "差别只在置信区间宽度\n",
                         "空心圈标记该方法下 95% 区间不跨零（即显著）"),
       x = "标准化回归系数（95% 置信区间）", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = .35),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11, hjust = 0),
        legend.position = "top", legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2))
ggsave(file.path(OUT, "fig_robust_se_tp8.png"), p, width = 11, height = 6.8, dpi = 300)
cat("已输出 fig_robust_se_tp8.png\n")
cat("\n各方法显著变量数:\n")
print(dcast(d[, .(n_sig = sum(sig)), by = .(grp, method)], grp ~ method, value.var = "n_sig"))
cat("\n标准误相对普通 OLS 的变化倍数(中位):\n")
w <- dcast(d, grp + var ~ method, value.var = "se")
setnames(w, 3:5, c("ols","hc3","con"))
print(w[, .(HC3倍数 = round(median(hc3/ols), 3),
            Conley倍数 = round(median(con/ols), 3)), by = grp])
