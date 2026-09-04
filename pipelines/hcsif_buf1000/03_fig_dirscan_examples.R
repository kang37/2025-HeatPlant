#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03_fig_dirscan_examples.R — 反向CCM诊断的示例站 rho-vs-tp 曲线
#
# 参考 Ye et al. 2015 (Sci Rep) 的画法：同一站把 fwd(VPD->SIF) 和
# rev(SIF->VPD) 两条 cross-map skill 曲线叠在一张图上，横轴是滞后 tp，
# 用来直观判断"方向是否干净"(只有一个方向在 tp=0 附近显著) vs
# "双向同时显著"(可能是强耦合/同步伪影，见 docs/02_ccm.md 2026-09-04 条目)。
# 误差棒 = rho_sd(50次随机子抽样在最大库长处的标准差)，数据来自
# 03b_ccm_buf1000_dirscan_examples_sd.R(02 主诊断脚本没存这个量)。
#
# 每类三个示例站(取自 03b 脚本对 9 个示例站的重跑结果)：
#   clean(清晰方向型，tp=0 处仅 fwd 显著)    : 57612, 57517, 58467
#   ambiguous(双向可疑型，tp=0 处两个方向都显著): 59238, 57273, 54818
#   null(无信号对照，两个方向基本都不显著)     : 57171, 59493, 57206
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
# 显式注册中文字体：showtext 在这台机器上不会自动探测到系统中文字体
# (环境问题，和 rEDM 丢失同一批 R 升级导致，见 docs/02_ccm.md)，
# 不指定的话中文全部渲染成占位符。
sysfonts::font_add("cjk", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/ccm_hcsif_buf1000_dirscan_examples"

r <- fread(file.path(OUT, "dirscan_buf1000_20260904_0149.csv"))
r[, sig := !is.na(p_surr) & p_surr < 0.1 & (rho - rho_min) > 0]

TYPE_LAB <- c(clean = "清晰方向型（tp=0 仅 fwd 显著）",
              ambiguous = "双向可疑型（tp=0 两个方向都显著）",
              null = "无信号对照（两个方向均不显著）")
r[, type_lab := factor(TYPE_LAB[example_type], levels = TYPE_LAB)]
# 每类内部按站号排一个顺序号，拼成"清晰方向型 · 57612"这样的分面标题
setorder(r, type_lab, meteo_stat)
r[, panel := factor(paste0(type_lab, " · ", meteo_stat),
                     levels = unique(paste0(type_lab, " · ", meteo_stat)))]
r[, dir_lab := fifelse(direction == "fwd", "fwd  VPD → SIF", "rev  SIF → VPD")]

# optimal tp = 显著tp里rho(基于最大库长的判据，不是60%库长的画图值)最高的那个，
# 每站每方向一条竖线；两个方向都不显著的站(null类)没有线可画。
opt_lines <- r[sig == TRUE, .SD[which.max(rho)], by = .(panel, dir_lab)][, .(panel, dir_lab, tp)]

PAL <- c("fwd  VPD → SIF" = "#2E7D74", "rev  SIF → VPD" = "#C2703A")
th <- theme_minimal(base_size = 10.5, base_family = "cjk") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9.5),
        legend.position = "top", legend.title = element_blank())

p <- ggplot(r, aes(tp, rho_ye, colour = dir_lab, group = dir_lab)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = .3, linetype = "22") +
  geom_vline(data = opt_lines, aes(xintercept = tp, colour = dir_lab),
             linetype = "42", linewidth = .55, alpha = .8, show.legend = FALSE) +
  geom_errorbar(aes(ymin = rho_ye - rho_ye_sd, ymax = rho_ye + rho_ye_sd),
                width = 0, alpha = .35, linewidth = .4) +
  geom_line(linewidth = .7, alpha = .5) +
  geom_point(aes(shape = sig, fill = ifelse(sig, dir_lab, NA)), size = 2.2, stroke = .9) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 21), guide = "none") +
  scale_colour_manual(values = PAL) +
  scale_fill_manual(values = PAL, guide = "none", na.value = "white") +
  scale_x_continuous(breaks = seq(-8, 8, 4)) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  labs(title = "反向 CCM 诊断示例：三种 rho-vs-tp 模式（每类 3 站）",
       subtitle = "实心点=替代检验显著(基于最大库长的 p_surr<0.1 且 Δρ>0，判据不变)；空心点=不显著；灰色竖虚线=tp=0(同期)\n彩色竖虚线=各方向的 optimal tp(显著tp中rho最高点，无显著tp则不画)\n点位置与误差棒改用 Ye et al. 2015 画法：固定库长=60%×最大库长，同一批50次抽样算均值±1SD",
       x = "滞后 tp（负值＝驱动变量取未来值，正值＝取过去值；每步 8 天）",
       y = expression(paste("跨映射技能 ", rho, "（固定库长=60%×最大库长）"))) + th

ggsave(file.path(OUT, "fig_dirscan_examples.png"), p, width = 10.5, height = 8.6, dpi = 300)
cat("已保存:", file.path(OUT, "fig_dirscan_examples.png"), "\n")
