#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 26_fig_rho_tp.R — 因果强度 rho 随滞后 tp 的变化
#
# rho = CCM 在最大库长度处的跨映射技能, 即"从 SIF 的重构吸引子还原 VPD 历史
# 状态"的准确度, 是本研究衡量 VPD->SIF 因果强度的量。
# trend = rho 随库长度增长的相关, 是收敛性判据(CCM 区别于普通相关之处)。
#
# 三个面板:
#   A  全体站点 rho 的中位数与四分位区间
#   B  按柯本气候带分组的 rho 中位数
#
# 原来还有一个 C 面板(收敛性 trend 中位数与"可信站点占比"), 已删除。
# 理由是它不含独立信息: trend = cor(LibSize, rho) 与 rho 本身 Spearman 相关 0.897,
# 而"trend 随 tp 下降"完全由"rho 随 tp 下降"带出来——
#   naive   cor(trend, tp)        = -0.211
#   partial cor(trend, tp | rho)  = +0.058
# 在 rho 的每个五分位箱内, trend 的中位数横跨 tp=0..8 基本不动
# (如 Q5 组 0.833 -> 0.822, Q3 组 0.537 -> 0.550)。
# 既然条件在 rho 上以后 trend 不再随 tp 变, C 面板就只是 A 面板的重画。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
old <- as.data.table(readRDS("data_proc/output_10y_built_up_05_01/station_covariates.rds"))
setnames(old, "meteo_stat_id", "meteo_stat"); old[, meteo_stat := as.integer(meteo_stat)]
cm <- merge(cm, unique(old, by = "meteo_stat")[, .(meteo_stat, koppen_group)],
            by = "meteo_stat", all.x = TRUE)

PAL <- c(B = "#C2703A", C = "#2E7D74", D = "#4B6BA8")
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey35", size = 9.5),
        legend.position = "top", legend.title = element_blank())

# --- A: 全体 ---
a <- cm[, .(med = median(rho, na.rm = TRUE),
            q25 = quantile(rho, .25, na.rm = TRUE),
            q75 = quantile(rho, .75, na.rm = TRUE)), by = tp][order(tp)]
pA <- ggplot(a, aes(tp, med)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2E7D74", alpha = .16) +
  geom_line(colour = "#2E7D74", linewidth = .9) +
  geom_point(colour = "#2E7D74", size = 2.1) +
  scale_x_continuous(breaks = 0:8, sec.axis = sec_axis(~ . * 8, name = "滞后天数",
                                                       breaks = seq(0, 64, 16))) +
  labs(title = "A  因果强度随滞后单调衰减",
       subtitle = "线为中位数，带为四分位区间；n = 904 站",
       x = "滞后 tp（每步 8 天）", y = expression(paste("跨映射技能 ", rho))) + th

# --- B: 分气候带 ---
b <- cm[koppen_group %in% c("B", "C", "D"),
        .(med = median(rho, na.rm = TRUE), n = .N), by = .(tp, koppen_group)]
lab <- c(B = "B 干旱带", C = "C 温带季风", D = "D 温带大陆性")
b[, kg := factor(lab[koppen_group], levels = lab)]
pB <- ggplot(b, aes(tp, med, colour = kg)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = setNames(PAL, lab)) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "B  干旱带的因果信号始终最强",
       subtitle = "各气候带 rho 中位数；干旱带在所有滞后上均高于季风区",
       x = "滞后 tp（每步 8 天）", y = expression(paste("中位 ", rho))) + th

p <- pA / pB + plot_annotation(
  title = "VPD → SIF 因果强度随滞后的变化",
  subtitle = "HCSIF 1000 m 缓冲区，2000–2022 年 5–9 月，8 天合成",
  theme = theme(plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(colour = "grey35")))
ggsave(file.path(OUT, "fig_rho_vs_tp.png"), p, width = 7.2, height = 7.4, dpi = 300)
cat("已保存 fig_rho_vs_tp.png\n")
print(a[, .(tp, rho中位 = round(med, 4), IQR = paste0(round(q25,3), "–", round(q75,3)))])
print(dcast(b, tp ~ kg, value.var = "med")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,4) else x)])
