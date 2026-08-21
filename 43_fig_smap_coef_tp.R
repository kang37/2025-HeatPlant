#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 43_fig_smap_coef_tp.R — S-map 系数(因果效应的方向与大小)随滞后 tp 的变化
#
# 与 26 号的分工:
#   26 号画 rho = 跨映射技能, 衡量"因果关系有多可靠"(证据强度)
#   本脚本画 mean_coef = S-map 同期偏导数 dSIF/dVPD, 衡量"效应本身有多大、朝哪边"
# 两者是不同的量: rho 只能为正且无方向, coef 带符号。
#
# 分组按 tp=0 处 mean_coef 的符号:
#   抑制组 coef(tp=0) < 0    促进组 coef(tp=0) > 0
# 这与 23_hazard_all.R 的 start_dir 定义一致(那里 start_dir 就是 tp=0 的符号),
# 区别是本脚本不剔除多次变号的站——描述性图不该先做筛选。
#
# 【必须声明的一件事】分组用的就是 tp=0 的符号, 所以 tp=0 处两组分开是
# 构造出来的, 不是发现。有信息的是: (a) 各组衰减回零的快慢,
# (b) 是否真的越过零(方向反转)而不只是收缩。C 面板用"负号占比"直接看这件事,
# 它对均值回归的敏感度比系数曲线低。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"
log_msg <- function(...) cat(..., "\n", sep = "")

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setnames(cm, "meteo_stat", "stat_id")
setorder(cm, stat_id, tp)

# 只保留 9 个滞后齐全的站, 否则各 tp 的样本不同, 曲线的升降会混进组成变化
full <- cm[, .(n = sum(!is.na(mean_coef))), by = stat_id][n == 9, stat_id]
cm <- cm[stat_id %in% full]
log_msg("滞后齐全的站: ", length(full), " / ", uniqueN(fread(
  "data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")$meteo_stat))

g0 <- cm[tp == 0, .(stat_id, grp = fifelse(mean_coef < 0, "抑制组", "促进组"))]
cm <- merge(cm, g0, by = "stat_id")
cm[, grp := factor(grp, levels = c("抑制组", "促进组"))]
NG <- g0[, .N, by = grp]
log_msg("分组(按 tp=0 符号): ", paste(sprintf("%s %d 站", NG$grp, NG$N), collapse = " | "))

qs <- function(x) as.list(setNames(c(median(x, na.rm = TRUE),
        quantile(x, c(.25, .75), na.rm = TRUE)), c("med", "q25", "q75")))
a  <- cm[, qs(mean_coef), by = tp][order(tp)]
b  <- cm[, qs(mean_coef), by = .(tp, grp)][order(grp, tp)]
fr <- cm[, .(neg = mean(mean_coef < 0, na.rm = TRUE), n = .N), by = .(tp, grp)]
fr_all <- cm[, .(neg = mean(mean_coef < 0, na.rm = TRUE), grp = "全体"), by = tp]

PAL <- c(抑制组 = "#1C6B66", 促进组 = "#C2703A")
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey35", size = 9.5, lineheight = 1.2),
        legend.position = "top", legend.title = element_blank())
XA <- scale_x_continuous(breaks = 0:8)

pA <- ggplot(a, aes(tp, med)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#7A4A78", alpha = .16) +
  geom_line(colour = "#7A4A78", linewidth = .9) +
  geom_point(colour = "#7A4A78", size = 2.1) +
  XA + labs(title = "A  全体：中位效应由零缓慢漂向正值",
       subtitle = sprintf(paste0("线为中位数，带为四分位区间；n = %d 站。",
                          "注意 tp>=5 时 rho 中位数已低于 0.03，\n系数近乎不可识别，",
                          "这段正漂应作存疑处理而非结论"), length(full)),
       x = "滞后 tp（每步 8 天）",
       y = expression(paste("S-map 系数  ", partialdiff, "SIF/", partialdiff, "VPD"))) + th

pB <- ggplot(b, aes(tp, med, colour = grp, fill = grp)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = .14, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
  XA + labs(title = "B  分组：抑制组越过零并转正，促进组只收缩不反号",
       subtitle = paste0("按 tp=0 的符号分组（抑制组 ", NG[grp == "抑制组", N],
                         " 站、促进组 ", NG[grp == "促进组", N],
                         " 站）；tp=0 处的分离是分组定义带来的，不是结果\n",
                         "有信息的是收缩速度不对称：抑制组在 tp≈3 越零，促进组到 tp=8 仍为正"),
       x = "滞后 tp（每步 8 天）", y = "S-map 系数中位数") + th

pC <- ggplot(fr, aes(tp, neg, colour = grp)) +
  geom_hline(yintercept = .5, linetype = "22", colour = "grey55", linewidth = .35) +
  geom_line(data = fr_all, aes(tp, neg), colour = "grey45", linewidth = .7) +
  geom_point(data = fr_all, aes(tp, neg), colour = "grey45", size = 1.8) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = PAL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  XA + labs(title = "C  负号占比：两组交汇在 30–45%，低于随机的 50%",
       subtitle = paste0("灰线为全体；虚线 50% = 方向完全随机。这是不依赖系数大小的看法。\n",
                         "两组并非收敛到 50%，而是双双落到 50% 以下——滞后不只抹平方向记忆，",
                         "还把整体推向正号"),
       x = "滞后 tp（每步 8 天）", y = "系数为负（抑制）的站点占比") + th

p <- pA / pB / pC + plot_annotation(
  title = "VPD → SIF 的 S-map 系数随滞后的变化",
  subtitle = "HCSIF 1000 m 缓冲区，2000–2022 年 5–9 月，8 天合成；系数为带符号的局部偏导数",
  theme = theme(plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(colour = "grey35")))
ggsave(file.path(OUT, "fig_smap_coef_tp.png"), p, width = 7.2, height = 10.5, dpi = 300)
log_msg("已保存 fig_smap_coef_tp.png")

fwrite(merge(b, fr, by = c("tp", "grp"))[order(grp, tp)], file.path(OUT, "smap_coef_by_tp.csv"))
log_msg("\n-- 全体 --")
print(a[, .(tp, 中位 = round(med, 4), IQR = sprintf("%.3f–%.3f", q25, q75))])
log_msg("\n-- 分组中位数与负号占比 --")
print(dcast(b, tp ~ grp, value.var = "med")[, lapply(.SD, function(x) round(x, 4))])
print(dcast(fr, tp ~ grp, value.var = "neg")[, lapply(.SD, function(x) round(x, 3))])

# 稳健性: 与 23 号严格判据的子集(剔除多次变号)对比, 看曲线是否被筛选改变
sg <- cm[, .(s = list(sign(mean_coef))), by = .(stat_id, grp)]
sg[, nflip := sapply(s, function(x) sum(diff(x) != 0))]
log_msg("\n-- 稳健性: 变号次数分布 --")
print(sg[, .N, by = nflip][order(nflip)])
st <- cm[stat_id %in% sg[nflip <= 1, stat_id]]
log_msg("严格判据子集(变号<=1 次) n = ", uniqueN(st$stat_id), " 站, 其分组中位数:")
print(dcast(st[, .(m = median(mean_coef, na.rm = TRUE)), by = .(tp, grp)],
            tp ~ grp, value.var = "m")[, lapply(.SD, function(x) round(x, 4))])
