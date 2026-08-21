#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 37_fig_survival_concept.R — 说明图: CCM 的 S-map 系数序列如何变成生存数据
#
# 左: 6 个真实站点的 mean_coef 随滞后的轨迹, 标出"首次变号"发生在哪一步
# 右: 同 6 站的泳道图(swimmer plot) —— 生存分析真正看到的就是这个:
#     一条从 tp0 开始的线段, 末端要么是事件(实心/空心标记), 要么是删失(箭头)
#
# 这张图用于方法部分, 也是给不熟悉生存分析的读者的入口。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
setorder(cm, meteo_stat, tp)
D <- readRDS(file.path(OUT, "surv_base.rds"))

# 每种 status 各取两个抑制组站点作例子(固定站号, 便于文中引用)
IDS <- c(54249, 54569, 53984, 54852, 59061, 58940)
SL  <- c("0" = "删失: 到 tp8 从未变号",
         "1" = "事件1: 持续翻转",
         "2" = "事件2: 瞬时翻转(竞争事件)")
CL  <- c("删失: 到 tp8 从未变号" = "#7F7F7F",
         "事件1: 持续翻转"       = "#1A7F37",
         "事件2: 瞬时翻转(竞争事件)" = "#B2182B")

meta <- D[stat_id %in% IDS, .(stat_id, time, status)]
meta[, lab_st := SL[as.character(status)]]
meta[, lab := sprintf("站 %d  (time=%d, status=%d)", stat_id, time, status)]
setorder(meta, status, stat_id)
meta[, lab := factor(lab, levels = lab)]

tr <- merge(cm[meteo_stat %in% IDS & tp <= 8, .(stat_id = meteo_stat, tp, mean_coef)],
            meta, by = "stat_id")
tr[, lab := factor(lab, levels = levels(meta$lab))]

th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), legend.key.size = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        strip.text = element_text(size = 8))

# ---- (a) 系数轨迹 ----------------------------------------------------------
pa <- ggplot(tr, aes(tp, mean_coef)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = .4) +
  geom_vline(data = meta[status > 0], aes(xintercept = time - .5, colour = lab_st),
             linetype = 2, linewidth = .5, show.legend = FALSE) +
  geom_line(colour = "grey45", linewidth = .5) +
  geom_point(aes(fill = mean_coef > 0), shape = 21, size = 2, stroke = .3,
             colour = "grey20") +
  scale_fill_manual(values = c("TRUE" = "#EF8A62", "FALSE" = "#67A9CF"),
                    labels = c("TRUE" = "系数为正(促进)", "FALSE" = "系数为负(抑制)")) +
  scale_colour_manual(values = CL) +
  scale_x_continuous(breaks = 0:8) +
  facet_wrap(~ lab, ncol = 2, scales = "free_y") +
  labs(title = "(a) 原始数据: 每站 9 个滞后上的 S-map 同期偏导数",
       subtitle = "虚线 = 首次变号发生的位置, 即生存分析里的事件时间 time",
       x = "滞后 tp(每步 8 天)", y = "mean_coef") + th

# ---- (b) 泳道图 ------------------------------------------------------------
sw <- copy(meta)
sw[, y := .I]
seg <- sw[, .(lab, y, x0 = 0, x1 = time, lab_st)]

pb <- ggplot(seg) +
  geom_segment(aes(x = x0, xend = x1, y = y, yend = y, colour = lab_st),
               linewidth = 1.6, lineend = "butt") +
  # 事件: 实心点; 删失: 右向箭头(表示"至少这么久, 之后不知道")
  geom_point(data = sw[status > 0], aes(time, y, colour = lab_st),
             size = 3.4, show.legend = FALSE) +
  geom_segment(data = sw[status == 0],
               aes(x = time, xend = time + .7, y = y, yend = y, colour = lab_st),
               arrow = arrow(length = unit(5, "pt"), type = "closed"),
               linewidth = 1.1, show.legend = FALSE) +
  geom_text(aes(x = x1 + ifelse(lab_st == SL[1], 1.5, .35), y = y,
                label = ifelse(lab_st == SL[1], "删失", paste0("事件 @ tp", x1))),
            hjust = 0, size = 2.7, colour = "grey25") +
  scale_colour_manual(values = CL) +
  scale_y_reverse(breaks = sw$y, labels = sw$lab) +
  scale_x_continuous(breaks = 0:8, limits = c(0, 12)) +
  labs(title = "(b) 生存分析看到的数据: (time, status) 二元因变量",
       subtitle = "线段长度 = 撑了多久; 末端标记 = 怎么结束的。删失者只贡献'至少撑到 tp8'这一信息",
       x = "滞后 tp(每步 8 天)", y = NULL) + th +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(OUT, "fig_survival_concept.png"), pa / pb + plot_layout(heights = c(1.5, 1)),
       width = 10, height = 9, dpi = 300)
cat("已写出", file.path(OUT, "fig_survival_concept.png"), "\n")
