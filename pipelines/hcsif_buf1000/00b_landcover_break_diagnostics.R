#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00b_landcover_break_diagnostics.R — 地表覆盖断点筛选的两张诊断图
#
# 图1(方法指标): magnitude x p值 散点, 全部924站做背景, 36个被排除站高亮+
#                标注站号, 虚线标出判据阈值(magnitude>0.08 且 p<0.01)。
# 图2(原始数据拼接): 36个被排除站的森林占比年际序列小多图, 按 magnitude 排序,
#                红色虚线=检出的断点年。
#
# 依赖 00_landcover_break_buf1000.R 已经跑过、landcover_break_stations.csv 已存在。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
  library(showtext); library(sysfonts)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 160)

GLC_DIR <- "data_raw/covariates_1km"; YEARS <- 2000:2022
FOREST_COLS <- sprintf("LC%02d", 4:13)
OUT_DIR <- "data_proc/output_hcsif_buf1000"
MAG_THRESH <- 0.08; P_THRESH <- 0.01
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

res <- fread(file.path(OUT_DIR, "landcover_break_stations.csv"))
flagged <- res[has_break == TRUE][order(-magnitude)]
log_msg("被排除站点: ", nrow(flagged))

# ===========================================================================
# 图1: 方法指标——magnitude x p值，全站背景 + 36站高亮
# ===========================================================================
res[, direction := fifelse(forest_frac_2022 > forest_frac_2000, "增(疑似造林)", "减(疑似砍伐)")]
res[, excluded := has_break == TRUE]

p1 <- ggplot(res[!is.na(pval)], aes(magnitude, pval)) +
  geom_point(data = res[excluded == FALSE | is.na(excluded)],
             color = "grey75", size = 1.2, alpha = 0.6) +
  geom_point(data = res[excluded == TRUE], aes(color = direction, size = direction)) +
  scale_size_manual(values = c("增(疑似造林)" = 5, "减(疑似砍伐)" = 2.6), guide = "none") +
  geom_text_repel(data = res[excluded == TRUE], aes(label = stat_id),
                   size = 2.6, max.overlaps = 40, segment.size = 0.2, seed = 1) +
  geom_vline(xintercept = MAG_THRESH, linetype = "dashed", color = "black") +
  geom_hline(yintercept = P_THRESH, linetype = "dashed", color = "black") +
  scale_y_log10() +
  scale_color_manual(values = c("增(疑似造林)" = "#1b9e77", "减(疑似砍伐)" = "#d95f02")) +
  labs(title = "地表覆盖断点筛选：方法判据指标(全部924站)",
       subtitle = sprintf("虚线=判据阈值(magnitude>%.2f 且 p<%.2f)；灰点=未排除(888站)；彩点=排除(%d站)",
                           MAG_THRESH, P_THRESH, nrow(flagged)),
       x = "magnitude(断点前后森林占比差, 绝对值)", y = "置换检验 p 值(对数轴)",
       color = "断点方向") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "heiti"),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "landcover_break_diag_metrics.png"), p1,
       width = 11, height = 8.5, dpi = 160)
log_msg("已写出 landcover_break_diag_metrics.png")

# ===========================================================================
# 图2: 原始数据拼接——36站森林占比年际序列小多图
# ===========================================================================
log_msg("读取 GLC 逐年数据...")
glc <- rbindlist(lapply(YEARS, function(y) {
  f <- file.path(GLC_DIR, sprintf("glc_station_%d.csv", y))
  d <- fread(f, select = c("stat_id", "year", FOREST_COLS))
  d[, forest_frac := rowSums(.SD), .SDcols = FOREST_COLS]
  d[, .(stat_id, year, forest_frac)]
}))

pd <- glc[stat_id %in% flagged$stat_id]
pd <- merge(pd, flagged[, .(stat_id, break_year, magnitude)], by = "stat_id")
# 按 magnitude 降序排面板顺序
pd[, stat_id_f := factor(stat_id, levels = flagged$stat_id)]
pd[, panel_lab := sprintf("%s (Δ=%.2f)", stat_id, magnitude)]
lab_order <- pd[order(-magnitude), unique(panel_lab)]
pd[, panel_lab := factor(panel_lab, levels = lab_order)]

p2 <- ggplot(pd, aes(year, forest_frac)) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.6) +
  geom_vline(aes(xintercept = break_year), color = "red", linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~panel_lab, scales = "free_y", ncol = 6) +
  labs(title = "地表覆盖断点筛选：36 个被排除站点的原始森林占比序列",
       subtitle = "按 Δ(断点幅度) 降序排列；红色虚线=检出的断点年；纵轴各站自适应尺度",
       x = "年份", y = "缓冲区森林占比") +
  theme_minimal(base_size = 10) +
  theme(text = element_text(family = "heiti"),
        strip.text = element_text(size = 7.5, face = "bold"),
        axis.text = element_text(size = 6.5),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "landcover_break_diag_rawseries_all36.png"), p2,
       width = 16, height = 13, dpi = 160)
log_msg("已写出 landcover_break_diag_rawseries_all36.png")
log_msg("完成")
