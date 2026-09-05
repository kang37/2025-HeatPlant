#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 11_fig_smap_yearly_heatmap.R — 134(实际127)站二变量S-map方向逐年热力图
#
# 用户要求: 横轴=年份(+最后一列Overall=06主线的pooled结果), 纵轴=站点,
# 不用格子用实心圆圈, 圆圈大小表示数值大小; 做两版——
#   fig3: 颜色/大小 = S-map系数中位数(median_coef, 连续量, 也是06主线判方向用的量)
#   fig4: 颜色/大小 = 系数为正的时间点占比(frac_positive, 06主线判方向的另一半信息
#         ——06是"中位数+主导符号占比>=75%"两个条件一起用, 这里让占比单独成一张图,
#         能看出"符号一致但中位数小"和"中位数不小但符号来回翻"这两种不同的不稳定)
#
# 数据来自 10_smap_yearly_134.R 的输出, 只画 n_pts>=5 的站x年格子(样本太少的
# 格子留空, 不画点), Overall列直接用06已发布的 median_coef_2v/frac_positive_2v。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(sysfonts); library(scales)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
font_add("cjk", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)

OUT_DIR <- file.path(PROJ, "data_proc/smap_bivar_134")
yr  <- readRDS(file.path(OUT_DIR, "smap_yearly_134.rds"))$yearly[valid == TRUE]
st  <- readRDS(file.path(OUT_DIR, "smap_bivar_134.rds"))[!is.na(direction_2v)]

YEARS <- 2000:2022
dir_colors <- c(Promote = "#C2703A", Inhibit = "#2E7D74", Ambiguous = "grey65")

# ---- 站点排序: 先按 Overall 方向分组(Promote/Ambiguous/Inhibit), 组内按
#      median_coef_2v 降序 ----------------------------------------------------
st[, direction_2v := factor(direction_2v, levels = c("Promote", "Ambiguous", "Inhibit"))]
setorder(st, direction_2v, -median_coef_2v)
st_order <- st$meteo_stat                       # 从上到下的顺序(第一个排最上面)
st[, meteo_stat_f := factor(meteo_stat, levels = rev(st_order))]

n_st <- nrow(st)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
log_msg("站数: ", n_st, "; 年份: ", min(YEARS), "-", max(YEARS))

# ---- 拼长表: 逐年格子 + Overall 列 -----------------------------------------
yearly_df <- yr[meteo_stat %in% st_order,
                .(meteo_stat, year_lab = as.character(year), median_coef, frac_positive, n_pts)]
overall_df <- st[, .(meteo_stat, year_lab = "Overall",
                      median_coef = median_coef_2v, frac_positive = frac_positive_2v,
                      n_pts = n_coef_pts)]
plot_df <- rbind(yearly_df, overall_df, fill = TRUE)
plot_df <- merge(plot_df, st[, .(meteo_stat, meteo_stat_f, direction_2v)], by = "meteo_stat")
plot_df[, year_lab := factor(year_lab, levels = c(as.character(YEARS), "Overall"))]

fig_height <- max(6, 1.2 + n_st * 0.085)

theme_heat <- theme_minimal(base_size = 11, base_family = "cjk") +
  theme(panel.grid = element_line(color = "grey93", linewidth = 0.2),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9.5),
        legend.position = "right")

gap_x <- length(YEARS) + 0.5   # Overall 列前面画一条竖线分隔

# ===========================================================================
# fig3: 颜色/大小 = median_coef (S-map系数中位数, 因果方向的强度+符号)
# ===========================================================================
lim3 <- 0.3
plot_df[, median_coef_clip := pmin(pmax(median_coef, -lim3), lim3)]   # 手动裁剪极值, 避免个别站拉爆色阶
p3 <- ggplot(plot_df, aes(x = year_lab, y = meteo_stat_f)) +
  geom_vline(xintercept = gap_x, linewidth = 0.4, color = "grey60") +
  geom_point(aes(color = median_coef_clip, size = abs(median_coef_clip)), shape = 16) +
  scale_color_gradient2(low = "#2E7D74", mid = "grey92", high = "#C2703A", midpoint = 0,
                         limits = c(-lim3, lim3),
                         name = expression(paste("系数中位数 ", partialdiff, "SIF/", partialdiff, "VPD"))) +
  scale_size_continuous(range = c(0.3, 3.2), limits = c(0, lim3),
                         name = "|系数中位数|") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "134站二变量S-map方向逐年变化(系数中位数版)",
       subtitle = "每格=某站某年S-map系数的中位数(n_pts>=5); 最后一列Overall=06主线pooled全部年份的结果; 站按Overall方向x强度排序",
       x = NULL, y = NULL) +
  theme_heat

ggsave(file.path(OUT_DIR, "fig3_heatmap_yearly_median.png"), p3,
       width = 11, height = fig_height, dpi = 300, limitsize = FALSE)
log_msg("已保存 fig3_heatmap_yearly_median.png (height=", round(fig_height, 1), "in)")

# ===========================================================================
# fig4: 颜色/大小 = frac_positive (系数为正的时间点占比, 主导符号法)
# ===========================================================================
plot_df[, dominance := abs(frac_positive - 0.5) * 2]   # 0=完全对半分, 1=符号完全一致

p4 <- ggplot(plot_df, aes(x = year_lab, y = meteo_stat_f)) +
  geom_vline(xintercept = gap_x, linewidth = 0.4, color = "grey60") +
  geom_point(aes(color = frac_positive, size = dominance), shape = 16) +
  scale_color_gradient2(low = "#2E7D74", mid = "grey92", high = "#C2703A", midpoint = 0.5,
                         limits = c(0, 1), name = "系数为正的\n时间点占比") +
  scale_size_continuous(range = c(0.3, 3.2), limits = c(0, 1), name = "符号一致度\n|占比-0.5|x2") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "134站二变量S-map方向逐年变化(主导符号占比版)",
       subtitle = "每格=某站某年S-map系数为正的时间点占比(n_pts>=5); >=75%记Promote/<=25%记Inhibit(06主线阈值); 最后一列Overall=06主线pooled结果; 站排序同上图",
       x = NULL, y = NULL) +
  theme_heat

ggsave(file.path(OUT_DIR, "fig4_heatmap_yearly_signshare.png"), p4,
       width = 11, height = fig_height, dpi = 300, limitsize = FALSE)
log_msg("已保存 fig4_heatmap_yearly_signshare.png (height=", round(fig_height, 1), "in)")

# ---- 控制台快速摘要, 方便写文档 --------------------------------------------
m <- merge(yr, st[, .(meteo_stat, direction_2v, median_coef_2v)], by = "meteo_stat")
m[, sign_year := sign(median_coef)]; m[, sign_overall := sign(median_coef_2v)]
m[, flip := sign_year != sign_overall & sign_year != 0 & sign_overall != 0]
cat("\n=== 摘要(写文档用) ===\n")
cat("站x年格子(n_pts>=5):", nrow(m), "; median符号与Overall相反的格子占比:",
    round(mean(m$flip) * 100, 1), "%\n")
byst <- m[, .(any_flip = any(flip)), by = .(meteo_stat, direction_2v)]
cat("Promote/Inhibit站里, 至少1年median符号与Overall相反的站数:",
    byst[direction_2v != "Ambiguous" & any_flip == TRUE, .N], "/",
    byst[direction_2v != "Ambiguous", .N], "\n")
