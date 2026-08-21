#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 38_fig_competing_explain.R — 拆解四联图的 (d) 面板: 两种风险到底差在哪
#
# (a) 两种风险集的大小对比。病因别风险集里, 发生过竞争事件的站被移走;
#     子分布风险集把它们留着。tp6 处 62 vs 238, 差近 4 倍 —— 分母不同,
#     算出来的"风险"当然不同, 这就是两套系数分岔的全部来源。
# (b) 抑制组逐变量的两种系数对照(全标签版, 原图太挤)
# (c) 验证: 偏离对角线的幅度 gap = beta_FG - beta_CS1 应当由"该变量对竞争
#     事件的作用 beta_CS2"驱动。若二者强负相关, 说明 (d) 的读法成立。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"

D  <- readRDS(file.path(OUT, "surv_base.rds"))
cp <- fread(file.path(OUT, "surv_competing_coef.csv"))
CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率",
        elev="海拔", ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
GL <- c(inhibit_first = "抑制组", promote_first = "促进组")
CL <- c(抑制组 = "#2166AC", 促进组 = "#B2182B")
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), legend.key.size = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

# ---- (a) 两种风险集 --------------------------------------------------------
d <- D[start_dir == "inhibit_first"]
rs <- rbindlist(lapply(1:8, function(t) data.table(
  tp = t,
  病因别 = d[time >= t, .N],
  已瞬时翻转 = d[status == 2 & time < t, .N])))
rs[, 子分布 := 病因别 + 已瞬时翻转]
rl <- melt(rs[, .(tp, `病因别风险集(还没变号的站)` = 病因别,
                  `子分布风险集(再加上已瞬时翻转的站)` = 子分布)],
           id.vars = "tp", variable.name = "type", value.name = "n")

pa <- ggplot(rl, aes(tp, n, colour = type, shape = type)) +
  geom_line(linewidth = .8) + geom_point(size = 2.4) +
  geom_segment(data = rs, aes(x = tp, xend = tp, y = 病因别, yend = 子分布),
               inherit.aes = FALSE, colour = "grey60", linetype = 3) +
  annotate("text", x = 6, y = 150, hjust = 0, size = 2.8, colour = "grey25",
           label = "tp6: 62 vs 238\n差的 176 站 = 已发生\n竞争事件、但被子分布\n风险集强行留下的站") +
  scale_colour_manual(values = c("#1A7F37", "#8C6BB1")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "(a) 两种风险的分母不一样",
       subtitle = "抑制组 406 站。删失只发生在 tp8, 故 Fine-Gray 权重在 tp1–7 全为 1",
       x = "滞后 tp", y = "风险集站数") + th +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

# ---- (b) 抑制组两种系数对照(全标签) ---------------------------------------
get <- function(g) {
  cs1 <- cp[tag == paste0("cause_specific_ev1_", g), .(var, cs1 = beta, p1 = p_conley)]
  cs2 <- cp[tag == paste0("cause_specific_ev2_", g), .(var, cs2 = beta, p2 = p_conley)]
  fg  <- cp[tag == paste0("fine_gray_ev1_",  g),     .(var, fg = beta,  pf = p_conley)]
  m <- Reduce(function(a, b) merge(a, b, by = "var"), list(cs1, cs2, fg))
  m[, `:=`(gap = fg - cs1, grp = GL[g], cn = CN[var])][]
}
mi <- get("inhibit_first"); mp <- get("promote_first"); mm <- rbind(mi, mp)
rr <- range(c(mi$cs1, mi$fg))

pb <- ggplot(mi, aes(cs1, fg)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey88") +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_segment(aes(xend = cs1, yend = cs1), colour = "grey75", linewidth = .3) +
  geom_point(aes(fill = p1 < .05 | pf < .05), shape = 21, size = 2.6,
             colour = "grey20", stroke = .3) +
  ggrepel::geom_text_repel(aes(label = cn), size = 2.6, max.overlaps = 30,
                           min.segment.length = 0, segment.colour = "grey70") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "white"),
                    labels = c(`TRUE` = "任一口径 p<0.05", `FALSE` = "均不显著")) +
  coord_equal(xlim = rr, ylim = rr) +
  labs(title = "(b) 抑制组: 病因别 β 与 Fine-Gray β",
       subtitle = "灰线 = 到对角线的垂距, 即偏离量 gap",
       x = "病因别 β(机制: 谁更快持续翻转)",
       y = "Fine-Gray β(占比: 谁最终更可能持续翻转)") + th

# ---- (c) 偏离量由什么驱动 --------------------------------------------------
rs_lab <- mm[, .(r = cor(gap, cs2)), by = grp]
rs_lab[, lab := sprintf("%s: r = %.2f", grp, r)]

pc <- ggplot(mm, aes(cs2, gap, colour = grp)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6, linetype = 2) +
  geom_point(size = 2.2, alpha = .85) +
  ggrepel::geom_text_repel(data = mm[abs(gap) > .12], aes(label = cn), size = 2.5,
                           show.legend = FALSE, max.overlaps = 20,
                           min.segment.length = 0, segment.colour = "grey70") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = c(1.6, 3.1), size = 2.9,
           colour = CL[rs_lab$grp], label = rs_lab$lab) +
  scale_colour_manual(values = CL) +
  labs(title = "(c) 偏离对角线的幅度, 由该变量对竞争事件的作用决定",
       subtitle = "横轴越大 = 越促进'瞬时翻转' -> 越多站被竞争事件用掉 -> Fine-Gray β 越被压低",
       x = "对竞争事件的病因别 β(事件2: 瞬时翻转)",
       y = "gap = Fine-Gray β − 病因别 β") + th

ggsave(file.path(OUT, "fig_competing_explain.png"),
       pa / (pb | pc) + plot_layout(heights = c(1, 1.15)),
       width = 11, height = 9.5, dpi = 300)
cat("已写出", file.path(OUT, "fig_competing_explain.png"), "\n")
