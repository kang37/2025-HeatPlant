#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 36_fig_survival.R — 35 号生存分析的四联图
#
#   (a) Kaplan-Meier: 到第 t 步仍未变号的比例, 两组对比 + log-rank
#   (b) 竞争风险累积发生率(Aalen-Johansen), 并叠上 1-KM 以显示"把竞争事件
#       当删失"会高估多少 —— 这是本轮方法改动最直观的一张图
#   (c) Cox 森林图: 18 协变量的 HR 与 Conley 200km 区间, 分组
#   (d) cause-specific 与 Fine-Gray 系数对照: 偏离 45 度线的变量, 其作用
#       有相当部分是通过竞争事件(瞬时翻转)实现的
#
# 依赖 35_survival_deep.R 的输出。先跑 35 再跑本脚本。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
  library(showtext); library(patchwork)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; CUT <- 200; MAXTP <- 8L
log_msg <- function(...) cat(..., "\n", sep = "")

D <- readRDS(file.path(OUT, "surv_base.rds"))
GL <- c(inhibit_first = "抑制组(tp0 为负)", promote_first = "促进组(tp0 为正)")
CL <- c("抑制组(tp0 为负)" = "#2166AC", "促进组(tp0 为正)" = "#B2182B")
CN <- c(imperv = "不透水面", grass = "草地", water = "水体",
        ever_needle = "常绿针叶林", deci_needle = "落叶针叶林",
        ever_broad = "常绿阔叶林", deci_broad = "落叶阔叶林", mixedleaf = "混交林",
        tavg = "平均气温", rh = "相对湿度", cloud = "云量", precip = "降水量",
        rsds_mean = "辐射均值", rsds_sd = "辐射年际变率", elev = "海拔",
        ntl = "夜间灯光", invest = "绿地投资", urban_rate = "城镇化率")
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(10, "pt"))

# ---- (a) Kaplan-Meier ------------------------------------------------------
km <- survfit(Surv(time, ev_any) ~ start_dir, data = D)
kd <- data.table(time = km$time, surv = km$surv, lo = km$lower, hi = km$upper,
                 grp = GL[sub("start_dir=", "", rep(names(km$strata), km$strata))])
# 曲线要从 (0, 1) 起画, KM 对象不含这一点
kd <- rbind(data.table(time = 0, surv = 1, lo = 1, hi = 1, grp = unique(kd$grp)), kd)
setorder(kd, grp, time)
sd0 <- survdiff(Surv(time, ev_any) ~ start_dir, data = D)
p_lr <- pchisq(sd0$chisq, length(sd0$n) - 1, lower.tail = FALSE)

# 空间稳健的 log-rank: log-rank 就是 Cox 中二值协变量的得分检验,
# 把方差换成 Conley 即可。与 35 号脚本同法, 此处重算以免副标题写死数字。
dist_km <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2) }
m_g <- coxph(Surv(time, ev_any) ~ start_dir, data = D, ties = "efron")
Dfb <- as.matrix(residuals(m_g, type = "dfbeta"))
DM <- dist_km(D$longitude, D$latitude); K <- DM; K[] <- pmax(0, 1 - DM / CUT)
se_c <- sqrt(diag(t(Dfb) %*% K %*% Dfb))
p_conley_lr <- 2 * pnorm(-abs(coef(m_g) / se_c))

pa <- ggplot(kd, aes(time, surv, colour = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_step(linewidth = .8) +
  geom_hline(yintercept = .5, linetype = 3, colour = "grey50") +
  scale_colour_manual(values = CL) + scale_fill_manual(values = CL) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "(a) 尚未变号的站点比例(Kaplan-Meier)",
       subtitle = sprintf("事件 = S-map 系数首次变号 | log-rank χ²=%.1f, p=%.2g | 空间稳健(Conley %dkm) p=%.2g",
                          sd0$chisq, p_lr, CUT, p_conley_lr),
       x = "滞后 tp(每步 8 天)", y = "S(t)") + th

# ---- (b) 竞争风险 CIF vs 1-KM ---------------------------------------------
aj <- survfit(Surv(time, status_f) ~ start_dir, data = D)
gvec <- sub("start_dir=", "", rep(names(aj$strata), aj$strata))
cif <- rbindlist(lapply(c("sustained", "transient"), function(s)
  data.table(time = aj$time, cif = aj$pstate[, s], state = s, grp = GL[gvec])))
kmc <- survfit(Surv(time, status == 1L) ~ start_dir, data = D)
gv2 <- sub("start_dir=", "", rep(names(kmc$strata), kmc$strata))
naive <- data.table(time = kmc$time, cif = 1 - kmc$surv, grp = GL[gv2])
z0 <- CJ(time = 0, grp = unique(cif$grp), state = c("sustained", "transient"))[, cif := 0]
cif <- rbind(z0, cif, use.names = TRUE)
naive <- rbind(data.table(time = 0, cif = 0, grp = unique(naive$grp)), naive)
setorder(cif, grp, state, time); setorder(naive, grp, time)
SL <- c(sustained = "持续翻转(事件1)", transient = "瞬时翻转(事件2)")
cif[, state := SL[state]]

pb <- ggplot(cif, aes(time, cif, colour = grp, linetype = state)) +
  geom_step(linewidth = .8) +
  geom_step(data = naive, aes(time, cif, colour = grp), linetype = 3,
            linewidth = .6, inherit.aes = FALSE, alpha = .8) +
  scale_colour_manual(values = CL) +
  scale_linetype_manual(values = c("持续翻转(事件1)" = 1, "瞬时翻转(事件2)" = 2)) +
  scale_x_continuous(breaks = 0:8) + scale_y_continuous(labels = scales::percent) +
  labs(title = "(b) 累积发生率(Aalen-Johansen)与 1−KM 的偏差",
       subtitle = sprintf("点线 = 1−KM(把竞争事件当删失): 抑制组 tp8 处高估至 %.0f%% vs 真实 %.0f%%",
                          100 * naive[grp == GL[1] & time == MAXTP, cif][1],
                          100 * cif[grp == GL[1] & state == SL[1] & time == MAXTP, cif][1]),
       x = "滞后 tp(每步 8 天)", y = "累积发生率") + th +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2))

# ---- (c) Cox 森林图 --------------------------------------------------------
cf <- fread(file.path(OUT, "surv_cox_coef.csv"))
cf <- cf[tag %in% c("cox_any_inhibit_first", "cox_any_promote_first")]
cf[, grp := GL[sub("cox_any_", "", tag)]]
cf[, lab := CN[var]]
ord <- cf[grp == GL[1]][order(HR), var]
cf[, lab := factor(CN[var], levels = CN[ord])]
cf[, sig := fifelse(p_conley < .05, "p<0.05", "n.s.")]

pc <- ggplot(cf, aes(HR, lab, colour = grp, shape = sig)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 position = position_dodge(width = .6), linewidth = .5) +
  geom_point(position = position_dodge(width = .6), size = 1.9) +
  scale_colour_manual(values = CL) +
  scale_shape_manual(values = c("p<0.05" = 16, "n.s." = 1)) +
  scale_x_log10() +
  labs(title = "(c) Cox 风险比(事件 = 首次变号)",
       subtitle = "每 1 个标准差; 区间为 Conley 200km 空间稳健; HR>1 = 变号更早",
       x = "风险比 HR(对数轴)", y = NULL) + th

# ---- (d) cause-specific vs Fine-Gray --------------------------------------
cp <- fread(file.path(OUT, "surv_competing_coef.csv"))
cs <- cp[grepl("^cause_specific_ev1_", tag), .(var, cs = beta,
          grp = GL[sub("cause_specific_ev1_", "", tag)], p_cs = p_conley)]
fg <- cp[grepl("^fine_gray_ev1_", tag), .(var, fg = beta,
          grp = GL[sub("fine_gray_ev1_", "", tag)])]
mg <- merge(cs, fg, by = c("var", "grp"))
mg[, lab := fifelse(p_cs < .05 | abs(cs - fg) > .15, CN[var], NA_character_)]
rng <- range(c(mg$cs, mg$fg))

pd <- ggplot(mg, aes(cs, fg, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 1.9, alpha = .85) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.6, na.rm = TRUE,
                           show.legend = FALSE, min.segment.length = 0,
                           segment.colour = "grey70", max.overlaps = 20) +
  scale_colour_manual(values = CL) + coord_equal(xlim = rng, ylim = rng) +
  labs(title = "(d) 病因别风险 vs 子分布风险(事件1)",
       subtitle = "偏离虚线 = 该变量的作用有一部分经由竞争事件(瞬时翻转)实现",
       x = "cause-specific β", y = "Fine-Gray β") + th

fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1.25))
ggsave(file.path(OUT, "fig_survival_4panel.png"), fig,
       width = 12, height = 10, dpi = 300)
log_msg("已写出 ", file.path(OUT, "fig_survival_4panel.png"))
