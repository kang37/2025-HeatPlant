#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 34_zone_forest.R — 分气候区的得分回归(OLS + Conley 200km)与森林图
#
# 与全体模型同法: OLS 估系数, Conley 空间 HAC(Bartlett 核, 200 km)估方差。
#
# 分区后样本大幅缩减, 18 个自变量放不下, 改用跨三维度的 5 变量精简集。
# 各格样本量差异很大(C 区抑制组 93 站, D 区促进组仅 17 站), 故:
#   - n < 25 的格不拟合;
#   - 每自变量观测数 < 10 的格照常拟合但标注"样本偏小", 由读者自行折扣;
#   - 同时给出"全体"(不分区)估计作为参照行, 这是森林图的惯例。
#
# 注意: 分区后站点在空间上更聚集, 同样的 200 km 截断覆盖的邻居比例更高,
# 故各区的有效独立样本量一并报出。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(grid) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; M <- 7L; CUT <- 200; MIN_N <- 25   # 该格站数下限
log_msg <- function(...) cat(..., "\n", sep = "")

env <- new.env(); L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
pat <- env$pat; dt <- env$dt; PRED <- env$PRED; winz <- env$winz
pat[, score := fifelse(start_dir == "inhibit_first",
                       fifelse(event == 1, (M + 2) - event_tp, 1),
                       fifelse(event == 1, as.numeric(event_tp), M + 1))]
st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
d0 <- merge(merge(pat[, .(stat_id, start_dir, score)],
                  dt[, c("stat_id", "koppen_group", PRED), with = FALSE], by = "stat_id"),
            st[, .(stat_id, longitude, latitude)], by = "stat_id")

# 与全体模型完全相同的 18 变量集。分区后每自变量观测数降至 3.5-4.9,
# 远低于常用的 10 的经验下限, 估计不稳、区间偏宽, 结果须据此折扣。
# 但保持设定一致才能与全体模型直接对照, 且 5 变量集经检验存在遗漏变量偏误
# (C 区抑制组 urban_rate 由 p=0.47 变为 p=0.0008)。
S5 <- PRED
CN <- c(imperv = "不透水面", grass = "草地", water = "水体",
        ever_needle = "常绿针叶林", deci_needle = "落叶针叶林",
        ever_broad = "常绿阔叶林", deci_broad = "落叶阔叶林", mixedleaf = "混交林",
        tavg = "平均气温", rh = "相对湿度", cloud = "云量", precip = "降水量",
        rsds_mean = "辐射均值", rsds_sd = "辐射年际变率", elev = "海拔",
        ntl = "夜间灯光", invest = "绿地投资", urban_rate = "城镇化率")

conley_se <- function(m, lon, lat, cut) {
  X <- model.matrix(m); u <- residuals(m); n <- nrow(X)
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  D <- sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2)
  br <- solve(crossprod(X))
  se <- sqrt(diag(br %*% crossprod(X, (pmax(0, 1 - D / cut) * outer(u, u)) %*% X) %*% br) *
             (n / (n - ncol(X))))
  list(se = se, eff_n = n / max(1, mean(rowSums(D <= cut)) - 1))
}

fit_cell <- function(d, zone, g) {
  d <- d[complete.cases(d[, c("score", S5), with = FALSE])]
  V <- S5[sapply(S5, function(v) sd(d[[v]], na.rm = TRUE) > 0)]   # 子集内无变异的列剔除
  # 需留出足够残差自由度, 否则 Conley 的小样本校正 n/(n-k) 会爆掉
  if (nrow(d) < MIN_N || uniqueN(d$score) < 3 || nrow(d) < length(V) + 10) return(NULL)
  z <- copy(d)[, c("score", V), with = FALSE]
  for (v in V) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  m <- lm(as.formula(paste("score ~", paste(V, collapse = "+"))), z)
  cs <- conley_se(m, d$longitude, d$latitude, CUT)
  b <- coef(m)[V]; se <- cs$se[V]
  data.table(zone = zone, grp = g, var = V, beta = b, se = se,
             p = 2 * pnorm(-abs(b / se)), n = nrow(d),
             obs_per_var = nrow(d) / length(V), eff_n = cs$eff_n,
             r2 = summary(m)$r.squared)
}

res <- list()
for (g in c("inhibit_first", "promote_first")) {
  res[[paste0("ALL_", g)]] <- fit_cell(d0[start_dir == g], "全体", g)
  for (kg in c("B", "C", "D"))
    res[[paste0(kg, "_", g)]] <- fit_cell(d0[start_dir == g & koppen_group == kg], kg, g)
}
r <- rbindlist(res)
ZL <- c(全体 = "全体（不分区）", B = "B 干旱带", C = "C 温带季风", D = "D 温带大陆性")



# 未拟合的格补空行, 使面板照常出现并标注原因(否则读者会误以为漏画)
cells <- CJ(zone = names(ZL), grp = c("inhibit_first", "promote_first"), unique = TRUE)
have <- unique(r[, .(zone, grp)])
miss <- cells[!have, on = .(zone, grp)]
if (nrow(miss)) {
  nz <- d0[complete.cases(d0[, c("score", S5), with = FALSE])][
    , .N, by = .(zone = koppen_group, grp = start_dir)]
  miss <- merge(miss, nz, by = c("zone", "grp"), all.x = TRUE)
  miss[is.na(N), N := 0L]
  r <- rbind(r, miss[, .(zone, grp, var = S5[1], beta = NA_real_, se = NA_real_,
                         p = NA_real_, n = N, obs_per_var = NA_real_, eff_n = NA_real_,
                         r2 = NA_real_)], fill = TRUE)
}
r[, `:=`(zone_cn = factor(ZL[zone], levels = ZL),
         cn = factor(CN[var], levels = rev(CN)),
         grp_cn = factor(fifelse(grp == "inhibit_first",
                                 "抑制组：得分高 = 越早脱离抑制",
                                 "促进组：得分高 = 促进维持越久"),
                levels = c("抑制组：得分高 = 越早脱离抑制", "促进组：得分高 = 促进维持越久")),
         lo = beta - 1.96 * se, hi = beta + 1.96 * se,
         small = obs_per_var < 10, sig = p < .05)]
log_msg("=== 各格拟合概况 ===")
print(unique(r[, .(气候区 = zone, 组 = fifelse(grp == "inhibit_first", "抑制", "促进"),
                   n, 每自变量 = round(obs_per_var, 1), 有效独立样本 = round(eff_n, 1),
                   R2 = round(r2, 3), 样本偏小 = fifelse(small, "是", ""))]))
log_msg("\n=== 显著结果 ===")
print(r[sig == TRUE][order(grp_cn, cn), .(组 = grp_cn, 气候区 = zone, 变量 = cn,
        系数 = round(beta, 3), CI = sprintf("[%.2f, %.2f]", lo, hi),
        p = signif(p, 3), 样本偏小 = fifelse(small, "是", ""))])
fwrite(r, file.path(OUT, sprintf("zone_conley%d_tp%d.csv", CUT, M)))

lab <- unique(r[!is.na(beta), .(zone_cn, grp_cn,
                 txt = sprintf("n = %d，R² = %.2f%s", n, r2,
                               fifelse(small, "，样本偏小", "")))])
gap <- unique(r[is.na(beta), .(zone_cn, grp_cn, txt = sprintf("样本不足（n = %d）", n))])

# B 区抑制组每自变量仅 2.1 个观测, 区间宽达 ±10, 会把其余面板压成一条线。
# 故显示范围截到 ±3, 超界的区间末端画箭头示意, 数据本身不删。
XL <- 3
r[, `:=`(lo_c = pmax(lo, -XL), hi_c = pmin(hi, XL),
         out_lo = lo < -XL, out_hi = hi > XL)]

PAL <- c("全体（不分区）" = "#3A3A3A", "B 干旱带" = "#C2703A",
         "C 温带季风" = "#2E7D74", "D 温带大陆性" = "#4B6BA8")
p <- ggplot(r[!is.na(beta)], aes(beta, cn, colour = zone_cn)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .35) +
  geom_errorbar(aes(xmin = lo_c, xmax = hi_c), orientation = "y", width = 0, linewidth = .7) +
  geom_segment(data = r[out_hi == TRUE], aes(x = XL - .25, xend = XL, y = cn, yend = cn),
               arrow = arrow(length = unit(.055, "in"), type = "closed"), linewidth = .5) +
  geom_segment(data = r[out_lo == TRUE], aes(x = -XL + .25, xend = -XL, y = cn, yend = cn),
               arrow = arrow(length = unit(.055, "in"), type = "closed"), linewidth = .5) +
  geom_point(aes(shape = sig, fill = zone_cn, size = small), stroke = .75) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.5, size = 2.9, colour = "grey40") +
  geom_text(data = gap, aes(x = 0, y = 3, label = txt), inherit.aes = FALSE,
            size = 3.2, colour = "grey55") +
  facet_grid(grp_cn ~ zone_cn) +
  coord_cartesian(xlim = c(-XL, XL)) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.7, `FALSE` = 2.5), guide = "none") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(title = "分气候区的转变时间得分回归",
       subtitle = paste0("OLS 系数 + Conley 空间 HAC 95% 置信区间（截断 200 km）；tp = 0–7（末步截断）\n",
                         "与全体模型相同的 18 变量集；实心 = p<0.05；",
                         "小点 = 每自变量观测数 < 10，结果需折扣；箭头 = 区间超出显示范围"),
       x = "标准化回归系数", y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = .3),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(.9, "lines"),
        strip.text.x = element_text(face = "bold", size = 10.5),
        strip.text.y = element_text(face = "bold", size = 9.5, angle = 0, hjust = 0),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9, lineheight = 1.2))
ggsave(file.path(OUT, sprintf("fig_zone_forest_tp%d.png", M)), p, width = 13, height = 11, dpi = 300)
log_msg("\n已输出 fig_zone_forest_tp", M, ".png")
