#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 31_conley_cutoff.R — Conley 空间 HAC 的截断距离该怎么定
#
# 三条依据, 本脚本把三条都算出来:
#   (1) 空间相关的实际作用距离: 残差相关图(按距离分带算 Moran's I),
#       看到哪个距离带上相关性落回零, 截断距离就该覆盖到那里。
#   (2) 标准误对截断距离的敏感曲线: 好的截断应落在曲线的平台段,
#       曲线还在爬说明截得太短, 开始剧烈抖动说明截得太长。
#   (3) 有效独立样本量: n / 平均邻居数。这个数太小(经验上 <30-40)
#       时 meat 矩阵不稳, 标准误会被人为压低。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; M <- 8L; set.seed(42)
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
                  dt[, c("stat_id", PRED), with = FALSE], by = "stat_id"),
            st[, .(stat_id, longitude, latitude)], by = "stat_id")

distmat <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2)
}
conley_se <- function(m, D, cut) {
  X <- model.matrix(m); u <- residuals(m); n <- nrow(X)
  K <- pmax(0, 1 - D / cut)
  br <- solve(crossprod(X))
  sqrt(diag(br %*% crossprod(X, (K * outer(u, u)) %*% X) %*% br) * (n / (n - ncol(X))))
}
# 距离分带的 Moran's I
band_I <- function(r, D, lo, hi, nperm = 499) {
  W <- (D > lo & D <= hi) * 1; diag(W) <- 0
  if (sum(W) < 20) return(c(NA, NA, sum(W)/2))
  z <- r - mean(r); n <- length(r)
  I <- (n / sum(W)) * sum(W * outer(z, z)) / sum(z^2)
  pm <- replicate(nperm, { zp <- sample(z)
    (n / sum(W)) * sum(W * outer(zp, zp)) / sum(zp^2) })
  c(I, (sum(abs(pm) >= abs(I)) + 1) / (nperm + 1), sum(W) / 2)
}

CUTS <- c(50, 75, 100, 150, 200, 300, 400, 600, 800)
BANDS <- c(0, 100, 200, 300, 400, 600, 800, 1200, 2000)
cor_all <- list(); se_all <- list()

for (g in c("inhibit_first", "promote_first")) {
  d <- d0[start_dir == g]; d <- d[complete.cases(d[, c("score", PRED), with = FALSE])]
  z <- copy(d)[, c("score", PRED), with = FALSE]
  for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  m <- lm(as.formula(paste("score ~", paste(PRED, collapse = "+"))), z)
  D <- distmat(d$longitude, d$latitude); r <- residuals(m)
  lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")

  log_msg("\n########## ", lab, " | n = ", nrow(z), " ##########")
  log_msg("站点两两距离: 中位 ", round(median(D[upper.tri(D)])), " km | 最大 ",
          round(max(D)), " km")

  log_msg("\n--- (1) 残差相关图: 各距离带的 Moran's I ---")
  cg <- rbindlist(lapply(seq_len(length(BANDS) - 1), function(i) {
    v <- band_I(r, D, BANDS[i], BANDS[i+1])
    data.table(grp = lab, lo = BANDS[i], hi = BANDS[i+1],
               I = v[1], p = v[2], npair = v[3]) }))
  print(cg[, .(距离带 = paste0(lo, "-", hi, "km"), I = round(I, 4),
               p = signif(p, 3), 站对数 = npair,
               显著 = fifelse(!is.na(p) & p < .05, "*", ""))])
  cor_all[[g]] <- cg

  log_msg("\n--- (2)(3) 标准误敏感曲线与有效样本量 ---")
  se_tab <- rbindlist(lapply(CUTS, function(k) {
    s <- conley_se(m, D, k)[PRED]
    nb <- mean(rowSums(D <= k)) - 1
    data.table(grp = lab, cut = k, var = PRED, se = s,
               nb = nb, eff_n = nrow(z) / max(1, nb)) }))
  se0 <- sqrt(diag(vcov(m)))[PRED]
  se_tab[, ratio := se / se0[var]]
  sm <- se_tab[, .(平均邻居 = round(nb[1], 1), 有效样本 = round(eff_n[1], 1),
                   SE比中位 = round(median(ratio), 3),
                   SE比范围 = sprintf("%.2f-%.2f", min(ratio), max(ratio))), by = cut]
  print(sm)
  se_all[[g]] <- se_tab
}

cg <- rbindlist(cor_all); se <- rbindlist(se_all)
fwrite(cg, file.path(OUT, "conley_correlogram.csv"))
fwrite(se, file.path(OUT, "conley_se_by_cutoff.csv"))

cg[, mid := (lo + hi) / 2]
p1 <- ggplot(cg, aes(mid, I, colour = grp)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_line(linewidth = .8) +
  geom_point(aes(shape = !is.na(p) & p < .05), size = 2.4) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = c(抑制组 = "#1C6B66", 促进组 = "#A2542C"), name = NULL) +
  labs(title = "A  残差相关图：空间相关在多远的距离上消失",
       subtitle = "实心点 = 该距离带的 Moran's I 显著（置换检验 p<0.05）",
       x = "距离带中点（km）", y = "Moran's I") +
  theme_minimal(base_size = 11) + theme(legend.position = "top",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

KEY <- c(urban_rate = "城镇化率", mixedleaf = "混交林",
         ever_broad = "常绿阔叶林", tavg = "平均气温")
p2 <- ggplot(se[var %in% names(KEY)], aes(cut, ratio, colour = grp)) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = .35) +
  geom_vline(xintercept = 200, linetype = "22", colour = "grey45") +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  facet_wrap(~ factor(KEY[var], levels = KEY), nrow = 1) +
  scale_x_continuous(breaks = c(50, 200, 400, 800)) +
  scale_colour_manual(values = c(抑制组 = "#1C6B66", 促进组 = "#A2542C"), name = NULL) +
  labs(title = "B  标准误对截断距离的敏感曲线",
       subtitle = "纵轴 = Conley 标准误 / 普通 OLS 标准误；虚线为本文采用的 200 km",
       x = "截断距离（km）", y = "标准误比值") +
  theme_minimal(base_size = 11) + theme(legend.position = "top",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"))

ggsave(file.path(OUT, "fig_conley_cutoff.png"), p1 / p2 + plot_layout(heights = c(1, 1.1)),
       width = 9.5, height = 8.2, dpi = 300)
log_msg("\n已输出 fig_conley_cutoff.png")
