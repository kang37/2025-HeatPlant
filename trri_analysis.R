# =============================================================================
# TRRI Analysis: Thermal Response Resilience Index
# 将全程促进、全程抑制、HTW（促→抑）、ITW（抑→促）统一到一个 12 级有序尺度
# 并行保留连续版本 mean_coef_all 作为鲁棒性验证
#
# 三层分析：
#   第一层（描述）：TRRI 分布 × 气候带
#   第二层（归因）：方差分解 — 投资 / 气候 / 地理对 TRRI 的独立贡献
#   第三层（阈值）：断点回归 — 绿化投资对 TRRI 的阈值效应
# =============================================================================

pacman::p_load(
  dplyr, ggplot2, tidyr, purrr, targets,
  segmented, vegan, MASS,
  scales, RColorBrewer, ggrepel
)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
if (!dir.exists("data_proc")) dir.create("data_proc")

# =============================================================================
# 0. 读取数据
# =============================================================================
cat("【0. 读取数据】\n")

results_weekly <- readRDS("data_proc/results_weekly_0_5.rds")
suppressMessages(tar_load(invest_metrics))
suppressMessages(tar_load(data_heat_sif_weekly))

cat("站点数:", length(unique(results_weekly$meteo_stat_id)), "\n")

# =============================================================================
# 1. 构建 TRRI（12 级有序）和 mean_coef_all（连续）
# =============================================================================
cat("\n【1. 构建 TRRI 变量】\n")

# 辅助：找首次稳定转变点（通用版，方向可指定）
find_transition <- function(coefs, to_negative = TRUE) {
  cond <- if (to_negative) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs)) {
    if (!is.na(coefs[i]) && cond(coefs[i])) {
      remaining <- coefs[i:length(coefs)]
      if (mean(sapply(remaining, cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1))
    }
  }
  NA_integer_
}

station_df <- results_weekly %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq     = list(mean_coef),
    coef_tp0     = mean_coef[tp == 0],
    coef_tp5     = mean_coef[tp == 5],
    rho_mean     = mean(rho,       na.rm = TRUE),
    mean_coef_all = mean(mean_coef, na.rm = TRUE),
    longitude    = first(longitude),
    latitude     = first(latitude),
    koppen_class = first(koppen_class),
    koppen_group = first(koppen_group),
    n_tp         = n(),
    .groups = "drop"
  ) %>%
  filter(n_tp == 6) %>%
  mutate(
    HTW = map_int(coef_seq, find_transition, to_negative = TRUE),
    ITW = map_int(coef_seq, find_transition, to_negative = FALSE),
    # 站点类型
    station_type = case_when(
      coef_tp0 > 0 & is.na(HTW)  ~ "全程促进",
      coef_tp0 < 0 & is.na(ITW)  ~ "全程抑制",
      coef_tp0 > 0 & !is.na(HTW) ~ "促进→抑制",
      coef_tp0 < 0 & !is.na(ITW) ~ "抑制→促进",
      TRUE                        ~ "全程抑制"
    ),
    # TRRI：12 级有序（1=最脆弱，12=最韧性）
    TRRI = case_when(
      station_type == "全程抑制"   ~ 1L,
      station_type == "抑制→促进"  ~ as.integer(1L + (6L - ITW)),  # ITW=5→2, ITW=1→6
      station_type == "促进→抑制"  ~ as.integer(6L + HTW),         # HTW=1→7, HTW=5→11
      station_type == "全程促进"   ~ 12L
    )
  )

# TRRI 标签（用于图例）
trri_labels <- c(
  "1"  = "全程抑制",
  "2"  = "抑→促 (ITW=5)",
  "3"  = "抑→促 (ITW=4)",
  "4"  = "抑→促 (ITW=3)",
  "5"  = "抑→促 (ITW=2)",
  "6"  = "抑→促 (ITW=1)",
  "7"  = "促→抑 (HTW=1)",
  "8"  = "促→抑 (HTW=2)",
  "9"  = "促→抑 (HTW=3)",
  "10" = "促→抑 (HTW=4)",
  "11" = "促→抑 (HTW=5)",
  "12" = "全程促进"
)

# 颜色：红→橙→黄→绿（连续渐变）
trri_colors <- colorRampPalette(c("#B2182B","#EF8A62","#FDDBC7",
                                   "#D9F0D3","#74C476","#1B7837"))(12)
names(trri_colors) <- 1:12

cat("TRRI 分布:\n")
print(table(station_df$TRRI))
cat("\n站点类型:\n")
print(count(station_df, station_type))

# =============================================================================
# 第一层：描述 — TRRI 分布 × 气候带
# =============================================================================
cat("\n【第一层：描述性分析】\n")

# 气候带标签
climate_labels <- c(
  Am="Am 热带季风", As="As 热带草原",
  BSh="BSh 热半干旱", BWh="BWh 热荒漠", BWk="BWk 冷荒漠",
  Csc="Csc 温带地中海", Cwa="Cwa 温带季风", Cwc="Cwc 温带高原",
  Dsd="Dsd 大陆干旱", Dwa="Dwa 大陆季风", Dwd="Dwd 大陆严寒"
)

# 气候带排序（A→B→C→D，组内按总站点降序）
climate_order <- station_df %>%
  count(koppen_class, koppen_group) %>%
  arrange(koppen_group, desc(n)) %>%
  pull(koppen_class)

# ── 1a. 热力图：气候带 × TRRI 等级 ──────────────────────────────────────────
hm_data <- station_df %>%
  mutate(
    climate = factor(climate_labels[koppen_class], levels = rev(climate_labels[climate_order])),
    TRRI_f  = factor(TRRI, levels = 1:12)
  ) %>%
  count(climate, koppen_group, TRRI_f, .drop = FALSE) %>%
  group_by(climate) %>%
  mutate(total = sum(n), pct = n / total * 100) %>%
  ungroup() %>%
  filter(!is.na(climate), total > 0)

p_heatmap <- ggplot(hm_data, aes(x = TRRI_f, y = climate, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 0, as.character(n), "")),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#F7FBFF","#9ECAE1","#2171B5","#08306B"),
    name = "气候带内\n占比 (%)", limits = c(0, 100)
  ) +
  scale_x_discrete(
    labels = c("1\n全程\n抑制",
               "2","3","4","5","6\n↑\nITW",
               "7\n↑\nHTW","8","9","10","11",
               "12\n全程\n促进"),
    position = "bottom"
  ) +
  annotate("rect", xmin = 0.5, xmax = 1.5,  ymin = -Inf, ymax = Inf,
           fill = "#B2182B", alpha = 0.08) +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
           fill = "#1B7837", alpha = 0.08) +
  annotate("rect", xmin = 1.5, xmax = 6.5,  ymin = -Inf, ymax = Inf,
           fill = "#FC8D59", alpha = 0.04) +
  annotate("rect", xmin = 6.5, xmax = 11.5, ymin = -Inf, ymax = Inf,
           fill = "#74C476", alpha = 0.04) +
  labs(
    title    = "热胁迫韧性（TRRI）× 柯本气候带分布热力图",
    subtitle = "格内数字为站点数；列从左到右代表韧性由弱到强",
    x = "TRRI 等级（←脆弱　　　　　　　　　韧性→）",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x     = element_text(size = 8.5, lineheight = 1),
    axis.text.y     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave("data_proc/trri_heatmap_climate.png", p_heatmap,
       width = 13, height = 7, dpi = 300)
cat("-> data_proc/trri_heatmap_climate.png\n")

# ── 1b. 箱线图：TRRI 按气候带 ───────────────────────────────────────────────
p_box_climate <- station_df %>%
  mutate(climate = factor(climate_labels[koppen_class],
                          levels = climate_labels[climate_order])) %>%
  ggplot(aes(x = climate, y = TRRI, fill = koppen_group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5, outlier.alpha = 0.5) +
  geom_hline(yintercept = 6.5, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  annotate("text", x = 0.6, y = 3.5, label = "抑制区", color = "#B2182B",
           size = 3.5, fontface = "italic", hjust = 0) +
  annotate("text", x = 0.6, y = 9.5, label = "促进区", color = "#1B7837",
           size = 3.5, fontface = "italic", hjust = 0) +
  scale_fill_manual(values = c(A="#E41A1C", B="#FF7F00", C="#4DAF4A", D="#377EB8"),
                    name = "柯本大类") +
  scale_y_continuous(breaks = c(1,6,7,12),
                     labels = c("1\n全程抑制","6\nITW=1","7\nHTW=1","12\n全程促进")) +
  labs(
    title = "各气候带的热胁迫韧性（TRRI）分布",
    subtitle = "虚线：TRRI=6.5 为抑制区/促进区分界",
    x = NULL, y = "TRRI（韧性等级）"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9.5),
    panel.grid.major.x = element_blank()
  )

ggsave("data_proc/trri_boxplot_climate.png", p_box_climate,
       width = 11, height = 6, dpi = 300)
cat("-> data_proc/trri_boxplot_climate.png\n")

# ── 1c. 散点图：coef_tp0 vs coef_tp5（四象限）──────────────────────────────
type_colors <- c("全程促进"="#1B7837","全程抑制"="#B2182B",
                 "促进→抑制"="#74C476","抑制→促进"="#FC8D59")

p_quad <- station_df %>%
  mutate(
    coef_tp0_clip = pmax(pmin(coef_tp0, 0.3), -0.3),
    coef_tp5_clip = pmax(pmin(coef_tp5, 0.3), -0.3)
  ) %>%
  ggplot(aes(x = coef_tp0_clip, y = coef_tp5_clip, color = station_type)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_point(alpha = 0.55, size = 1.8) +
  annotate("text", x =  0.22, y =  0.22, label = "全程促进", color = "#1B7837",
           fontface = "bold", size = 4) +
  annotate("text", x = -0.22, y = -0.22, label = "全程抑制", color = "#B2182B",
           fontface = "bold", size = 4) +
  annotate("text", x =  0.22, y = -0.22, label = "促进→抑制", color = "#4DAF4A",
           fontface = "bold", size = 4) +
  annotate("text", x = -0.22, y =  0.22, label = "抑制→促进", color = "#E6550D",
           fontface = "bold", size = 4) +
  scale_color_manual(values = type_colors, guide = "none") +
  labs(
    title    = "热胁迫因果系数轨迹图（tp=0 vs tp=5）",
    subtitle = "每点为一个站点；截断于 ±0.3 以改善可读性",
    x = "即时效应系数（tp=0）",
    y = "延迟效应系数（tp=5）"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("data_proc/trri_quadrant_scatter.png", p_quad,
       width = 7, height = 7, dpi = 300)
cat("-> data_proc/trri_quadrant_scatter.png\n")

# =============================================================================
# 第二层：归因 — 方差分解（TRRI + mean_coef_all）
# =============================================================================
cat("\n【第二层：方差分解】\n")

climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    vpd_mean_summer     = mean(vpd_mean,            na.rm = TRUE),
    heat_freq_mean      = mean(heat_event_freq,      na.rm = TRUE),
    heat_intensity_mean = mean(heat_index_composite, na.rm = TRUE),
    .groups = "drop"
  )

df_full <- station_df %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  left_join(climate_bg,    by = "meteo_stat_id") %>%
  filter(!is.na(invest_pa_tot), invest_pa_tot > 0,
         !is.na(vpd_mean_summer), !is.na(latitude)) %>%
  mutate(
    TRRI_num   = as.numeric(TRRI),
    log_invest = log10(invest_pa_tot),
    log_invest_r = log1p(coalesce(invest_ratio, 0))
  )

cat("方差分解样本数:", nrow(df_full), "\n")

run_varpart <- function(Y_vec, label) {
  X1 <- df_full %>%
    dplyr::select(log_invest, log_invest_r) %>%
    mutate(across(everything(), ~scale(.)[,1])) %>% as.matrix()
  X2 <- df_full %>%
    dplyr::select(vpd_mean_summer, heat_freq_mean, heat_intensity_mean) %>%
    mutate(across(everything(), ~scale(.)[,1])) %>% as.matrix()
  X3 <- df_full %>%
    dplyr::select(latitude, longitude) %>%
    mutate(across(everything(), ~scale(.)[,1])) %>% as.matrix()

  vp <- varpart(Y_vec, X1, X2, X3)
  r2 <- vp$part$indfract$Adj.R.square

  cat(sprintf("\n  [%s] 绿化投资=%.1f%% | 气候=%.1f%% | 地理=%.1f%% | 残差=%.1f%%\n",
              label, r2[1]*100, r2[2]*100, r2[3]*100, r2[8]*100))
  list(vp = vp, r2 = r2, label = label)
}

vp_trri  <- run_varpart(df_full$TRRI_num,    "TRRI（12级）")
vp_coef  <- run_varpart(df_full$mean_coef_all, "mean_coef_all（连续）")

# ── 并排对比柱状图 ────────────────────────────────────────────────────────────
comp_labels <- c("绿化投资\n独立","气候背景\n独立","地理位置\n独立",
                 "投资×气候","投资×地理","气候×地理","三者共享")
comp_colors <- c("#4DAF4A","#FF7F00","#377EB8",
                 "#A6D96A","#ABD9E9","#FDAE61","#D9EF8B")

vp_cmp_df <- bind_rows(
  data.frame(Y = "TRRI（12级）",       component = comp_labels,
             r2 = pmax(vp_trri$r2[1:7], 0)),
  data.frame(Y = "mean_coef_all（连续）", component = comp_labels,
             r2 = pmax(vp_coef$r2[1:7], 0))
) %>%
  mutate(
    component = factor(component, levels = comp_labels),
    label     = ifelse(r2 > 0.005, sprintf("%.1f%%", r2*100), "")
  )

p_vp <- ggplot(vp_cmp_df, aes(x = component, y = r2, fill = component)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.8) +
  facet_wrap(~Y, ncol = 2) +
  scale_fill_manual(values = comp_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "TRRI 与 mean_coef_all 方差分解对比",
    subtitle = sprintf("N = %d 站点 | 两种 Y 变量的独立贡献对比", nrow(df_full)),
    x = NULL, y = "解释方差（Adj. R²）"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9, lineheight = 1.1),
    strip.text  = element_text(face = "bold", size = 11)
  )

ggsave("data_proc/trri_varpart_compare.png", p_vp,
       width = 12, height = 5, dpi = 300)
cat("-> data_proc/trri_varpart_compare.png\n")

# =============================================================================
# 第三层：阈值 — 断点回归（TRRI + mean_coef_all）
# =============================================================================
cat("\n【第三层：断点回归】\n")

run_segmented_plot <- function(y_var, y_label, out_file, line_color = "steelblue") {
  df_s <- df_full %>% mutate(Y = .data[[y_var]])
  if (nrow(df_s) < 20) { cat("  样本不足\n"); return(invisible(NULL)) }

  lm0 <- lm(Y ~ log_invest, data = df_s)

  seg_res <- tryCatch({
    seg <- segmented(lm0, seg.Z = ~log_invest,
                     npsi = 1, control = seg.control(it.max = 200))
    bp  <- seg$psi[, "Est."]
    dv  <- davies.test(lm0, ~log_invest, k = 10)
    cat(sprintf("  [%s] 阈值 T=%.3f (投资=%.2f), Davies p=%.4f\n",
                y_label, bp, 10^bp, dv$p.value))
    list(seg = seg, bp = bp, davies_p = dv$p.value)
  }, error = function(e) {
    cat("  断点回归失败:", conditionMessage(e), "\n"); NULL
  })

  # 散点 + 拟合
  x_seq <- seq(min(df_s$log_invest, na.rm=TRUE),
               max(df_s$log_invest, na.rm=TRUE), length.out = 300)
  model_for_pred <- if (!is.null(seg_res)) seg_res$seg else lm0
  fit_df <- data.frame(x = x_seq,
                       y = predict(model_for_pred,
                                   newdata = data.frame(log_invest = x_seq)))

  subtitle <- if (!is.null(seg_res))
    sprintf("断点回归 | Davies p = %.3f | N = %d", seg_res$davies_p, nrow(df_s))
  else sprintf("线性拟合（断点不收敛）| N = %d", nrow(df_s))

  p <- ggplot(df_s, aes(x = log_invest, y = Y)) +
    geom_point(aes(color = station_type), alpha = 0.5, size = 1.8) +
    geom_line(data = fit_df, aes(x = x, y = y),
              color = "black", linewidth = 1.2, inherit.aes = FALSE) +
    { if (!is.null(seg_res))
        list(
          geom_vline(xintercept = seg_res$bp, linetype = "dashed",
                     color = "red", linewidth = 1),
          annotate("label",
                   x = seg_res$bp + 0.08,
                   y = quantile(df_s$Y, 0.95, na.rm=TRUE),
                   label = sprintf("阈值\n投资=%.1f", 10^seg_res$bp),
                   hjust = 0, color = "red", size = 3.5,
                   fill = "white", alpha = 0.85)
        )
    } +
    { if (y_var == "TRRI_num")
        list(
          geom_hline(yintercept = 6.5, linetype = "dotted",
                     color = "grey40", linewidth = 0.7),
          annotate("text", x = min(df_s$log_invest, na.rm=TRUE) + 0.05,
                   y = 7.2, label = "促进区 ↑", color = "grey40",
                   size = 3.2, hjust = 0)
        )
    } +
    scale_color_manual(values = c("全程促进"="#1B7837","全程抑制"="#B2182B",
                                   "促进→抑制"="#74C476","抑制→促进"="#FC8D59"),
                       name = NULL) +
    scale_x_continuous(
      sec.axis = sec_axis(~10^., name = "绿化投资强度（原始尺度）",
                          breaks = c(0.1, 1, 10, 100, 1000))
    ) +
    labs(
      title    = sprintf("绿化投资对%s的阈值效应", y_label),
      subtitle = subtitle,
      x = "log10(绿化投资强度 invest_pa_tot)",
      y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.text     = element_text(size = 9)
    )

  ggsave(out_file, p, width = 9, height = 6, dpi = 300)
  cat(sprintf("  -> %s\n", out_file))
  invisible(seg_res)
}

seg_trri <- run_segmented_plot("TRRI_num",    "TRRI（韧性等级）",
                                "data_proc/trri_seg_trri.png")
seg_coef <- run_segmented_plot("mean_coef_all", "mean_coef_all（均值系数）",
                                "data_proc/trri_seg_coef.png")

# ── 有序逻辑回归（TRRI 作为有序因变量）──────────────────────────────────────
cat("\n【补充：有序逻辑回归（TRRI）】\n")

df_olr <- df_full %>%
  mutate(
    TRRI_ord     = factor(TRRI, levels = 1:12, ordered = TRUE),
    log_invest_s = scale(log_invest)[,1],
    vpd_s        = scale(vpd_mean_summer)[,1],
    heat_s       = scale(heat_freq_mean)[,1],
    lat_s        = scale(latitude)[,1]
  )

olr <- tryCatch(
  polr(TRRI_ord ~ log_invest_s + vpd_s + heat_s + lat_s,
       data = df_olr, Hess = TRUE),
  error = function(e) { cat("有序回归失败:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(olr)) {
  ct   <- coef(summary(olr))
  pv   <- pnorm(abs(ct[,"t value"]), lower.tail = FALSE) * 2
  tbl  <- round(cbind(as.data.frame(ct), p_value = pv)[1:4,], 4)
  cat("\n有序逻辑回归结果（标准化系数）:\n")
  print(tbl)
  beta_inv <- coef(olr)["log_invest_s"]
  cat(sprintf("\n绿化投资系数 = %.3f → %s\n", beta_inv,
              ifelse(beta_inv > 0,
                     "投资越高，TRRI 越大（韧性越强）✓",
                     "投资越高，TRRI 越小（方向不符）✗")))
}

# =============================================================================
# 汇总
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("【TRRI 分析汇总】\n")
cat(sprintf("总站点数: %d\n", nrow(station_df)))
cat(sprintf("  全程促进: %d | 促→抑(HTW): %d | 抑→促(ITW): %d | 全程抑制: %d\n",
            sum(station_df$station_type=="全程促进"),
            sum(station_df$station_type=="促进→抑制"),
            sum(station_df$station_type=="抑制→促进"),
            sum(station_df$station_type=="全程抑制")))
cat(sprintf("方差分解（N=%d）:\n", nrow(df_full)))
cat(sprintf("  TRRI:          投资=%.1f%% | 气候=%.1f%% | 地理=%.1f%%\n",
            vp_trri$r2[1]*100, vp_trri$r2[2]*100, vp_trri$r2[3]*100))
cat(sprintf("  mean_coef_all: 投资=%.1f%% | 气候=%.1f%% | 地理=%.1f%%\n",
            vp_coef$r2[1]*100, vp_coef$r2[2]*100, vp_coef$r2[3]*100))
if (!is.null(seg_trri))
  cat(sprintf("断点回归 TRRI: 阈值=%.2f（投资=%.1f），Davies p=%.4f\n",
              seg_trri$bp, 10^seg_trri$bp, seg_trri$davies_p))
cat(strrep("=", 60), "\n")
