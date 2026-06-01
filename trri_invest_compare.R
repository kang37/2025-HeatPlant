# =============================================================================
# TRRI 分析扩展：
#   Part A — 四种投资变量 × 三个气候大类的分层回归对比
#   Part B — TRRI 空间分布地图（中国）
# =============================================================================

pacman::p_load(dplyr, ggplot2, tidyr, purrr, targets,
               MASS, scales, sf, rnaturalearth, rnaturalearthdata,
               showtext, sysfonts)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
if (!dir.exists("data_proc")) dir.create("data_proc")

# ── 中文字体设置 ──────────────────────────────────────────────────────────────
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto()
showtext_opts(dpi = 300)   # 与 ggsave dpi 保持一致

# ── 全局 ggplot 主题（字体 + 字号）──────────────────────────────────────────
BASE_SIZE <- 18
theme_cn <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size) +
    theme(text = element_text(family = "heiti"),
          plot.title    = element_text(face = "bold", hjust = 0.5,
                                       size = base_size + 4),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40",
                                       size = base_size),
          strip.text    = element_text(face = "bold", size = base_size),
          axis.text     = element_text(size = base_size - 2),
          axis.title    = element_text(size = base_size),
          legend.text   = element_text(size = base_size - 2),
          legend.title  = element_text(face = "bold", size = base_size))
}
theme_set(theme_cn())

# =============================================================================
# 0. 数据重建（与 trri_stratified.R 相同基础）
# =============================================================================
results_weekly <- readRDS("data_proc/results_weekly_0_5.rds")
suppressMessages(tar_load(invest_metrics))
suppressMessages(tar_load(data_heat_sif_weekly))

climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    vpd_mean_summer = mean(vpd_mean,           na.rm = TRUE),
    heat_freq_mean  = mean(heat_event_freq,     na.rm = TRUE),
    .groups = "drop"
  )

find_transition <- function(coefs, to_negative = TRUE) {
  cond <- if (to_negative) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1))
  NA_integer_
}

df_base <- results_weekly %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq     = list(mean_coef),
    coef_tp0     = mean_coef[tp == 0],
    n_tp         = n(),
    longitude    = first(longitude),
    latitude     = first(latitude),
    koppen_group = first(koppen_group),
    .groups = "drop"
  ) %>%
  filter(n_tp == 6) %>%
  mutate(
    HTW = map_int(coef_seq, find_transition, TRUE),
    ITW = map_int(coef_seq, find_transition, FALSE),
    station_type = case_when(
      coef_tp0 > 0 & is.na(HTW)  ~ "全程促进",
      coef_tp0 < 0 & is.na(ITW)  ~ "全程抑制",
      coef_tp0 > 0 & !is.na(HTW) ~ "促进→抑制",
      coef_tp0 < 0 & !is.na(ITW) ~ "抑制→促进",
      TRUE ~ "全程抑制"
    ),
    TRRI = case_when(
      station_type == "全程抑制"  ~ 1L,
      station_type == "抑制→促进" ~ as.integer(1L + (6L - ITW)),
      station_type == "促进→抑制" ~ as.integer(6L + HTW),
      station_type == "全程促进"  ~ 12L
    ),
    TRRI_num = as.numeric(TRRI)
  ) %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  left_join(climate_bg,    by = "meteo_stat_id") %>%
  filter(!is.na(TRRI), !is.na(vpd_mean_summer), !is.na(latitude))

# =============================================================================
# Part A — 四种投资变量 × 气候大类分层回归对比
# =============================================================================
cat("【Part A：四种投资变量分层回归对比】\n")

invest_vars <- c(
  invest_ratio    = "invest_ratio",
  invest_pa_tot   = "invest_pa_tot",
  invest_pa_built = "invest_pa_built",
  invest_pa_park  = "invest_pa_park"
)
invest_labels <- c(
  invest_ratio    = "占GDP比例\n(invest_ratio)",
  invest_pa_tot   = "÷总绿地面积\n(pa_tot)",
  invest_pa_built = "÷建成区绿地\n(pa_built)",
  invest_pa_park  = "÷公园绿地\n(pa_park)"
)
group_levels <- c("B", "C", "D")
group_labels <- c(B = "B类(热带干旱)", C = "C类(温带)", D = "D类(大陆)")
group_colors <- c("B类(热带干旱)" = "#FF7F00",
                  "C类(温带)"     = "#4DAF4A",
                  "D类(大陆)"     = "#377EB8")

# 对每个投资变量 × 气候大类跑 lm（标准化系数）
run_one <- function(df, inv_var, group) {
  df_sub <- df %>%
    filter(koppen_group == group, !is.na(.data[[inv_var]]),
           .data[[inv_var]] > 0) %>%
    mutate(
      inv_s   = scale(log10(.data[[inv_var]]))[, 1],
      vpd_s   = scale(vpd_mean_summer)[, 1],
      heat_s  = scale(heat_freq_mean)[, 1],
      lat_s   = scale(latitude)[, 1]
    )

  if (nrow(df_sub) < 20) return(NULL)

  fit <- lm(TRRI_num ~ inv_s + vpd_s + heat_s + lat_s, data = df_sub)
  ct  <- coef(summary(fit))

  data.frame(
    invest_var   = inv_var,
    invest_label = invest_labels[inv_var],
    group        = group,
    group_label  = group_labels[group],
    n            = nrow(df_sub),
    beta         = ct["inv_s", "Estimate"],
    se           = ct["inv_s", "Std. Error"],
    p            = ct["inv_s", "Pr(>|t|)"],
    r2_total     = summary(fit)$r.squared,
    stringsAsFactors = FALSE
  )
}

coef_df <- map_dfr(invest_vars, function(iv) {
  map_dfr(group_levels, function(g) run_one(df_base, iv, g))
}) %>%
  mutate(
    ci_lo       = beta - 1.96 * se,
    ci_hi       = beta + 1.96 * se,
    sig         = case_when(
      p < 0.001 ~ "***", p < 0.01 ~ "**",
      p < 0.05  ~ "*",   p < 0.10 ~ ".",
      TRUE      ~ ""
    ),
    invest_label = factor(invest_label, levels = invest_labels),
    group_label  = factor(group_label,  levels = group_labels)
  )

cat("\n投资系数汇总：\n")
coef_df %>%
  dplyr::select(invest_var, group, n, beta, se, p, sig) %>%
  mutate(across(c(beta, se, p), ~round(., 4))) %>%
  print()

# ── 森林图：行 = 投资变量，列 = 气候大类 ──────────────────────────────────
p_forest <- ggplot(coef_df,
                   aes(x = beta, y = invest_label, color = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.3, linewidth = 1,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(shape = p < 0.05),
             size = 3.5, position = position_dodge(width = 0.7)) +
  geom_text(aes(x = ci_hi + 0.02, label = sig),
            position = position_dodge(width = 0.7),
            hjust = 0, size = 6, fontface = "bold",
            family = "heiti", show.legend = FALSE) +
  scale_color_manual(values = group_colors, name = "气候大类") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("不显著 (p≥0.05)", "显著 (p<0.05)"),
                     name = "显著性") +
  facet_wrap(~group_label, ncol = 3) +
  labs(
    title    = "四种绿化投资变量对 TRRI 的标准化回归系数（分气候大类）",
    subtitle = "控制了 VPD、热事件频率、纬度 | 横线为 95% CI",
    x = "标准化系数 (beta) | β>0 表示投资越高韧性越强",
    y = "投资变量"
  ) +
  theme_cn() +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("data_proc/trri_invest4_forest.png", p_forest,
       width = 13, height = 6, dpi = 300)
cat("\n-> data_proc/trri_invest4_forest.png\n")

# ── 气泡图：β 大小 × 显著性（更紧凑的总览）──────────────────────────────
p_bubble <- ggplot(coef_df,
                   aes(x = invest_label, y = group_label,
                       size = abs(beta), color = beta,
                       shape = p < 0.05)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = ifelse(p < 0.1, sprintf("β=%.2f%s", beta, sig), "")),
            vjust = -1.1, size = 6, color = "grey20", fontface = "bold",
            family = "heiti") +
  scale_color_gradient2(low = "#B2182B", mid = "white", high = "#1B7837",
                        midpoint = 0, name = "β 值",
                        limits = c(-1, 1), oob = scales::squish) +
  scale_size_continuous(range = c(3, 12), name = "|β|") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("p≥0.05", "p<0.05"), name = "显著性") +
  labs(
    title    = "四种投资变量对 TRRI 的效应强度与方向总览",
    subtitle = "圆圈大小 = |β|；颜色：绿=正效应，红=负效应；实心点 = p<0.05",
    x = "投资变量", y = "气候大类"
  ) +
  theme_cn() +
  theme(axis.text.x      = element_text(angle = 15, hjust = 0.8),
        panel.grid.major = element_line(color = "grey90"),
        legend.position  = "right")

ggsave("data_proc/trri_invest4_bubble.png", p_bubble,
       width = 10, height = 5, dpi = 300)
cat("-> data_proc/trri_invest4_bubble.png\n")

# =============================================================================
# Part B — TRRI 空间分布地图
# =============================================================================
cat("\n【Part B：TRRI 空间分布地图】\n")

# 中国地图底图
china      <- ne_countries(country = "China", scale = "medium", returnclass = "sf")
china_prov <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent == "Asia") %>%
  slice(0)  # 省级用国界代替，避免依赖 rnaturalearthhires

# 用 medium scale 的 10° 格网省级近似（直接用国界+站点分布）
china_border <- ne_countries(country = "China", scale = "medium", returnclass = "sf")

# TRRI 颜色方案（1=深红→12=深绿，12级渐变）
trri_palette <- colorRampPalette(
  c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7",
    "#D9F0D3", "#74C476", "#1B7837")
)(12)

# 站点数据（含坐标和TRRI）
map_df <- df_base %>%
  filter(!is.na(longitude), !is.na(latitude), !is.na(TRRI)) %>%
  mutate(
    TRRI_f = factor(TRRI, levels = 1:12),
    # 简化标签：分5段
    TRRI_group = case_when(
      TRRI == 1           ~ "1 全程抑制",
      TRRI %in% 2:4       ~ "2-4 抑→促(慢)",
      TRRI %in% 5:6       ~ "5-6 抑→促(快)",
      TRRI %in% 7:8       ~ "7-8 促→抑(快)",
      TRRI %in% 9:11      ~ "9-11 促→抑(慢)",
      TRRI == 12          ~ "12 全程促进"
    ),
    TRRI_group = factor(TRRI_group, levels = c(
      "1 全程抑制","2-4 抑→促(慢)","5-6 抑→促(快)",
      "7-8 促→抑(快)","9-11 促→抑(慢)","12 全程促进"
    ))
  )

group_palette <- c(
  "1 全程抑制"     = "#B2182B",
  "2-4 抑→促(慢)"  = "#EF8A62",
  "5-6 抑→促(快)"  = "#FDBF6F",
  "7-8 促→抑(快)"  = "#A6D96A",
  "9-11 促→抑(慢)" = "#4DAF4A",
  "12 全程促进"    = "#1B7837"
)

cat("各 TRRI 组站点数：\n")
print(count(map_df, TRRI_group))

# ── 主地图：6色分组 ────────────────────────────────────────────────────────
p_map <- ggplot() +
  geom_sf(data = china_border, fill = "grey95", color = "grey40",
          linewidth = 0.5) +
  geom_point(
    data = map_df,
    aes(x = longitude, y = latitude, color = TRRI_group),
    size = 2.2, alpha = 0.8, shape = 16
  ) +
  scale_color_manual(values = group_palette,
                     name = "TRRI 等级（热胁迫韧性）",
                     guide = guide_legend(
                       override.aes = list(size = 4),
                       title.position = "top"
                     )) +
  coord_sf(xlim = c(73, 136), ylim = c(17, 54)) +
  labs(
    title    = "中国气象站点热胁迫响应韧性（TRRI）空间分布",
    subtitle = sprintf("N = %d 个站点 | 颜色：红→绿 = 脆弱→韧性", nrow(map_df)),
    x = "经度", y = "纬度"
  ) +
  theme_cn() +
  theme(legend.position = "right",
        panel.grid      = element_line(color = "grey85", linetype = "dotted"))

ggsave("data_proc/trri_spatial_map.png", p_map,
       width = 12, height = 8, dpi = 300)
cat("-> data_proc/trri_spatial_map.png\n")

# ── 分面地图：每个 TRRI 组单独高亮 ──────────────────────────────────────────
p_map_facet <- ggplot() +
  geom_sf(data = china_border, fill = "grey95", color = "grey60",
          linewidth = 0.3) +
  # 背景：所有其他站点（灰色）
  geom_point(
    data = map_df,
    aes(x = longitude, y = latitude),
    color = "grey80", size = 1, alpha = 0.4
  ) +
  # 前景：当前分面的站点
  geom_point(
    data = map_df,
    aes(x = longitude, y = latitude, color = TRRI_group),
    size = 2, alpha = 0.9
  ) +
  scale_color_manual(values = group_palette, guide = "none") +
  facet_wrap(~TRRI_group, ncol = 3) +
  coord_sf(xlim = c(73, 136), ylim = c(17, 54)) +
  labs(
    title    = "TRRI 各等级的空间分布（分面）",
    subtitle = "灰色为其余站点背景",
    x = NULL, y = NULL
  ) +
  theme_cn() +
  theme(axis.text  = element_blank(),
        panel.grid = element_line(color = "grey88", linetype = "dotted"))

ggsave("data_proc/trri_spatial_facet.png", p_map_facet,
       width = 13, height = 9, dpi = 300)
cat("-> data_proc/trri_spatial_facet.png\n")

# =============================================================================
# 汇总
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("【汇总：最强投资关联变量（按 |beta| 排序）】\n")
coef_df %>%
  dplyr::select(group, invest_var, beta, p, sig) %>%
  arrange(group, desc(abs(beta))) %>%
  mutate(across(c(beta, p), ~round(., 4))) %>%
  print()
cat(strrep("=", 60), "\n")
