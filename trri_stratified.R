# =============================================================================
# TRRI 分层回归分析
# 在气候大类 B / C / D 内分别跑含投资的多元回归
# 输出：系数方向、显著性、投资贡献大小（标准化系数 + 方差分解）
# =============================================================================

pacman::p_load(dplyr, ggplot2, tidyr, purrr, targets,
               MASS, vegan, scales, ggrepel, showtext, sysfonts)

font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto()
showtext_opts(dpi = 300)

BASE_SIZE <- 18
theme_cn <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size) +
    theme(text          = element_text(family = "heiti"),
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

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
if (!dir.exists("data_proc")) dir.create("data_proc")

# =============================================================================
# 0. 数据准备（复用 trri_analysis.R 的构建逻辑）
# =============================================================================
results_weekly <- readRDS("data_proc/results_weekly_0_5.rds")
suppressMessages(tar_load(invest_metrics))
suppressMessages(tar_load(data_heat_sif_weekly))

climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    vpd_mean_summer     = mean(vpd_mean,            na.rm = TRUE),
    heat_freq_mean      = mean(heat_event_freq,      na.rm = TRUE),
    heat_intensity_mean = mean(heat_index_composite, na.rm = TRUE),
    .groups = "drop"
  )

find_transition <- function(coefs, to_negative = TRUE) {
  cond <- if (to_negative) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs)) {
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1))
  }
  NA_integer_
}

df_all <- results_weekly %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq    = list(mean_coef),
    coef_tp0    = mean_coef[tp == 0],
    n_tp        = n(),
    longitude   = first(longitude),
    latitude    = first(latitude),
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
      TRUE                        ~ "全程抑制"   # 兜底（极少数边缘案例）
    ),
    TRRI = case_when(
      station_type == "全程抑制"  ~ 1L,
      station_type == "抑制→促进" ~ as.integer(1L + (6L - ITW)),
      station_type == "促进→抑制" ~ as.integer(6L + HTW),
      station_type == "全程促进"  ~ 12L
    )
  ) %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  left_join(climate_bg,    by = "meteo_stat_id") %>%
  filter(
    !is.na(TRRI),
    !is.na(invest_pa_tot), invest_pa_tot > 0,
    !is.na(vpd_mean_summer),
    !is.na(latitude)
  ) %>%
  mutate(
    TRRI_num     = as.numeric(TRRI),
    log_invest   = log10(invest_pa_tot),
    # 所有自变量标准化（组内，在后续 map 里做）
  )

cat("各气候大类样本量：\n")
print(count(df_all, koppen_group))

# =============================================================================
# 1. 分层回归函数
#    对每个气候大类：
#    (a) 线性回归 lm（系数 + p值，易于解读）
#    (b) 有序逻辑回归 polr（更严格，TRRI 是有序变量）
#    (c) 方差分解 varpart（投资 / 气候 / 地理的独立贡献）
# =============================================================================
run_stratum <- function(df, group_label) {
  cat("\n", strrep("-", 50), "\n")
  cat(sprintf("【气候大类 %s】  N = %d\n", group_label, nrow(df)))
  cat(strrep("-", 50), "\n")

  if (nrow(df) < 20) {
    cat("  样本量不足，跳过\n")
    return(NULL)
  }

  # 标准化（组内）
  df <- df %>%
    mutate(
      invest_s = scale(log_invest)[, 1],
      vpd_s    = scale(vpd_mean_summer)[, 1],
      heat_s   = scale(heat_freq_mean)[, 1],
      lat_s    = scale(latitude)[, 1]
    )

  # ── (a) 线性回归 ────────────────────────────────────────────────────────────
  lm_fit <- lm(TRRI_num ~ invest_s + vpd_s + heat_s + lat_s, data = df)
  lm_tbl <- as.data.frame(coef(summary(lm_fit))) %>%
    setNames(c("Estimate","SE","t","p")) %>%
    mutate(
      var    = rownames(.),
      sig    = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05 ~ "*",   p < 0.10 ~ ".",  TRUE ~ ""),
      group  = group_label,
      method = "lm"
    ) %>%
    filter(var != "(Intercept)")
  cat("\n线性回归（标准化系数）：\n")
  print(lm_tbl[, c("var","Estimate","SE","p","sig")])

  # R² 分解（各预测变量的半偏 R²，用 lmg 法近似）
  r2_total <- summary(lm_fit)$r.squared
  cat(sprintf("  总 R² = %.3f\n", r2_total))

  # ── (b) 有序逻辑回归 ────────────────────────────────────────────────────────
  df <- df %>%
    mutate(TRRI_ord = factor(TRRI_num, levels = sort(unique(TRRI_num)), ordered = TRUE))

  olr_fit <- tryCatch(
    polr(TRRI_ord ~ invest_s + vpd_s + heat_s + lat_s,
         data = df, Hess = TRUE),
    error = function(e) { cat("  有序回归失败:", conditionMessage(e), "\n"); NULL }
  )

  olr_tbl <- NULL
  if (!is.null(olr_fit)) {
    ct <- coef(summary(olr_fit))
    pv <- pnorm(abs(ct[, "t value"]), lower.tail = FALSE) * 2
    olr_tbl <- as.data.frame(ct[1:4, ]) %>%
      setNames(c("Estimate","SE","t")) %>%
      mutate(
        p      = pv[1:4],
        var    = rownames(.),
        sig    = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                           p < 0.05 ~ "*",   p < 0.10 ~ ".", TRUE ~ ""),
        group  = group_label,
        method = "polr"
      )
    cat("\n有序逻辑回归（标准化系数）：\n")
    print(olr_tbl[, c("var","Estimate","SE","p","sig")])
  }

  # ── (c) 方差分解 ────────────────────────────────────────────────────────────
  X1 <- df %>% dplyr::select(invest_s) %>% as.matrix()
  X2 <- df %>% dplyr::select(vpd_s, heat_s) %>% as.matrix()
  X3 <- df %>% dplyr::select(lat_s) %>% as.matrix()

  vp <- tryCatch(
    varpart(df$TRRI_num, X1, X2, X3),
    error = function(e) NULL
  )

  vp_r2 <- if (!is.null(vp)) vp$part$indfract$Adj.R.square else rep(NA, 8)
  cat(sprintf("\n方差分解：投资独立=%.1f%% | 气候独立=%.1f%% | 地理独立=%.1f%% | 残差=%.1f%%\n",
              vp_r2[1]*100, vp_r2[2]*100, vp_r2[3]*100, vp_r2[8]*100))

  list(
    group   = group_label,
    n       = nrow(df),
    lm_tbl  = lm_tbl,
    olr_tbl = olr_tbl,
    vp_r2   = vp_r2,
    r2_lm   = r2_total
  )
}

# 只跑 B / C / D（A 组 n=2，不可行）
groups_to_run <- c("B", "C", "D")
results_list  <- map(groups_to_run, function(g) {
  run_stratum(df_all %>% filter(koppen_group == g), g)
})
names(results_list) <- groups_to_run

# =============================================================================
# 2. 森林图：各气候大类的标准化回归系数（lm）
# =============================================================================
cat("\n【绘制森林图】\n")

var_labels <- c(
  invest_s = "绿化投资",
  vpd_s    = "VPD（干燥度）",
  heat_s   = "热事件频率",
  lat_s    = "纬度"
)
group_labels <- c(B = "B类(热带干旱)", C = "C类(温带)", D = "D类(大陆)")
group_colors <- c("B类(热带干旱)" = "#FF7F00", "C类(温带)" = "#4DAF4A", "D类(大陆)" = "#377EB8")

# 合并 lm 结果
lm_all <- map_dfr(results_list, ~ .x$lm_tbl) %>%
  mutate(
    var_label   = var_labels[var],
    group_label = group_labels[group],
    group_label = factor(group_label, levels = group_labels),
    # 95% CI
    ci_lo = Estimate - 1.96 * SE,
    ci_hi = Estimate + 1.96 * SE
  )

p_forest <- ggplot(lm_all,
                   aes(x = Estimate, y = group_label, color = group_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.25, linewidth = 0.9) +
  geom_point(aes(shape = sig %in% c("*","**","***")), size = 4) +
  geom_text(aes(x = ci_hi + 0.02, label = sig),
            hjust = 0, size = 6, fontface = "bold", family = "heiti") +
  scale_color_manual(values = group_colors, name = "气候大类") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("不显著","p<0.05"), name = "显著性") +
  facet_wrap(~var_label, scales = "free_x", ncol = 2) +
  labs(
    title    = "各气候大类中各自变量对 TRRI 的标准化回归系数",
    subtitle = "点为系数估计值，横线为 95% 置信区间；*** p<0.001  ** p<0.01  * p<0.05",
    x = "标准化系数 (beta)",
    y = NULL
  ) +
  theme_cn() +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("data_proc/trri_strat_forest.png", p_forest,
       width = 10, height = 8, dpi = 300)
cat("-> data_proc/trri_strat_forest.png\n")

# =============================================================================
# 3. 方差分解对比图（三气候大类并排）
# =============================================================================
vp_labels <- c("投资\n独立","气候\n独立","地理\n独立",
                "投资×气候","投资×地理","气候×地理","三者共享")
vp_colors <- c("#4DAF4A","#FF7F00","#377EB8",
               "#A6D96A","#ABD9E9","#FDAE61","#D9EF8B")

vp_df <- map_dfr(results_list, function(r) {
  data.frame(
    group     = r$group,
    component = vp_labels,
    r2        = pmax(r$vp_r2[1:7], 0)
  )
}) %>%
  mutate(
    group_label = group_labels[group],
    group_label = factor(group_label, levels = group_labels),
    component   = factor(component, levels = vp_labels),
    label       = ifelse(r2 > 0.005, sprintf("%.1f%%", r2 * 100), "")
  )

p_vp <- ggplot(vp_df, aes(x = component, y = r2, fill = component)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.4, size = 5.5, family = "heiti") +
  facet_wrap(~group_label, ncol = 3) +
  scale_fill_manual(values = vp_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    title    = "各气候大类内 TRRI 方差分解",
    subtitle = "投资 / 气候 / 地理对韧性等级变异的独立解释量",
    x = NULL, y = "解释方差（Adj. R²）"
  ) +
  theme_cn() +
  theme(panel.grid.major.x = element_blank())

ggsave("data_proc/trri_strat_varpart.png", p_vp,
       width = 12, height = 5, dpi = 300)
cat("-> data_proc/trri_strat_varpart.png\n")

# =============================================================================
# 4. 汇总表：投资系数 × 气候大类（lm + polr 并排）
# =============================================================================
cat("\n【汇总：投资系数对比】\n")

invest_summary <- map_dfr(results_list, function(r) {
  lm_row  <- r$lm_tbl  %>% filter(var == "invest_s")
  olr_row <- r$olr_tbl %>% filter(var == "invest_s")
  data.frame(
    group        = r$group,
    n            = r$n,
    lm_beta      = round(lm_row$Estimate, 3),
    lm_p         = round(lm_row$p, 4),
    lm_sig       = lm_row$sig,
    polr_beta    = round(olr_row$Estimate, 3),
    polr_p       = round(olr_row$p, 4),
    polr_sig     = olr_row$sig,
    invest_indR2 = sprintf("%.1f%%", r$vp_r2[1] * 100),
    total_R2_lm  = sprintf("%.1f%%", r$r2_lm  * 100)
  )
})
print(invest_summary)

cat("\n解读说明：\n")
cat("  β > 0：投资↑ → TRRI↑（韧性增强）\n")
cat("  β < 0：投资↑ → TRRI↓（韧性减弱，需检查）\n")
cat("  投资独立R²：控制气候和地理后，投资单独解释的方差比例\n")
