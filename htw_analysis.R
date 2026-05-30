# =============================================================================
# HTW Analysis: Heat Tolerance Window
# 因变量：站点从"促进（mean_coef>0）"转为"抑制（mean_coef<0）"的最小滞后周数
# 三步走：
#   步骤1：统计"促→抑"模式站点数量，构建 HTW 变量
#   步骤2：断点回归，找绿化投资对 HTW 的阈值效应
#   步骤3：方差分解，量化投资/气候/地理的独立贡献百分比
# 注意：用 mean_coef 符号判断促进/抑制，规避中文编码比较问题
# =============================================================================

pacman::p_load(
  dplyr, ggplot2, tidyr, purrr, targets,
  segmented, vegan, MASS,
  ggrepel, scales, viridis, RColorBrewer
)

# =============================================================================
# 0. 读取数据
# =============================================================================
cat("【0. 读取数据】\n")

results_weekly <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  mutate(
    is_promote = (mean_coef > 0),   # 用系数符号判断，规避编码问题
    is_inhibit = (mean_coef < 0)
  )

suppressMessages(tar_load(invest_metrics))
suppressMessages(tar_load(data_heat_sif_weekly))

cat("CCM 结果:", nrow(results_weekly), "行,", length(unique(results_weekly$meteo_stat_id)), "个站点\n")

# =============================================================================
# 步骤1：构建 HTW 变量
# 定义：
#   - tp=0 时 mean_coef > 0（促进）
#   - 存在某个 tp* >= 1，使 mean_coef < 0（抑制），
#     且此后（tp* 到 5）超过半数滞后仍为抑制
#   HTW = 该 tp*（首次稳定转抑的滞后周数）
# =============================================================================
cat("\n【步骤1：构建 HTW 变量】\n")

find_htw <- function(coefs) {
  # coefs: 按 tp=0,1,2,3,4,5 排序的 mean_coef 向量
  if (length(coefs) < 2) return(NA_integer_)
  if (is.na(coefs[1]) || coefs[1] <= 0) return(NA_integer_)  # tp=0 须为促进

  for (i in 2:length(coefs)) {
    if (!is.na(coefs[i]) && coefs[i] < 0) {
      remaining <- coefs[i:length(coefs)]
      if (mean(remaining < 0, na.rm = TRUE) >= 0.5) {
        return(as.integer(i - 1))  # i=2 → tp=1
      }
    }
  }
  return(NA_integer_)
}

htw_df <- results_weekly %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    n_tp          = n(),
    coef_seq      = list(mean_coef),
    effect_seq    = list(is_promote),
    coef_tp0      = mean_coef[tp == 0],
    rho_mean      = mean(rho, na.rm = TRUE),
    longitude     = first(longitude),
    latitude      = first(latitude),
    koppen_group  = first(koppen_group),
    .groups = "drop"
  ) %>%
  filter(n_tp == 6) %>%
  mutate(
    HTW          = map_int(coef_seq, find_htw),
    promote_tp0  = (coef_tp0 > 0)
  )

n_total     <- nrow(htw_df)
n_promote0  <- sum(htw_df$promote_tp0, na.rm = TRUE)
n_htw_valid <- sum(!is.na(htw_df$HTW))

cat("总站点数（有完整6阶结果）:", n_total, "\n")
cat(sprintf("tp=0 为促进的站点: %d (%.1f%%)\n", n_promote0, n_promote0 / n_total * 100))
cat(sprintf("有效 HTW（促→抑转变）站点: %d (%.1f%% of 促进站点)\n",
            n_htw_valid, n_htw_valid / max(n_promote0, 1) * 100))
cat("\nHTW 分布:\n")
print(table(htw_df$HTW, useNA = "ifany"))

# 可视化 HTW 分布
p_htw_dist <- htw_df %>%
  filter(!is.na(HTW)) %>%
  ggplot(aes(x = factor(HTW))) +
  geom_bar(fill = "steelblue", color = "white", width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 5) +
  labs(
    title = "热胁迫缓冲窗口（HTW）分布",
    subtitle = sprintf("N = %d 个'促→抑'站点 | HTW = 首次稳定转抑的滞后周数", n_htw_valid),
    x = "HTW（周）",
    y = "站点数量"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("data_proc/htw_step1_distribution.png", p_htw_dist,
       width = 7, height = 5, dpi = 300)
cat("-> 图已保存: data_proc/htw_step1_distribution.png\n")

if (n_htw_valid < 30) {
  cat("\n[!] 有效站点数 <30，请查看下方按柯本气候带的分解，考虑调整 HTW 定义。\n")
  cat("    各气候带的促进站点数:\n")
  print(htw_df %>% filter(promote_tp0) %>% count(koppen_group))
}

# =============================================================================
# 合并投资与气候背景变量
# =============================================================================

climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    vpd_mean_summer     = mean(vpd_mean,             na.rm = TRUE),
    heat_freq_mean      = mean(heat_event_freq,       na.rm = TRUE),
    heat_intensity_mean = mean(heat_index_composite,  na.rm = TRUE),
    .groups = "drop"
  )

htw_full <- htw_df %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  left_join(climate_bg,    by = "meteo_stat_id") %>%
  filter(!is.na(HTW))

cat("\n合并后有效样本数:", nrow(htw_full), "\n")
cat("其中有投资数据:", sum(!is.na(htw_full$invest_pa_tot)), "\n")

# =============================================================================
# 步骤2：断点回归 —— 绿化投资对 HTW 的阈值效应
# =============================================================================
cat("\n【步骤2：断点回归 —— 投资阈值对 HTW】\n")

df_seg <- htw_full %>%
  filter(!is.na(invest_pa_tot), invest_pa_tot > 0) %>%
  mutate(
    log_invest = log10(invest_pa_tot),
    HTW_num    = as.numeric(HTW)
  )

cat("断点回归样本数:", nrow(df_seg), "\n")

# ── 断点回归函数 ──────────────────────────────────────────────────────────────
run_segmented_htw <- function(df) {
  if (nrow(df) < 20) {
    cat("  [!] 样本量不足（<20），跳过断点回归\n")
    return(invisible(NULL))
  }

  lm_base <- lm(HTW_num ~ log_invest, data = df)

  tryCatch({
    seg <- segmented(lm_base, seg.Z = ~log_invest,
                     npsi = 1, control = seg.control(it.max = 200))
    bp    <- seg$psi[, "Est."]
    bp_se <- seg$psi[, "St.Err"]
    dv    <- davies.test(lm_base, ~log_invest, k = 10)

    cat(sprintf("  阈值 T = %.3f (投资强度 = %.2f)，SE = %.3f，Davies p = %.4f\n",
                bp, 10^bp, bp_se, dv$p.value))

    list(seg = seg, bp = bp, bp_se = bp_se, davies_p = dv$p.value, df = df)
  }, error = function(e) {
    cat("  断点回归失败:", conditionMessage(e), "\n"); NULL
  })
}

seg_res <- run_segmented_htw(df_seg)

# ── 断点回归可视化 ────────────────────────────────────────────────────────────
if (!is.null(seg_res)) {
  df_plot <- seg_res$df
  x_rng   <- range(df_plot$log_invest, na.rm = TRUE)
  x_seq   <- seq(x_rng[1], x_rng[2], length.out = 300)
  y_fit   <- predict(seg_res$seg, newdata = data.frame(log_invest = x_seq))
  fit_df  <- data.frame(x = x_seq, y = y_fit)

  p_seg <- ggplot(df_plot, aes(x = log_invest, y = HTW_num)) +
    geom_jitter(aes(color = koppen_group), width = 0.03, height = 0.12,
                alpha = 0.65, size = 2.2) +
    geom_line(data = fit_df, aes(x = x, y = y),
              color = "black", linewidth = 1.3, inherit.aes = FALSE) +
    geom_vline(xintercept = seg_res$bp, linetype = "dashed",
               color = "red", linewidth = 1) +
    annotate("label", x = seg_res$bp + 0.05,
             y = max(df_plot$HTW_num, na.rm = TRUE) - 0.5,
             label = sprintf("阈值 T = %.2f\n投资强度 = %.1f",
                             seg_res$bp, 10^seg_res$bp),
             hjust = 0, color = "red", size = 3.8, fill = "white", alpha = 0.8) +
    scale_color_brewer(palette = "Set1", name = "柯本气候带") +
    scale_x_continuous(
      sec.axis = sec_axis(~10^., name = "绿化投资强度（原始尺度）",
                          breaks = c(1, 5, 10, 50, 100, 500))
    ) +
    labs(
      title = "绿化投资对热胁迫缓冲窗口（HTW）的阈值效应",
      subtitle = sprintf("断点回归 | Davies p = %.3f | N = %d 站点",
                         seg_res$davies_p, nrow(df_plot)),
      x = "log10(绿化投资强度 invest_pa_tot)",
      y = "HTW（首次转抑的滞后周数）"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("data_proc/htw_step2_segmented.png", p_seg,
         width = 9, height = 6, dpi = 300)
  cat("-> 图已保存: data_proc/htw_step2_segmented.png\n")
}

# ── 箱线图：按投资四分位看 HTW ─────────────────────────────────────────────
if (nrow(df_seg) >= 10) {
  q_breaks <- unique(quantile(df_seg$log_invest, probs = 0:4/4, na.rm = TRUE))
  n_q      <- length(q_breaks) - 1
  q_labels <- c("Q1\n（低投资）", "Q2", "Q3", "Q4\n（高投资）")[1:n_q]

  df_box <- df_seg %>%
    mutate(invest_q = cut(log_invest, breaks = q_breaks,
                          labels = q_labels, include.lowest = TRUE)) %>%
    filter(!is.na(invest_q))

  p_box <- ggplot(df_box, aes(x = invest_q, y = HTW_num, fill = invest_q)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1.8) +
    stat_summary(fun = mean, geom = "point", shape = 23,
                 size = 4, fill = "white", color = "black") +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    labs(
      title = "绿化投资分位数 vs 热胁迫缓冲窗口（HTW）",
      subtitle = "◆ 为均值；若投资越高 HTW 越大，则绿化能延缓热胁迫的负效应",
      x = "绿化投资强度分位数",
      y = "HTW（周）"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("data_proc/htw_step2_boxplot.png", p_box,
         width = 7, height = 5, dpi = 300)
  cat("-> 图已保存: data_proc/htw_step2_boxplot.png\n")

  # ANOVA
  aov_res <- aov(HTW_num ~ invest_q, data = df_box)
  cat("\nANOVA（投资四分位 vs HTW）:\n")
  print(summary(aov_res))
}

# =============================================================================
# 步骤3：方差分解 —— 投资 / 气候 / 地理对 HTW 的独立贡献
# =============================================================================
cat("\n【步骤3：方差分解（vegan::varpart）】\n")

df_vp <- htw_full %>%
  filter(
    !is.na(invest_pa_tot), invest_pa_tot > 0,
    !is.na(vpd_mean_summer),
    !is.na(latitude)
  ) %>%
  mutate(
    HTW_num        = as.numeric(HTW),
    log_invest     = log10(invest_pa_tot),
    log_invest_r   = log1p(coalesce(invest_ratio, 0))
  )

cat("方差分解样本数:", nrow(df_vp), "\n")

if (nrow(df_vp) >= 20) {

  Y_vp <- as.matrix(df_vp$HTW_num)

  X1_invest  <- df_vp %>%
    dplyr::select(log_invest, log_invest_r) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>%
    as.matrix()

  X2_climate <- df_vp %>%
    dplyr::select(vpd_mean_summer, heat_freq_mean, heat_intensity_mean) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>%
    as.matrix()

  X3_geo <- df_vp %>%
    dplyr::select(latitude, longitude) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>%
    as.matrix()

  vp <- varpart(Y_vp, X1_invest, X2_climate, X3_geo)

  fracs <- vp$part$indfract
  r2    <- fracs$Adj.R.square   # vegan 用 Adj.R.square（无 d）
  cat("\n方差分解独立贡献（Adjusted R²）:\n")
  cat(sprintf("  绿化投资 (X1): %.3f (%.1f%%)\n", r2[1], r2[1] * 100))
  cat(sprintf("  气候背景 (X2): %.3f (%.1f%%)\n", r2[2], r2[2] * 100))
  cat(sprintf("  地理位置 (X3): %.3f (%.1f%%)\n", r2[3], r2[3] * 100))
  cat(sprintf("  残差（未解释）: %.3f (%.1f%%)\n", r2[8], r2[8] * 100))

  # ── 维恩图 ─────────────────────────────────────────────────────────────────
  png("data_proc/htw_step3_varpart_venn.png", width = 900, height = 800, res = 150)
  plot(vp,
       Xnames = c("绿化投资", "气候背景", "地理位置"),
       bg     = c("#4DAF4A80", "#FF7F0080", "#377EB880"),
       digits = 3, cex = 1.3)
  title(main = sprintf("HTW 方差分解  (N = %d)", nrow(df_vp)),
        cex.main = 1.2)
  dev.off()
  cat("-> 图已保存: data_proc/htw_step3_varpart_venn.png\n")

  # ── 柱状图 ─────────────────────────────────────────────────────────────────
  labels_all <- c("绿化投资\n独立", "气候背景\n独立", "地理位置\n独立",
                  "投资×气候", "投资×地理", "气候×地理", "三者共享", "残差")
  colors_all <- c("#4DAF4A","#FF7F00","#377EB8",
                  "#A6D96A","#ABD9E9","#FDAE61","#D9EF8B","#CCCCCC")

  vp_tidy <- data.frame(
    component = factor(labels_all, levels = labels_all),
    r2        = r2,
    color     = colors_all
  ) %>%
    mutate(
      r2_show = pmax(r2, 0),
      label   = ifelse(abs(r2) > 0.005, sprintf("%.1f%%", r2 * 100), "")
    )

  p_bar <- vp_tidy %>%
    filter(component != "残差") %>%
    ggplot(aes(x = component, y = r2_show, fill = component)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = label), vjust = -0.4, size = 4.2) +
    scale_fill_manual(values = colors_all[-8]) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "HTW 方差分解：各成分对缓冲窗口变异的贡献",
      subtitle = sprintf("N = %d 站点 | 总解释量 = %.1f%%",
                         nrow(df_vp),
                         sum(pmax(r2[-8], 0)) * 100),
      x = NULL, y = "解释方差（Adjusted R²）"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, lineheight = 1.1)
    )

  ggsave("data_proc/htw_step3_varpart_bar.png", p_bar,
         width = 9, height = 5, dpi = 300)
  cat("-> 图已保存: data_proc/htw_step3_varpart_bar.png\n")

} else {
  cat("[!] 样本量不足（<20），跳过方差分解。\n")
}

# =============================================================================
# 补充：有序逻辑回归（更适合 HTW 是离散有序变量的情况）
# =============================================================================
cat("\n【补充：有序逻辑回归（Ordinal Logistic Regression）】\n")

if (nrow(df_vp) >= 20) {
  df_olr <- df_vp %>%
    mutate(
      HTW_ord      = factor(HTW, levels = sort(unique(HTW)), ordered = TRUE),
      log_invest_s = scale(log_invest)[, 1],
      vpd_s        = scale(vpd_mean_summer)[, 1],
      heat_s       = scale(heat_freq_mean)[, 1],
      lat_s        = scale(latitude)[, 1]
    )

  olr <- tryCatch(
    polr(HTW_ord ~ log_invest_s + vpd_s + heat_s + lat_s,
         data = df_olr, Hess = TRUE),
    error = function(e) { cat("有序回归失败:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(olr)) {
    coef_tbl <- coef(summary(olr))
    pvals    <- pnorm(abs(coef_tbl[, "t value"]), lower.tail = FALSE) * 2
    result_tbl <- cbind(as.data.frame(coef_tbl), p_value = round(pvals, 4))

    cat("\n有序逻辑回归结果（标准化系数）:\n")
    print(result_tbl[1:4, ])  # 只显示预测因子行

    beta_inv <- coef(olr)["log_invest_s"]
    cat(sprintf("\n绿化投资系数 = %.3f → %s\n", beta_inv,
                ifelse(beta_inv > 0,
                       "投资越高，HTW 越大（缓冲窗口越宽）✓",
                       "投资越高，HTW 越小（方向不符）✗")))
  }
}

# =============================================================================
# ▶▶ ITW 分析：Inhibition-to-Promotion Window（抑→促恢复窗口）
# 对称于 HTW，分析 tp=0 为抑制、后来转为促进的站点
# 管理意义：投资是否能缩短恢复窗口（ITW 越小 = 恢复越快）？
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("【ITW 分析：抑制→促进恢复窗口】\n")
cat(strrep("=", 60), "\n")

# ── ITW 构建函数（与 HTW 镜像对称）─────────────────────────────────────────
find_itw <- function(coefs) {
  if (length(coefs) < 2 || is.na(coefs[1]) || coefs[1] >= 0) return(NA_integer_)
  for (i in 2:length(coefs)) {
    if (!is.na(coefs[i]) && coefs[i] > 0) {
      remaining <- coefs[i:length(coefs)]
      if (mean(remaining > 0, na.rm = TRUE) >= 0.5) return(as.integer(i - 1))
    }
  }
  NA_integer_
}

# ── 步骤 ITW-1：构建 ITW ──────────────────────────────────────────────────
cat("\n【ITW 步骤1：构建 ITW 变量】\n")

itw_df <- results_weekly %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    n_tp         = n(),
    coef_seq     = list(mean_coef),
    coef_tp0     = mean_coef[tp == 0],
    rho_mean     = mean(rho, na.rm = TRUE),
    longitude    = first(longitude),
    latitude     = first(latitude),
    koppen_group = first(koppen_group),
    .groups = "drop"
  ) %>%
  filter(n_tp == 6) %>%
  mutate(
    ITW          = map_int(coef_seq, find_itw),
    inhibit_tp0  = (coef_tp0 < 0)
  )

n_inhibit0    <- sum(itw_df$inhibit_tp0, na.rm = TRUE)
n_itw_valid   <- sum(!is.na(itw_df$ITW))

cat(sprintf("tp=0 为抑制的站点: %d (%.1f%%)\n", n_inhibit0, n_inhibit0 / n_total * 100))
cat(sprintf("有效 ITW（抑→促转变）站点: %d (%.1f%% of 抑制站点)\n",
            n_itw_valid, n_itw_valid / max(n_inhibit0, 1) * 100))
cat("\nITW 分布:\n")
print(table(itw_df$ITW, useNA = "ifany"))

# 对比 HTW vs ITW 分布图
p_compare <- bind_rows(
  htw_df %>% filter(!is.na(HTW)) %>% transmute(window = as.integer(HTW), type = "HTW（促→抑）"),
  itw_df %>% filter(!is.na(ITW)) %>% transmute(window = as.integer(ITW), type = "ITW（抑→促）")
) %>%
  count(type, window) %>%
  ggplot(aes(x = factor(window), y = n, fill = type)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = n), position = position_dodge(0.65),
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("HTW（促→抑）" = "steelblue",
                                "ITW（抑→促）" = "tomato"),
                    name = NULL) +
  labs(
    title = "HTW vs ITW 分布对比",
    subtitle = sprintf("HTW N=%d（促→抑）vs ITW N=%d（抑→促）",
                       n_htw_valid, n_itw_valid),
    x = "转变发生的滞后周数", y = "站点数量"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top")

ggsave("data_proc/itw_step1_htw_itw_compare.png", p_compare,
       width = 8, height = 5, dpi = 300)
cat("-> 图已保存: data_proc/itw_step1_htw_itw_compare.png\n")

# ── 合并投资与气候 ───────────────────────────────────────────────────────────
itw_full <- itw_df %>%
  left_join(invest_metrics, by = "meteo_stat_id") %>%
  left_join(climate_bg,    by = "meteo_stat_id") %>%
  filter(!is.na(ITW))

cat("\n合并后有效样本数:", nrow(itw_full), "\n")
cat("其中有投资数据:", sum(!is.na(itw_full$invest_pa_tot)), "\n")

# ── 步骤 ITW-2：断点回归 ─────────────────────────────────────────────────────
cat("\n【ITW 步骤2：断点回归 —— 投资阈值对 ITW】\n")

df_seg_itw <- itw_full %>%
  filter(!is.na(invest_pa_tot), invest_pa_tot > 0) %>%
  mutate(log_invest = log10(invest_pa_tot), ITW_num = as.numeric(ITW))

cat("断点回归样本数:", nrow(df_seg_itw), "\n")

seg_res_itw <- tryCatch({
  lm_itw <- lm(ITW_num ~ log_invest, data = df_seg_itw)
  seg    <- segmented(lm_itw, seg.Z = ~log_invest,
                      npsi = 1, control = seg.control(it.max = 200))
  bp     <- seg$psi[, "Est."]
  bp_se  <- seg$psi[, "St.Err"]
  dv     <- davies.test(lm_itw, ~log_invest, k = 10)
  cat(sprintf("  阈值 T = %.3f (投资强度 = %.2f)，SE = %.3f，Davies p = %.4f\n",
              bp, 10^bp, bp_se, dv$p.value))
  list(seg = seg, bp = bp, bp_se = bp_se, davies_p = dv$p.value)
}, error = function(e) { cat("  断点回归失败:", conditionMessage(e), "\n"); NULL })

# 散点图 + 断点拟合
{
  x_rng <- range(df_seg_itw$log_invest, na.rm = TRUE)
  x_seq <- seq(x_rng[1], x_rng[2], length.out = 300)

  if (!is.null(seg_res_itw)) {
    y_fit  <- predict(seg_res_itw$seg, newdata = data.frame(log_invest = x_seq))
    fit_df <- data.frame(x = x_seq, y = y_fit)
    bp_val <- seg_res_itw$bp
  } else {
    # 退而用 LOESS
    lm_simple <- lm(ITW_num ~ log_invest, data = df_seg_itw)
    y_fit  <- predict(lm_simple, newdata = data.frame(log_invest = x_seq))
    fit_df <- data.frame(x = x_seq, y = y_fit)
    bp_val <- NULL
  }

  p_seg_itw <- ggplot(df_seg_itw, aes(x = log_invest, y = ITW_num)) +
    geom_jitter(aes(color = koppen_group), width = 0.03, height = 0.12,
                alpha = 0.6, size = 2) +
    geom_line(data = fit_df, aes(x = x, y = y),
              color = "black", linewidth = 1.3, inherit.aes = FALSE) +
    { if (!is.null(bp_val))
        list(
          geom_vline(xintercept = bp_val, linetype = "dashed",
                     color = "tomato", linewidth = 1),
          annotate("label", x = bp_val + 0.05,
                   y = max(df_seg_itw$ITW_num, na.rm = TRUE) - 0.3,
                   label = sprintf("阈值 T = %.2f\n投资强度 = %.1f",
                                   bp_val, 10^bp_val),
                   hjust = 0, color = "tomato", size = 3.8,
                   fill = "white", alpha = 0.85)
        )
    } +
    scale_color_brewer(palette = "Set1", name = "柯本气候带") +
    scale_x_continuous(
      sec.axis = sec_axis(~10^., name = "绿化投资强度（原始尺度）",
                          breaks = c(1, 5, 10, 50, 100, 500))
    ) +
    labs(
      title = "绿化投资对抑→促恢复窗口（ITW）的阈值效应",
      subtitle = if (!is.null(seg_res_itw))
        sprintf("断点回归 | Davies p = %.3f | N = %d 站点",
                seg_res_itw$davies_p, nrow(df_seg_itw))
      else
        sprintf("线性拟合（断点不显著）| N = %d 站点", nrow(df_seg_itw)),
      x = "log10(绿化投资强度 invest_pa_tot)",
      y = "ITW（首次稳定转促的滞后周数）"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("data_proc/itw_step2_segmented.png", p_seg_itw,
         width = 9, height = 6, dpi = 300)
  cat("-> 图已保存: data_proc/itw_step2_segmented.png\n")
}

# 箱线图
if (nrow(df_seg_itw) >= 10) {
  q_breaks <- unique(quantile(df_seg_itw$log_invest, probs = 0:4/4, na.rm = TRUE))
  n_q      <- length(q_breaks) - 1
  q_labels <- c("Q1\n（低投资）","Q2","Q3","Q4\n（高投资）")[1:n_q]

  df_box_itw <- df_seg_itw %>%
    mutate(invest_q = cut(log_invest, breaks = q_breaks,
                          labels = q_labels, include.lowest = TRUE)) %>%
    filter(!is.na(invest_q))

  p_box_itw <- ggplot(df_box_itw, aes(x = invest_q, y = ITW_num, fill = invest_q)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 23,
                 size = 4, fill = "white", color = "black") +
    scale_fill_brewer(palette = "Reds", guide = "none") +
    labs(
      title = "绿化投资分位数 vs 抑→促恢复窗口（ITW）",
      subtitle = "◆ 为均值；若投资越高 ITW 越小，则绿化加速了热胁迫后的恢复",
      x = "绿化投资强度分位数",
      y = "ITW（周）"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("data_proc/itw_step2_boxplot.png", p_box_itw,
         width = 7, height = 5, dpi = 300)
  cat("-> 图已保存: data_proc/itw_step2_boxplot.png\n")

  aov_itw <- aov(ITW_num ~ invest_q, data = df_box_itw)
  cat("\nANOVA（投资四分位 vs ITW）:\n")
  print(summary(aov_itw))
}

# ── 步骤 ITW-3：方差分解 ─────────────────────────────────────────────────────
cat("\n【ITW 步骤3：方差分解（vegan::varpart）】\n")

df_vp_itw <- itw_full %>%
  filter(!is.na(invest_pa_tot), invest_pa_tot > 0,
         !is.na(vpd_mean_summer), !is.na(latitude)) %>%
  mutate(
    ITW_num      = as.numeric(ITW),
    log_invest   = log10(invest_pa_tot),
    log_invest_r = log1p(coalesce(invest_ratio, 0))
  )

cat("方差分解样本数:", nrow(df_vp_itw), "\n")

if (nrow(df_vp_itw) >= 20) {
  Y_itw <- as.matrix(df_vp_itw$ITW_num)

  X1i <- df_vp_itw %>%
    dplyr::select(log_invest, log_invest_r) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>% as.matrix()
  X2i <- df_vp_itw %>%
    dplyr::select(vpd_mean_summer, heat_freq_mean, heat_intensity_mean) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>% as.matrix()
  X3i <- df_vp_itw %>%
    dplyr::select(latitude, longitude) %>%
    mutate(across(everything(), ~scale(.)[, 1])) %>% as.matrix()

  vp_itw <- varpart(Y_itw, X1i, X2i, X3i)
  r2i    <- vp_itw$part$indfract$Adj.R.square

  cat("\n方差分解独立贡献（Adjusted R²）:\n")
  cat(sprintf("  绿化投资 (X1): %.3f (%.1f%%)\n", r2i[1], r2i[1] * 100))
  cat(sprintf("  气候背景 (X2): %.3f (%.1f%%)\n", r2i[2], r2i[2] * 100))
  cat(sprintf("  地理位置 (X3): %.3f (%.1f%%)\n", r2i[3], r2i[3] * 100))
  cat(sprintf("  残差（未解释）: %.3f (%.1f%%)\n", r2i[8], r2i[8] * 100))

  # 并排对比两张方差分解柱状图（HTW vs ITW）
  labels_comp <- c("绿化投资\n独立","气候背景\n独立","地理位置\n独立",
                   "投资×气候","投资×地理","气候×地理","三者共享")
  colors_comp <- c("#4DAF4A","#FF7F00","#377EB8",
                   "#A6D96A","#ABD9E9","#FDAE61","#D9EF8B")

  vp_cmp <- bind_rows(
    data.frame(window = "HTW（促→抑）", component = labels_comp,
               r2 = pmax(r2[1:7], 0), color = colors_comp),
    data.frame(window = "ITW（抑→促）", component = labels_comp,
               r2 = pmax(r2i[1:7], 0), color = colors_comp)
  ) %>%
    mutate(
      component = factor(component, levels = labels_comp),
      label     = ifelse(r2 > 0.005, sprintf("%.1f%%", r2 * 100), "")
    )

  p_vp_cmp <- ggplot(vp_cmp, aes(x = component, y = r2, fill = component)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.8) +
    facet_wrap(~window, ncol = 2) +
    scale_fill_manual(values = colors_comp) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "HTW vs ITW 方差分解对比",
      subtitle = "各成分对转变窗口变异的独立贡献（Adjusted R²）",
      x = NULL, y = "解释方差"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 9, lineheight = 1.1),
      strip.text  = element_text(face = "bold", size = 12)
    )

  ggsave("data_proc/itw_step3_varpart_compare.png", p_vp_cmp,
         width = 12, height = 5, dpi = 300)
  cat("-> 图已保存: data_proc/itw_step3_varpart_compare.png\n")
}

# ── 补充：有序逻辑回归（ITW）──────────────────────────────────────────────
cat("\n【ITW 补充：有序逻辑回归】\n")

if (nrow(df_vp_itw) >= 20) {
  df_olr_itw <- df_vp_itw %>%
    mutate(
      ITW_ord      = factor(ITW, levels = sort(unique(ITW)), ordered = TRUE),
      log_invest_s = scale(log_invest)[, 1],
      vpd_s        = scale(vpd_mean_summer)[, 1],
      heat_s       = scale(heat_freq_mean)[, 1],
      lat_s        = scale(latitude)[, 1]
    )

  olr_itw <- tryCatch(
    polr(ITW_ord ~ log_invest_s + vpd_s + heat_s + lat_s,
         data = df_olr_itw, Hess = TRUE),
    error = function(e) { cat("有序回归失败:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(olr_itw)) {
    ct   <- coef(summary(olr_itw))
    pv   <- pnorm(abs(ct[, "t value"]), lower.tail = FALSE) * 2
    cat("\n有序逻辑回归结果（标准化系数）:\n")
    print(round(cbind(as.data.frame(ct), p_value = pv)[1:4, ], 4))

    beta_inv_itw <- coef(olr_itw)["log_invest_s"]
    cat(sprintf("\n绿化投资系数 = %.3f → %s\n", beta_inv_itw,
                ifelse(beta_inv_itw < 0,
                       "投资越高，ITW 越小（恢复更快）✓",
                       "投资越高，ITW 越大（恢复更慢，方向不符）✗")))
  }
}

# =============================================================================
# 最终汇总（HTW + ITW）
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("【完整分析汇总：HTW vs ITW】\n")
cat(strrep("=", 60), "\n")
cat(sprintf("HTW（促→抑）有效样本: %d 站点\n", n_htw_valid))
cat(sprintf("ITW（抑→促）有效样本: %d 站点\n", n_itw_valid))
if (exists("seg_res") && !is.null(seg_res))
  cat(sprintf("HTW 断点: log10(T)=%.3f，Davies p=%.4f\n",
              seg_res$bp, seg_res$davies_p))
if (exists("seg_res_itw") && !is.null(seg_res_itw))
  cat(sprintf("ITW 断点: log10(T)=%.3f，Davies p=%.4f\n",
              seg_res_itw$bp, seg_res_itw$davies_p))
if (exists("r2") && exists("r2i"))
  cat(sprintf("方差分解（投资独立贡献）—— HTW: %.1f%% | ITW: %.1f%%\n",
              r2[1]*100, r2i[1]*100))
cat(strrep("=", 60), "\n")
