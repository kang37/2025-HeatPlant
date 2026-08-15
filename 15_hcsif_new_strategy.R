#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 15_hcsif_new_strategy.R
#
# 用新的 HCSIF (500m, 8天) VPD->SIF 因果结果，按"新策略"重做类似 05_1 的分析：
#   问题1（属于哪种模式）：两阶段 —— 起始方向二元 logit + 随机森林分类
#   问题2（反转快慢）      ：生存分析 —— person-period 展开 + 离散时间风险模型
#                            + Kaplan-Meier / 竞争风险描述
#
# 时间轴 = CCM 滞后 tp = 0..8（每步 8 天）；事件 = 因果系数符号翻转；
# always_ 类 = 右删失（8 步内未翻转）。
#
# ★ 切换到全量结果：只改下面 CCM_RDS 一行即可（其余自动适配）。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(stringr)
  library(survival)
  library(ranger)
  library(showtext)
})
showtext_auto(); showtext_opts(dpi = 300)

# ===========================================================================
# 0. 配置（切换全量数据时改这里）
# ===========================================================================
PROJ    <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
CCM_RDS <- file.path(PROJ, "data_proc/ccm_hcsif_500m/ccm_hcsif_vpd_20260726_2018.rds")  # 全量(901站)
COV_RDS <- file.path(PROJ, "data_proc/output_10y_built_up_05_01/station_covariates.rds")
OUT     <- file.path(PROJ, "data_proc/output_hcsif_new")
SIF_Y   <- "SIF_buf750_dt"     # 主 SIF 定义：750m 邻域均值
VPD_X   <- "vpd_mean_dt"       # 驱动：VPD

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_cn <- function() theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
winsorize <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

cat("=== 15_hcsif_new_strategy ===\n")
cat("CCM :", basename(CCM_RDS), "\n")
cat("关系:", VPD_X, "->", SIF_Y, "\n\n")

# ===========================================================================
# 1. 读取新 CCM，抽取 VPD->SIF(buf750) 的系数序列
# ===========================================================================
ccm <- readRDS(CCM_RDS) %>% as.data.frame()
ccm <- ccm %>% filter(y_var == SIF_Y, x_var == VPD_X)
N_TP     <- length(unique(ccm$tp))
TRRI_MAX <- 2L * N_TP
cat(sprintf("站点数: %d | tp: %d..%d (共%d步) | 行: %d\n",
            length(unique(ccm$meteo_stat)),
            min(ccm$tp), max(ccm$tp), N_TP, nrow(ccm)))

# ===========================================================================
# 2. 模式分类（复用 05_1 的 find_transition 逻辑，仅用符号）
# ===========================================================================
find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

pat <- ccm %>%
  arrange(meteo_stat, tp) %>%
  group_by(meteo_stat) %>%
  summarise(coef_seq  = list(mean_coef),
            longitude = first(longitude),
            latitude  = first(latitude),
            .groups   = "drop") %>%
  filter(map_lgl(coef_seq, ~length(.x) == N_TP)) %>%
  mutate(
    tp0_pos = map_lgl(coef_seq, ~!is.na(.x[[1]]) && .x[[1]] > 0),
    HTW     = map_int(coef_seq, ~find_transition(.x, to_neg = TRUE)),   # 促->抑 翻转tp
    ITW     = map_int(coef_seq, ~find_transition(.x, to_neg = FALSE)),  # 抑->促 翻转tp
    stype   = case_when(
      tp0_pos  & is.na(HTW)  ~ "always_promote",
      !tp0_pos & is.na(ITW)  ~ "always_inhibit",
      tp0_pos  & !is.na(HTW) ~ "promote_inhibit",
      !tp0_pos & !is.na(ITW) ~ "inhibit_promote",
      TRUE                   ~ "always_inhibit"
    ),
    # 起始方向（问题1第一阶段的二元结果）
    start_dir = ifelse(tp0_pos, "promote_first", "inhibit_first"),
    # 生存框架：事件是否发生 + 事件发生的 tp（未发生则 NA）
    flipped   = stype %in% c("inhibit_promote", "promote_inhibit"),
    event_tp  = case_when(
      stype == "inhibit_promote" ~ ITW,
      stype == "promote_inhibit" ~ HTW,
      TRUE                       ~ NA_integer_
    )
  )

cat("\n模式分布:\n"); print(count(pat, stype))
cat("\n起始方向分布:\n"); print(count(pat, start_dir))

# ===========================================================================
# 3. 合并协变量（按 meteo_stat 复用 05_1 站点级协变量）
# ===========================================================================
cov <- readRDS(COV_RDS) %>% rename(meteo_stat = meteo_stat_id) %>%
  mutate(meteo_stat = as.character(meteo_stat))
dat <- pat %>%
  mutate(meteo_stat = as.character(meteo_stat)) %>%
  left_join(cov, by = "meteo_stat", suffix = c("", ".cov")) %>%
  mutate(
    koppen_group = ifelse(is.na(koppen_group), NA_character_, koppen_group),
    # 截尾管理/局地/背景变量
    inv_w    = winsorize(pa_built_10y),
    cgi_w    = winsorize(cgi_score),
    precip_w = winsorize(precip_mean),
    pop_w    = winsorize(pop_10y),
    road_w   = winsorize(road_density),
    bvol_w   = winsorize(building_vol_density),
    lon      = longitude,
    lat      = latitude
  )
cat(sprintf("\n合并协变量后：%d 站，其中有投资 %d，有降水 %d，有气候区 %d\n",
            nrow(dat), sum(!is.na(dat$inv_w)),
            sum(!is.na(dat$precip_w)), sum(!is.na(dat$koppen_group))))
write_csv(dat %>% select(-coef_seq), file.path(OUT, "hcsif_station_patterns.csv"))
cat("-> hcsif_station_patterns.csv\n")

# 预测变量集合（管理 + 局地 + 背景 + 地理）
pred_vars <- c("inv_w", "cgi_w", "precip_w", "pop_w", "road_w", "bvol_w", "lon", "lat")
var_labels <- c(inv_w="Investment", cgi_w="CGI", precip_w="Precip", pop_w="Population",
                road_w="Road", bvol_w="Bldg Vol.Density", lon="Longitude", lat="Latitude")

# 小样本安全检查
enough <- function(df, min_n = 12) nrow(df) >= min_n && n_distinct(df$y) >= 2

# ===========================================================================
# 4. 问题1 —— 两阶段模式分析
# ===========================================================================
cat("\n=== 问题1：模式的影响因素 ===\n")

# 4a. 第一阶段：起始方向二元 logit（inhibit_first=1）
s1 <- dat %>%
  mutate(y = as.integer(start_dir == "inhibit_first")) %>%
  select(y, all_of(pred_vars)) %>% drop_na()

if (enough(s1)) {
  f1 <- glm(y ~ ., data = s1, family = binomial)
  s1_tab <- broom_like <- {
    co <- summary(f1)$coefficients
    tibble(variable = rownames(co), coef = co[,1], se = co[,2], p_val = co[,4]) %>%
      filter(variable != "(Intercept)")
  }
  write_csv(s1_tab, file.path(OUT, "stage1_startdir_logit.csv"))
  cat("-> stage1_startdir_logit.csv （起始方向 logit）\n")

  p_s1 <- s1_tab %>%
    mutate(var_label = factor(var_labels[variable], levels = rev(var_labels[pred_vars])),
           sig = p_val < 0.05) %>%
    ggplot(aes(x = coef, y = var_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se), height = 0.25) +
    geom_point(aes(shape = sig), size = 3) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(title = "Stage 1: What drives starting direction (inhibit-first vs promote-first)",
         subtitle = "Logit coef > 0 → more likely inhibit-first; filled = p<0.05",
         x = "Logit coefficient (95% CI)", y = NULL) +
    theme_cn()
  ggsave(file.path(OUT, "stage1_startdir_logit.png"), p_s1, width = 9, height = 5, dpi = 300)
  cat("-> stage1_startdir_logit.png\n")
} else {
  cat("[跳过] 起始方向 logit：样本不足（试跑数据正常）\n")
}

# 4b. 随机森林：4 模式分类 + 变量重要性
rf_df <- dat %>% select(stype, all_of(pred_vars)) %>% drop_na() %>%
  mutate(stype = factor(stype))
if (nrow(rf_df) >= 12 && n_distinct(rf_df$stype) >= 2) {
  rf <- ranger(stype ~ ., data = rf_df, num.trees = 500,
               importance = "permutation", probability = TRUE, seed = 42)
  imp <- tibble(variable = names(rf$variable.importance),
                importance = as.numeric(rf$variable.importance))
  write_csv(imp, file.path(OUT, "pattern_rf_importance.csv"))
  cat("-> pattern_rf_importance.csv （模式RF重要性, OOB brier:",
      round(rf$prediction.error, 3), ")\n")

  p_rf <- imp %>%
    mutate(var_label = reorder(var_labels[variable], importance)) %>%
    ggplot(aes(x = importance, y = var_label)) +
    geom_col(fill = "#4575B4", alpha = 0.85) +
    labs(title = "Pattern classification (RF): permutation importance",
         subtitle = sprintf("4 patterns; %d stations; OOB Brier = %.3f",
                            nrow(rf_df), rf$prediction.error),
         x = "Permutation importance", y = NULL) +
    theme_cn()
  ggsave(file.path(OUT, "pattern_rf_importance.png"), p_rf, width = 9, height = 5, dpi = 300)
  cat("-> pattern_rf_importance.png\n")
} else {
  cat("[跳过] 模式RF：样本不足（试跑数据正常）\n")
}

# ===========================================================================
# 5. 问题2 —— 生存分析（反转快慢）
# ===========================================================================
cat("\n=== 问题2：反转时机的生存分析 ===\n")

# 5a. 构造生存数据：按起始方向分层
#     time  = 事件发生的 tp（删失则 = 最大 tp）
#     event = 是否观测到翻转
surv_df <- dat %>%
  mutate(
    time  = ifelse(flipped, event_tp, max(ccm$tp)),
    event = as.integer(flipped)
  )

# 5b. person-period 长表展开（离散时间风险模型的关键结构）
#     每个站点在 tp=1..(事件tp或删失tp) 各一行，event 指示该 tp 是否翻转
expand_person_period <- function(df) {
  df %>%
    filter(!is.na(time), time >= 1) %>%
    rowwise() %>%
    do({
      row <- .
      tps <- 1:row$time
      tibble(meteo_stat = row$meteo_stat,
             start_dir  = row$start_dir,
             tp         = tps,
             event      = as.integer(tps == row$time & row$event == 1)) %>%
        bind_cols(row[pred_vars])
    }) %>%
    ungroup()
}
pp <- expand_person_period(surv_df)
write_csv(pp, file.path(OUT, "person_period.csv"))
cat(sprintf("-> person_period.csv （%d 行 = 站点×滞后步）\n", nrow(pp)))

# 5c. Kaplan-Meier / 累积翻转曲线（分起始方向，描述性）
km_df <- surv_df %>% filter(!is.na(time), time >= 1)
if (nrow(km_df) >= 8) {
  km <- survfit(Surv(time, event) ~ start_dir, data = km_df)
  km_tidy <- tibble(
    time     = km$time,
    surv     = km$surv,
    cum_flip = 1 - km$surv,
    strata   = rep(names(km$strata), km$strata)
  )
  write_csv(km_tidy, file.path(OUT, "km_transition.csv"))

  p_km <- km_tidy %>%
    ggplot(aes(x = time, y = cum_flip, color = strata)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = c("start_dir=inhibit_first" = "#3182BD",
                                  "start_dir=promote_first" = "#E6550D"),
                       labels = c("Inhibit-first (→flip to +)",
                                  "Promote-first (→flip to −)"),
                       name = "Starting direction") +
    labs(title = "Cumulative sign-flip probability along CCM lag",
         subtitle = "Survival view: 'time' = lag step (8-day); event = coefficient sign flip",
         x = "CCM lag tp (each = 8 days)", y = "Cumulative fraction flipped") +
    theme_cn()
  ggsave(file.path(OUT, "km_transition.png"), p_km, width = 9, height = 5, dpi = 300)
  cat("-> km_transition.png / km_transition.csv\n")
} else {
  cat("[跳过] KM 曲线：样本不足\n")
}

# 5d. 离散时间风险模型：分起始方向各拟合一个（cause-specific hazard）
fit_discrete_hazard <- function(pp_sub, label, tag) {
  d <- pp_sub %>% select(event, tp, all_of(pred_vars)) %>% drop_na()
  if (nrow(d) < 15 || sum(d$event) < 3) {
    cat(sprintf("[跳过] 离散时间风险(%s)：事件/样本不足（试跑数据正常）\n", label))
    return(invisible(NULL))
  }
  # tp 作为基线风险（因子；样本极少时退化为线性）
  form <- if (n_distinct(d$tp) >= 3 && nrow(d) > 40)
    as.formula(paste("event ~ factor(tp) +", paste(pred_vars, collapse = " + ")))
  else
    as.formula(paste("event ~ tp +", paste(pred_vars, collapse = " + ")))
  fit <- tryCatch(glm(form, data = d, family = binomial), error = function(e) NULL)
  if (is.null(fit)) { cat(sprintf("[跳过] 离散时间风险(%s)：拟合失败\n", label)); return(invisible(NULL)) }

  co <- summary(fit)$coefficients
  tab <- tibble(variable = rownames(co), coef = co[,1], se = co[,2], p_val = co[,4]) %>%
    filter(variable %in% pred_vars)
  write_csv(tab, file.path(OUT, sprintf("hazard_%s.csv", tag)))
  cat(sprintf("-> hazard_%s.csv （%s；事件=%d/%d）\n", tag, label, sum(d$event), nrow(d)))

  p <- tab %>%
    mutate(var_label = factor(var_labels[variable], levels = rev(var_labels[pred_vars])),
           sig = p_val < 0.05) %>%
    ggplot(aes(x = coef, y = var_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef - 1.96*se, xmax = coef + 1.96*se), height = 0.25) +
    geom_point(aes(shape = sig), size = 3) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(title = sprintf("Discrete-time hazard of transition: %s", label),
         subtitle = "coef > 0 → higher per-step hazard → EARLIER flip; filled = p<0.05",
         x = "Log-hazard coefficient (95% CI)", y = NULL) +
    theme_cn()
  ggsave(file.path(OUT, sprintf("hazard_%s.png", tag)), p, width = 9, height = 5, dpi = 300)
  cat(sprintf("-> hazard_%s.png\n", tag))
}

fit_discrete_hazard(filter(pp, start_dir == "inhibit_first"),
                    "Inhibit-first → flip to promote (ITW)", "inhibit_first")
fit_discrete_hazard(filter(pp, start_dir == "promote_first"),
                    "Promote-first → flip to inhibit (HTW)", "promote_first")

cat("\n", strrep("=", 60), "\n")
cat("15_hcsif_new_strategy.R 完成\n")
cat("输出目录:", OUT, "\n")
for (f in sort(list.files(OUT, pattern = "\\.(csv|png)$"))) cat("  ", f, "\n")
