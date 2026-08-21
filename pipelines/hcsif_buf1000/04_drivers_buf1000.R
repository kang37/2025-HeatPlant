#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04_drivers_buf1000.R  （仅 1000m 方案）
#   Q1: 什么因素决定一个站点被分入哪种模式
#   Q2: 哪个因素最重要
#   Q3: 2分类下(抑制先行 = 全抑制+先抑后促; 促进先行 = 全促进+先促后抑)，
#       什么因素决定其转变速度快慢、如何决定
#   —— 全部用「城市聚类稳健标准误」得到诚实 p 值（投资/CGI 为城市级变量）
# 全英文标签，避免非交互 Rscript 中文缺字形。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(ranger); library(sandwich); library(lmtest)
})

PROJ    <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
CCM_DIR <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000")
COV_RDS <- file.path(PROJ, "data_proc/output_10y_built_up_05_01/station_covariates.rds")
OUT     <- file.path(PROJ, "data_proc/output_hcsif_buf1000")
SIF_Y   <- "SIF_buf1000_dt"; VPD_X <- "vpd_mean_dt"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
theme_cn <- function() theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
winsorize <- function(x, p = 0.01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

CCM_RDS <- {
  fs <- list.files(CCM_DIR, pattern = "^ccm_hcsif_buf1000_vpd_.*[0-9]\\.rds$", full.names = TRUE)
  tail(sort(fs[!grepl("_raw\\.rds$", fs)]), 1)
}
cat("CCM:", basename(CCM_RDS), "\n\n")

# ===========================================================================
# 1. 模式分类（宽松"多数持续"规则：每个站点都归入4模式之一，无 other）
# ===========================================================================
ccm <- readRDS(CCM_RDS) %>% as.data.frame() %>% filter(y_var == SIF_Y, x_var == VPD_X)
N_TP <- length(unique(ccm$tp))

find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

pat <- ccm %>% arrange(meteo_stat, tp) %>% group_by(meteo_stat) %>%
  summarise(coef_seq = list(mean_coef), .groups = "drop") %>%
  filter(map_lgl(coef_seq, ~length(.x) == N_TP)) %>%
  mutate(
    tp0_pos = map_lgl(coef_seq, ~!is.na(.x[[1]]) && .x[[1]] > 0),
    HTW = map_int(coef_seq, ~find_transition(.x, TRUE)),
    ITW = map_int(coef_seq, ~find_transition(.x, FALSE)),
    stype = case_when(
      tp0_pos  & is.na(HTW)  ~ "always_promote",
      !tp0_pos & is.na(ITW)  ~ "always_inhibit",
      tp0_pos  & !is.na(HTW) ~ "promote_inhibit",
      !tp0_pos & !is.na(ITW) ~ "inhibit_promote",
      TRUE                   ~ "always_inhibit"),
    # 2分类（用户定义）
    class2 = ifelse(stype %in% c("always_inhibit","inhibit_promote"),
                    "Inhibit-first", "Promote-first"),
    flipped  = stype %in% c("inhibit_promote","promote_inhibit"),
    event_tp = case_when(stype=="inhibit_promote" ~ ITW,
                         stype=="promote_inhibit" ~ HTW, TRUE ~ NA_integer_)
  )

# ===========================================================================
# 2. 合并协变量
# ===========================================================================
cov <- readRDS(COV_RDS) %>% rename(meteo_stat = meteo_stat_id) %>%
  mutate(meteo_stat = as.character(meteo_stat))
dat <- pat %>% mutate(meteo_stat = as.character(meteo_stat)) %>%
  left_join(cov, by = "meteo_stat") %>%
  mutate(inv_w=winsorize(pa_built_10y), cgi_w=winsorize(cgi_score),
         precip_w=winsorize(precip_mean), pop_w=winsorize(pop_10y),
         road_w=winsorize(road_density), bvol_w=winsorize(building_vol_density),
         lon=longitude, lat=latitude, city=city_name)

pred_vars  <- c("inv_w","cgi_w","precip_w","pop_w","road_w","bvol_w","lon","lat")
var_labels <- c(inv_w="Investment", cgi_w="CGI governance", precip_w="Precipitation",
                pop_w="Population", road_w="Road density", bvol_w="Building vol.density",
                lon="Longitude", lat="Latitude")

# 城市聚类稳健系数表
clustered_coef <- function(fit, keep = pred_vars) {
  ct <- lmtest::coeftest(fit, vcov = sandwich::vcovCL, cluster = ~ city)
  tibble(variable = rownames(ct), coef = ct[,1], se = ct[,2], p_val = ct[,4]) %>%
    filter(variable %in% keep)
}
forest <- function(tab, title, sub, xlab) {
  tab %>% mutate(var_label = factor(var_labels[variable], levels = rev(var_labels[pred_vars])),
                 sig = p_val < 0.05) %>%
    ggplot(aes(coef, var_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = coef-1.96*se, xmax = coef+1.96*se), height = 0.25) +
    geom_point(aes(shape = sig), size = 3) +
    scale_shape_manual(values = c(`TRUE`=16,`FALSE`=1), guide = "none") +
    labs(title = title, subtitle = sub, x = xlab, y = NULL) + theme_cn()
}

# ===========================================================================
# Q1 + Q2: 什么决定模式 / 最重要因素
# ===========================================================================
cat("=== Q1/Q2: 模式的影响因素 ===\n")
# (a) 随机森林 4 模式分类 → 排列重要性（回答"最重要"）
rf_df <- dat %>% select(stype, all_of(pred_vars)) %>% drop_na() %>% mutate(stype = factor(stype))
rf <- ranger(stype ~ ., data = rf_df, num.trees = 1000, importance = "permutation",
             probability = TRUE, seed = 42)
imp <- tibble(variable = names(rf$variable.importance),
              importance = as.numeric(rf$variable.importance)) %>% arrange(desc(importance))
write_csv(imp, file.path(OUT, "drivers_pattern_rf_importance.csv"))
cat("模式RF重要性(降序), OOB Brier =", round(rf$prediction.error,3), ":\n"); print(as.data.frame(imp))

p_imp <- imp %>% mutate(vl = reorder(var_labels[variable], importance)) %>%
  ggplot(aes(importance, vl)) + geom_col(fill = "#4575B4", alpha = .85) +
  labs(title = "Q1/Q2: What determines the pattern (RF permutation importance)",
       subtitle = sprintf("4 modes; n=%d; OOB Brier=%.3f; top factor = %s",
                          nrow(rf_df), rf$prediction.error, var_labels[imp$variable[1]]),
       x = "Permutation importance", y = NULL) + theme_cn()
ggsave(file.path(OUT, "drivers_pattern_rf_importance.png"), p_imp, width = 9, height = 5, dpi = 300)

# (b) 起始方向(=2分类) logit，城市聚类稳健SE → 给方向与显著性
s1 <- dat %>% mutate(y = as.integer(class2 == "Inhibit-first")) %>%
  select(y, city, all_of(pred_vars)) %>% drop_na()
f1 <- glm(y ~ . -city, data = s1, family = binomial)
s1_tab <- clustered_coef(f1)
write_csv(s1_tab, file.path(OUT, "drivers_startdir_logit_clustered.csv"))
cat("\n起始方向(2分类) logit[城市聚类]，coef>0 → 更可能 Inhibit-first:\n"); print(as.data.frame(s1_tab))
ggsave(file.path(OUT, "drivers_startdir_logit_clustered.png"),
       forest(s1_tab, "Q1: Drivers of class membership (Inhibit-first vs Promote-first)",
              "Logit + city-clustered SE; coef>0 -> Inhibit-first; filled = p<0.05",
              "Logit coefficient (95% CI)"), width = 9, height = 5, dpi = 300)

# ===========================================================================
# Q3: 2分类下，什么决定转变速度（城市聚类稳健离散时间风险模型）
# ===========================================================================
cat("\n=== Q3: 转变速度的影响因素（城市聚类）===\n")
surv <- dat %>% mutate(time = ifelse(flipped, event_tp, max(ccm$tp)), event = as.integer(flipped))
pp <- surv %>% filter(!is.na(time), time >= 1) %>%
  rowwise() %>% do({ r <- .; tibble(meteo_stat=r$meteo_stat, class2=r$class2, city=r$city,
      tp=1:r$time, event=as.integer((1:r$time)==r$time & r$event==1)) %>%
      bind_cols(r[pred_vars]) }) %>% ungroup()

fit_haz <- function(cl, label, tag) {
  d <- pp %>% filter(class2 == cl) %>% select(event, tp, city, all_of(pred_vars)) %>% drop_na()
  form <- if (n_distinct(d$tp) >= 3)
    as.formula(paste("event ~ factor(tp) +", paste(pred_vars, collapse=" + ")))
  else as.formula(paste("event ~ tp +", paste(pred_vars, collapse=" + ")))
  fit <- glm(form, data = d, family = binomial)
  tab <- clustered_coef(fit)
  write_csv(tab, file.path(OUT, sprintf("drivers_hazard_%s_clustered.csv", tag)))
  cat(sprintf("\n[%s] 事件=%d/%d 站行; coef>0 → 每步翻转风险更高 → 转变更早 [城市聚类]:\n",
              label, sum(d$event), nrow(d))); print(as.data.frame(tab))
  ggsave(file.path(OUT, sprintf("drivers_hazard_%s_clustered.png", tag)),
         forest(tab, sprintf("Q3: Transition-speed drivers — %s", label),
                "Discrete-time hazard + city-clustered SE; coef>0 -> earlier transition; filled = p<0.05",
                "Log-hazard coefficient (95% CI)"), width = 9, height = 5, dpi = 300)
  tab %>% mutate(class = label)
}
t_inh <- fit_haz("Inhibit-first", "Inhibit-first -> flip to Promote", "inhibit_first")
t_pro <- fit_haz("Promote-first", "Promote-first -> flip to Inhibit", "promote_first")

# 合并的双面板对比图
bind_rows(t_inh, t_pro) %>%
  mutate(var_label = factor(var_labels[variable], levels = rev(var_labels[pred_vars])),
         sig = p_val < 0.05) %>%
  ggplot(aes(coef, var_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = coef-1.96*se, xmax = coef+1.96*se), height = .25) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c(`TRUE`=16,`FALSE`=1), guide = "none") +
  facet_wrap(~class) +
  labs(title = "Q3: Transition-speed drivers by 2-class (city-clustered)",
       subtitle = "coef>0 -> higher per-step hazard -> EARLIER transition; filled = p<0.05",
       x = "Log-hazard coefficient (95% CI)", y = NULL) + theme_cn()
ggsave(file.path(OUT, "drivers_hazard_2class_clustered.png"), width = 12, height = 5, dpi = 300)

cat("\n完成. 输出前缀 drivers_ 于", OUT, "\n")
