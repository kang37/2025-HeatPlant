# =============================================================================
# 02_trri_invest.R  (更新版 — tp自动适配, 18级TRRI, 分组分析)
# 输入: ccm_vpd_results.rds  (根目录，VPD均值 & VPD超量累积, tp=0..8)
# 流程:
#   A. TRRI 分类（自动检测N_TP，动态生成级别）
#   B. 投资变量（多年均值、斜率；行位置对齐）
#   C. OLS回归 + 方差分解（全部/B/C/D × 7投资变量 × 2X变量）
#   D. 分组分析（always_inhibit/always_promote=逻辑回归；
#                inhibit_promote→ITW/promote_inhibit→HTW=OLS）
#   E. 输出图表 + CSV
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr,
               vegan, tibble, scales, showtext, sysfonts, targets)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_vpd"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 字体 -------------------------------------------------------------------------
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 16
theme_cn <- function(bs = BS) {
  theme_minimal(base_size = bs) +
    theme(text          = element_text(family = "heiti"),
          plot.title    = element_text(face = "bold", hjust = 0.5, size = bs + 4),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = bs),
          strip.text    = element_text(face = "bold", size = bs),
          axis.text     = element_text(size = bs - 2),
          axis.title    = element_text(size = bs),
          legend.text   = element_text(size = bs - 2),
          legend.title  = element_text(face = "bold", size = bs))
}
theme_set(theme_cn())

# =============================================================================
# A.1  读取 CCM 结果，自动检测 N_TP
# =============================================================================
ccm <- readRDS("ccm_vpd_results.rds")
cat("CCM行数:", nrow(ccm),
    "| 站点数:", n_distinct(ccm$meteo_stat_id),
    "| X变量:", paste(unique(ccm$x_var), collapse = ", "),
    "| tp范围:", min(ccm$tp), "~", max(ccm$tp), "\n")

N_TP   <- n_distinct(ccm$tp)          # 例：tp=0..8 → N_TP=9
MAX_TP <- max(ccm$tp)
TRRI_MAX <- 2L * N_TP                 # 例：18
cat(sprintf("N_TP=%d, TRRI范围: 1..%d\n", N_TP, TRRI_MAX))

# =============================================================================
# A.2  TRRI 分类函数（动态N_TP）
# =============================================================================
find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

make_trri <- function(df_xvar, n_tp) {
  df_xvar %>%
    arrange(meteo_stat_id, tp) %>%
    group_by(meteo_stat_id) %>%
    summarise(
      coef_seq     = list(mean_coef),
      longitude    = first(longitude),
      latitude     = first(latitude),
      koppen_class = first(koppen_class),
      koppen_group = first(koppen_group),
      .groups      = "drop"
    ) %>%
    filter(map_lgl(coef_seq, ~length(.x) == n_tp)) %>%   # 必须有完整序列
    mutate(
      tp0_pos = map_lgl(coef_seq, ~!is.na(.x[[1]]) && .x[[1]] > 0),
      HTW     = map_int(coef_seq, ~find_transition(.x, to_neg = TRUE)),
      ITW     = map_int(coef_seq, ~find_transition(.x, to_neg = FALSE)),
      stype   = case_when(
        tp0_pos  & is.na(HTW)  ~ "always_promote",
        !tp0_pos & is.na(ITW)  ~ "always_inhibit",
        tp0_pos  & !is.na(HTW) ~ "promote_inhibit",
        !tp0_pos & !is.na(ITW) ~ "inhibit_promote",
        TRUE                   ~ "always_inhibit"
      ),
      TRRI = case_when(
        stype == "always_inhibit"  ~ 1L,
        stype == "inhibit_promote" ~ as.integer(1L + (n_tp - ITW)),
        stype == "promote_inhibit" ~ as.integer(n_tp + HTW),
        stype == "always_promote"  ~ as.integer(2L * n_tp)
      )
    )
}

# 对两个 X 变量分别分类
trri_list <- ccm %>%
  group_by(x_var) %>%
  group_split() %>%
  setNames(sort(unique(ccm$x_var))) %>%
  map(~make_trri(.x, N_TP))

cat("\n=== 站点分类分布 ===\n")
imap(trri_list, ~{ cat("[", .y, "]\n"); print(count(.x, stype)); cat("\n") })

# =============================================================================
# B.  投资变量（行位置对齐）
# =============================================================================
tar_load(green_invest_2020)
tar_load(gdp_data_2020)
tar_load(green_area_2020)
tar_load(station_city_map)

invest_raw0 <- read.csv("data_raw/green_invest/city_invest_data.csv",
                        check.names = FALSE)
nc_inv <- ncol(invest_raw0)
names(invest_raw0) <- c("city_raw", paste0("inv_", 2002:(2001 + nc_inv - 1)))

invest_raw <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  dplyr::select(invest_raw0, -city_raw)
)

year_cols <- paste0("inv_", 2002:2024)
year_cols <- year_cols[year_cols %in% names(invest_raw)]
years_num <- as.integer(str_extract(year_cols, "\\d{4}$"))

get_mean_inv <- function(y1, y2) {
  cols <- year_cols[years_num >= y1 & years_num <= y2]
  rowMeans(dplyr::select(invest_raw, all_of(cols)) %>%
             mutate(across(everything(), as.numeric)), na.rm = TRUE)
}

slope_mat <- dplyr::select(invest_raw,
                            all_of(year_cols[years_num %in% 2005:2020])) %>%
  mutate(across(everything(), as.numeric))
yr_sub <- years_num[years_num %in% 2005:2020]
inv_slope <- apply(slope_mat, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 5) return(NA_real_)
  coef(lm(x[ok] ~ yr_sub[ok]))[2]
})

invest_multi <- invest_raw %>%
  mutate(
    inv_2020  = as.numeric(inv_2020),
    inv_5y    = get_mean_inv(2016, 2020),
    inv_10y   = get_mean_inv(2011, 2020),
    inv_lag10 = get_mean_inv(2005, 2015),
    inv_slope = inv_slope
  ) %>%
  dplyr::select(city_name, inv_2020, inv_5y, inv_10y, inv_lag10, inv_slope)

gdp <- gdp_data_2020 %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))
green <- green_area_2020 %>%
  dplyr::select(city_name,
                area_built = area_green_built,
                area_park  = area_green_park) %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))

invest_full <- station_city_map %>%
  left_join(invest_multi, by = "city_name") %>%
  left_join(gdp,          by = "city_name") %>%
  left_join(green,        by = "city_name") %>%
  mutate(
    ratio_2020   = inv_2020  / gdp * 100,
    ratio_5y     = inv_5y    / gdp * 100,
    ratio_10y    = inv_10y   / gdp * 100,
    ratio_lag10  = inv_lag10 / gdp * 100,
    pa_built_10y = inv_10y   / area_built,
    pa_park_10y  = inv_10y   / area_park
  )

# =============================================================================
# 气候背景变量
# =============================================================================
tar_load(data_heat_sif_weekly)
climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(vpd_bg       = mean(vpd_mean,       na.rm = TRUE),
            heat_freq_bg = mean(heat_event_freq, na.rm = TRUE),
            .groups = "drop")

# 合并
anal_list <- map(trri_list, ~{
  .x %>%
    left_join(invest_full, by = "meteo_stat_id") %>%
    left_join(climate_bg,  by = "meteo_stat_id") %>%
    filter(!is.na(TRRI), !is.na(koppen_group))
})

cat("\n各 X 变量分析数据框行数:\n")
map_int(anal_list, nrow) %>% print()

# =============================================================================
# C.  OLS 回归 + 方差分解：7投资变量 × 2X变量 × 4气候组
# =============================================================================
inv_vars <- c("ratio_2020", "ratio_5y", "ratio_10y", "ratio_lag10",
              "inv_slope",  "pa_built_10y", "pa_park_10y")
inv_labels <- c(
  ratio_2020   = "2020单年/GDP(%)",
  ratio_5y     = "近5年均值/GDP(%)",
  ratio_10y    = "近10年均值/GDP(%)",
  ratio_lag10  = "滞后10年均值/GDP(%)",
  inv_slope    = "投资增长斜率",
  pa_built_10y = "10年均值/建成区绿地",
  pa_park_10y  = "10年均值/公园绿地"
)
ctrl_vars  <- c("vpd_bg", "heat_freq_bg", "latitude")
groups_ols <- c("ALL", "B", "C", "D")

run_ols <- function(df, inv, outcome = "TRRI") {
  d <- df %>%
    filter(!is.na(.data[[inv]]),
           if_all(all_of(ctrl_vars), ~!is.na(.x)),
           !is.na(.data[[outcome]]),
           .data[[inv]] > 0,
           .data[[inv]] < quantile(.data[[inv]], 0.99, na.rm = TRUE))
  if (nrow(d) < 15) return(NULL)
  fml <- as.formula(paste(outcome, "~", inv, "+",
                          paste(ctrl_vars, collapse = "+")))
  m   <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s  <- summary(m)
  cr <- s$coefficients[inv, , drop = FALSE]
  tibble(inv_var = inv, outcome = outcome, n = nrow(d),
         beta = cr[1, 1], se = cr[1, 2], t_val = cr[1, 3],
         p_val = cr[1, 4], r2_adj = s$adj.r.squared)
}

reg_grid <- imap_dfr(anal_list, function(df, xv) {
  map_dfr(groups_ols, function(grp) {
    df_g <- if (grp == "ALL") df else filter(df, koppen_group == grp)
    res  <- map_dfr(inv_vars, ~run_ols(df_g, .x, "TRRI"))
    if (nrow(res) == 0) return(NULL)
    mutate(res, group = grp, x_var = xv)
  })
}) %>%
  mutate(sig = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
                         p_val < 0.05  ~ "*",   p_val < 0.1  ~ ".", TRUE ~ "ns"),
         inv_label = inv_labels[inv_var])

write_csv(reg_grid, file.path(OUT, "reg_grid_vpd.csv"))

cat("\n=== 回归结果（|t|前15）===\n")
reg_grid %>% arrange(desc(abs(t_val))) %>% head(15) %>%
  mutate(across(c(beta, r2_adj), ~round(.x, 4)), p_val = round(p_val, 4)) %>%
  dplyr::select(x_var, group, inv_var, n, beta, p_val, sig, r2_adj) %>%
  print()

# ---- 方差分解 ----------------------------------------------------------------
run_vp <- function(df, inv) {
  d <- df %>%
    filter(!is.na(.data[[inv]]), !is.na(longitude),
           if_all(all_of(ctrl_vars), ~!is.na(.x)),
           .data[[inv]] > 0,
           .data[[inv]] < quantile(.data[[inv]], 0.99, na.rm = TRUE))
  if (nrow(d) < 20) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI,
                   dplyr::select(d, all_of(inv)),
                   dplyr::select(d, vpd_bg, heat_freq_bg),
                   dplyr::select(d, latitude, longitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(inv_var    = inv,
         n          = nrow(d),
         invest_R2  = fr$Adj.R.square[1],
         climate_R2 = fr$Adj.R.square[2],
         geo_R2     = fr$Adj.R.square[3])
}

vp_results <- imap_dfr(anal_list, function(df, xv) {
  map_dfr(groups_ols, function(grp) {
    df_g <- if (grp == "ALL") df else filter(df, koppen_group == grp)
    res  <- map_dfr(inv_vars, ~run_vp(df_g, .x))
    if (nrow(res) == 0) return(NULL)
    mutate(res, group = grp, x_var = xv)
  })
})

write_csv(vp_results, file.path(OUT, "varpart_vpd.csv"))

cat("\n=== 方差分解（全部站点，按投资贡献排序）===\n")
vp_results %>% filter(group == "ALL") %>%
  arrange(x_var, desc(invest_R2)) %>%
  mutate(across(c(invest_R2, climate_R2, geo_R2), ~sprintf("%.2f%%", .x * 100))) %>%
  dplyr::select(x_var, inv_var, n, invest_R2, climate_R2, geo_R2) %>%
  print(n = 20)

# =============================================================================
# D.  分组分析 D1–D4
# =============================================================================
# D1: always_inhibit vs 其他（逻辑回归）
# D2: always_promote  vs 其他（逻辑回归）
# D3: inhibit_promote → ITW（OLS；越小=恢复越快）
# D4: promote_inhibit → HTW（OLS；越大=促进期越长）

run_logistic <- function(df, inv, outcome_col) {
  d <- df %>%
    filter(!is.na(.data[[inv]]),
           if_all(all_of(ctrl_vars), ~!is.na(.x)),
           !is.na(.data[[outcome_col]]),
           .data[[inv]] > 0,
           .data[[inv]] < quantile(.data[[inv]], 0.99, na.rm = TRUE))
  if (nrow(d) < 15 || length(unique(d[[outcome_col]])) < 2) return(NULL)
  fml <- as.formula(paste(outcome_col, "~", inv, "+",
                          paste(ctrl_vars, collapse = "+")))
  m   <- tryCatch(glm(fml, data = d, family = binomial()),
                  error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s  <- summary(m)
  cr <- s$coefficients[inv, , drop = FALSE]
  tibble(inv_var  = inv, outcome = outcome_col,
         n        = nrow(d),
         n_pos    = sum(d[[outcome_col]] == 1),
         beta     = cr[1, 1], se = cr[1, 2],
         z_val    = cr[1, 3], p_val = cr[1, 4])
}

grouped_res <- imap_dfr(anal_list, function(df, xv) {
  # 准备数据：添加二元标记列
  df_g <- df %>%
    mutate(
      is_always_inhibit = as.integer(stype == "always_inhibit"),
      is_always_promote = as.integer(stype == "always_promote")
    )

  results <- map_dfr(inv_vars, function(inv) {
    # D1: always_inhibit
    r1 <- run_logistic(df_g, inv, "is_always_inhibit") %>%
      { if (!is.null(.)) mutate(., analysis = "D1_always_inhibit") else NULL }

    # D2: always_promote
    r2 <- run_logistic(df_g, inv, "is_always_promote") %>%
      { if (!is.null(.)) mutate(., analysis = "D2_always_promote") else NULL }

    # D3: inhibit_promote → ITW
    df_itw <- df_g %>% filter(stype == "inhibit_promote", !is.na(ITW))
    r3 <- run_ols(df_itw, inv, "ITW") %>%
      { if (!is.null(.)) mutate(., analysis = "D3_inhibit_promote_ITW",
                                 n_pos = nrow(df_itw)) else NULL }

    # D4: promote_inhibit → HTW
    df_htw <- df_g %>% filter(stype == "promote_inhibit", !is.na(HTW))
    r4 <- run_ols(df_htw, inv, "HTW") %>%
      { if (!is.null(.)) mutate(., analysis = "D4_promote_inhibit_HTW",
                                 n_pos = nrow(df_htw)) else NULL }

    bind_rows(r1, r2, r3, r4)
  })

  if (nrow(results) == 0) return(NULL)
  mutate(results, x_var = xv)
}) %>%
  mutate(sig = case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
                         p_val < 0.05  ~ "*",   p_val < 0.1  ~ ".", TRUE ~ "ns"),
         inv_label = inv_labels[inv_var])

write_csv(grouped_res, file.path(OUT, "grouped_analysis.csv"))

cat("\n=== 分组分析结果（|t|或|z|前20）===\n")
grouped_res %>%
  mutate(stat = coalesce(t_val, z_val)) %>%
  arrange(analysis, desc(abs(stat))) %>%
  mutate(across(c(beta), ~round(.x, 4)),
         p_val = round(p_val, 4)) %>%
  dplyr::select(x_var, analysis, inv_var, n, beta, p_val, sig) %>%
  print(n = 40)

# =============================================================================
# E.  图表输出
# =============================================================================
xvar_labels <- c(vpd_mean_dt = "X=VPD均值", heat_over_dt = "X=VPD超量")
grp_labels  <- c(ALL = "全部", B = "B(干旱)", C = "C(温带)", D = "D(大陆)")

# 图1：OLS 森林图（TRRI）
p_forest <- reg_grid %>%
  filter(!is.na(beta)) %>%
  mutate(
    ci_lo    = beta - 1.96 * se,
    ci_hi    = beta + 1.96 * se,
    sig_flag = p_val < 0.05,
    group    = factor(group, levels = c("ALL", "B", "C", "D")),
    x_label  = xvar_labels[x_var]
  ) %>%
  ggplot(aes(x = beta, y = inv_label, color = sig_flag)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D73027"),
                     labels = c("ns", "p<0.05"), name = "") +
  facet_grid(x_label ~ group) +
  labs(title    = "投资强度对 TRRI 的回归系数（VPD变量版）",
       subtitle = "控制VPD背景、热频率、纬度；红色=显著(p<0.05)",
       x = "回归系数 β", y = NULL) +
  theme_cn() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

ggsave(file.path(OUT, "forest_trri.png"), p_forest, width = 16, height = 8, dpi = 300)
cat("-> forest_trri.png\n")

# 图2：分组分析森林图（D1–D4）
analysis_labels <- c(
  D1_always_inhibit      = "D1: 全程抑制（逻辑回归）",
  D2_always_promote      = "D2: 全程促进（逻辑回归）",
  D3_inhibit_promote_ITW = "D3: 抑制转促进→ITW（OLS）",
  D4_promote_inhibit_HTW = "D4: 促进转抑制→HTW（OLS）"
)

p_grouped <- grouped_res %>%
  filter(!is.na(beta)) %>%
  mutate(
    ci_lo     = beta - 1.96 * se,
    ci_hi     = beta + 1.96 * se,
    sig_flag  = p_val < 0.05,
    x_label   = xvar_labels[x_var],
    anal_label = factor(analysis_labels[analysis], levels = analysis_labels)
  ) %>%
  ggplot(aes(x = beta, y = inv_label, color = sig_flag)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D73027"),
                     labels = c("ns", "p<0.05"), name = "") +
  facet_grid(x_label ~ anal_label) +
  labs(title    = "投资强度分组分析：热响应恢复力各指标",
       subtitle = "D1/D2=发生概率(log-odds)；D3=ITW越小恢复越快；D4=HTW越大促进越长",
       x = "回归系数 β", y = NULL) +
  theme_cn() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        strip.text.x = element_text(size = BS - 4))

ggsave(file.path(OUT, "forest_grouped.png"), p_grouped, width = 20, height = 8, dpi = 300)
cat("-> forest_grouped.png\n")

# 图3：方差分解气泡图
p_vp <- vp_results %>%
  pivot_longer(c(invest_R2, climate_R2, geo_R2),
               names_to = "comp", values_to = "r2") %>%
  mutate(
    r2         = pmax(r2 * 100, 0),
    comp_label = case_when(comp == "invest_R2"  ~ "投资",
                           comp == "climate_R2" ~ "气候",
                           comp == "geo_R2"     ~ "地理"),
    inv_label  = inv_labels[inv_var],
    grp_label  = factor(grp_labels[group], levels = grp_labels),
    x_label    = xvar_labels[x_var]
  ) %>%
  filter(!is.na(grp_label)) %>%
  ggplot(aes(x = comp_label, y = inv_label, size = r2, color = comp)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = ifelse(r2 > 0.1, sprintf("%.1f%%", r2), "")),
            size = 3.2, color = "black", vjust = -1.4, family = "heiti") +
  scale_size_continuous(range = c(1, 12), name = "独立贡献(%)") +
  scale_color_manual(values = c(invest_R2  = "#D73027",
                                 climate_R2 = "#4575B4",
                                 geo_R2     = "#1A9850"),
                     guide = "none") +
  facet_grid(x_label ~ grp_label) +
  labs(title = "方差分解：投资/气候/地理对TRRI独立贡献（VPD变量版）",
       x = "方差来源", y = "投资变量") +
  theme_cn() +
  theme(panel.grid.major = element_line(color = "grey92"))

ggsave(file.path(OUT, "varpart_vpd.png"), p_vp, width = 16, height = 8, dpi = 300)
cat("-> varpart_vpd.png\n")

# 图4：站点类型分布
p_stype <- imap_dfr(trri_list, ~mutate(.x, x_var = .y)) %>%
  mutate(x_label = xvar_labels[x_var]) %>%
  count(x_label, stype) %>%
  group_by(x_label) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(stype = factor(stype,
                        levels = c("always_inhibit", "inhibit_promote",
                                   "promote_inhibit", "always_promote"),
                        labels = c("全程抑制", "抑制→促进",
                                   "促进→抑制", "全程促进"))) %>%
  ggplot(aes(x = stype, y = pct, fill = stype)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
            vjust = -0.3, size = 4, family = "heiti") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"),
                    guide = "none") +
  facet_wrap(~x_label) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title    = "热响应类型分布（VPD变量版）",
       subtitle = sprintf("共%d个站点", n_distinct(ccm$meteo_stat_id)),
       x = NULL, y = "比例 (%)") +
  theme_cn()

ggsave(file.path(OUT, "stype_dist.png"), p_stype, width = 12, height = 6, dpi = 300)
cat("-> stype_dist.png\n")

# =============================================================================
# 汇总表
# =============================================================================
summary_vpd <- reg_grid %>%
  dplyr::select(x_var, group, inv_var, inv_label, n, beta, se,
                t_val, p_val, sig, r2_adj) %>%
  left_join(vp_results %>%
              dplyr::select(x_var, group, inv_var,
                            invest_R2_pct = invest_R2,
                            climate_R2_pct = climate_R2,
                            geo_R2_pct = geo_R2),
            by = c("x_var", "group", "inv_var")) %>%
  mutate(across(c(invest_R2_pct, climate_R2_pct, geo_R2_pct),
                ~round(.x * 100, 2))) %>%
  arrange(x_var, group, desc(abs(t_val)))

write_csv(summary_vpd, file.path(OUT, "summary_vpd.csv"))
cat("-> summary_vpd.csv\n")

# =============================================================================
# 最终摘要
# =============================================================================
cat("\n", strrep("=", 65), "\n")
cat(sprintf("TRRI量表：1..%d（N_TP=%d, tp=0..%d）\n", TRRI_MAX, N_TP, MAX_TP))
cat(strrep("=", 65), "\n\n")

best_r <- reg_grid %>% filter(!is.na(p_val)) %>% arrange(desc(abs(t_val))) %>% head(1)
cat(sprintf("【TRRI回归最强】\n  X=%s | 投资=%s | 组=%s\n  β=%.4f (p=%.4f%s) R²=%.3f\n\n",
            best_r$x_var, best_r$inv_var, best_r$group,
            best_r$beta, best_r$p_val, best_r$sig, best_r$r2_adj))

best_vp <- vp_results %>% filter(group == "ALL") %>% arrange(desc(invest_R2)) %>% head(1)
cat(sprintf("【方差分解最强（全部）】\n  X=%s | 投资=%s\n  投资=%.2f%% 气候=%.2f%% 地理=%.2f%%\n\n",
            best_vp$x_var, best_vp$inv_var,
            best_vp$invest_R2 * 100, best_vp$climate_R2 * 100, best_vp$geo_R2 * 100))

best_g <- grouped_res %>%
  filter(!is.na(p_val)) %>%
  mutate(stat = abs(coalesce(t_val, z_val))) %>%
  arrange(desc(stat)) %>% head(1)
cat(sprintf("【分组分析最强】\n  X=%s | 分析=%s | 投资=%s\n  β=%.4f (p=%.4f%s)\n",
            best_g$x_var, best_g$analysis, best_g$inv_var,
            best_g$beta, best_g$p_val, best_g$sig))

cat("\n输出目录:", OUT, "\n")
cat("文件列表:\n")
cat(paste(" -", list.files(OUT)), sep = "\n")
