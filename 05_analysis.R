# =============================================================================
# 05_analysis.R  ——  完整独立分析脚本
#
# 与 04_analysis.R 的唯一区别：
#   CCM 自变量改为 vpd_mean_dt（去趋势后的周均VPD）
#   而非 heat_over_dt（VPD超量累积）
#
# X变量   ：vpd_mean_dt（周均VPD，去趋势）
# 投资变量：近10年均值 / 建成区绿地（pa_built_10y）
# 方差分解：投资（1组） vs 地理（经纬度 + Köppen气候区，1组）
#
# 04_analysis.R CCM结果备份路径：
#   data_proc/output_10y_built_up_04_20260601/ccm_04_results.rds
#
# 运行顺序：
#   A. 加载数据 & 去趋势
#   B. SIF结构突变检测（CCM前过滤）& 可视化
#   C. CCM 分析（vpd_mean_dt × tp=0..8，全量站点）
#   D. TRRI 分类
#   E. 投资变量构建（pa_built_10y）
#   F. 合并分析数据框
#   G. OLS 回归
#   H. 方差分解（投资 vs 地理）
#   I. 图表 & 汇总输出
#   J. 分层方差分解 1：按响应类型分层
#   K. 分层方差分解 2：按气候区分层
#   L. 城市级聚合 + 方差分解
#   M. 城市内部站点异质性可视化
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr,
               vegan, tibble, scales, showtext, sysfonts,
               targets, rEDM, strucchange, MASS)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT     <- "data_proc/output_10y_built_up_05"
CCM_RDS <- "data_proc/output_10y_built_up_05/ccm_05_results.rds"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

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
          axis.title    = element_text(size = bs))
}
theme_set(theme_cn())


# A. 加载数据 & 去趋势
# -----------------------------------------------------------------------------

tar_load(data_heat_sif_weekly)

safe_detrend <- function(x, t) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch(
    residuals(lm(x ~ t, na.action = na.exclude)),
    error = function(e) rep(NA_real_, length(x))
  )
}

cat("去趋势 vpd_mean 和 SIF...\n")
df_weekly <- data_heat_sif_weekly %>%
  group_by(meteo_stat_id) %>%
  arrange(year, week) %>%
  mutate(
    time_idx      = row_number(),
    vpd_mean_dt   = safe_detrend(vpd_mean,    time_idx),   # ← 05版主X变量
    heat_over_dt  = safe_detrend(heat_over_sum, time_idx), # 保留备用
    sif_detrended = safe_detrend(sif_interp,  time_idx)
  ) %>%
  ungroup()

# 站点元数据（Köppen等）
meta <- readRDS("data_proc/results_weekly_0_5.rds") %>%
  distinct(meteo_stat_id, longitude, latitude, koppen_class, koppen_group)


# B. SIF 结构突变检测 & 可视化
# -----------------------------------------------------------------------------

cat("\nSIF突变检测（p<0.01 且幅度>1SD）...\n")

detect_sif_break <- function(sid, df) {
  d <- df %>%
    filter(meteo_stat_id == sid, week %in% 20:39) %>%
    arrange(year, week) %>%
    filter(!is.na(sif_interp))
  if (nrow(d) < 30)
    return(tibble(meteo_stat_id = sid, has_break = NA_integer_,
                  break_pval = NA_real_, break_magnitude = NA_real_,
                  break_year = NA_real_))
  tryCatch({
    fs  <- strucchange::Fstats(sif_interp ~ 1, data = d,
                               from = 0.15, to = 0.85)
    pv  <- strucchange::sctest(fs)$p.value
    bp  <- strucchange::breakpoints(sif_interp ~ 1, data = d, breaks = 1)
    idx <- bp$breakpoints
    if (!is.na(idx) && idx > 1 && idx < nrow(d)) {
      mu1 <- mean(d$sif_interp[1:idx],           na.rm = TRUE)
      mu2 <- mean(d$sif_interp[(idx+1):nrow(d)], na.rm = TRUE)
      mag <- abs(mu2 - mu1) / sd(d$sif_interp, na.rm = TRUE)
      yr  <- d$year[idx]
    } else { mag <- 0; yr <- NA_real_ }
    tibble(meteo_stat_id   = sid,
           has_break       = as.integer(pv < 0.01 & mag > 1.0),
           break_pval      = pv,
           break_magnitude = mag,
           break_year      = yr)
  }, error = function(e)
    tibble(meteo_stat_id = sid, has_break = 0L,
           break_pval = NA_real_, break_magnitude = NA_real_,
           break_year = NA_real_))
}

all_sids <- df_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(sif_interp)) >= 30) %>%
  summarise(.groups = "drop") %>%
  pull(meteo_stat_id)

# 若04分析中已有突变结果则复用（SIF突变检测与X变量无关）
break_rds_04 <- "data_proc/output_10y_built_up_04_20260601/sif_break_results.rds"
break_rds_05 <- file.path(OUT, "sif_break_results.rds")

if (file.exists(break_rds_04)) {
  cat("复用04分析的SIF突变结果:", break_rds_04, "\n")
  break_res <- readRDS(break_rds_04)
} else if (file.exists(break_rds_05)) {
  cat("读取已有突变结果:", break_rds_05, "\n")
  break_res <- readRDS(break_rds_05)
} else {
  break_res <- map_dfr(all_sids, ~detect_sif_break(.x, df_weekly))
  saveRDS(break_res, break_rds_05)
}

n_break <- sum(break_res$has_break == 1L, na.rm = TRUE)
n_clean <- sum(break_res$has_break == 0L, na.rm = TRUE)
cat(sprintf("  突变站点: %d (%.1f%%)  |  正常站点: %d (%.1f%%)\n",
            n_break, n_break / length(all_sids) * 100,
            n_clean, n_clean / length(all_sids) * 100))

# 可视化：各抽4个站点对比
set.seed(7)
n_break_avail <- sum(break_res$has_break == 1L, na.rm = TRUE)
n_clean_avail <- sum(break_res$has_break == 0L, na.rm = TRUE)
sids_break <- break_res %>% filter(has_break == 1L) %>%
  slice_sample(n = min(4L, n_break_avail)) %>% pull(meteo_stat_id)
sids_clean <- break_res %>% filter(has_break == 0L) %>%
  slice_sample(n = min(4L, n_clean_avail)) %>% pull(meteo_stat_id)

plot_sids  <- c(sids_break, sids_clean)
plot_label <- tibble(
  meteo_stat_id = plot_sids,
  group_label   = c(rep("SIF突变站点（已排除）", length(sids_break)),
                    rep("正常站点（纳入CCM）",   length(sids_clean)))
)

sif_ts <- df_weekly %>%
  filter(meteo_stat_id %in% plot_sids, week %in% 20:39) %>%
  arrange(meteo_stat_id, year, week) %>%
  left_join(plot_label, by = "meteo_stat_id") %>%
  left_join(dplyr::select(break_res, meteo_stat_id, break_year),
            by = "meteo_stat_id")

p_sif <- ggplot(sif_ts,
                aes(x = year + (week - 20) / 20, y = sif_interp)) +
  geom_line(color = "steelblue", linewidth = 0.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  geom_vline(aes(xintercept = break_year), color = "#D73027",
             linetype = "solid", linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ meteo_stat_id + group_label, scales = "free_y", ncol = 4) +
  labs(title    = "SIF 时间序列：突变站点 vs 正常站点",
       subtitle = "红色竖线 = 检测到的突变年份（p<0.01 且幅度>1SD）",
       x = "年份", y = "SIF（插值后）") +
  theme_cn() +
  theme(strip.text = element_text(size = BS - 4))

ggsave(file.path(OUT, "sif_break_vs_clean.png"),
       p_sif, width = 16, height = 8, dpi = 300)
cat("-> sif_break_vs_clean.png\n")


# C. CCM 分析（vpd_mean_dt × tp=0..8，全量站点）
# -----------------------------------------------------------------------------
# !! 与04_analysis.R的核心区别：X变量为 vpd_mean_dt（周均VPD，去趋势）

stations_clean <- break_res %>%
  filter(has_break == 0L) %>%
  pull(meteo_stat_id)

stations_ok <- df_weekly %>%
  filter(week %in% 20:39, meteo_stat_id %in% stations_clean) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(vpd_mean_dt))  >= 40,    # ← 检查 vpd_mean_dt 可用性
         sum(!is.na(sif_detrended)) >= 40) %>%
  summarise(.groups = "drop") %>%
  pull(meteo_stat_id)

cat(sprintf("\n突变过滤后有效站点: %d\n", length(stations_ok)))

if (file.exists(CCM_RDS)) {
  cat("发现已有CCM结果，直接读取:", CCM_RDS, "\n")
  ccm_results <- readRDS(CCM_RDS)
} else {

  perform_ccm <- function(station_id, df, tp_x = 0, min_pts = 40) {
    d <- df %>%
      filter(meteo_stat_id == station_id, week %in% 20:39) %>%
      arrange(year, week) %>%
      mutate(x_var = dplyr::lag(vpd_mean_dt, n = tp_x)) %>%   # ← vpd_mean_dt
      dplyr::select(sif = sif_detrended, heat_index = x_var) %>%
      filter(!is.na(sif), !is.na(heat_index)) %>%
      mutate(time = row_number(), .before = 1) %>%
      as.data.frame()
    if (nrow(d) < min_pts) return(NULL)
    n <- nrow(d)
    tryCatch({
      best_E <- rEDM::EmbedDimension(
        dataFrame = d, columns = "sif", target = "sif",
        lib = paste("1", n), pred = paste("1", n),
        maxE = min(8, floor(n / 10)), showPlot = FALSE
      ) %>% { .$E[which.max(.$rho)] }

      ccm_res <- rEDM::CCM(
        dataFrame = d, E = best_E, Tp = 0,
        columns = "heat_index", target = "sif",
        libSizes = paste(best_E + 2, n - best_E,
                         max(2, floor((n - 2 * best_E - 2) / 10))),
        sample = 50, random = TRUE, showPlot = FALSE
      )
      ccm_sum <- ccm_res %>%
        group_by(LibSize) %>%
        summarise(rho_m = mean(`heat_index:sif`, na.rm = TRUE),
                  .groups = "drop")

      sm <- rEDM::SMap(
        dataFrame = d, E = best_E, theta = 2,
        lib = paste("1", n), pred = paste("1", n),
        columns = "heat_index", target = "sif",
        embedded = FALSE
      )
      coef_col  <- which(grepl("heat", colnames(sm$coefficients),
                               ignore.case = TRUE))[1]
      if (is.na(coef_col)) coef_col <- 2
      mean_coef <- mean(sm$coefficients[, coef_col], na.rm = TRUE)

      tibble(meteo_stat_id = station_id,
             tp            = tp_x,
             rho           = ccm_sum$rho_m[nrow(ccm_sum)],
             trend         = cor(ccm_sum$LibSize, ccm_sum$rho_m),
             mean_coef     = mean_coef,
             n_obs         = n)
    }, error = function(e) NULL)
  }

  tp_seq      <- 0:8
  n_stations  <- length(stations_ok)
  total_calls <- length(tp_seq) * n_stations
  cat(sprintf("开始CCM（X=vpd_mean_dt）：%d站 × %d tp = %d次\n",
              n_stations, length(tp_seq), total_calls))

  ccm_results <- map_dfr(tp_seq, function(tp_now) {
    cat(sprintf("  tp=%d ...\n", tp_now))
    map_dfr(stations_ok, ~perform_ccm(.x, df_weekly, tp_now))
  })

  ccm_results <- ccm_results %>%
    left_join(meta, by = "meteo_stat_id")

  saveRDS(ccm_results, CCM_RDS)
  cat(sprintf("CCM结果已保存 (%d行): %s\n", nrow(ccm_results), CCM_RDS))
}

cat(sprintf("站点数: %d | tp: %d..%d\n",
            length(unique(ccm_results$meteo_stat_id)),
            min(ccm_results$tp), max(ccm_results$tp)))


# D. TRRI 分类
# -----------------------------------------------------------------------------

N_TP     <- length(unique(ccm_results$tp))   # 应为 9
MAX_TP   <- max(ccm_results$tp)
TRRI_MAX <- 2L * N_TP
cat(sprintf("\nN_TP=%d, TRRI范围: 1..%d\n", N_TP, TRRI_MAX))

find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm = TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

trri_df <- ccm_results %>%
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
  filter(map_lgl(coef_seq, ~length(.x) == N_TP)) %>%
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
      stype == "inhibit_promote" ~ as.integer(1L + (N_TP - ITW)),
      stype == "promote_inhibit" ~ as.integer(N_TP + HTW),
      stype == "always_promote"  ~ as.integer(2L * N_TP)
    )
  )

cat("\n站点类型分布:\n")
print(count(trri_df, stype))


# E. 投资与绿地结构变量
# -----------------------------------------------------------------------------
# 投资变量（3个，分别单独进入模型）：
#   pa_built_10y   : 2011-2020年均值 / 建成区绿地面积
#   pa_built_5y    : 近5年均值（可用最新5年）/ 建成区绿地面积
#   pa_built_latest: 最新一年 / 建成区绿地面积

tar_load(green_invest_2020)
tar_load(green_area_2020)
tar_load(station_city_map)

# --- E1: 读取投资原始数据 ---
invest_raw0 <- read.csv("data_raw/green_invest/city_invest_data.csv",
                        check.names = FALSE)
nc_inv <- ncol(invest_raw0)
names(invest_raw0) <- c("city_raw", paste0("inv_", 2002:(2001 + nc_inv - 1)))

invest_raw <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  dplyr::select(invest_raw0, -city_raw)
)

year_cols <- paste0("inv_", 2002:2030)
year_cols <- year_cols[year_cols %in% names(invest_raw)]
years_num <- as.integer(str_extract(year_cols, "\\d{4}$"))

# 10年均值：2011-2020
cols_10y <- year_cols[years_num >= 2011 & years_num <= 2020]
inv_10y  <- rowMeans(
  dplyr::select(invest_raw, all_of(cols_10y)) %>%
    mutate(across(everything(), as.numeric)),
  na.rm = TRUE
)

# 近5年均值：可用年份中最新5年
cols_5y <- tail(year_cols, 5)
inv_5y  <- rowMeans(
  dplyr::select(invest_raw, all_of(cols_5y)) %>%
    mutate(across(everything(), as.numeric)),
  na.rm = TRUE
)

# 最新一年
col_latest <- tail(year_cols, 1)
inv_latest <- as.numeric(invest_raw[[col_latest]])

cat(sprintf("\n投资数据年份范围: %d - %d\n", min(years_num), max(years_num)))
cat(sprintf("10年均值年份: %s\n", paste(years_num[years_num >= 2011 & years_num <= 2020], collapse = " ")))
cat(sprintf("近5年均值年份: %s\n", paste(str_extract(cols_5y, "\\d{4}$"), collapse = " ")))
cat(sprintf("最新年份: %s\n", str_extract(col_latest, "\\d{4}$")))

# --- E2: 建成区绿地面积（用于计算单位面积投资）---
green_2020 <- green_area_2020 %>%
  dplyr::select(city_name, area_built = area_green_built) %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))

invest_tbl <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  tibble(inv_10y = inv_10y, inv_5y = inv_5y, inv_latest = inv_latest)
) %>%
  left_join(green_2020, by = "city_name") %>%
  mutate(
    pa_built_10y    = inv_10y    / area_built,
    pa_built_5y     = inv_5y     / area_built,
    pa_built_latest = inv_latest / area_built
  )

city_vars <- invest_tbl

invest_station <- station_city_map %>%
  left_join(city_vars, by = "city_name")

cat(sprintf("\npa_built_10y    非NA站点: %d\n", sum(!is.na(invest_station$pa_built_10y))))
cat(sprintf("pa_built_5y     非NA站点: %d\n", sum(!is.na(invest_station$pa_built_5y))))
cat(sprintf("pa_built_latest 非NA站点: %d\n", sum(!is.na(invest_station$pa_built_latest))))


# F. 合并分析数据框
# -----------------------------------------------------------------------------

winsorize <- function(x, p = 0.99) {
  q <- quantile(x, p, na.rm = TRUE)
  ifelse(x > q, q, x)
}

anal_df <- trri_df %>%
  left_join(invest_station, by = "meteo_stat_id") %>%
  filter(!is.na(TRRI), !is.na(longitude), !is.na(latitude),
         !is.na(koppen_group)) %>%
  mutate(
    koppen_B           = as.integer(koppen_group == "B"),
    koppen_C           = as.integer(koppen_group == "C"),
    koppen_D           = as.integer(koppen_group == "D"),
    pa_built_10y_w     = winsorize(pa_built_10y),
    pa_built_5y_w      = winsorize(pa_built_5y),
    pa_built_latest_w  = winsorize(pa_built_latest)
  )

cat(sprintf("\n站点总数（TRRI有效）: %d\n", nrow(anal_df)))
cat(sprintf("pa_built_10y    非NA: %d\n", sum(!is.na(anal_df$pa_built_10y))))
cat(sprintf("pa_built_5y     非NA: %d\n", sum(!is.na(anal_df$pa_built_5y))))
cat(sprintf("pa_built_latest 非NA: %d\n", sum(!is.na(anal_df$pa_built_latest))))
cat("气候组分布:\n"); print(count(anal_df, koppen_group))


# G. OLS 回归（3个投资变量分别单独入模）
# -----------------------------------------------------------------------------
geo_vars <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

# 3个投资变量定义（原始变量名 + winsorize后变量名 + 标签）
invest_vars <- list(
  list(raw = "pa_built_10y",    w = "pa_built_10y_w",    label = "近10年均值（2011-2020）"),
  list(raw = "pa_built_5y",     w = "pa_built_5y_w",     label = "近5年均值"),
  list(raw = "pa_built_latest", w = "pa_built_latest_w", label = "最新一年")
)

run_ols_var <- function(df, var, label) {
  d <- df %>% filter(!is.na(.data[[var]]), .data[[var]] > 0)
  if (nrow(d) < 20) return(NULL)
  fml <- as.formula(paste("TRRI ~", var, "+",
                          paste(geo_vars, collapse = "+")))
  m   <- lm(fml, data = d)
  s   <- summary(m)
  cr  <- s$coefficients[var, ]
  sig <- case_when(cr[4]<0.001~"***", cr[4]<0.01~"**",
                   cr[4]<0.05~"*",   cr[4]<0.1~".", TRUE~"ns")
  cat(sprintf("\n[%s] n=%d  β=%.4f  SE=%.4f  t=%.3f  p=%.4f%s  R²=%.3f\n",
              label, nrow(d), cr[1], cr[2], cr[3], cr[4], sig, s$adj.r.squared))
  list(model=m, summary=s, coef=cr, sig=sig, n=nrow(d), var=var, label=label)
}

cat("\n=== G. OLS 回归（因变量：TRRI，3个投资变量分别入模）===\n")
ols_list <- map(invest_vars, function(v)
  run_ols_var(anal_df, v$w, v$label))
names(ols_list) <- sapply(invest_vars, `[[`, "raw")

# 取pa_built_10y作为后续代表性汇总变量
ols_main <- ols_list[["pa_built_10y"]]
s_ols    <- ols_main$summary
cr_inv   <- ols_main$coef
sig_inv  <- ols_main$sig
main_var <- "pa_built_10y_w"

# 回归系数汇总表
coef_tbl <- map_dfr(invest_vars, function(v) {
  obj <- ols_list[[v$raw]]
  if (is.null(obj)) return(NULL)
  cr <- obj$coef
  tibble(
    variable    = v$raw,
    label       = v$label,
    n           = obj$n,
    beta        = round(cr[1], 6),
    se          = round(cr[2], 6),
    t_val       = round(cr[3], 3),
    p_val       = round(cr[4], 4),
    sig         = obj$sig,
    r2_adj_full = round(obj$summary$adj.r.squared, 4),
    direction   = ifelse(cr[1] > 0,
                         "↑ 投资增加→TRRI上升（韧性增强）",
                         "↓ 投资增加→TRRI下降（韧性减弱）")
  )
})

cat("\n=== G. 回归系数汇总（站点级，3个投资变量）===\n")
print(coef_tbl %>% dplyr::select(label, n, beta, se, p_val, sig, direction, r2_adj_full))
write_csv(coef_tbl, file.path(OUT, "ols_coef_invest.csv"))
cat("-> ols_coef_invest.csv\n")


# G2. 有序 Logit 回归（MASS::polr，TRRI视为有序因子）
# -----------------------------------------------------------------------------
# TRRI取值1..18，天然是有序离散变量，有序Logit比OLS更符合数据结构。
# polr() 不直接提供p值，用 z = coef/SE（大样本正态近似）计算。
# 系数解释：正值表示投资增加 → 倾向于更高TRRI等级（韧性增强）。
# 伪R²：McFadden's R² = 1 - logLik(full)/logLik(null)
# -----------------------------------------------------------------------------

run_olr_var <- function(df, var, label, outcome = "TRRI") {
  d <- df %>% filter(!is.na(.data[[var]]), .data[[var]] > 0,
                     !is.na(.data[[outcome]]))
  if (nrow(d) < 20) return(NULL)
  # 因变量转有序因子（若为连续则四舍五入后转）
  ord_vals <- sort(unique(round(d[[outcome]])))
  d$Y_ord  <- factor(round(d[[outcome]]), levels = ord_vals, ordered = TRUE)
  fml <- as.formula(paste("Y_ord ~", var, "+",
                           paste(geo_vars, collapse = "+")))
  m <- tryCatch(
    MASS::polr(fml, data = d, Hess = TRUE, method = "logistic"),
    error = function(e) { cat(sprintf("  [polr error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s  <- summary(m)
  cr <- s$coefficients[var, ]           # Estimate / Std. Error / t value
  z  <- cr["t value"]
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  sig  <- case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**",
                    pval < 0.05  ~ "*",   pval < 0.1  ~ ".",  TRUE ~ "ns")
  # McFadden 伪 R²
  m0 <- tryCatch(
    MASS::polr(Y_ord ~ 1, data = d, Hess = FALSE, method = "logistic"),
    error = function(e) NULL)
  mcf <- if (!is.null(m0))
    round(1 - as.numeric(logLik(m)) / as.numeric(logLik(m0)), 4)
  else NA_real_
  cat(sprintf("\n[%s] n=%d  coef=%.4f  SE=%.4f  z=%.3f  p=%.4f%s  McF_R²=%.4f\n",
              label, nrow(d), cr[1], cr[2], z, pval, sig, mcf))
  list(model = m, coef = cr, pval = pval, sig = sig,
       n = nrow(d), var = var, label = label, mcfadden = mcf)
}

cat("\n=== G2. 有序Logit（站点级，因变量=TRRI有序因子）===\n")
olr_list <- map(invest_vars, function(v)
  run_olr_var(anal_df, v$w, v$label, outcome = "TRRI"))
names(olr_list) <- sapply(invest_vars, `[[`, "raw")

olr_coef_tbl <- map_dfr(invest_vars, function(v) {
  obj <- olr_list[[v$raw]]
  if (is.null(obj)) return(NULL)
  cr <- obj$coef
  tibble(
    model       = "有序Logit",
    variable    = v$raw,
    label       = v$label,
    n           = obj$n,
    coef        = round(cr[1], 6),   # log-odds scale
    se          = round(cr[2], 6),
    z_val       = round(cr[3], 3),
    p_val       = round(obj$pval, 4),
    sig         = obj$sig,
    mcfadden_r2 = obj$mcfadden,
    direction   = ifelse(cr[1] > 0,
                         "↑ 投资增加→倾向更高TRRI（韧性增强）",
                         "↓ 投资增加→倾向更低TRRI（韧性减弱）")
  )
})

# OLS与有序Logit并排对比
olr_compare <- bind_rows(
  coef_tbl %>%
    mutate(model = "OLS") %>%
    dplyr::rename(coef = beta, z_val = t_val, mcfadden_r2 = r2_adj_full) %>%
    dplyr::select(model, variable, label, n, coef, se, z_val, p_val, sig,
                  mcfadden_r2, direction),
  olr_coef_tbl
) %>% arrange(variable, model)

cat("\n=== G2. OLS vs 有序Logit 系数对比（站点级）===\n")
print(olr_compare %>% dplyr::select(model, label, n, coef, se, p_val, sig, direction))
write_csv(olr_coef_tbl,  file.path(OUT, "olr_coef_invest.csv"))
write_csv(olr_compare,   file.path(OUT, "olr_vs_ols_compare.csv"))
cat("-> olr_coef_invest.csv\n-> olr_vs_ols_compare.csv\n")


# H. 方差分解（3个投资变量分别 vs 地理）
# -----------------------------------------------------------------------------
run_varpart <- function(df, inv_var, label) {
  d <- df %>% filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0)
  if (nrow(d) < 20) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI,
                   dplyr::select(d, all_of(inv_var)),
                   dplyr::select(d, all_of(geo_vars))),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  cat(sprintf("[%s] n=%d  投资=%.2f%%  地理=%.2f%%  共享=%.2f%%\n",
              label, nrow(d),
              fr$Adj.R.square[1]*100, fr$Adj.R.square[2]*100,
              fr$Adj.R.square[3]*100))
  tibble(variable=inv_var, label=label, n=nrow(d),
         invest_R2 = round(fr$Adj.R.square[1]*100, 2),
         geo_R2    = round(fr$Adj.R.square[2]*100, 2),
         shared_R2 = round(fr$Adj.R.square[3]*100, 2))
}

cat("\n=== H. 方差分解（TRRI ~ 投资变量 | 地理，3个变量分别）===\n")
vp_compare <- map_dfr(invest_vars, function(v)
  run_varpart(anal_df, v$w, v$label))
print(vp_compare)
write_csv(vp_compare, file.path(OUT, "varpart_invest_compare.csv"))
cat("-> varpart_invest_compare.csv\n")

# 主方差分解对象取pa_built_10y（供 L 段站点级 vs 城市级对比）
d_main <- anal_df %>% filter(!is.na(pa_built_10y_w), pa_built_10y_w > 0)
vp_main <- tryCatch(
  vegan::varpart(d_main$TRRI,
                 dplyr::select(d_main, pa_built_10y_w),
                 dplyr::select(d_main, all_of(geo_vars))),
  error = function(e) NULL)
fr <- if (!is.null(vp_main)) vp_main$part$indfract else NULL

vp_tbl <- if (!is.null(fr)) tibble(
  component  = c("近10年均值（独立）", "地理（独立）", "共享部分", "未解释"),
  adj_R2_pct = round(fr$Adj.R.square * 100, 2)
)
cat("\n"); print(vp_tbl)
write_csv(vp_tbl, file.path(OUT, "varpart_05.csv"))


# I. 图表 & 汇总输出
# -----------------------------------------------------------------------------

# 图1：散点图（3个投资变量 vs TRRI，三面板）
scatter_data <- map_dfr(invest_vars, function(v) {
  anal_df %>%
    filter(!is.na(.data[[v$raw]]), .data[[v$raw]] > 0) %>%
    mutate(x_val = .data[[v$w]], panel = v$label)
})

p_scatter <- ggplot(scatter_data,
                    aes(x = x_val, y = TRRI, color = koppen_group)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linetype = "dashed", linewidth = 0.9) +
  scale_color_brewer(palette = "Set1", name = "Köppen气候区") +
  scale_y_continuous(breaks = seq(1, TRRI_MAX, by = 2)) +
  facet_wrap(~panel, scales = "free_x", nrow = 1) +
  labs(
    title    = "单位面积投资（不同时间窗口）与热响应韧性指数（TRRI）",
    subtitle = "X=VPD均值（去趋势）；三个时间窗口分别展示（99%分位截尾）",
    x = NULL, y = sprintf("TRRI（1–%d）", TRRI_MAX)
  ) +
  theme_cn()

ggsave(file.path(OUT, "scatter_trri_invest.png"),
       p_scatter, width = 18, height = 6, dpi = 300)
cat("-> scatter_trri_invest.png\n")

# 图2：方差分解对比条形图（3个投资变量并排）
p_vp <- vp_compare %>%
  pivot_longer(c(invest_R2, geo_R2, shared_R2),
               names_to = "comp", values_to = "r2") %>%
  mutate(
    r2_show    = pmax(r2, 0),
    comp_label = factor(comp,
                        levels = c("invest_R2", "shared_R2", "geo_R2"),
                        labels = c("投资变量", "共享", "地理")),
    label      = factor(label, levels = sapply(invest_vars, `[[`, "label"))
  ) %>%
  ggplot(aes(x = comp_label, y = r2_show, fill = comp_label)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f%%", r2_show)),
            vjust = -0.4, size = 4.5, family = "heiti") +
  scale_fill_manual(
    values = c("投资变量" = "#D73027", "共享" = "#FDAE61", "地理" = "#4575B4"),
    guide  = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~label, nrow = 1) +
  labs(
    title    = "方差分解：3个投资变量 vs 地理的独立贡献对比",
    subtitle = "X=VPD均值（去趋势）；地理 = 经纬度 + Köppen气候区哑变量",
    x = NULL, y = "调整R²（%）"
  ) +
  theme_cn()

ggsave(file.path(OUT, "varpart_bar_05.png"),
       p_vp, width = 14, height = 6, dpi = 300)
cat("-> varpart_bar_05.png\n")

# 图3：S-map系数轨迹（4类站点均值±SE）
coef_ts <- ccm_results %>%
  left_join(dplyr::select(trri_df, meteo_stat_id, stype), by = "meteo_stat_id") %>%
  filter(!is.na(stype)) %>%
  group_by(stype, tp) %>%
  summarise(
    mean_c = mean(mean_coef, na.rm = TRUE),
    se_c   = sd(mean_coef, na.rm = TRUE) / sqrt(sum(!is.na(mean_coef))),
    .groups = "drop"
  ) %>%
  mutate(stype_label = factor(stype,
    levels = c("always_inhibit", "inhibit_promote",
               "promote_inhibit", "always_promote"),
    labels = c("全程抑制", "抑制→促进", "促进→抑制", "全程促进")))

p_coef <- ggplot(coef_ts,
                 aes(x = tp, y = mean_c,
                     color = stype_label, fill = stype_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = mean_c - se_c, ymax = mean_c + se_c),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("全程抑制"  = "#4575B4", "抑制→促进" = "#74ADD1",
               "促进→抑制" = "#FDAE61", "全程促进"  = "#D73027"),
    name = "响应类型") +
  scale_fill_manual(
    values = c("全程抑制"  = "#4575B4", "抑制→促进" = "#74ADD1",
               "促进→抑制" = "#FDAE61", "全程促进"  = "#D73027"),
    name = "响应类型") +
  scale_x_continuous(breaks = 0:MAX_TP) +
  labs(
    title    = "各类型站点的热响应轨迹（S-map系数均值）",
    subtitle = "X=VPD均值（去趋势）；色带=±1SE；tp=滞后周数",
    x        = "时间滞后 tp（周）",
    y        = "S-map系数（VPD→SIF）"
  ) +
  theme_cn() +
  theme(legend.position = "right")

ggsave(file.path(OUT, "coef_trajectory.png"),
       p_coef, width = 11, height = 7, dpi = 300)
cat("-> coef_trajectory.png\n")

# 图4：站点类型分布条形
p_stype <- trri_df %>%
  count(stype) %>%
  mutate(
    pct        = n / sum(n) * 100,
    stype_label = factor(stype,
      levels = c("always_inhibit", "inhibit_promote",
                 "promote_inhibit", "always_promote"),
      labels = c("全程抑制", "抑制→促进", "促进→抑制", "全程促进"))
  ) %>%
  ggplot(aes(x = stype_label, y = pct, fill = stype_label)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
            vjust = -0.3, size = 4.5, family = "heiti") +
  scale_fill_manual(
    values = c("全程抑制"  = "#4575B4", "抑制→促进" = "#74ADD1",
               "促进→抑制" = "#FDAE61", "全程促进"  = "#D73027"),
    guide = "none") +
  scale_y_continuous(limits = c(0, 65),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title    = "热响应类型分布",
       subtitle = sprintf("共%d站（X=VPD均值去趋势，tp=0..%d）",
                          nrow(trri_df), MAX_TP),
       x = NULL, y = "比例 (%)") +
  theme_cn()

ggsave(file.path(OUT, "stype_dist.png"),
       p_stype, width = 9, height = 6, dpi = 300)
cat("-> stype_dist.png\n")

# 汇总CSV（以pa_built_10y为代表）
result_summary <- tibble(
  x_var         = "vpd_mean_dt",
  inv_var       = "pa_built_10y（近10年均值，代表）",
  n_ccm         = length(unique(ccm_results$meteo_stat_id)),
  n_regression  = nrow(anal_df),
  beta          = cr_inv[1],
  se            = cr_inv[2],
  t_val         = cr_inv[3],
  p_val         = cr_inv[4],
  sig           = sig_inv,
  r2_adj_full   = s_ols$adj.r.squared,
  invest_R2_pct = if (!is.null(fr)) round(fr$Adj.R.square[1] * 100, 2) else NA,
  geo_R2_pct    = if (!is.null(fr)) round(fr$Adj.R.square[2] * 100, 2) else NA,
  shared_R2_pct = if (!is.null(fr)) round(fr$Adj.R.square[3] * 100, 2) else NA
)
write_csv(result_summary, file.path(OUT, "summary_05.csv"))
cat("-> summary_05.csv\n")

cat("\n", strrep("=", 60), "\n")
cat("分析摘要（05版：X=vpd_mean_dt）\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("X变量      : vpd_mean_dt（周均VPD，去趋势）\n"))
cat(sprintf("投资变量   : pa_built_10y / pa_built_5y / pa_built_latest（分别入模）\n"))
cat(sprintf("TRRI量表   : 1..%d（N_TP=%d, tp=0..%d）\n",
            TRRI_MAX, N_TP, MAX_TP))
cat(sprintf("CCM站点数  : %d\n", length(unique(ccm_results$meteo_stat_id))))
cat(sprintf("回归站点数 : %d\n", nrow(anal_df)))
if (!is.null(fr)) {
  cat(sprintf("\n【方差分解（pa_built_10y）】\n  投资独立贡献 = %.2f%%\n  地理独立贡献 = %.2f%%\n  共享部分     = %.2f%%\n",
              fr$Adj.R.square[1] * 100,
              fr$Adj.R.square[2] * 100,
              fr$Adj.R.square[3] * 100))
}
cat("\n输出目录:", OUT, "\n")


# J. 分层方差分解 1：按响应类型分层
# -----------------------------------------------------------------------------

stype_levels <- c("always_promote", "promote_inhibit",
                  "inhibit_promote", "always_inhibit")
stype_labels <- c("全程促进", "促进→抑制", "抑制→促进", "全程抑制")

stype_outcome <- list(
  inhibit_promote  = list(var = "ITW",              label = "ITW（转折周数，小=恢复快）"),
  promote_inhibit  = list(var = "HTW",              label = "HTW（促进周数，大=持续长）"),
  always_inhibit   = list(var = "inhibit_strength", label = "抑制强度（|coef@tp0|）"),
  always_promote   = list(var = "promote_strength", label = "促进强度（coef@tp0）")
)

tp0_coef <- ccm_results %>%
  filter(tp == 0) %>%
  dplyr::select(meteo_stat_id, coef_tp0 = mean_coef)

anal_df2 <- anal_df %>%
  left_join(tp0_coef, by = "meteo_stat_id") %>%
  mutate(
    inhibit_strength = -coef_tp0,
    promote_strength =  coef_tp0
  )

run_vp2 <- function(df, outcome_var, label) {
  d <- df %>%
    filter(!is.na(.data[[outcome_var]]),
           !is.na(pa_built_10y), !is.na(longitude), !is.na(latitude),
           pa_built_10y > 0)
  if (nrow(d) < 15) {
    cat(sprintf("    [跳过] n=%d < 15\n", nrow(d)))
    return(NULL)
  }
  kopp_cols <- c("koppen_B","koppen_C","koppen_D")
  kopp_use  <- kopp_cols[sapply(kopp_cols, function(v) var(d[[v]], na.rm=TRUE) > 0)]
  geo_df    <- dplyr::select(d, longitude, latitude, all_of(kopp_use))

  vp <- tryCatch(
    vegan::varpart(d[[outcome_var]],
                   dplyr::select(d, pa_built_10y),
                   geo_df),
    error = function(e) { cat("    varpart error:", e$message, "\n"); NULL }
  )
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(
    outcome     = outcome_var,
    label       = label,
    n           = nrow(d),
    invest_R2   = round(fr$Adj.R.square[1] * 100, 2),
    geo_R2      = round(fr$Adj.R.square[2] * 100, 2),
    shared_R2   = round(fr$Adj.R.square[3] * 100, 2),
    unexplained = round(fr$Adj.R.square[4] * 100, 2)
  )
}

cat("\n=== J. 按响应类型分层的方差分解 ===\n")
vp_stype <- map_dfr(names(stype_outcome), function(st) {
  info <- stype_outcome[[st]]
  df_s <- anal_df2 %>% filter(stype == st)
  cat(sprintf("\n[%s] n=%d, 因变量=%s\n", st, nrow(df_s), info$label))
  res <- run_vp2(df_s, info$var, info$label)
  if (!is.null(res)) mutate(res, stype = st)
})

if (nrow(vp_stype) > 0) {
  vp_stype <- vp_stype %>%
    mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))
  write_csv(vp_stype, file.path(OUT, "varpart_by_stype.csv"))
  cat("\n"); print(vp_stype %>% dplyr::select(stype, label, n, invest_R2, geo_R2))

  p_vp_stype <- vp_stype %>%
    pivot_longer(c(invest_R2, geo_R2, shared_R2),
                 names_to = "comp", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(comp,
                          levels = c("invest_R2", "shared_R2", "geo_R2"),
                          labels = c("投资", "共享", "地理"))
    ) %>%
    ggplot(aes(x = comp_label, y = stype_label)) +
    geom_point(aes(size = r2_show, color = comp), alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", r2_show)),
              size = 4, vjust = -1.5, family = "heiti") +
    scale_size_continuous(range = c(2, 14), name = "独立解释力(%)") +
    scale_color_manual(
      values = c(invest_R2 = "#D73027", shared_R2 = "#FDAE61", geo_R2 = "#4575B4"),
      guide  = "none") +
    labs(
      title    = "分层方差分解：响应类型分层",
      subtitle = "投资=pa_built_10y；地理=经纬度+Köppen哑变量\nX=VPD均值（去趋势）",
      x = "方差来源", y = "响应类型"
    ) +
    theme_cn() +
    theme(panel.grid.major = element_line(color = "grey92"))

  ggsave(file.path(OUT, "varpart_by_stype.png"),
         p_vp_stype, width = 9, height = 7, dpi = 300)
  cat("-> varpart_by_stype.png\n")
}


# K. 分层方差分解 2：按气候区分层（3个投资变量分别）
# -----------------------------------------------------------------------------

cat("\n=== K. 按气候区分层的方差分解 ===\n")

run_vp_geo <- function(df, inv_var, grp_label) {
  d <- df %>%
    filter(!is.na(TRRI), !is.na(.data[[inv_var]]),
           !is.na(longitude), !is.na(latitude),
           .data[[inv_var]] > 0)
  if (nrow(d) < 15) {
    cat(sprintf("    [%s/%s 跳过] n=%d < 15\n", inv_var, grp_label, nrow(d)))
    return(NULL)
  }
  vp <- tryCatch(
    vegan::varpart(d$TRRI,
                   dplyr::select(d, all_of(inv_var)),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) { cat("    varpart error:", e$message, "\n"); NULL }
  )
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(
    inv_var      = inv_var,
    koppen_group = grp_label,
    n            = nrow(d),
    invest_R2    = round(fr$Adj.R.square[1] * 100, 2),
    geo_R2       = round(fr$Adj.R.square[2] * 100, 2),
    shared_R2    = round(fr$Adj.R.square[3] * 100, 2),
    unexplained  = round(fr$Adj.R.square[4] * 100, 2)
  )
}

koppen_groups <- sort(unique(anal_df$koppen_group))

# 3个变量分别跑
vp_koppen_all3 <- map_dfr(invest_vars, function(v) {
  cat(sprintf("\n--- %s ---\n", v$label))
  bind_rows(
    run_vp_geo(anal_df, v$raw, "ALL"),
    map_dfr(koppen_groups, function(g) {
      cat(sprintf("  [%s] n=%d\n", g, sum(anal_df$koppen_group == g)))
      run_vp_geo(filter(anal_df, koppen_group == g), v$raw, g)
    })
  ) %>% mutate(inv_label = v$label)
})

write_csv(vp_koppen_all3, file.path(OUT, "varpart_by_koppen_all3.csv"))
cat("\n"); print(vp_koppen_all3 %>% dplyr::select(inv_label, koppen_group, n, invest_R2, geo_R2))
cat("-> varpart_by_koppen_all3.csv\n")

# 保留pa_built_10y单独的对象供K2/L使用
vp_koppen_full <- vp_koppen_all3 %>% filter(inv_var == "pa_built_10y")

write_csv(vp_koppen_full, file.path(OUT, "varpart_by_koppen.csv"))
cat("-> varpart_by_koppen.csv（pa_built_10y）\n")

if (nrow(vp_koppen_all3) > 0) {
  koppen_order <- intersect(c("ALL","A","B","C","D"), vp_koppen_all3$koppen_group)
  koppen_name  <- c(ALL="全部", A="A(热带)", B="B(干旱)",
                    C="C(温带)", D="D(大陆)")
  inv_label_levels <- sapply(invest_vars, `[[`, "label")

  # 图：各气候区×各投资变量的投资独立R²对比（面板=气候区，x=投资变量）
  p_vp_koppen <- vp_koppen_all3 %>%
    mutate(
      r2_show   = pmax(invest_R2, 0),
      grp_label = factor(
        dplyr::recode(koppen_group, !!!koppen_name),
        levels = koppen_name[koppen_order]
      ),
      inv_label = factor(inv_label, levels = inv_label_levels)
    ) %>%
    filter(!is.na(grp_label)) %>%
    ggplot(aes(x = inv_label, y = r2_show, fill = inv_label)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.2f%%", invest_R2)),
              vjust = -0.4, size = 3.8, family = "heiti") +
    scale_fill_manual(
      values = c("#D73027", "#FC8D59", "#FEE090")[seq_along(inv_label_levels)],
      name = "投资变量") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    facet_wrap(~grp_label, nrow = 1) +
    labs(
      title    = "气候区分层方差分解：3个投资变量投资独立解释力对比",
      subtitle = "X=VPD均值（去趋势）；因变量=TRRI；地理=经纬度（不含气候区哑变量）",
      x = NULL, y = "投资独立调整R²（%）"
    ) +
    theme_cn() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(OUT, "varpart_by_koppen.png"),
         p_vp_koppen, width = 14, height = 6, dpi = 300)
  cat("-> varpart_by_koppen.png\n")
}


# K2. 气候区分层：站点级 vs 城市级对比（3个投资变量）
# -----------------------------------------------------------------------------

cat("\n=== K2. 气候区分层：站点级 vs 城市级对比 ===\n")

# 城市级数据（L段会再次使用）
city_agg_k <- anal_df %>%
  group_by(city_name) %>%
  summarise(
    TRRI_mean       = mean(TRRI,      na.rm = TRUE),
    pa_built_10y    = first(pa_built_10y),
    pa_built_5y     = first(pa_built_5y),
    pa_built_latest = first(pa_built_latest),
    longitude       = mean(longitude, na.rm = TRUE),
    latitude        = mean(latitude,  na.rm = TRUE),
    koppen_group    = first(koppen_group),
    .groups         = "drop"
  ) %>%
  filter(!is.na(TRRI_mean), !is.na(longitude), !is.na(latitude))

run_vp_city_var <- function(df, inv_var, grp_label) {
  d <- df %>% filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0)
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, all_of(inv_var)),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(inv_var = inv_var, koppen_group = grp_label, n = nrow(d),
         invest_R2 = round(fr2$Adj.R.square[1]*100, 2),
         geo_R2    = round(fr2$Adj.R.square[2]*100, 2),
         shared_R2 = round(fr2$Adj.R.square[3]*100, 2),
         unexplained = round(fr2$Adj.R.square[4]*100, 2))
}

# 3个变量分别跑城市级气候区分层
vp_city_koppen_all3 <- map_dfr(invest_vars, function(v) {
  bind_rows(
    run_vp_city_var(city_agg_k, v$raw, "ALL"),
    map_dfr(sort(unique(city_agg_k$koppen_group)), ~run_vp_city_var(
      filter(city_agg_k, koppen_group == .x), v$raw, .x))
  ) %>% mutate(inv_label = v$label)
})

cat("\n城市级气候区方差分解（3变量）：\n")
print(vp_city_koppen_all3 %>% dplyr::select(inv_label, koppen_group, n, invest_R2, geo_R2))

# pa_built_10y城市级单独保留供L段
vp_city_koppen_full <- vp_city_koppen_all3 %>% filter(inv_var == "pa_built_10y")

# 合并站点级与城市级（以pa_built_10y为代表）
koppen_compare <- bind_rows(
  vp_koppen_full %>%
    mutate(level = "站点级", data_n_label = paste0("n=", n, "站")),
  vp_city_koppen_full %>%
    mutate(level = "城市级（聚合）", data_n_label = paste0("n=", n, "城市"))
) %>%
  mutate(level = factor(level, levels = c("站点级", "城市级（聚合）")))

write_csv(koppen_compare, file.path(OUT, "varpart_koppen_compare.csv"))
cat("-> varpart_koppen_compare.csv（pa_built_10y代表）\n")

# 全3变量的站点级+城市级对比表
koppen_compare_all3 <- bind_rows(
  vp_koppen_all3 %>% mutate(level = "站点级"),
  vp_city_koppen_all3 %>% mutate(level = "城市级（聚合）")
) %>% mutate(level = factor(level, levels = c("站点级", "城市级（聚合）")))
write_csv(koppen_compare_all3, file.path(OUT, "varpart_koppen_compare_all3.csv"))
cat("-> varpart_koppen_compare_all3.csv（3变量完整对比）\n")

cat("\n站点级 vs 城市级气候区对比（pa_built_10y）：\n")
print(koppen_compare %>% dplyr::select(koppen_group, level, n, invest_R2, geo_R2, unexplained))

# 对比图：面板=气候区，x=投资变量，颜色=站点/城市级
if (nrow(koppen_compare_all3) > 0) {
  koppen_order2  <- intersect(c("ALL","A","B","C","D"), koppen_compare_all3$koppen_group)
  koppen_name2   <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  inv_label_levels <- sapply(invest_vars, `[[`, "label")

  p_koppen_cmp <- koppen_compare_all3 %>%
    mutate(
      r2_show   = pmax(invest_R2, 0),
      r2_label  = sprintf("%.2f%%", invest_R2),
      grp_label = factor(dplyr::recode(koppen_group, !!!koppen_name2),
                         levels = koppen_name2[koppen_order2]),
      inv_label = factor(inv_label, levels = inv_label_levels)
    ) %>%
    filter(!is.na(grp_label)) %>%
    ggplot(aes(x = inv_label, y = r2_show, fill = level)) +
    geom_col(width = 0.65, alpha = 0.85, position = position_dodge(width = 0.75)) +
    geom_text(aes(label = r2_label),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 3.2, family = "heiti") +
    scale_fill_manual(
      values = c("站点级" = "#D73027", "城市级（聚合）" = "#4575B4"),
      name = "分析层级") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, NA)) +
    facet_wrap(~grp_label, nrow = 1) +
    labs(
      title    = "气候区分层方差分解：3个投资变量 × 站点级/城市级对比",
      subtitle = "X=VPD均值（去趋势）；地理=经纬度；标注值含负数（显示为0）",
      x = NULL, y = "投资独立调整R²（%）"
    ) +
    theme_cn() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(OUT, "varpart_koppen_compare.png"),
         p_koppen_cmp, width = 16, height = 6, dpi = 300)
  cat("-> varpart_koppen_compare.png\n")
}


# L. 城市级聚合 + 方差分解（消除空间伪重复）
# -----------------------------------------------------------------------------

cat("\n=== L. 城市级聚合方差分解 ===\n")

city_agg <- anal_df %>%
  group_by(city_name) %>%
  summarise(
    n_stations   = n(),
    TRRI_mean    = mean(TRRI,        na.rm = TRUE),
    TRRI_sd      = sd(TRRI,          na.rm = TRUE),
    TRRI_range   = max(TRRI) - min(TRRI),
    pa_built_10y    = first(pa_built_10y),
    pa_built_5y     = first(pa_built_5y),
    pa_built_latest = first(pa_built_latest),
    longitude    = mean(longitude,   na.rm = TRUE),
    latitude     = mean(latitude,    na.rm = TRUE),
    koppen_group = first(koppen_group),
    koppen_B     = first(koppen_B),
    koppen_C     = first(koppen_C),
    koppen_D     = first(koppen_D),
    .groups      = "drop"
  ) %>%
  filter(!is.na(TRRI_mean), !is.na(pa_built_10y),
         !is.na(longitude), !is.na(latitude))

cat(sprintf("城市数: %d（站点数中位数: %.0f）\n",
            nrow(city_agg), median(city_agg$n_stations)))

# --- 城市级 OLS 回归系数（投资对TRRI均值的效应）---
run_ols_city <- function(df, inv_var, var_label) {
  d <- df %>% filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0,
                     !is.na(TRRI_mean))
  if (nrow(d) < 10) return(NULL)
  d$inv_w <- as.numeric(scale(d[[inv_var]]))
  m <- lm(TRRI_mean ~ inv_w + longitude + latitude + koppen_B + koppen_C + koppen_D, data = d)
  s <- summary(m)
  cr <- s$coefficients["inv_w", ]
  sig <- case_when(cr[4]<0.001~"***", cr[4]<0.01~"**", cr[4]<0.05~"*",
                   cr[4]<0.1~".", TRUE~"ns")
  tibble(
    level       = "城市级（聚合）",
    variable    = inv_var,
    label       = var_label,
    n           = nrow(d),
    beta        = round(cr[1], 6),
    se          = round(cr[2], 6),
    t_val       = round(cr[3], 3),
    p_val       = round(cr[4], 4),
    sig         = sig,
    r2_adj_full = round(s$adj.r.squared, 4),
    direction   = ifelse(cr[1] > 0,
                         "↑ 投资增加→TRRI上升（韧性增强）",
                         "↓ 投资增加→TRRI下降（韧性减弱）")
  )
}

ols_city_tbl <- map_dfr(invest_vars, function(v)
  run_ols_city(city_agg, v$raw, paste0(v$label, "（城市级）")))

cat("\n=== L. 城市级OLS回归系数 ===\n")
print(ols_city_tbl %>% dplyr::select(label, n, beta, se, p_val, sig, direction, r2_adj_full))
write_csv(ols_city_tbl, file.path(OUT, "ols_coef_invest_city.csv"))
cat("-> ols_coef_invest_city.csv\n")

# 与站点级系数合并输出一张对比表
if (file.exists(file.path(OUT, "ols_coef_invest.csv"))) {
  station_coef <- read_csv(file.path(OUT, "ols_coef_invest.csv"), show_col_types = FALSE) %>%
    mutate(level = "站点级")
  ols_compare <- bind_rows(station_coef, ols_city_tbl) %>%
    dplyr::select(level, label, n, beta, se, p_val, sig, direction, r2_adj_full)
  write_csv(ols_compare, file.path(OUT, "ols_coef_invest_compare.csv"))
  cat("-> ols_coef_invest_compare.csv（站点级+城市级对比）\n")
  cat("\n=== 站点级 vs 城市级 OLS系数对比 ===\n")
  print(ols_compare)
}

# G2b. 有序Logit：城市级（TRRI_mean四舍五入后视为有序因子）
# 注：城市级因变量为站点TRRI均值（连续），四舍五入后作为近似有序处理，
#     结果仅供参考，OLS仍为城市级主要回归方法。
cat("\n=== G2b. 有序Logit（城市级，TRRI_mean四舍五入）===\n")

run_olr_city <- function(df, inv_var, var_label) {
  d <- df %>% filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0,
                     !is.na(TRRI_mean))
  if (nrow(d) < 15) return(NULL)
  d$inv_w <- as.numeric(scale(d[[inv_var]]))
  d$Y_ord <- factor(round(d$TRRI_mean),
                    levels = sort(unique(round(d$TRRI_mean))),
                    ordered = TRUE)
  m <- tryCatch(
    MASS::polr(Y_ord ~ inv_w + longitude + latitude + koppen_B + koppen_C + koppen_D,
               data = d, Hess = TRUE, method = "logistic"),
    error = function(e) { cat(sprintf("  [polr error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s    <- summary(m)
  cr   <- s$coefficients["inv_w", ]
  z    <- cr["t value"]
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  sig  <- case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**",
                    pval < 0.05  ~ "*",   pval < 0.1  ~ ".", TRUE ~ "ns")
  m0 <- tryCatch(
    MASS::polr(Y_ord ~ 1, data = d, Hess = FALSE, method = "logistic"),
    error = function(e) NULL)
  mcf <- if (!is.null(m0))
    round(1 - as.numeric(logLik(m)) / as.numeric(logLik(m0)), 4)
  else NA_real_
  cat(sprintf("\n[%s] n=%d  coef=%.4f  SE=%.4f  z=%.3f  p=%.4f%s  McF_R²=%.4f\n",
              var_label, nrow(d), cr[1], cr[2], z, pval, sig, mcf))
  tibble(
    model       = "有序Logit",
    level       = "城市级（聚合，TRRI_mean四舍五入）",
    variable    = inv_var,
    label       = var_label,
    n           = nrow(d),
    coef        = round(cr[1], 6),
    se          = round(cr[2], 6),
    z_val       = round(z, 3),
    p_val       = round(pval, 4),
    sig         = sig,
    mcfadden_r2 = mcf,
    direction   = ifelse(cr[1] > 0,
                         "↑ 投资增加→倾向更高TRRI（韧性增强）",
                         "↓ 投资增加→倾向更低TRRI（韧性减弱）")
  )
}

olr_city_tbl <- map_dfr(invest_vars, function(v)
  run_olr_city(city_agg, v$raw, v$label))

print(olr_city_tbl %>% dplyr::select(label, n, coef, se, p_val, sig, direction))
write_csv(olr_city_tbl, file.path(OUT, "olr_coef_invest_city.csv"))
cat("-> olr_coef_invest_city.csv\n")


geo_vars_city <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

# 城市级varpart：3个变量分别跑
run_vp_city_full <- function(df, inv_var, label) {
  d <- df %>% filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0, !is.na(TRRI_mean))
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, all_of(inv_var)),
                   dplyr::select(d, all_of(geo_vars_city))),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(level="城市级（聚合）", variable=inv_var, label=label, n=nrow(d),
         invest_R2 = round(fr2$Adj.R.square[1]*100, 2),
         geo_R2    = round(fr2$Adj.R.square[2]*100, 2),
         shared_R2 = round(fr2$Adj.R.square[3]*100, 2),
         unexplained = round(fr2$Adj.R.square[4]*100, 2))
}

vp_city_all3 <- map_dfr(invest_vars, function(v)
  run_vp_city_full(city_agg, v$raw, v$label))
cat("\n城市级方差分解（3变量）：\n")
print(vp_city_all3 %>% dplyr::select(label, n, invest_R2, geo_R2, shared_R2))

# 站点级+城市级汇总对比（pa_built_10y为代表）
vp_city_10y <- vp_city_all3 %>% filter(variable == "pa_built_10y")
fr_city_10y <- vp_city_10y$invest_R2  # 仅作参考

compare_tbl <- bind_rows(
  tibble(level="站点级（原始）", variable="pa_built_10y", label="近10年均值（2011-2020）",
         invest_R2 = if (!is.null(fr)) fr$Adj.R.square[1]*100 else NA,
         geo_R2    = if (!is.null(fr)) fr$Adj.R.square[2]*100 else NA,
         shared_R2 = if (!is.null(fr)) fr$Adj.R.square[3]*100 else NA,
         n         = nrow(anal_df)),
  vp_city_all3 %>% filter(variable == "pa_built_10y")
) %>% mutate(across(c(invest_R2, geo_R2, shared_R2), ~round(.x, 2)))

# 全3变量站点+城市合并
compare_tbl_all3 <- bind_rows(
  vp_compare %>% mutate(level = "站点级（原始）"),
  vp_city_all3
) %>% mutate(across(c(invest_R2, geo_R2, shared_R2), ~round(.x, 2)))

cat("\n--- 站点级 vs 城市级对比（3变量）---\n")
print(compare_tbl_all3 %>% dplyr::select(level, label, n, invest_R2, geo_R2))
write_csv(compare_tbl_all3, file.path(OUT, "varpart_city_vs_station.csv"))
cat("-> varpart_city_vs_station.csv\n")

# 图：3变量 × 站点/城市级投资解释力对比
p_compare <- compare_tbl_all3 %>%
  mutate(
    r2_show   = pmax(invest_R2, 0),
    r2_label  = sprintf("%.2f%%", invest_R2),
    label     = factor(label, levels = sapply(invest_vars, `[[`, "label")),
    level     = factor(level, levels = c("站点级（原始）", "城市级（聚合）"))
  ) %>%
  ggplot(aes(x = label, y = r2_show, fill = level)) +
  geom_col(width = 0.6, alpha = 0.85, position = position_dodge(width = 0.7)) +
  geom_text(aes(label = r2_label),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4, family = "heiti") +
  scale_fill_manual(
    values = c("站点级（原始）" = "#D73027", "城市级（聚合）" = "#4575B4"),
    name = "分析层级") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title    = "方差分解：3个投资变量 × 站点级 vs 城市级",
       subtitle = "X=VPD均值（去趋势）；消除城市内伪重复后投资解释力的变化",
       x = NULL, y = "投资独立调整R²（%）") +
  theme_cn() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(OUT, "varpart_city_compare.png"),
       p_compare, width = 10, height = 6, dpi = 300)
cat("-> varpart_city_compare.png\n")


# M. 城市内部站点异质性可视化（分面板输出）
# -----------------------------------------------------------------------------

cat("\n=== M. 城市内部站点异质性可视化 ===\n")

stype_colors <- c(
  "always_inhibit"  = "#4575B4",
  "inhibit_promote" = "#74ADD1",
  "promote_inhibit" = "#FDAE61",
  "always_promote"  = "#D73027"
)
stype_labels_map <- c(
  always_inhibit  = "全程抑制",
  inhibit_promote = "抑制→促进",
  promote_inhibit = "促进→抑制",
  always_promote  = "全程促进"
)

city_stype <- anal_df %>%
  dplyr::select(city_name, stype, koppen_group) %>%
  mutate(stype_label = factor(stype_labels_map[stype],
                              levels = stype_labels_map))

city_order <- city_stype %>%
  group_by(city_name) %>%
  summarise(
    n_total      = n(),
    n_types      = n_distinct(stype),
    dominant     = names(sort(table(stype), decreasing=TRUE))[1],
    koppen_group = first(koppen_group),
    .groups      = "drop"
  ) %>%
  arrange(koppen_group, dominant, desc(n_total)) %>%
  mutate(city_rank = row_number())

N_PER   <- 50
n_city  <- nrow(city_order)
n_panel <- ceiling(n_city / N_PER)

city_order <- city_order %>%
  mutate(panel = ceiling(city_rank / N_PER))

cat(sprintf("城市总数: %d，分 %d 个面板（每面板≤%d城市）\n",
            n_city, n_panel, N_PER))

walk(1:n_panel, function(pid) {
  cities_p   <- city_order %>% filter(panel == pid) %>% arrange(city_rank)
  city_levels <- cities_p$city_name

  df_p <- city_stype %>%
    filter(city_name %in% city_levels) %>%
    group_by(city_name, stype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(city_name) %>%
    mutate(pct = n / sum(n) * 100, n_total = sum(n)) %>%
    ungroup() %>%
    left_join(dplyr::select(city_order, city_name, koppen_group, n_types),
              by = "city_name") %>%
    mutate(
      city_name  = factor(city_name, levels = city_levels),
      city_label = ifelse(n_types > 1,
                          paste0(as.character(city_name), " *"),
                          as.character(city_name))
    )

  label_order <- df_p %>%
    distinct(city_name, city_label) %>%
    arrange(match(city_name, city_levels)) %>%
    pull(city_label)
  df_p <- df_p %>%
    mutate(city_label = factor(city_label, levels = label_order))

  koppen_breaks <- cities_p %>%
    group_by(koppen_group) %>%
    summarise(first_rank = min(city_rank), .groups = "drop") %>%
    rowwise() %>%
    mutate(x_pos = match(
      cities_p$city_name[cities_p$city_rank == first_rank][1],
      city_levels
    ) - 0.5) %>%
    ungroup() %>%
    filter(!is.na(x_pos), x_pos > 0.5)

  p <- ggplot(df_p, aes(x = city_label, y = pct, fill = stype_label)) +
    geom_col(width = 0.8, alpha = 0.9) +
    geom_vline(xintercept = koppen_breaks$x_pos,
               color = "grey30", linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = setNames(stype_colors, stype_labels_map),
                      name = "响应类型") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = sprintf("城市内部站点热响应类型分布（面板 %d/%d）", pid, n_panel),
      subtitle = "X=VPD均值（去趋势）；* = 城市内存在多种类型；虚线 = 气候区边界",
      x = NULL, y = "站点比例"
    ) +
    geom_text(data = df_p %>% distinct(city_label, n_total),
              aes(x = city_label, y = 102, label = n_total, fill = NULL),
              size = 2.8, family = "heiti", color = "grey30") +
    theme_cn(bs = 11) +
    theme(
      axis.text.x     = element_text(angle = 55, hjust = 1, size = 7.5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )

  fname <- sprintf("city_stype_panel%02d.png", pid)
  ggsave(file.path(OUT, fname), p, width = 16, height = 7, dpi = 300)
  cat(sprintf("-> %s\n", fname))
})

# 城市异质性汇总
city_diversity <- city_order %>%
  left_join(
    city_stype %>% group_by(city_name, stype) %>% summarise(n=n(),.groups="drop") %>%
      pivot_wider(names_from=stype, values_from=n, values_fill=0),
    by = "city_name"
  )
write_csv(city_diversity, file.path(OUT, "city_stype_diversity.csv"))
cat("-> city_stype_diversity.csv\n")

cat(sprintf("\n城市内类型数分布:\n"))
print(table(city_order$n_types))
cat(sprintf("异质城市（>1种类型）: %d / %d (%.1f%%)\n",
            sum(city_order$n_types > 1), nrow(city_order),
            mean(city_order$n_types > 1)*100))

cat("\n全部输出文件:\n")
cat(paste(" -", list.files(OUT)), sep = "\n")
