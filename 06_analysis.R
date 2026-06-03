# =============================================================================
# 06_analysis.R  ——  完整独立分析脚本
#
# 与 04/05_analysis.R 的区别：
#   CCM 自变量改为 heat_index_norm_dt（归一化热事件指数，去趋势）
#   该指数在 A 段构建：
#     1. heat_over_sum（VPD超量累积）和 heat_event_freq（热事件频率）
#        各自在全局99%分位截尾，再 min-max 归一化到 [0, 1]
#     2. 等权平均：heat_index_norm = 0.5 × norm(heat_over) + 0.5 × norm(heat_freq)
#     3. 按站点去线性趋势 → heat_index_norm_dt
#   指数范围：[0, 1]（截尾前极端值>99%分位被压缩至1）
#   M段城市异质性图改为单张多行拼图输出（不再分多个文件）。
#
# X变量   ：heat_index_norm_dt（归一化热事件指数，去趋势）
# 投资变量：近10年均值 / 建成区绿地（pa_built_10y）+ 绿地存量（green_built_10y）
# 方差分解：投资（1组） vs 地理（经纬度 + Köppen气候区，1组）
#
# 04_analysis.R CCM结果备份路径：
#   data_proc/output_10y_built_up_04_20260601/ccm_04_results.rds
#
# 运行顺序：
#   A. 加载数据 & 构建归一化热事件指数
#   B. SIF结构突变检测（复用已有结果）
#   C. CCM 分析（heat_index_norm_dt × tp=0..8，全量站点）
#   D. TRRI 分类
#   E. 投资变量构建（pa_built_10y + green_built_10y）
#   F. 合并分析数据框
#   G. OLS 回归 + 系数汇总CSV
#   H. 方差分解（投资 vs 地理）
#   I. 图表 & 汇总输出
#   J. 分层方差分解 1：按响应类型分层
#   K. 分层方差分解 2：按气候区分层（站点级）
#   K2. 气候区分层：站点级 vs 城市级对比
#   L. 城市级聚合 + 方差分解
#   M. 城市内部站点异质性可视化（单张多行拼图）
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr,
               vegan, tibble, scales, showtext, sysfonts,
               targets, rEDM, strucchange, patchwork)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT     <- "data_proc/output_10y_built_up_06"
CCM_RDS <- "data_proc/output_10y_built_up_06/ccm_06_results.rds"
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


# A. 加载数据 & 构建归一化热事件指数
# -----------------------------------------------------------------------------
# 归一化热事件指数（heat_index_norm）构建方法：
#   1. heat_over_sum（VPD超量累积，kPa/周）：全局99%分位截尾 → min-max → [0,1]
#   2. heat_event_freq（热事件天占比，0~1/周）：全局99%分位截尾 → min-max → [0,1]
#   3. heat_index_norm = 0.5 × norm(heat_over) + 0.5 × norm(heat_freq)  → [0,1]
#   4. 按站点去线性趋势 → heat_index_norm_dt（CCM X变量）
# 先做全局截尾再归一化，避免极端值（原始最大值>23σ）压缩正常区间的区分度

tar_load(data_heat_sif_weekly)

safe_detrend <- function(x, t) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  tryCatch(
    residuals(lm(x ~ t, na.action = na.exclude)),
    error = function(e) rep(NA_real_, length(x))
  )
}

# 全局截尾+归一化辅助函数（在函数外计算分位数，保证全局一致）
winsorize_global <- function(x, p = 0.99) {
  q <- quantile(x, p, na.rm = TRUE)
  pmin(x, q)
}
minmax_norm <- function(x) {
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  if (is.na(mn) || is.na(mx) || mx == mn) return(rep(0, length(x)))
  (x - mn) / (mx - mn)
}

cat("构建归一化热事件指数（heat_index_norm）...\n")

# 全局截尾（在 mutate 前统一处理，确保同一阈值）
heat_over_w <- winsorize_global(data_heat_sif_weekly$heat_over_sum)
heat_freq_w <- winsorize_global(data_heat_sif_weekly$heat_event_freq)

cat(sprintf("  heat_over_sum：99%%分位截尾阈值 = %.3f kPa\n",
            quantile(data_heat_sif_weekly$heat_over_sum, 0.99, na.rm=TRUE)))
cat(sprintf("  heat_event_freq：99%%分位截尾阈值 = %.3f\n",
            quantile(data_heat_sif_weekly$heat_event_freq, 0.99, na.rm=TRUE)))

heat_over_norm_global <- minmax_norm(heat_over_w)
heat_freq_norm_global <- minmax_norm(heat_freq_w)

df_weekly <- data_heat_sif_weekly %>%
  mutate(
    heat_over_norm  = heat_over_norm_global,
    heat_freq_norm  = heat_freq_norm_global,
    heat_index_norm = 0.5 * heat_over_norm + 0.5 * heat_freq_norm
  ) %>%
  group_by(meteo_stat_id) %>%
  arrange(year, week) %>%
  mutate(
    time_idx           = row_number(),
    heat_index_norm_dt = safe_detrend(heat_index_norm, time_idx),
    sif_detrended      = if (all(is.na(sif_detrended)))
                           safe_detrend(sif_interp, time_idx)
                         else sif_detrended
  ) %>%
  ungroup()

cat(sprintf("  heat_index_norm 范围：[%.4f, %.4f]，均值=%.4f\n",
            min(df_weekly$heat_index_norm, na.rm=TRUE),
            max(df_weekly$heat_index_norm, na.rm=TRUE),
            mean(df_weekly$heat_index_norm, na.rm=TRUE)))
cat(sprintf("  heat_index_norm_dt 非NA行：%d / %d\n",
            sum(!is.na(df_weekly$heat_index_norm_dt)),
            nrow(df_weekly)))

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

# 优先复用已有突变结果（SIF突变检测与X变量无关）
break_rds_prior <- "data_proc/output_10y_built_up_04_20260601/sif_break_results.rds"
break_rds_06    <- file.path(OUT, "sif_break_results.rds")

if (file.exists(break_rds_prior)) {
  cat("复用已有SIF突变结果:", break_rds_prior, "\n")
  break_res <- readRDS(break_rds_prior)
} else if (file.exists(break_rds_06)) {
  cat("读取本地突变结果:", break_rds_06, "\n")
  break_res <- readRDS(break_rds_06)
} else {
  break_res <- map_dfr(all_sids, ~detect_sif_break(.x, df_weekly))
  saveRDS(break_res, break_rds_06)
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


# C. CCM 分析（heat_index_norm_dt × tp=0..8，全量站点）
# -----------------------------------------------------------------------------
# !! 与04/05的核心区别：X变量为归一化热事件指数（A段构建，去趋势）

stations_clean <- break_res %>%
  filter(has_break == 0L) %>%
  pull(meteo_stat_id)

stations_ok <- df_weekly %>%
  filter(week %in% 20:39, meteo_stat_id %in% stations_clean) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(heat_index_norm_dt)) >= 40,
         sum(!is.na(sif_detrended))      >= 40) %>%
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
      # 使用归一化热事件指数（已去趋势）
      mutate(x_var = dplyr::lag(heat_index_norm_dt, n = tp_x)) %>%
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
  cat(sprintf("开始CCM（X=heat_index_norm_dt）：%d站 × %d tp = %d次\n",
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

tar_load(green_invest_2020)
tar_load(green_area_2020)
tar_load(station_city_map)

# --- E1: 投资变量 ---
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
cols_10y  <- year_cols[years_num >= 2011 & years_num <= 2020]
inv_10y   <- rowMeans(
  dplyr::select(invest_raw, all_of(cols_10y)) %>%
    mutate(across(everything(), as.numeric)),
  na.rm = TRUE
)

green_2020 <- green_area_2020 %>%
  dplyr::select(city_name, area_built = area_green_built) %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))

invest_tbl <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  tibble(inv_10y = inv_10y)
) %>%
  left_join(green_2020, by = "city_name") %>%
  mutate(pa_built_10y = inv_10y / area_built)

# --- E2: 多年建成区绿地面积 ---
ga_raw <- readxl::read_excel("data_raw/green_area_data.xlsx")
ga_data <- ga_raw[4:nrow(ga_raw), ]

built_idx   <- seq(3, ncol(ga_raw), by = 3)
built_yrs   <- as.integer(str_extract(names(ga_raw)[built_idx], "^\\d{4}"))
cols_ga_10y <- names(ga_raw)[built_idx[!is.na(built_yrs) &
                                         built_yrs >= 2011 & built_yrs <= 2020]]
cat(sprintf("\n绿地面积：2011-2020可用年份（%d年）: %s\n",
            length(cols_ga_10y),
            paste(str_extract(cols_ga_10y, "^\\d{4}"), collapse = " ")))

ga_tbl <- ga_data %>%
  rename(city_raw = 1) %>%
  dplyr::select(city_raw, all_of(cols_ga_10y)) %>%
  mutate(
    across(-city_raw, as.numeric),
    green_built_10y = rowMeans(across(-city_raw), na.rm = TRUE),
    city_name = paste0(enc2utf8(as.character(city_raw)), "市")
  ) %>%
  dplyr::select(city_name, green_built_10y)

city_vars <- invest_tbl %>%
  left_join(ga_tbl, by = "city_name")

invest_station <- station_city_map %>%
  left_join(city_vars, by = "city_name")

cat(sprintf("\npa_built_10y    非NA站点: %d\n",
            sum(!is.na(invest_station$pa_built_10y))))
cat(sprintf("green_built_10y 非NA站点: %d\n",
            sum(!is.na(invest_station$green_built_10y))))


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
    koppen_B          = as.integer(koppen_group == "B"),
    koppen_C          = as.integer(koppen_group == "C"),
    koppen_D          = as.integer(koppen_group == "D"),
    pa_built_10y_w    = winsorize(pa_built_10y),
    green_built_10y_w = winsorize(green_built_10y)
  )

cat(sprintf("\n站点总数（TRRI有效）: %d\n", nrow(anal_df)))
cat(sprintf("pa_built_10y    非NA: %d\n", sum(!is.na(anal_df$pa_built_10y))))
cat(sprintf("green_built_10y 非NA: %d\n", sum(!is.na(anal_df$green_built_10y))))
cat("气候组分布:\n"); print(count(anal_df, koppen_group))


# G. OLS 回归 + 系数汇总CSV
# -----------------------------------------------------------------------------
geo_vars <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

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

cat("\n=== G. OLS 回归（因变量：TRRI）===\n")
ols_inv   <- run_ols_var(anal_df, "pa_built_10y_w",    "投资/绿地面积（pa_built_10y）")
ols_green <- run_ols_var(anal_df, "green_built_10y_w", "建成区绿地面积（green_built_10y）")

d_both <- anal_df %>%
  filter(!is.na(pa_built_10y), pa_built_10y > 0,
         !is.na(green_built_10y), green_built_10y > 0)
if (nrow(d_both) >= 20) {
  m_both <- lm(TRRI ~ pa_built_10y_w + green_built_10y_w +
                 longitude + latitude + koppen_B + koppen_C + koppen_D,
               data = d_both)
  cat("\n=== 两变量同时入模（n=", nrow(d_both), "）===\n")
  print(summary(m_both))
}

s_ols   <- if (!is.null(ols_green)) ols_green$summary else ols_inv$summary
cr_inv  <- if (!is.null(ols_green)) ols_green$coef    else ols_inv$coef
sig_inv <- if (!is.null(ols_green)) ols_green$sig     else ols_inv$sig
main_var <- if (!is.null(ols_green)) "green_built_10y_w" else "pa_built_10y_w"

# 回归系数汇总表（含效应方向）
extract_coef_row <- function(ols_obj, var_label) {
  if (is.null(ols_obj)) return(NULL)
  cr <- ols_obj$coef
  tibble(
    variable    = ols_obj$var,
    label       = var_label,
    n           = ols_obj$n,
    beta        = round(cr[1], 6),
    se          = round(cr[2], 6),
    t_val       = round(cr[3], 3),
    p_val       = round(cr[4], 4),
    sig         = ols_obj$sig,
    r2_adj_full = round(ols_obj$summary$adj.r.squared, 4),
    direction   = ifelse(cr[1] > 0,
                         "↑ 投资增加→TRRI上升（韧性增强）",
                         "↓ 投资增加→TRRI下降（韧性减弱）")
  )
}

coef_tbl <- bind_rows(
  extract_coef_row(ols_inv,   "投资强度（pa_built_10y，单变量）"),
  extract_coef_row(ols_green, "绿地存量（green_built_10y，单变量）")
)

if (nrow(d_both) >= 20 && exists("m_both")) {
  s_both <- summary(m_both)
  for (v in c("pa_built_10y_w", "green_built_10y_w")) {
    if (v %in% rownames(s_both$coefficients)) {
      cr2 <- s_both$coefficients[v, ]
      coef_tbl <- bind_rows(coef_tbl, tibble(
        variable    = v,
        label       = paste0(sub("_w$", "", v), "（双变量同时入模）"),
        n           = nrow(d_both),
        beta        = round(cr2[1], 6),
        se          = round(cr2[2], 6),
        t_val       = round(cr2[3], 3),
        p_val       = round(cr2[4], 4),
        sig         = case_when(cr2[4]<0.001~"***", cr2[4]<0.01~"**",
                                cr2[4]<0.05~"*",    cr2[4]<0.1~".", TRUE~"ns"),
        r2_adj_full = round(s_both$adj.r.squared, 4),
        direction   = ifelse(cr2[1] > 0,
                             "↑ 投资增加→TRRI上升（韧性增强）",
                             "↓ 投资增加→TRRI下降（韧性减弱）")
      ))
    }
  }
}

cat("\n=== G. 回归系数汇总（投资对TRRI的效应方向）===\n")
print(coef_tbl %>% dplyr::select(label, n, beta, se, p_val, sig, direction, r2_adj_full))
write_csv(coef_tbl, file.path(OUT, "ols_coef_invest.csv"))
cat("-> ols_coef_invest.csv\n")


# H. 方差分解（两个绿地变量分别 vs 地理）
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

cat("\n=== H. 方差分解（TRRI ~ 绿地变量 | 地理）===\n")
vp_compare <- bind_rows(
  run_varpart(anal_df, "pa_built_10y_w",    "投资/绿地面积（pa_built_10y）"),
  run_varpart(anal_df, "green_built_10y_w", "建成区绿地面积（green_built_10y）")
)
print(vp_compare)
write_csv(vp_compare, file.path(OUT, "varpart_green_compare.csv"))

d_main  <- anal_df %>% filter(!is.na(green_built_10y_w), green_built_10y_w > 0)
vp_main <- tryCatch(
  vegan::varpart(d_main$TRRI,
                 dplyr::select(d_main, green_built_10y_w),
                 dplyr::select(d_main, all_of(geo_vars))),
  error = function(e) NULL)
fr <- if (!is.null(vp_main)) vp_main$part$indfract else NULL

vp_tbl <- if (!is.null(fr)) tibble(
  component  = c("绿地存量（独立）", "地理（独立）", "共享部分", "未解释"),
  adj_R2_pct = round(fr$Adj.R.square * 100, 2)
)
cat("\n"); print(vp_tbl)
write_csv(vp_tbl, file.path(OUT, "varpart_06.csv"))


# I. 图表 & 汇总输出
# -----------------------------------------------------------------------------

# 图1：散点图（两个绿地变量 vs TRRI，双面板）
scatter_data <- bind_rows(
  anal_df %>% filter(!is.na(pa_built_10y), pa_built_10y > 0) %>%
    mutate(x_val = pa_built_10y_w,
           panel = "投资强度（投资额/建成区绿地，pa_built_10y）"),
  anal_df %>% filter(!is.na(green_built_10y), green_built_10y > 0) %>%
    mutate(x_val = green_built_10y_w,
           panel = "绿地存量（建成区绿地面积均值，green_built_10y，公顷）")
)

p_scatter <- ggplot(scatter_data,
                    aes(x = x_val, y = TRRI, color = koppen_group)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linetype = "dashed", linewidth = 0.9) +
  scale_color_brewer(palette = "Set1", name = "Köppen气候区") +
  scale_y_continuous(breaks = seq(1, TRRI_MAX, by = 2)) +
  facet_wrap(~panel, scales = "free_x") +
  labs(
    title    = "绿地变量与热响应韧性指数（TRRI）",
    subtitle = "X=归一化热事件指数[0,1]（去趋势）；左=投资强度，右=绿地存量（99%分位截尾）",
    x = NULL, y = sprintf("TRRI（1–%d）", TRRI_MAX)
  ) +
  theme_cn()

ggsave(file.path(OUT, "scatter_trri_green.png"),
       p_scatter, width = 14, height = 7, dpi = 300)
cat("-> scatter_trri_green.png\n")

# 图2：方差分解对比条形图
p_vp <- vp_compare %>%
  pivot_longer(c(invest_R2, geo_R2, shared_R2),
               names_to = "comp", values_to = "r2") %>%
  mutate(
    r2_show    = pmax(r2, 0),
    comp_label = factor(comp,
                        levels = c("invest_R2", "shared_R2", "geo_R2"),
                        labels = c("绿地变量", "共享", "地理"))
  ) %>%
  ggplot(aes(x = comp_label, y = r2_show, fill = comp_label)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f%%", r2_show)),
            vjust = -0.4, size = 4.5, family = "heiti") +
  scale_fill_manual(
    values = c("绿地变量" = "#D73027", "共享" = "#FDAE61", "地理" = "#4575B4"),
    guide  = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~label) +
  labs(
    title    = "方差分解：两个绿地变量 vs 地理的独立贡献对比",
    subtitle = "X=归一化热事件指数[0,1]（去趋势）；地理 = 经纬度 + Köppen气候区哑变量",
    x = NULL, y = "调整R²（%）"
  ) +
  theme_cn()

ggsave(file.path(OUT, "varpart_bar_06.png"),
       p_vp, width = 8, height = 6, dpi = 300)
cat("-> varpart_bar_06.png\n")

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
    subtitle = "X=归一化热事件指数[0,1]（去趋势）；色带=±1SE；tp=滞后周数",
    x        = "时间滞后 tp（周）",
    y        = "S-map系数（热事件指数→SIF）"
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
       subtitle = sprintf("共%d站（X=归一化热事件指数[0,1]，tp=0..%d）",
                          nrow(trri_df), MAX_TP),
       x = NULL, y = "比例 (%)") +
  theme_cn()

ggsave(file.path(OUT, "stype_dist.png"),
       p_stype, width = 9, height = 6, dpi = 300)
cat("-> stype_dist.png\n")

# 汇总CSV
result_summary <- tibble(
  x_var         = "heat_index_norm_dt",
  inv_var       = "green_built_10y",
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
write_csv(result_summary, file.path(OUT, "summary_06.csv"))
cat("-> summary_06.csv\n")

cat("\n", strrep("=", 60), "\n")
cat("分析摘要（06版：X=heat_index_norm_dt，归一化热事件指数）\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("X变量      : heat_index_norm_dt（归一化热事件指数，VPD超量+频率各50%%，[0,1]截尾归一化后去趋势）\n"))
cat(sprintf("投资变量   : pa_built_10y + green_built_10y\n"))
cat(sprintf("TRRI量表   : 1..%d（N_TP=%d, tp=0..%d）\n", TRRI_MAX, N_TP, MAX_TP))
cat(sprintf("CCM站点数  : %d\n", length(unique(ccm_results$meteo_stat_id))))
cat(sprintf("回归站点数 : %d\n", nrow(anal_df)))
if (!is.null(fr)) {
  cat(sprintf("\n【方差分解】\n  绿地存量独立贡献 = %.2f%%\n  地理独立贡献     = %.2f%%\n  共享部分         = %.2f%%\n",
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
      subtitle = "投资=pa_built_10y；地理=经纬度+Köppen哑变量\nX=归一化热事件指数[0,1]（去趋势）",
      x = "方差来源", y = "响应类型"
    ) +
    theme_cn() +
    theme(panel.grid.major = element_line(color = "grey92"))

  ggsave(file.path(OUT, "varpart_by_stype.png"),
         p_vp_stype, width = 9, height = 7, dpi = 300)
  cat("-> varpart_by_stype.png\n")
}


# K. 分层方差分解 2：按气候区分层（站点级）
# -----------------------------------------------------------------------------

cat("\n=== K. 按气候区分层的方差分解 ===\n")

run_vp_geo <- function(df_input, grp_label, outcome_col = "TRRI") {
  d <- df_input %>%
    filter(!is.na(.data[[outcome_col]]), !is.na(pa_built_10y),
           !is.na(longitude), !is.na(latitude), pa_built_10y > 0)
  if (nrow(d) < 15) {
    cat(sprintf("    [%s 跳过] n=%d < 15\n", grp_label, nrow(d)))
    return(NULL)
  }
  vp <- tryCatch(
    vegan::varpart(d[[outcome_col]],
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) { cat("    varpart error:", e$message, "\n"); NULL }
  )
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(
    koppen_group = grp_label, n = nrow(d),
    invest_R2    = round(fr$Adj.R.square[1] * 100, 2),
    geo_R2       = round(fr$Adj.R.square[2] * 100, 2),
    shared_R2    = round(fr$Adj.R.square[3] * 100, 2),
    unexplained  = round(fr$Adj.R.square[4] * 100, 2)
  )
}

koppen_groups <- sort(unique(anal_df$koppen_group))
vp_koppen <- map_dfr(koppen_groups, function(g) {
  df_g <- filter(anal_df, koppen_group == g)
  cat(sprintf("\n[%s] n=%d\n", g, nrow(df_g)))
  run_vp_geo(df_g, g)
})

vp_all_latlon  <- run_vp_geo(anal_df, "ALL")
vp_koppen_full <- bind_rows(
  if (!is.null(vp_all_latlon)) mutate(vp_all_latlon, koppen_group = "ALL"),
  vp_koppen
)

write_csv(vp_koppen_full, file.path(OUT, "varpart_by_koppen.csv"))
cat("\n"); print(vp_koppen_full)

if (nrow(vp_koppen_full) > 0) {
  koppen_order <- intersect(c("ALL","A","B","C","D"), vp_koppen_full$koppen_group)
  koppen_name  <- c(ALL="全部", A="A(热带)", B="B(干旱)",
                    C="C(温带)", D="D(大陆)")

  p_vp_koppen <- vp_koppen_full %>%
    pivot_longer(c(invest_R2, geo_R2, shared_R2),
                 names_to = "comp", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(comp,
                          levels = c("invest_R2", "shared_R2", "geo_R2"),
                          labels = c("投资", "共享", "地理")),
      grp_label  = factor(
        dplyr::recode(koppen_group, !!!koppen_name),
        levels = koppen_name[koppen_order]
      )
    ) %>%
    filter(!is.na(grp_label)) %>%
    ggplot(aes(x = comp_label, y = r2_show, fill = comp_label)) +
    geom_col(width = 0.55, alpha = 0.85, position = "dodge") +
    geom_text(aes(label = sprintf("%.1f%%", r2_show)),
              vjust = -0.4, size = 3.8, family = "heiti") +
    scale_fill_manual(
      values = c("投资" = "#D73027", "共享" = "#FDAE61", "地理" = "#4575B4"),
      name = "方差来源") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~grp_label, nrow = 1) +
    labs(
      title    = "分层方差分解：气候区分层（站点级）",
      subtitle = "X=归一化热事件指数[0,1]（去趋势）；因变量=TRRI；地理=经纬度",
      x = NULL, y = "调整R²（%）"
    ) +
    theme_cn() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(OUT, "varpart_by_koppen.png"),
         p_vp_koppen, width = 12, height = 6, dpi = 300)
  cat("-> varpart_by_koppen.png\n")
}


# K2. 气候区分层：站点级 vs 城市级对比
# -----------------------------------------------------------------------------

cat("\n=== K2. 气候区分层：站点级 vs 城市级对比 ===\n")

city_agg_k <- anal_df %>%
  group_by(city_name) %>%
  summarise(
    TRRI_mean    = mean(TRRI,        na.rm = TRUE),
    pa_built_10y = first(pa_built_10y),
    longitude    = mean(longitude,   na.rm = TRUE),
    latitude     = mean(latitude,    na.rm = TRUE),
    koppen_group = first(koppen_group),
    .groups      = "drop"
  ) %>%
  filter(!is.na(TRRI_mean), !is.na(pa_built_10y),
         !is.na(longitude),  !is.na(latitude))

run_vp_city <- function(df, grp_label) {
  d <- df %>% filter(!is.na(pa_built_10y), pa_built_10y > 0)
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(koppen_group = grp_label, n = nrow(d),
         invest_R2   = round(fr2$Adj.R.square[1]*100, 2),
         geo_R2      = round(fr2$Adj.R.square[2]*100, 2),
         shared_R2   = round(fr2$Adj.R.square[3]*100, 2),
         unexplained = round(fr2$Adj.R.square[4]*100, 2))
}

vp_city_koppen_full <- bind_rows(
  run_vp_city(city_agg_k, "ALL"),
  map_dfr(sort(unique(city_agg_k$koppen_group)),
          ~run_vp_city(filter(city_agg_k, koppen_group == .x), .x))
)

koppen_compare <- bind_rows(
  vp_koppen_full %>% mutate(level = "站点级"),
  vp_city_koppen_full %>% mutate(level = "城市级（聚合）")
) %>%
  mutate(level = factor(level, levels = c("站点级", "城市级（聚合）")))

write_csv(koppen_compare, file.path(OUT, "varpart_koppen_compare.csv"))
cat("\n站点级 vs 城市级气候区对比：\n")
print(koppen_compare %>%
        dplyr::select(koppen_group, level, n, invest_R2, geo_R2, unexplained))

if (nrow(koppen_compare) > 0) {
  koppen_order2 <- intersect(c("ALL","A","B","C","D"), koppen_compare$koppen_group)
  koppen_name2  <- c(ALL="全部", A="A(热带)", B="B(干旱)",
                     C="C(温带)", D="D(大陆)")

  p_koppen_cmp <- koppen_compare %>%
    pivot_longer(c(invest_R2, geo_R2),
                 names_to = "comp", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      r2_label   = sprintf("%.2f%%", r2),
      comp_label = factor(comp,
                          levels = c("invest_R2", "geo_R2"),
                          labels = c("投资（独立）", "地理（独立）")),
      grp_label  = factor(
        dplyr::recode(koppen_group, !!!koppen_name2),
        levels = koppen_name2[koppen_order2]
      )
    ) %>%
    filter(!is.na(grp_label)) %>%
    ggplot(aes(x = level, y = r2_show, fill = comp_label)) +
    geom_col(width = 0.6, alpha = 0.85,
             position = position_dodge(width = 0.7)) +
    geom_text(aes(label = r2_label),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3.5, family = "heiti") +
    scale_fill_manual(
      values = c("投资（独立）" = "#D73027", "地理（独立）" = "#4575B4"),
      name = "方差来源") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, NA)) +
    facet_wrap(~grp_label, nrow = 1) +
    labs(
      title    = "气候区分层方差分解：站点级 vs 城市级聚合",
      subtitle = "X=归一化热事件指数[0,1]；投资=pa_built_10y；地理=经纬度；标注值含负数",
      x = NULL, y = "调整R²（%）"
    ) +
    theme_cn() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 15, hjust = 1))

  ggsave(file.path(OUT, "varpart_koppen_compare.png"),
         p_koppen_cmp, width = 13, height = 6, dpi = 300)
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
    pa_built_10y = first(pa_built_10y),
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

geo_vars_city <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

vp_city <- vegan::varpart(
  city_agg$TRRI_mean,
  dplyr::select(city_agg, pa_built_10y),
  dplyr::select(city_agg, all_of(geo_vars_city))
)

fr_city <- vp_city$part$indfract
vp_city_tbl <- tibble(
  level      = "城市级（聚合）",
  component  = c("投资（独立）", "地理（独立）", "共享部分", "未解释"),
  adj_R2_pct = round(fr_city$Adj.R.square * 100, 2)
)
cat("\n"); print(vp_city_tbl)

compare_tbl <- bind_rows(
  tibble(level="站点级（原始）",
         invest_R2 = if (!is.null(fr)) fr$Adj.R.square[1]*100 else NA,
         geo_R2    = if (!is.null(fr)) fr$Adj.R.square[2]*100 else NA,
         shared_R2 = if (!is.null(fr)) fr$Adj.R.square[3]*100 else NA,
         n         = nrow(anal_df)),
  tibble(level="城市级（聚合）",
         invest_R2 = fr_city$Adj.R.square[1]*100,
         geo_R2    = fr_city$Adj.R.square[2]*100,
         shared_R2 = fr_city$Adj.R.square[3]*100,
         n         = nrow(city_agg))
) %>% mutate(across(c(invest_R2, geo_R2, shared_R2), ~round(.x, 2)))

cat("\n--- 站点级 vs 城市级对比 ---\n")
print(compare_tbl)
write_csv(compare_tbl, file.path(OUT, "varpart_city_vs_station.csv"))

vp_city_koppen <- map_dfr(sort(unique(city_agg$koppen_group)), function(g) {
  d <- filter(city_agg, koppen_group == g)
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(koppen=g, n=nrow(d),
         invest_R2 = round(fr2$Adj.R.square[1]*100, 2),
         geo_R2    = round(fr2$Adj.R.square[2]*100, 2),
         shared_R2 = round(fr2$Adj.R.square[3]*100, 2))
})
print(vp_city_koppen)
write_csv(vp_city_koppen, file.path(OUT, "varpart_city_koppen.csv"))

p_compare <- compare_tbl %>%
  pivot_longer(c(invest_R2, geo_R2, shared_R2),
               names_to = "comp", values_to = "r2") %>%
  mutate(
    r2_show    = pmax(r2, 0),
    comp_label = factor(comp,
                        levels = c("invest_R2","shared_R2","geo_R2"),
                        labels = c("投资","共享","地理")),
    level      = factor(level, levels = c("站点级（原始）","城市级（聚合）"))
  ) %>%
  ggplot(aes(x = comp_label, y = r2_show, fill = comp_label)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f%%", r2_show)),
            vjust = -0.4, size = 5, family = "heiti") +
  scale_fill_manual(
    values = c("投资"="#D73027","共享"="#FDAE61","地理"="#4575B4"),
    guide  = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~level) +
  labs(title    = "方差分解：站点级 vs 城市级聚合",
       subtitle = "X=归一化热事件指数[0,1]；消除城市内伪重复后投资解释力的变化",
       x = NULL, y = "调整R²（%）") +
  theme_cn()

ggsave(file.path(OUT, "varpart_city_compare.png"),
       p_compare, width = 10, height = 6, dpi = 300)
cat("-> varpart_city_compare.png\n")


# M. 城市内部站点异质性可视化（单张多行拼图）
# -----------------------------------------------------------------------------
# 所有城市排成多行，输出为单张图（city_stype_all.png）
# 每行最多 N_ROW_WIDTH 个城市；颜色=响应类型；*=城市内有多种类型

cat("\n=== M. 城市内部站点异质性可视化（单张多行拼图）===\n")

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

n_city       <- nrow(city_order)
N_ROW_WIDTH  <- 50   # 每行城市数
n_rows       <- ceiling(n_city / N_ROW_WIDTH)

cat(sprintf("城市总数: %d，排成 %d 行（每行≤%d城市）\n",
            n_city, n_rows, N_ROW_WIDTH))

# 为每个城市计算所在行和行内位置
city_order <- city_order %>%
  mutate(
    row_id  = ceiling(city_rank / N_ROW_WIDTH),
    col_pos = city_rank - (row_id - 1) * N_ROW_WIDTH
  )

# 构建绘图数据：每城市×每类型的比例
df_plot <- city_stype %>%
  group_by(city_name, stype_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(city_name) %>%
  mutate(pct = n / sum(n) * 100, n_total = sum(n)) %>%
  ungroup() %>%
  left_join(dplyr::select(city_order, city_name, koppen_group,
                           n_types, city_rank, row_id, col_pos),
            by = "city_name") %>%
  mutate(
    city_label = ifelse(n_types > 1,
                        paste0(as.character(city_name), "*"),
                        as.character(city_name))
  )

# 每行独立绘图，然后用 patchwork 垂直拼合
plots_m <- map(1:n_rows, function(rid) {
  df_r <- df_plot %>% filter(row_id == rid)

  # 该行城市的顺序（按col_pos）
  city_lv <- city_order %>%
    filter(row_id == rid) %>%
    arrange(col_pos) %>%
    pull(city_name)

  # city_label 顺序
  label_lv <- df_r %>%
    distinct(city_name, city_label, col_pos) %>%
    arrange(col_pos) %>%
    pull(city_label)

  df_r <- df_r %>%
    mutate(
      city_name  = factor(city_name,  levels = city_lv),
      city_label = factor(city_label, levels = unique(label_lv))
    )

  # 气候区分隔线位置（该行内）
  kopp_meta <- city_order %>%
    filter(row_id == rid) %>%
    group_by(koppen_group) %>%
    summarise(first_col = min(col_pos), .groups = "drop") %>%
    filter(first_col > 1) %>%   # 第一组不需要左侧线
    mutate(x_pos = first_col - 0.5)

  # 是否显示y轴标签（每行都有）
  p_row <- ggplot(df_r,
                  aes(x = city_label, y = pct, fill = stype_label)) +
    geom_col(width = 0.85, alpha = 0.9) +
    geom_vline(data = kopp_meta,
               aes(xintercept = x_pos),
               color = "grey40", linewidth = 0.6,
               linetype = "dashed", inherit.aes = FALSE) +
    # 气候区标签（标在分隔线右上方）
    geom_text(
      data = kopp_meta %>%
        mutate(koppen_label = paste0(koppen_group, "区")),
      aes(x = x_pos + 0.6, y = 97, label = koppen_label),
      inherit.aes = FALSE,
      size = 3.2, family = "heiti", color = "grey30", hjust = 0
    ) +
    # 站点总数标注（顶部）
    geom_text(
      data = df_r %>% distinct(city_label, n_total, col_pos),
      aes(x = city_label, y = 104, label = n_total, fill = NULL),
      size = 2.5, family = "heiti", color = "grey40"
    ) +
    scale_fill_manual(
      values = setNames(stype_colors, stype_labels_map),
      name   = "响应类型"
    ) +
    scale_y_continuous(
      limits = c(0, 108),
      breaks = c(0, 50, 100),
      labels = c("0%", "50%", "100%"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL,
         y = if (rid == ceiling(n_rows / 2)) "站点比例" else NULL) +
    theme_cn(bs = 10) +
    theme(
      axis.text.x        = element_text(angle = 55, hjust = 1, size = 6.5),
      axis.text.y        = element_text(size = 8),
      legend.position    = "none",
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(2, 4, 2, 4)
    )

  p_row
})

# 用 patchwork 垂直拼合，共享图例
legend_plot <- ggplot(
  tibble(x=1, y=1,
         stype_label=factor("全程抑制", levels=stype_labels_map)),
  aes(x=x, y=y, fill=stype_label)
) +
  geom_col() +
  scale_fill_manual(values = setNames(stype_colors, stype_labels_map),
                    name = "响应类型") +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(family = "heiti", size = 11),
    legend.title     = element_text(family = "heiti", size = 11, face = "bold"),
    legend.key.size  = unit(0.5, "cm")
  )
shared_legend <- cowplot::get_legend(legend_plot)

# patchwork 拼图（加载 cowplot 用于提取图例）
pacman::p_load(cowplot)

combined <- wrap_plots(plots_m, ncol = 1) +
  plot_annotation(
    title    = "城市内部站点热响应类型分布",
    subtitle = sprintf(
      "X=归一化热事件指数[0,1]；共%d城市；* = 城市内存在多种类型；虚线 = 气候区边界；顶部数字 = 站点数",
      n_city),
    theme = theme(
      plot.title    = element_text(family = "heiti", face = "bold",
                                   hjust = 0.5, size = 18),
      plot.subtitle = element_text(family = "heiti", hjust = 0.5,
                                   color = "grey40", size = 13)
    )
  )

# 带图例的最终图
final_m <- plot_grid(
  combined,
  shared_legend,
  ncol    = 1,
  rel_heights = c(1, 0.03)
)

# 自动计算高度：每行约3英寸
fig_height <- n_rows * 3.2 + 1.5
ggsave(file.path(OUT, "city_stype_all.png"),
       final_m,
       width  = 18,
       height = fig_height,
       dpi    = 300,
       limitsize = FALSE)
cat(sprintf("-> city_stype_all.png  (%d行，%.0f×%.0f英寸)\n",
            n_rows, 18, fig_height))

# 城市异质性汇总CSV
city_diversity <- city_order %>%
  left_join(
    city_stype %>%
      group_by(city_name, stype) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = stype, values_from = n, values_fill = 0),
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
