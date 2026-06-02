# =============================================================================
# 04_analysis.R  ——  完整独立分析脚本
#
# X变量   ：VPD超量累积（heat_over_dt，去趋势）
# 投资变量：近10年均值 / 建成区绿地（pa_built_10y）
# 方差分解：投资（1组） vs 地理（经纬度 + Köppen气候区，1组）
#
# 运行顺序：
#   A. 加载数据 & 去趋势
#   B. SIF结构突变检测（CCM前过滤）& 可视化
#   C. CCM 分析（heat_over_dt × tp=0..8，全量站点）
#   D. TRRI 分类
#   E. 投资变量构建（pa_built_10y）
#   F. 合并分析数据框
#   G. OLS 回归
#   H. 方差分解（投资 vs 地理）
#   I. 图表 & 汇总输出
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr,
               vegan, tibble, scales, showtext, sysfonts,
               targets, rEDM, strucchange)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT     <- "data_proc/output_10y_built_up_04"
CCM_RDS <- "data_proc/output_10y_built_up_04/ccm_04_results.rds"
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

cat("去趋势 heat_over_sum 和 SIF...\n")
df_weekly <- data_heat_sif_weekly %>%
  group_by(meteo_stat_id) %>%
  arrange(year, week) %>%
  mutate(
    time_idx     = row_number(),
    heat_over_dt = safe_detrend(heat_over_sum, time_idx),
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

break_res <- map_dfr(all_sids, ~detect_sif_break(.x, df_weekly))
saveRDS(break_res, file.path(OUT, "sif_break_results.rds"))

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


# C. CCM 分析（heat_over_dt × tp=0..8，全量站点）
# -----------------------------------------------------------------------------

# 过滤突变站点，筛选有效站点
stations_clean <- break_res %>%
  filter(has_break == 0L) %>%
  pull(meteo_stat_id)

stations_ok <- df_weekly %>%
  filter(week %in% 20:39, meteo_stat_id %in% stations_clean) %>%
  group_by(meteo_stat_id) %>%
  filter(sum(!is.na(heat_over_dt))  >= 40,
         sum(!is.na(sif_detrended)) >= 40) %>%
  summarise(.groups = "drop") %>%
  pull(meteo_stat_id)

cat(sprintf("\n突变过滤后有效站点: %d\n", length(stations_ok)))

# 若CCM结果已存在则跳过（节省重复运行时间）
if (file.exists(CCM_RDS)) {
  cat("发现已有CCM结果，直接读取:", CCM_RDS, "\n")
  ccm_results <- readRDS(CCM_RDS)
} else {
  perform_ccm <- function(station_id, df, tp_x = 0, min_pts = 40) {
    d <- df %>%
      filter(meteo_stat_id == station_id, week %in% 20:39) %>%
      arrange(year, week) %>%
      mutate(x_var = dplyr::lag(heat_over_dt, n = tp_x)) %>%
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
  cat(sprintf("开始CCM：%d站 × %d tp = %d次\n",
              n_stations, length(tp_seq), total_calls))

  ccm_results <- map_dfr(tp_seq, function(tp_now) {
    cat(sprintf("  tp=%d ...\n", tp_now))
    map_dfr(stations_ok, ~perform_ccm(.x, df_weekly, tp_now))
  })

  ccm_results <- ccm_results %>%
    left_join(meta, by = "meteo_stat_id")

  saveRDS(ccm_results, CCM_RDS)
  cat(sprintf("CCM结果已保存 (%d行)\n", nrow(ccm_results)))
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


# E. 投资变量：pa_built_10y（近10年均值 / 建成区绿地）
# -----------------------------------------------------------------------------

tar_load(green_invest_2020)
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

cols_10y <- year_cols[years_num >= 2011 & years_num <= 2020]
inv_10y  <- rowMeans(
  dplyr::select(invest_raw, all_of(cols_10y)) %>%
    mutate(across(everything(), as.numeric)),
  na.rm = TRUE
)

green <- green_area_2020 %>%
  dplyr::select(city_name, area_built = area_green_built) %>%
  mutate(city_name = ifelse(str_detect(city_name, "市$"), city_name,
                            paste0(city_name, "市")))

invest_tbl <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  tibble(inv_10y = inv_10y)
) %>%
  left_join(green, by = "city_name") %>%
  mutate(pa_built_10y = inv_10y / area_built)

invest_station <- station_city_map %>%
  left_join(invest_tbl, by = "city_name")

cat(sprintf("\npa_built_10y 非NA站点: %d\n",
            sum(!is.na(invest_station$pa_built_10y))))


# F. 合并分析数据框
# -----------------------------------------------------------------------------

anal_df <- trri_df %>%
  left_join(invest_station, by = "meteo_stat_id") %>%
  filter(!is.na(TRRI), !is.na(pa_built_10y),
         !is.na(longitude), !is.na(latitude),
         !is.na(koppen_group),
         pa_built_10y > 0,
         pa_built_10y < quantile(pa_built_10y, 0.99, na.rm = TRUE)) %>%
  mutate(
    koppen_B = as.integer(koppen_group == "B"),
    koppen_C = as.integer(koppen_group == "C"),
    koppen_D = as.integer(koppen_group == "D")
  )

cat(sprintf("\n进入回归/方差分解的站点数: %d\n", nrow(anal_df)))
cat("气候组分布:\n"); print(count(anal_df, koppen_group))


# G. OLS 回归
# -----------------------------------------------------------------------------

m_ols <- lm(TRRI ~ pa_built_10y + longitude + latitude +
              koppen_B + koppen_C + koppen_D,
            data = anal_df)

cat("\n=== OLS 回归（因变量：TRRI）===\n")
print(summary(m_ols))

s_ols  <- summary(m_ols)
cr_inv <- s_ols$coefficients["pa_built_10y", ]
sig_inv <- case_when(cr_inv[4] < 0.001 ~ "***", cr_inv[4] < 0.01 ~ "**",
                     cr_inv[4] < 0.05  ~ "*",   cr_inv[4] < 0.1  ~ ".", TRUE ~ "ns")
cat(sprintf("\n【投资系数】β=%.4f  SE=%.4f  t=%.3f  p=%.4f%s\n",
            cr_inv[1], cr_inv[2], cr_inv[3], cr_inv[4], sig_inv))


# H. 方差分解（投资 vs 地理）
# -----------------------------------------------------------------------------
# 组1：投资  → pa_built_10y
# 组2：地理  → longitude + latitude + koppen_B + koppen_C + koppen_D

geo_vars <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

vp <- vegan::varpart(
  anal_df$TRRI,
  dplyr::select(anal_df, pa_built_10y),
  dplyr::select(anal_df, all_of(geo_vars))
)

cat("\n=== 方差分解（投资 vs 地理）===\n")
print(vp)

fr     <- vp$part$indfract
vp_tbl <- tibble(
  component  = c("投资（独立）", "地理（独立）", "共享部分", "未解释"),
  adj_R2_pct = round(fr$Adj.R.square * 100, 2)
)
cat("\n"); print(vp_tbl)
write_csv(vp_tbl, file.path(OUT, "varpart_04.csv"))


# I. 图表 & 汇总输出
# -----------------------------------------------------------------------------

# 图1：SIF突变 vs 正常（已在B节生成）

# 图2：散点图（pa_built_10y vs TRRI，按气候区着色）
p_scatter <- ggplot(anal_df,
                    aes(x = pa_built_10y, y = TRRI, color = koppen_group)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linetype = "dashed", linewidth = 0.9) +
  scale_color_brewer(palette = "Set1", name = "Köppen气候区") +
  scale_y_continuous(breaks = seq(1, TRRI_MAX, by = 2)) +
  labs(
    title    = "绿地投资强度与热响应韧性指数（TRRI）",
    subtitle = sprintf("X=VPD超量累积 | 投资=近10年均值/建成区绿地 | n=%d站",
                       nrow(anal_df)),
    x        = "单位面积绿地投资强度（pa_built_10y）",
    y        = sprintf("TRRI（1–%d）", TRRI_MAX)
  ) +
  theme_cn()

ggsave(file.path(OUT, "scatter_trri_invest.png"),
       p_scatter, width = 10, height = 7, dpi = 300)
cat("-> scatter_trri_invest.png\n")

# 图3：方差分解条形图
p_vp <- vp_tbl %>%
  filter(component != "未解释") %>%
  mutate(
    r2_show   = pmax(adj_R2_pct, 0),
    component = factor(component,
                       levels = c("投资（独立）", "共享部分", "地理（独立）"))
  ) %>%
  ggplot(aes(x = component, y = r2_show, fill = component)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f%%", r2_show)),
            vjust = -0.4, size = 5.5, family = "heiti") +
  scale_fill_manual(
    values = c("投资（独立）" = "#D73027",
               "共享部分"    = "#FDAE61",
               "地理（独立）" = "#4575B4"),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "方差分解：投资 vs 地理对TRRI的独立贡献",
    subtitle = "地理 = 经纬度 + Köppen气候区哑变量（B/C/D，以A为参照）",
    x        = NULL,
    y        = "调整R²（%）"
  ) +
  theme_cn()

ggsave(file.path(OUT, "varpart_bar_04.png"),
       p_vp, width = 8, height = 6, dpi = 300)
cat("-> varpart_bar_04.png\n")

# 图4：S-map系数轨迹（4类站点均值±SE）
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
    subtitle = "X=VPD超量累积；色带=±1SE；tp=滞后周数",
    x        = "时间滞后 tp（周）",
    y        = "S-map系数（热→SIF）"
  ) +
  theme_cn() +
  theme(legend.position = "right")

ggsave(file.path(OUT, "coef_trajectory.png"),
       p_coef, width = 11, height = 7, dpi = 300)
cat("-> coef_trajectory.png\n")

# 图5：站点类型分布饼/条形
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
  scale_y_continuous(limits = c(0, 60),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title    = "热响应类型分布",
       subtitle = sprintf("共%d站（X=VPD超量累积，tp=0..%d）",
                          nrow(trri_df), MAX_TP),
       x = NULL, y = "比例 (%)") +
  theme_cn()

ggsave(file.path(OUT, "stype_dist.png"),
       p_stype, width = 9, height = 6, dpi = 300)
cat("-> stype_dist.png\n")

# 汇总CSV
result_summary <- tibble(
  x_var         = "heat_over_dt",
  inv_var       = "pa_built_10y",
  n_ccm         = length(unique(ccm_results$meteo_stat_id)),
  n_regression  = nrow(anal_df),
  beta          = cr_inv[1],
  se            = cr_inv[2],
  t_val         = cr_inv[3],
  p_val         = cr_inv[4],
  sig           = sig_inv,
  r2_adj_full   = s_ols$adj.r.squared,
  invest_R2_pct = round(fr$Adj.R.square[1] * 100, 2),
  geo_R2_pct    = round(fr$Adj.R.square[2] * 100, 2),
  shared_R2_pct = round(fr$Adj.R.square[3] * 100, 2)
)
write_csv(result_summary, file.path(OUT, "summary_04.csv"))
cat("-> summary_04.csv\n")

# 终端摘要
cat("\n", strrep("=", 60), "\n")
cat("分析摘要\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("X变量      : heat_over_dt（VPD超量累积，去趋势）\n"))
cat(sprintf("投资变量   : pa_built_10y（近10年均值/建成区绿地）\n"))
cat(sprintf("TRRI量表   : 1..%d（N_TP=%d, tp=0..%d）\n",
            TRRI_MAX, N_TP, MAX_TP))
cat(sprintf("CCM站点数  : %d\n", length(unique(ccm_results$meteo_stat_id))))
cat(sprintf("回归站点数 : %d\n\n", nrow(anal_df)))
cat(sprintf("【OLS回归】 β=%.4f  p=%.4f%s  全模型R²=%.3f\n\n",
            cr_inv[1], cr_inv[4], sig_inv, s_ols$adj.r.squared))
cat(sprintf("【方差分解】\n  投资独立贡献 = %.2f%%\n  地理独立贡献 = %.2f%%\n  共享部分     = %.2f%%\n",
            fr$Adj.R.square[1] * 100,
            fr$Adj.R.square[2] * 100,
            fr$Adj.R.square[3] * 100))
cat("\n输出目录:", OUT, "\n")
cat("文件列表:\n")
cat(paste(" -", list.files(OUT)), sep = "\n")


# J. 分层方差分解 1：按响应类型分层
# -----------------------------------------------------------------------------
# 按站点类型分层；各层因变量不同（因为 always_* 组内TRRI无方差）：
#   inhibit_promote  → ITW（转折时间，越小越快恢复）
#   promote_inhibit  → HTW（促进持续时间，越大越好）
#   always_inhibit   → tp=0 时 mean_coef 的绝对值（抑制强度）
#   always_promote   → tp=0 时 mean_coef（促进强度）
# 投资组：pa_built_10y；地理组：longitude + latitude + Köppen哑变量

stype_levels <- c("always_promote", "promote_inhibit",
                  "inhibit_promote", "always_inhibit")
stype_labels <- c("全程促进", "促进→抑制", "抑制→促进", "全程抑制")

# 各类型对应的因变量名称和含义
stype_outcome <- list(
  inhibit_promote  = list(var = "ITW",       label = "ITW（转折周数，小=恢复快）"),
  promote_inhibit  = list(var = "HTW",       label = "HTW（促进周数，大=持续长）"),
  always_inhibit   = list(var = "inhibit_strength", label = "抑制强度（|coef@tp0|）"),
  always_promote   = list(var = "promote_strength", label = "促进强度（coef@tp0）")
)

# 补充额外因变量到 anal_df
tp0_coef <- ccm_results %>%
  filter(tp == 0) %>%
  dplyr::select(meteo_stat_id, coef_tp0 = mean_coef)

anal_df2 <- anal_df %>%
  left_join(tp0_coef, by = "meteo_stat_id") %>%
  mutate(
    inhibit_strength = -coef_tp0,   # 正值=抑制越强
    promote_strength =  coef_tp0
  )

run_vp2 <- function(df, outcome_var, label) {
  d <- df %>%
    filter(!is.na(.data[[outcome_var]]),
           !is.na(pa_built_10y), !is.na(longitude), !is.na(latitude),
           pa_built_10y > 0) %>%
    droplevels()
  if (nrow(d) < 15) {
    cat(sprintf("    [跳过] n=%d < 15\n", nrow(d)))
    return(NULL)
  }
  # 地理组：经纬度 + Köppen哑变量（若组内有多个气候区）
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
  cat(sprintf("\n[%s] n=%d, 因变量=%s\n",
              st, nrow(df_s), info$label))
  res <- run_vp2(df_s, info$var, info$label)
  if (!is.null(res)) mutate(res, stype = st)
})

if (nrow(vp_stype) > 0) {
  vp_stype <- vp_stype %>%
    mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))
  write_csv(vp_stype, file.path(OUT, "varpart_by_stype.csv"))
  cat("\n"); print(vp_stype %>% dplyr::select(stype, label, n, invest_R2, geo_R2, shared_R2))

  # 图：分层方差分解气泡图（响应类型）
  p_vp_stype <- vp_stype %>%
    pivot_longer(c(invest_R2, geo_R2, shared_R2),
                 names_to = "comp", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(comp,
                          levels = c("invest_R2", "shared_R2", "geo_R2"),
                          labels = c("投资", "共享", "地理")),
      fill_col   = case_when(comp == "invest_R2" ~ "#D73027",
                             comp == "shared_R2" ~ "#FDAE61",
                             comp == "geo_R2"    ~ "#4575B4")
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
      subtitle = "投资=pa_built_10y；地理=经纬度+Köppen哑变量\n因变量：ITW / HTW / 抑制强度 / 促进强度",
      x = "方差来源", y = "响应类型"
    ) +
    theme_cn() +
    theme(panel.grid.major = element_line(color = "grey92"))

  ggsave(file.path(OUT, "varpart_by_stype.png"),
         p_vp_stype, width = 9, height = 7, dpi = 300)
  cat("-> varpart_by_stype.png\n")
}


# K. 分层方差分解 2：按气候区分层
# -----------------------------------------------------------------------------
# 因变量：TRRI；投资组：pa_built_10y；地理组：longitude + latitude（仅经纬度）

cat("\n=== K. 按气候区分层的方差分解 ===\n")

run_vp_geo <- function(df, grp_label) {
  d <- df %>%
    filter(!is.na(TRRI), !is.na(pa_built_10y),
           !is.na(longitude), !is.na(latitude),
           pa_built_10y > 0)
  if (nrow(d) < 15) {
    cat(sprintf("    [%s 跳过] n=%d < 15\n", grp_label, nrow(d)))
    return(NULL)
  }
  vp <- tryCatch(
    vegan::varpart(d$TRRI,
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) { cat("    varpart error:", e$message, "\n"); NULL }
  )
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(
    koppen_group = grp_label,
    n            = nrow(d),
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

# 补上全样本行
vp_all_latlon <- run_vp_geo(anal_df, "ALL")
vp_koppen_full <- bind_rows(
  if (!is.null(vp_all_latlon)) mutate(vp_all_latlon, koppen_group = "ALL"),
  vp_koppen
)

write_csv(vp_koppen_full, file.path(OUT, "varpart_by_koppen.csv"))
cat("\n"); print(vp_koppen_full)

if (nrow(vp_koppen_full) == 0) {
  cat("警告：所有气候区 varpart 均失败，跳过绘图。\n")
} else {

# 图：分层方差分解条形图（气候区）
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
    title    = "分层方差分解：气候区分层",
    subtitle = "因变量=TRRI；投资=pa_built_10y；地理=经纬度（不含气候区哑变量）",
    x = NULL, y = "调整R²（%）"
  ) +
  theme_cn() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT, "varpart_by_koppen.png"),
       p_vp_koppen, width = 12, height = 6, dpi = 300)
cat("-> varpart_by_koppen.png\n")

} # end if nrow(vp_koppen_full) > 0


# 终端汇总补充
cat("\n=== J+K 分层解释力汇总 ===\n")
cat("\n[响应类型分层]\n")
if (exists("vp_stype") && nrow(vp_stype) > 0)
  print(vp_stype %>% dplyr::select(stype, n, invest_R2, geo_R2, shared_R2))
cat("\n[气候区分层]\n")
print(vp_koppen_full)

cat("\n全部输出文件:\n")
cat(paste(" -", list.files(OUT)), sep = "\n")


# L. 城市级聚合 + 方差分解（消除空间伪重复）
# -----------------------------------------------------------------------------
# 每城市取 TRRI 均值，投资变量已是城市级，地理取站点均值

cat("\n=== L. 城市级聚合方差分解 ===\n")

# 城市级汇总（city_name 已在 anal_df 中，来自 E 段投资连接）
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
cat("每城市站点数分布:\n")
print(table(city_agg$n_stations))

# 城市级方差分解：投资 vs 地理（经纬度 + Köppen哑变量）
geo_vars_city <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")

vp_city <- vegan::varpart(
  city_agg$TRRI_mean,
  dplyr::select(city_agg, pa_built_10y),
  dplyr::select(city_agg, all_of(geo_vars_city))
)

cat("\n--- 城市级方差分解（投资 vs 地理）---\n")
print(vp_city)

fr_city  <- vp_city$part$indfract
vp_city_tbl <- tibble(
  level      = "城市级（聚合）",
  component  = c("投资（独立）", "地理（独立）", "共享部分", "未解释"),
  adj_R2_pct = round(fr_city$Adj.R.square * 100, 2)
)
print(vp_city_tbl)

# 与站点级对比
compare_tbl <- bind_rows(
  tibble(level="站点级（原始）",
         invest_R2 = fr$Adj.R.square[1]*100,
         geo_R2    = fr$Adj.R.square[2]*100,
         shared_R2 = fr$Adj.R.square[3]*100,
         n         = nrow(anal_df)),
  tibble(level="城市级（聚合）",
         invest_R2 = fr_city$Adj.R.square[1]*100,
         geo_R2    = fr_city$Adj.R.square[2]*100,
         shared_R2 = fr_city$Adj.R.square[3]*100,
         n         = nrow(city_agg))
) %>% mutate(across(c(invest_R2,geo_R2,shared_R2), ~round(.x,2)))

cat("\n--- 站点级 vs 城市级对比 ---\n")
print(compare_tbl)
write_csv(compare_tbl, file.path(OUT, "varpart_city_vs_station.csv"))

# 同样按气候区分层（城市级）
cat("\n--- 城市级 × 气候区分层 ---\n")
vp_city_koppen <- map_dfr(c("B","C","D"), function(g) {
  d <- filter(city_agg, koppen_group == g)
  cat(sprintf("  [%s] n=%d 城市\n", g, nrow(d)))
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(koppen=g, n=nrow(d),
         invest_R2 = round(fr2$Adj.R.square[1]*100,2),
         geo_R2    = round(fr2$Adj.R.square[2]*100,2),
         shared_R2 = round(fr2$Adj.R.square[3]*100,2))
})
print(vp_city_koppen)
write_csv(vp_city_koppen, file.path(OUT, "varpart_city_koppen.csv"))

# 城市级方差分解条形图（站点级 vs 城市级）
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
       subtitle = "消除城市内伪重复后投资解释力的变化",
       x = NULL, y = "调整R²（%）") +
  theme_cn()

ggsave(file.path(OUT, "varpart_city_compare.png"),
       p_compare, width = 10, height = 6, dpi = 300)
cat("-> varpart_city_compare.png\n")


# M. 城市内部站点异质性可视化（分面板输出）
# -----------------------------------------------------------------------------
# 每城市：各站点的响应类型（stype）堆叠条形图
# 城市按气候区+站点数排序，每面板约50个城市

cat("\n=== M. 城市内部站点异质性可视化 ===\n")

# 准备数据：站点级 stype，附城市信息
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

# 城市排序：先按气候区，再按主导类型，再按站点数
city_order <- city_stype %>%
  group_by(city_name) %>%
  summarise(
    n_total     = n(),
    n_types     = n_distinct(stype),        # 城市内类型数（多样性）
    dominant    = names(sort(table(stype), decreasing=TRUE))[1],
    koppen_group = first(koppen_group),
    .groups     = "drop"
  ) %>%
  arrange(koppen_group, dominant, desc(n_total)) %>%
  mutate(city_rank = row_number())

# 划分面板（每面板 N_PER 个城市）
N_PER  <- 50
n_city <- nrow(city_order)
n_panel <- ceiling(n_city / N_PER)

city_order <- city_order %>%
  mutate(panel = ceiling(city_rank / N_PER))

cat(sprintf("城市总数: %d，分 %d 个面板（每面板≤%d城市）\n",
            n_city, n_panel, N_PER))

# 各面板分别绘图
walk(1:n_panel, function(pid) {
  cities_p <- city_order %>% filter(panel == pid) %>%
    arrange(city_rank)
  city_levels <- cities_p$city_name

  df_p <- city_stype %>%
    filter(city_name %in% city_levels) %>%
    group_by(city_name, stype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(city_name) %>%
    mutate(pct = n / sum(n) * 100,
           n_total = sum(n)) %>%
    ungroup() %>%
    left_join(dplyr::select(city_order, city_name, koppen_group, n_types),
              by = "city_name") %>%
    mutate(
      city_name = factor(city_name, levels = city_levels),
      # 城市标签：加星号表示内部有多种类型
      city_label = ifelse(n_types > 1,
                          paste0(as.character(city_name), " *"),
                          as.character(city_name))
    )

  # 重建 city_label 的 factor 顺序
  label_order <- df_p %>%
    distinct(city_name, city_label) %>%
    arrange(match(city_name, city_levels)) %>%
    pull(city_label)
  df_p <- df_p %>%
    mutate(city_label = factor(city_label, levels = label_order))

  # 气候区分隔线位置
  koppen_breaks <- cities_p %>%
    group_by(koppen_group) %>%
    summarise(first_rank = min(city_rank), .groups="drop") %>%
    mutate(x_pos = match(
      cities_p$city_name[cities_p$city_rank == first_rank],
      city_levels
    ) - 0.5) %>%
    filter(x_pos > 0.5)

  p <- ggplot(df_p, aes(x = city_label, y = pct, fill = stype_label)) +
    geom_col(width = 0.8, alpha = 0.9) +
    geom_vline(xintercept = koppen_breaks$x_pos,
               color = "grey30", linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = setNames(stype_colors, stype_labels_map),
                      name = "响应类型") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = sprintf("城市内部站点热响应类型分布（面板 %d/%d）",
                         pid, n_panel),
      subtitle = "* = 城市内存在多种类型；虚线 = 气候区边界；各城市站点数标注于顶部",
      x = NULL, y = "站点比例"
    ) +
    # 顶部标注站点总数
    geom_text(data = df_p %>% distinct(city_label, n_total),
              aes(x = city_label, y = 102, label = n_total, fill = NULL),
              size = 2.8, family = "heiti", color = "grey30") +
    theme_cn(bs = 11) +
    theme(
      axis.text.x    = element_text(angle = 55, hjust = 1, size = 7.5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )

  fname <- sprintf("city_stype_panel%02d.png", pid)
  ggsave(file.path(OUT, fname), p, width = 16, height = 7, dpi = 300)
  cat(sprintf("-> %s\n", fname))
})

# 汇总：各城市异质性统计
city_diversity <- city_order %>%
  left_join(
    city_stype %>% group_by(city_name, stype) %>% summarise(n=n(),.groups="drop") %>%
      pivot_wider(names_from=stype, values_from=n, values_fill=0),
    by = "city_name"
  )
write_csv(city_diversity, file.path(OUT, "city_stype_diversity.csv"))
cat("-> city_stype_diversity.csv\n")

cat(sprintf("\n城市内类型数分布（n_types=1为同质，>1为异质）:\n"))
print(table(city_order$n_types))
cat(sprintf("异质城市（>1种类型）: %d / %d (%.1f%%)\n",
            sum(city_order$n_types > 1), nrow(city_order),
            mean(city_order$n_types > 1)*100))

cat("\n全部输出文件:\n")
cat(paste(" -", list.files(OUT)), sep = "\n")

# 城市级 × 气候区分层方差分解（地理仅含经纬度）
vp_city_koppen <- map_dfr(c("B","C","D"), function(g) {
  d <- filter(city_agg, koppen_group == g)
  cat(sprintf("[%s] n=%d 城市\n", g, nrow(d)))
  if (nrow(d) < 15) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI_mean,
                   dplyr::select(d, pa_built_10y),
                   dplyr::select(d, longitude, latitude)),
    error = function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr2 <- vp$part$indfract
  tibble(koppen = g, n = nrow(d),
         invest_R2 = round(fr2$Adj.R.square[1]*100, 2),
         geo_R2    = round(fr2$Adj.R.square[2]*100, 2),
         shared_R2 = round(fr2$Adj.R.square[3]*100, 2))
})
print(vp_city_koppen)
write_csv(vp_city_koppen, file.path(OUT, "varpart_city_koppen.csv"))