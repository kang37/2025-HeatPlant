# =============================================================================
# 03_run_full.R
# 使用已有CCM结果(results_weekly_0_5.rds) → TRRI分类
# 新增：多年投资均值 + SIF结构突变检测（修复突然种树的问题）
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr, readxl,
               vegan, tibble, scales, showtext, sysfonts, targets,
               strucchange)   # 结构突变检测

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "analysis_vpd2/output"
dir.create(OUT, showWarnings = FALSE)

# 字体
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 16
theme_cn <- function(bs = BS) {
  theme_minimal(base_size = bs) +
    theme(text = element_text(family = "heiti"),
          plot.title   = element_text(face="bold", hjust=0.5, size=bs+4),
          plot.subtitle= element_text(hjust=0.5, color="grey40", size=bs),
          strip.text   = element_text(face="bold", size=bs),
          axis.text    = element_text(size=bs-2),
          axis.title   = element_text(size=bs),
          legend.text  = element_text(size=bs-2),
          legend.title = element_text(face="bold", size=bs))
}
theme_set(theme_cn())

# =============================================================================
# 1. 读取 CCM 结果
# =============================================================================
ccm <- readRDS("data_proc/results_weekly_0_5.rds")
cat("CCM结果:", nrow(ccm), "行，站点数:", n_distinct(ccm$meteo_stat_id), "\n")

# =============================================================================
# 2. SIF结构突变检测（修复"突然种树"问题）
# =============================================================================
# 问题：城市突然新增绿地会使SIF产生阶跃，线性去趋势无法消除这个阶跃，
#       导致CCM的S-map系数在跳跃点前后出现系统性偏差，混淆热事件的因果效应。
# 修复：用strucchange::Fstats检测每个站点SIF的结构突变；排除有显著突变的站点。
# =============================================================================

cat("\n检测SIF结构突变...\n")
tar_load(data_heat_sif_weekly)

detect_sif_break <- function(sid, df) {
  d <- df %>%
    filter(meteo_stat_id == sid, week %in% 20:39) %>%
    arrange(year, week) %>%
    filter(!is.na(sif_interp))
  if (nrow(d) < 30) return(tibble(meteo_stat_id = sid, has_break = NA,
                                   break_fstat = NA_real_, break_pval = NA_real_,
                                   break_magnitude = NA_real_))
  tryCatch({
    fs <- strucchange::Fstats(sif_interp ~ 1, data = d,
                              from = 0.15, to = 0.85)
    pv <- strucchange::sctest(fs)$p.value
    # 用 breakpoints 找突变点，计算突变幅度（相对于SIF标准差）
    bp <- strucchange::breakpoints(sif_interp ~ 1, data = d, breaks = 1)
    bp_idx <- bp$breakpoints
    if (!is.na(bp_idx) && bp_idx > 1 && bp_idx < nrow(d)) {
      mu1 <- mean(d$sif_interp[1:bp_idx],         na.rm=TRUE)
      mu2 <- mean(d$sif_interp[(bp_idx+1):nrow(d)], na.rm=TRUE)
      magnitude <- abs(mu2 - mu1) / sd(d$sif_interp, na.rm=TRUE)
    } else {
      magnitude <- 0
    }
    # 双重标准：显著 AND 幅度大（>1个SD），避免过度排除
    tibble(meteo_stat_id  = sid,
           has_break      = (pv < 0.01 & magnitude > 1.0),
           break_fstat    = max(fs$Fstats),
           break_pval     = pv,
           break_magnitude = magnitude)
  }, error = function(e)
    tibble(meteo_stat_id = sid, has_break = FALSE,
           break_fstat = NA_real_, break_pval = NA_real_,
           break_magnitude = NA_real_))
}

stations_all <- unique(ccm$meteo_stat_id)
break_res <- map_dfr(stations_all, ~detect_sif_break(.x, data_heat_sif_weekly))

n_break <- sum(break_res$has_break, na.rm = TRUE)
n_nobrk <- sum(!break_res$has_break, na.rm = TRUE)
cat(sprintf("  有显著SIF突变站点: %d (%.1f%%)\n", n_break,
            n_break/length(stations_all)*100))
cat(sprintf("  无突变站点（用于分析）: %d\n", n_nobrk))

# 排除有突变的站点
stations_clean <- break_res %>% filter(!is.na(has_break), !has_break) %>%
  pull(meteo_stat_id)

ccm_clean <- ccm %>% filter(meteo_stat_id %in% stations_clean)
cat("过滤后CCM行数:", nrow(ccm_clean), "| 站点数:", n_distinct(ccm_clean$meteo_stat_id), "\n")

# =============================================================================
# 3. TRRI 分类
# =============================================================================
find_transition <- function(coefs, to_neg = TRUE) {
  cond <- if (to_neg) function(x) x < 0 else function(x) x > 0
  if (length(coefs) < 2) return(NA_integer_)
  for (i in 2:length(coefs))
    if (!is.na(coefs[i]) && cond(coefs[i]))
      if (mean(sapply(coefs[i:length(coefs)], cond), na.rm=TRUE) >= 0.5)
        return(as.integer(i - 1L))
  NA_integer_
}

trri_df <- ccm_clean %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq     = list(mean_coef),
    longitude    = first(longitude),
    latitude     = first(latitude),
    koppen_class = first(koppen_class),
    koppen_group = first(koppen_group),
    .groups = "drop"
  ) %>%
  filter(map_lgl(coef_seq, ~length(.x) == 6)) %>%
  mutate(
    tp0_pos = map_lgl(coef_seq, ~.x[[1]] > 0),
    HTW     = map_int(coef_seq, ~find_transition(.x, to_neg=TRUE)),
    ITW     = map_int(coef_seq, ~find_transition(.x, to_neg=FALSE)),
    stype = case_when(
      tp0_pos  & is.na(HTW) ~ "always_promote",
      !tp0_pos & is.na(ITW) ~ "always_inhibit",
      tp0_pos  & !is.na(HTW) ~ "promote_inhibit",
      !tp0_pos & !is.na(ITW) ~ "inhibit_promote",
      TRUE ~ "always_inhibit"
    ),
    TRRI = case_when(
      stype == "always_inhibit"  ~ 1L,
      stype == "inhibit_promote" ~ as.integer(1L + (6L - ITW)),
      stype == "promote_inhibit" ~ as.integer(6L + HTW),
      stype == "always_promote"  ~ 12L
    )
  )

cat("\n站点类型分布（去除突变站点后）:\n")
print(count(trri_df, stype))

# =============================================================================
# 4. 构建多年投资变量
# =============================================================================
# 读取多年投资数据，并通过行位置与targets中的green_invest_2020对齐（同一CSV，行序一致）
tar_load(green_invest_2020)
invest_raw0 <- read.csv("data_raw/green_invest/city_invest_data.csv",
                        check.names=FALSE)
nc_inv <- ncol(invest_raw0)
names(invest_raw0) <- c("city_raw", paste0("inv_", 2002:(2001 + nc_inv - 1)))
# 按行绑定：city_name 来自 targets（UTF-8），年份数据来自 invest_raw0
invest_raw <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  dplyr::select(invest_raw0, -city_raw)
)

year_cols <- grep("^inv_\\d{4}$", names(invest_raw), value=TRUE)
years_num <- as.integer(str_extract(year_cols, "\\d{4}$"))

get_mean_inv <- function(df, y1, y2) {
  cols <- year_cols[years_num >= y1 & years_num <= y2]
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  rowMeans(
    mutate(dplyr::select(df, all_of(cols)), across(everything(), as.numeric)),
    na.rm = TRUE
  )
}

# 增长斜率（2005-2020）
slope_mat <- invest_raw %>%
  dplyr::select(all_of(year_cols[years_num %in% 2005:2020])) %>%
  mutate(across(everything(), as.numeric))
yr_sub <- years_num[years_num %in% 2005:2020]
inv_slope <- apply(slope_mat, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 5) return(NA_real_)
  coef(lm(x[ok] ~ yr_sub[ok]))[2]
})

invest_multi <- invest_raw %>%
  mutate(
    inv_2020     = as.numeric(inv_2020),
    inv_5y       = get_mean_inv(., 2016, 2020),
    inv_10y      = get_mean_inv(., 2011, 2020),
    inv_lag10    = get_mean_inv(., 2005, 2015),
    inv_slope    = inv_slope
  ) %>%
  dplyr::select(city_name, inv_2020, inv_5y, inv_10y, inv_lag10, inv_slope)

# GDP
gdp_raw <- read_excel("data_raw/china_city_db.xlsx")
# 列名是中文，用位置选取：col1=城市, 找年份列和GDP列
names(gdp_raw) <- iconv(names(gdp_raw), "UTF-8", "UTF-8")
# 直接用targets缓存的gdp对象
tar_load(gdp_data_2020)
gdp <- gdp_data_2020 %>%
  mutate(city_name=ifelse(str_detect(city_name,"市$"),city_name,paste0(city_name,"市")))

# 绿地面积 - 用targets缓存
tar_load(green_area_2020)
green <- green_area_2020 %>%
  dplyr::select(city_name, area_built=area_green_built, area_park=area_green_park) %>%
  mutate(city_name=ifelse(str_detect(city_name,"市$"),city_name,paste0(city_name,"市")))

# 站点→城市
tar_load(station_city_map)

invest_full <- station_city_map %>%
  left_join(invest_multi, by="city_name") %>%
  left_join(gdp,          by="city_name") %>%
  left_join(green,        by="city_name") %>%
  mutate(
    ratio_2020  = inv_2020  / gdp * 100,
    ratio_5y    = inv_5y    / gdp * 100,
    ratio_10y   = inv_10y   / gdp * 100,
    ratio_lag10 = inv_lag10 / gdp * 100,
    pa_built_10y = inv_10y  / area_built,
    pa_park_10y  = inv_10y  / area_park
  )

# =============================================================================
# 5. 气候背景变量（站点级多年均值）
# =============================================================================
climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  filter(meteo_stat_id %in% stations_clean) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    vpd_bg      = mean(vpd_mean,       na.rm=TRUE),
    heat_freq_bg = mean(heat_event_freq, na.rm=TRUE),
    .groups="drop"
  )

# 合并
anal_df <- trri_df %>%
  left_join(invest_full,  by="meteo_stat_id") %>%
  left_join(climate_bg,   by="meteo_stat_id") %>%
  filter(!is.na(TRRI), !is.na(koppen_group))

cat("\n分析数据框: ", nrow(anal_df), "行\n")

# =============================================================================
# 6. 回归网格：7个投资变量 × 全部/B/C/D气候组
# =============================================================================
inv_vars <- c("ratio_2020","ratio_5y","ratio_10y","ratio_lag10",
              "inv_slope","pa_built_10y","pa_park_10y")
inv_labels <- c(
  ratio_2020   = "2020单年/GDP(%)",
  ratio_5y     = "近5年均值/GDP(%)",
  ratio_10y    = "近10年均值/GDP(%)",
  ratio_lag10  = "滞后10年均值/GDP(%)",
  inv_slope    = "投资增长斜率",
  pa_built_10y = "10年均值/建成区绿地",
  pa_park_10y  = "10年均值/公园绿地"
)
ctrl_vars <- c("vpd_bg","heat_freq_bg","latitude")

run_ols <- function(df, inv) {
  d <- df %>%
    filter(!is.na(.data[[inv]]), if_all(all_of(ctrl_vars), ~!is.na(.x))) %>%
    filter(.data[[inv]] < quantile(.data[[inv]], 0.99, na.rm=TRUE),
           .data[[inv]] > 0)
  if (nrow(d) < 20) return(NULL)
  fm <- as.formula(paste("TRRI ~", inv, "+", paste(ctrl_vars, collapse="+")))
  m  <- tryCatch(lm(fm, data=d), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)
  cr <- s$coefficients[inv,,drop=FALSE]
  tibble(inv_var=inv, n=nrow(d), beta=cr[1,1], se=cr[1,2],
         t_val=cr[1,3], p_val=cr[1,4], r2_adj=s$adj.r.squared)
}

groups <- c("ALL","B","C","D")
reg_grid <- map_dfr(groups, function(grp) {
  df_g <- if (grp=="ALL") anal_df else filter(anal_df, koppen_group==grp)
  res <- map_dfr(inv_vars, ~run_ols(df_g, .x))
  if (nrow(res) == 0) return(NULL)
  res %>% mutate(group=grp)
})

if (nrow(reg_grid) > 0) {
  reg_grid <- reg_grid %>%
    mutate(
      sig = case_when(p_val<0.001~"***", p_val<0.01~"**",
                      p_val<0.05~"*", p_val<0.1~".", TRUE~"ns"),
      inv_label = inv_labels[inv_var]
    )
}

write_csv(reg_grid, file.path(OUT,"reg_grid.csv"))

cat("\n=== 回归结果（按|t值|排序，前15名）===\n")
reg_grid %>%
  filter(!is.na(p_val)) %>%
  arrange(desc(abs(t_val))) %>%
  head(15) %>%
  mutate(across(c(beta,se,r2_adj),~round(.x,4))) %>%
  dplyr::select(group, inv_var, n, beta, p_val, sig, r2_adj) %>%
  print()

# =============================================================================
# 7. 方差分解（varpart）——全部投资变量
# =============================================================================
run_vp <- function(df, inv) {
  d <- df %>%
    filter(!is.na(.data[[inv]]), if_all(all_of(ctrl_vars), ~!is.na(.x)),
           !is.na(longitude)) %>%
    filter(.data[[inv]] < quantile(.data[[inv]], 0.99, na.rm=TRUE),
           .data[[inv]] > 0)
  if (nrow(d) < 30) return(NULL)
  vp <- tryCatch(
    vegan::varpart(d$TRRI,
                   dplyr::select(d, all_of(inv)),
                   dplyr::select(d, vpd_bg, heat_freq_bg),
                   dplyr::select(d, latitude, longitude)),
    error=function(e) NULL)
  if (is.null(vp)) return(NULL)
  fr <- vp$part$indfract
  tibble(inv_var=inv, n=nrow(d),
         invest_R2  = fr$Adj.R.square[1],
         climate_R2 = fr$Adj.R.square[2],
         geo_R2     = fr$Adj.R.square[3])
}

vp_all <- map_dfr(c("ALL","B","C","D"), function(grp) {
  df_g <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  res  <- map_dfr(inv_vars, ~run_vp(df_g, .x))
  if (nrow(res) == 0) return(NULL)
  mutate(res, group = grp)
})

write_csv(vp_all, file.path(OUT,"varpart.csv"))

for (grp in c("ALL","B","C","D")) {
  sub <- filter(vp_all, group == grp)
  if (nrow(sub) == 0) { cat("\n=== 方差分解（", grp, "组）：样本不足，跳过 ===\n"); next }
  cat(sprintf("\n=== 方差分解（%s组，n站点≈%d）===\n", grp, sub$n[1]))
  sub %>%
    arrange(desc(invest_R2)) %>%
    mutate(across(c(invest_R2,climate_R2,geo_R2), ~sprintf("%.2f%%", .x*100))) %>%
    print()
}

# =============================================================================
# 8. 图 1：森林图（投资系数）
# =============================================================================
p_forest <- reg_grid %>%
  filter(!is.na(beta)) %>%
  mutate(ci_lo=beta-1.96*se, ci_hi=beta+1.96*se,
         sig_flag = p_val < 0.05,
         group = factor(group, levels=c("ALL","B","C","D"))) %>%
  ggplot(aes(x=beta, y=inv_label, color=sig_flag)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_errorbar(aes(xmin=ci_lo,xmax=ci_hi), width=0.25, linewidth=0.7) +
  geom_point(size=3) +
  scale_color_manual(values=c("FALSE"="grey60","TRUE"="#D73027"),
                     labels=c("ns","p<0.05"), name="") +
  facet_wrap(~group, ncol=4) +
  labs(title="投资强度变量对TRRI的回归系数（控制VPD、热频率、纬度）",
       subtitle="红色=显著(p<0.05)；误差线为95% CI；已排除SIF结构突变站点",
       x="回归系数 β", y=NULL) +
  theme_cn() +
  theme(legend.position="bottom",
        panel.grid.major.y=element_blank())

ggsave(file.path(OUT,"forest_invest.png"), p_forest, width=14, height=7, dpi=300)
cat("-> forest_invest.png\n")

# =============================================================================
# 图 2：方差分解气泡图
# =============================================================================
grp_labels <- c(ALL="全部站点", B="B类(干旱)", C="C类(温带)", D="D类(大陆)")
p_vp <- vp_all %>%
  pivot_longer(c(invest_R2,climate_R2,geo_R2),
               names_to="comp", values_to="r2") %>%
  mutate(
    r2 = pmax(r2*100, 0),
    comp_label = case_when(comp=="invest_R2"~"投资", comp=="climate_R2"~"气候",
                           comp=="geo_R2"~"地理"),
    inv_label  = inv_labels[inv_var],
    grp_label  = factor(grp_labels[group], levels=grp_labels)
  ) %>%
  filter(!is.na(grp_label)) %>%
  ggplot(aes(x=comp_label, y=inv_label, size=r2, color=comp)) +
  geom_point(alpha=0.8) +
  geom_text(aes(label=ifelse(r2>0.05, sprintf("%.1f%%",r2), "")),
            size=3.5, color="black", vjust=-1.5, family="heiti") +
  scale_size_continuous(range=c(1,14), name="独立贡献(%)") +
  scale_color_manual(values=c(invest_R2="#D73027",climate_R2="#4575B4",
                               geo_R2="#1A9850"), guide="none") +
  facet_wrap(~grp_label, nrow=1) +
  labs(title="方差分解：投资、气候、地理对TRRI的独立贡献（按气候区）",
       x="方差来源", y="投资变量") +
  theme_cn() +
  theme(panel.grid.major=element_line(color="grey90"))

ggsave(file.path(OUT,"varpart_bubble.png"), p_vp, width=12, height=6, dpi=300)
cat("-> varpart_bubble.png\n")

# =============================================================================
# 图 3：突变站点分布（地图）
# =============================================================================
break_map_df <- break_res %>%
  left_join(ccm %>% distinct(meteo_stat_id,longitude,latitude), by="meteo_stat_id") %>%
  filter(!is.na(longitude))

china_border <- rnaturalearth::ne_countries(country="China", scale="medium",
                                              returnclass="sf")
p_break_map <- ggplot() +
  geom_sf(data=china_border, fill="grey95", color="grey70") +
  geom_point(data=break_map_df, aes(x=longitude, y=latitude,
                                     color=has_break, size=has_break), alpha=0.7) +
  scale_color_manual(values=c("FALSE"="#2166AC","TRUE"="#D73027"),
                     labels=c("no break","structural break"), name="SIF break") +
  scale_size_manual(values=c("FALSE"=1.2,"TRUE"=2.5), guide="none") +
  coord_sf(xlim=c(73,135), ylim=c(18,54)) +
  labs(title="SIF结构突变站点空间分布",
       subtitle=sprintf("共%d站点有显著突变（p<0.05），已从TRRI分析中排除",n_break),
       x="经度", y="纬度") +
  theme_cn() +
  theme(legend.position="right", panel.grid=element_line(color="grey90"))

ggsave(file.path(OUT,"break_stations_map.png"), p_break_map,
       width=10, height=6, dpi=300)
cat("-> break_stations_map.png\n")

# =============================================================================
# 9. 最终汇总
# =============================================================================
best_reg <- reg_grid %>% filter(!is.na(p_val)) %>%
  arrange(desc(abs(t_val))) %>% head(1)
best_vp  <- vp_all %>% filter(group=="ALL") %>%
  arrange(desc(invest_R2)) %>% head(1)

cat("\n", strrep("=",65), "\n")
cat("最终结果摘要\n")
cat(strrep("=",65), "\n\n")
cat(sprintf("SIF突变站点: %d个（%.1f%%）已排除\n", n_break,
            n_break/length(stations_all)*100))
cat(sprintf("剩余分析站点: %d\n\n", nrow(trri_df)))

cat("【回归：投资效应最强的组合】\n")
cat(sprintf("  气候组: %s\n  投资变量: %s\n  β = %.4f (SE=%.4f, t=%.2f, p=%.4f%s)\n  调整R²= %.4f\n\n",
            best_reg$group, best_reg$inv_var, best_reg$beta,
            best_reg$se, best_reg$t_val, best_reg$p_val, best_reg$sig,
            best_reg$r2_adj))

cat("【方差分解：投资独立贡献最大的组合（全部站点）】\n")
cat(sprintf("  投资变量: %s\n  投资独立贡献: %.2f%%\n  气候独立贡献: %.2f%%\n  地理独立贡献: %.2f%%\n",
            best_vp$inv_var, best_vp$invest_R2*100,
            best_vp$climate_R2*100, best_vp$geo_R2*100))

cat("\n所有结果保存至:", OUT, "\n")

# =============================================================================
# 图 4：最优组合散点图（inv_slope vs TRRI，按气候组着色）
# =============================================================================
best_inv  <- best_reg$inv_var   # "inv_slope"
best_grp  <- best_reg$group     # "ALL"

plot_df <- anal_df %>%
  filter(!is.na(.data[[best_inv]]), !is.na(TRRI)) %>%
  filter(.data[[best_inv]] < quantile(.data[[best_inv]], 0.99, na.rm=TRUE),
         .data[[best_inv]] > 0)

p_scatter <- ggplot(plot_df,
                    aes(x = .data[[best_inv]], y = TRRI, color = koppen_group)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  scale_color_brewer(palette = "Set1", name = "Köppen") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title  = "最优组合：投资增长斜率 vs TRRI",
    subtitle = sprintf("β=%.4f, p=%.4f**, 调整R²=%.3f | 控制VPD、热频率、纬度",
                       best_reg$beta, best_reg$p_val, best_reg$r2_adj),
    x = "园林绿化投资增长斜率（万元/年，2005-2020）",
    y = "TRRI (1-12)"
  ) +
  theme_cn() +
  theme(legend.position = "right")

ggsave(file.path(OUT, "scatter_best_combo.png"),
       p_scatter, width = 10, height = 6, dpi = 300)
cat("-> scatter_best_combo.png\n")

# =============================================================================
# 图 5：TRRI 分布箱线图（按气候组）
# =============================================================================
p_trri_box <- anal_df %>%
  mutate(
    stype_label = case_when(
      stype == "always_inhibit"  ~ "全程抑制",
      stype == "always_promote"  ~ "全程促进",
      stype == "promote_inhibit" ~ "促进→抑制",
      stype == "inhibit_promote" ~ "抑制→促进"
    ),
    stype_label = factor(stype_label,
                         levels = c("全程促进","促进→抑制","抑制→促进","全程抑制"))
  ) %>%
  ggplot(aes(x = koppen_group, y = TRRI, fill = koppen_group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.8) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(title = "TRRI 分布（按Köppen气候大类）",
       subtitle = "已排除SIF结构突变站点（n=679）",
       x = "Köppen气候大类", y = "TRRI (1=全程抑制, 12=全程促进)") +
  theme_cn()

ggsave(file.path(OUT, "trri_by_climate.png"),
       p_trri_box, width = 8, height = 6, dpi = 300)
cat("-> trri_by_climate.png\n")

# =============================================================================
# 表格：完整汇总（reg + varpart 合并）导出 CSV
# =============================================================================
summary_table <- reg_grid %>%
  dplyr::select(group, inv_var, n, beta, se, t_val, p_val, sig, r2_adj) %>%
  left_join(
    vp_all %>% dplyr::select(group, inv_var, invest_R2, climate_R2, geo_R2),
    by = c("group", "inv_var")
  ) %>%
  mutate(
    inv_label    = inv_labels[inv_var],
    beta         = round(beta,        6),
    se           = round(se,          6),
    t_val        = round(t_val,       3),
    p_val        = round(p_val,       4),
    r2_adj       = round(r2_adj,      4),
    invest_R2_pct = round(invest_R2 * 100, 2),
    climate_R2_pct= round(climate_R2 * 100, 2),
    geo_R2_pct   = round(geo_R2 * 100, 2)
  ) %>%
  dplyr::select(group, inv_label, inv_var, n, beta, se, t_val, p_val, sig,
                r2_adj, invest_R2_pct, climate_R2_pct, geo_R2_pct) %>%
  arrange(group, desc(abs(t_val)))

write_csv(summary_table, file.path(OUT, "summary_table.csv"))
cat("-> summary_table.csv\n")

cat("\n全部", length(list.files(OUT)), "个文件已保存至:", OUT, "\n")
cat("文件列表:\n")
cat(paste(" -", list.files(OUT)), sep="\n")

