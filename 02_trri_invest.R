# =============================================================================
# 02_trri_invest.R  (重写版)
# 输入: analysis_vpd2/ccm_vpd_results.rds  (VPD均值 & VPD超量累积)
# 流程: TRRI分类 → 多年投资变量 → OLS回归 × 4组 → varpart × 4组 → 图表
# =============================================================================

pacman::p_load(dplyr, tidyr, purrr, ggplot2, stringr, readr,
               vegan, tibble, scales, showtext, sysfonts, targets)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "analysis_vpd2/output_vpd"
dir.create(OUT, showWarnings = FALSE)

# 字体 -------------------------------------------------------------------------
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 16
theme_cn <- function(bs = BS) {
  theme_minimal(base_size = bs) +
    theme(text         = element_text(family = "heiti"),
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
# 1. 读取 CCM 结果（VPD 变量版）
# =============================================================================
ccm <- readRDS("analysis_vpd2/ccm_vpd_results.rds")
cat("CCM行数:", nrow(ccm), "| 站点数:", n_distinct(ccm$meteo_stat_id),
    "| X变量:", paste(unique(ccm$x_var), collapse=", "), "\n")

# =============================================================================
# 2. TRRI 分类（4类 → 12级）
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

make_trri <- function(df_xvar) {
  df_xvar %>%
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
}

# 对两个 X 变量分别分类
trri_list <- ccm %>%
  group_by(x_var) %>%
  group_split(x_var) %>%
  setNames(sort(unique(ccm$x_var))) %>%
  map(make_trri)

cat("\n站点分类分布:\n")
imap(trri_list, ~{ cat("[", .y, "]\n"); print(count(.x, stype)) })

# =============================================================================
# 3. 投资变量（与 03_run_full.R 相同方法：行位置对齐）
# =============================================================================
tar_load(green_invest_2020)   # 来自 targets，UTF-8，有市后缀
tar_load(gdp_data_2020)
tar_load(green_area_2020)
tar_load(station_city_map)

# 读取多年数据并按行对齐
invest_raw0 <- read.csv("data_raw/green_invest/city_invest_data.csv",
                        check.names=FALSE)
nc_inv <- ncol(invest_raw0)
names(invest_raw0) <- c("city_raw", paste0("inv_", 2002:(2001+nc_inv-1)))

invest_raw <- bind_cols(
  dplyr::select(green_invest_2020, city_name),
  dplyr::select(invest_raw0, -city_raw)
)

year_cols <- paste0("inv_", 2002:2024)
year_cols <- year_cols[year_cols %in% names(invest_raw)]
years_num <- as.integer(str_extract(year_cols, "\\d{4}$"))

get_mean_inv <- function(y1, y2) {
  cols <- year_cols[years_num >= y1 & years_num <= y2]
  rowMeans(mutate(dplyr::select(invest_raw, all_of(cols)),
                  across(everything(), as.numeric)), na.rm=TRUE)
}

slope_mat <- mutate(dplyr::select(invest_raw,
                                   all_of(year_cols[years_num %in% 2005:2020])),
                    across(everything(), as.numeric))
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
  mutate(city_name = ifelse(str_detect(city_name,"市$"), city_name,
                            paste0(city_name,"市")))
green <- green_area_2020 %>%
  dplyr::select(city_name, area_built=area_green_built, area_park=area_green_park) %>%
  mutate(city_name = ifelse(str_detect(city_name,"市$"), city_name,
                            paste0(city_name,"市")))

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
# 4. 气候背景变量
# =============================================================================
tar_load(data_heat_sif_weekly)
climate_bg <- data_heat_sif_weekly %>%
  filter(week %in% 20:39) %>%
  group_by(meteo_stat_id) %>%
  summarise(vpd_bg       = mean(vpd_mean,        na.rm=TRUE),
            heat_freq_bg = mean(heat_event_freq,  na.rm=TRUE),
            .groups="drop")

# 合并：每个 X 变量对应一个分析数据框
anal_list <- map(trri_list, ~{
  .x %>%
    left_join(invest_full, by="meteo_stat_id") %>%
    left_join(climate_bg,  by="meteo_stat_id") %>%
    filter(!is.na(TRRI), !is.na(koppen_group))
})
cat("\n各 X 变量分析数据框行数:\n")
map_int(anal_list, nrow) %>% print()

# =============================================================================
# 5. OLS 回归网格：7 投资变量 × 2 X变量 × 4 气候组
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
    filter(!is.na(.data[[inv]]), if_all(all_of(ctrl_vars), ~!is.na(.x)),
           .data[[inv]] > 0,
           .data[[inv]] < quantile(.data[[inv]], 0.99, na.rm=TRUE))
  if (nrow(d) < 20) return(NULL)
  m <- tryCatch(lm(as.formula(paste("TRRI ~", inv, "+",
                                     paste(ctrl_vars, collapse="+"))), data=d),
                error=function(e) NULL)
  if (is.null(m)) return(NULL)
  s  <- summary(m)
  cr <- s$coefficients[inv,,drop=FALSE]
  tibble(inv_var=inv, n=nrow(d), beta=cr[1,1], se=cr[1,2],
         t_val=cr[1,3], p_val=cr[1,4], r2_adj=s$adj.r.squared)
}

groups <- c("ALL","B","C","D")
reg_grid <- imap_dfr(anal_list, function(df, xv) {
  map_dfr(groups, function(grp) {
    df_g <- if (grp=="ALL") df else filter(df, koppen_group==grp)
    res  <- map_dfr(inv_vars, ~run_ols(df_g, .x))
    if (nrow(res)==0) return(NULL)
    mutate(res, group=grp, x_var=xv)
  })
}) %>%
  mutate(sig = case_when(p_val<0.001~"***", p_val<0.01~"**",
                         p_val<0.05~"*",    p_val<0.1~".",  TRUE~"ns"),
         inv_label = inv_labels[inv_var])

write_csv(reg_grid, file.path(OUT,"reg_grid_vpd.csv"))

cat("\n=== 回归结果（|t值|前15，两个X变量合并）===\n")
reg_grid %>% arrange(desc(abs(t_val))) %>% head(15) %>%
  mutate(across(c(beta,r2_adj),~round(.x,4)), p_val=round(p_val,4)) %>%
  dplyr::select(x_var, group, inv_var, n, beta, p_val, sig, r2_adj) %>%
  print()

# =============================================================================
# 6. 方差分解：7 投资变量 × 2 X变量 × 4 气候组
# =============================================================================
run_vp <- function(df, inv) {
  d <- df %>%
    filter(!is.na(.data[[inv]]), !is.na(longitude),
           if_all(all_of(ctrl_vars), ~!is.na(.x)),
           .data[[inv]] > 0,
           .data[[inv]] < quantile(.data[[inv]], 0.99, na.rm=TRUE))
  if (nrow(d) < 20) return(NULL)
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

vp_results <- imap_dfr(anal_list, function(df, xv) {
  map_dfr(groups, function(grp) {
    df_g <- if (grp=="ALL") df else filter(df, koppen_group==grp)
    res  <- map_dfr(inv_vars, ~run_vp(df_g, .x))
    if (nrow(res)==0) return(NULL)
    mutate(res, group=grp, x_var=xv)
  })
})

write_csv(vp_results, file.path(OUT,"varpart_vpd.csv"))

cat("\n=== 方差分解（投资独立贡献，按大小排序，全部站点）===\n")
vp_results %>% filter(group=="ALL") %>%
  arrange(x_var, desc(invest_R2)) %>%
  mutate(across(c(invest_R2,climate_R2,geo_R2),~sprintf("%.2f%%",.x*100))) %>%
  dplyr::select(x_var, inv_var, n, invest_R2, climate_R2, geo_R2) %>%
  print(n=20)

# =============================================================================
# 7. 图1：森林图（投资系数 × 气候组 × X变量分面）
# =============================================================================
p_forest <- reg_grid %>%
  filter(!is.na(beta)) %>%
  mutate(
    ci_lo     = beta - 1.96*se,
    ci_hi     = beta + 1.96*se,
    sig_flag  = p_val < 0.05,
    group     = factor(group, levels=c("ALL","B","C","D")),
    x_label   = ifelse(x_var=="vpd_mean_dt","X=VPD均值","X=VPD超量累积")
  ) %>%
  ggplot(aes(x=beta, y=inv_label, color=sig_flag)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_errorbar(aes(xmin=ci_lo, xmax=ci_hi), width=0.25, linewidth=0.7) +
  geom_point(size=3) +
  scale_color_manual(values=c("FALSE"="grey60","TRUE"="#D73027"),
                     labels=c("ns","p<0.05"), name="") +
  facet_grid(x_label ~ group) +
  labs(title="VPD变量版：投资强度对TRRI的回归系数",
       subtitle="控制VPD背景、热频率、纬度；红色=显著(p<0.05)",
       x="回归系数 β", y=NULL) +
  theme_cn() +
  theme(legend.position="bottom", panel.grid.major.y=element_blank())

ggsave(file.path(OUT,"forest_vpd.png"), p_forest, width=16, height=8, dpi=300)
cat("-> forest_vpd.png\n")

# =============================================================================
# 图2：方差分解气泡图（2 X变量 × 4 气候组）
# =============================================================================
grp_labels <- c(ALL="全部",B="B(干旱)",C="C(温带)",D="D(大陆)")
xvar_labels <- c(vpd_mean_dt="X=VPD均值", heat_over_dt="X=VPD超量")

p_vp <- vp_results %>%
  pivot_longer(c(invest_R2,climate_R2,geo_R2),
               names_to="comp", values_to="r2") %>%
  mutate(
    r2         = pmax(r2*100, 0),
    comp_label = case_when(comp=="invest_R2"~"投资",
                           comp=="climate_R2"~"气候",
                           comp=="geo_R2"~"地理"),
    inv_label  = inv_labels[inv_var],
    grp_label  = factor(grp_labels[group], levels=grp_labels),
    x_label    = xvar_labels[x_var]
  ) %>%
  filter(!is.na(grp_label)) %>%
  ggplot(aes(x=comp_label, y=inv_label, size=r2, color=comp)) +
  geom_point(alpha=0.8) +
  geom_text(aes(label=ifelse(r2>0.1, sprintf("%.1f%%",r2), "")),
            size=3.2, color="black", vjust=-1.4, family="heiti") +
  scale_size_continuous(range=c(1,12), name="独立贡献(%)") +
  scale_color_manual(values=c(invest_R2="#D73027",climate_R2="#4575B4",
                               geo_R2="#1A9850"), guide="none") +
  facet_grid(x_label ~ grp_label) +
  labs(title="方差分解：投资/气候/地理对TRRI独立贡献（VPD变量版）",
       x="方差来源", y="投资变量") +
  theme_cn() +
  theme(panel.grid.major=element_line(color="grey92"))

ggsave(file.path(OUT,"varpart_vpd.png"), p_vp, width=16, height=8, dpi=300)
cat("-> varpart_vpd.png\n")

# =============================================================================
# 图3：TRRI箱线图（两个X变量对比）
# =============================================================================
p_box <- imap_dfr(trri_list, ~mutate(.x, x_var=.y)) %>%
  mutate(x_label = xvar_labels[x_var]) %>%
  ggplot(aes(x=koppen_group, y=TRRI, fill=koppen_group)) +
  geom_boxplot(alpha=0.7, outlier.size=1) +
  geom_jitter(width=0.2, alpha=0.2, size=0.8) +
  scale_fill_brewer(palette="Set1", guide="none") +
  facet_wrap(~x_label) +
  labs(title="TRRI分布（VPD变量版，按气候区）",
       subtitle=sprintf("共%d个站点（各X变量）", n_distinct(ccm$meteo_stat_id)),
       x="Köppen气候大类", y="TRRI (1-12)") +
  theme_cn()

ggsave(file.path(OUT,"trri_box_vpd.png"), p_box, width=12, height=6, dpi=300)
cat("-> trri_box_vpd.png\n")

# =============================================================================
# 汇总表
# =============================================================================
summary_vpd <- reg_grid %>%
  dplyr::select(x_var, group, inv_var, inv_label, n, beta, se,
                t_val, p_val, sig, r2_adj) %>%
  left_join(vp_results %>%
              dplyr::select(x_var, group, inv_var,
                            invest_R2_pct=invest_R2,
                            climate_R2_pct=climate_R2,
                            geo_R2_pct=geo_R2),
            by=c("x_var","group","inv_var")) %>%
  mutate(across(c(invest_R2_pct,climate_R2_pct,geo_R2_pct),
                ~round(.x*100, 2))) %>%
  arrange(x_var, group, desc(abs(t_val)))

write_csv(summary_vpd, file.path(OUT,"summary_vpd.csv"))
cat("-> summary_vpd.csv\n")

# 最终汇总
cat("\n", strrep("=",65), "\n")
cat("最终结果摘要（VPD变量版）\n")
cat(strrep("=",65), "\n\n")

best <- reg_grid %>% filter(!is.na(p_val)) %>% arrange(desc(abs(t_val))) %>% head(1)
cat(sprintf("【回归最强组合】\n  X变量=%s | 投资=%s | 气候组=%s\n  β=%.4f (p=%.4f%s) | R²=%.3f\n\n",
            best$x_var, best$inv_var, best$group,
            best$beta, best$p_val, best$sig, best$r2_adj))

best_vp <- vp_results %>% filter(group=="ALL") %>%
  arrange(desc(invest_R2)) %>% head(1)
cat(sprintf("【方差分解最强（全部站点）】\n  X变量=%s | 投资=%s\n  投资=%.2f%% 气候=%.2f%% 地理=%.2f%%\n",
            best_vp$x_var, best_vp$inv_var,
            best_vp$invest_R2*100, best_vp$climate_R2*100, best_vp$geo_R2*100))

cat("\n输出目录:", OUT, "\n")
cat("文件列表:\n")
cat(paste(" -", list.files(OUT)), sep="\n")
