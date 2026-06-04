# =============================================================================
# 05_1_analysis.R
#
# CCM  ：VPD均值（vpd_mean_dt），读取已有结果，不重新跑
# 投资 ：pa_built_10y（单位建成区绿地面积投资，2011-2020均值）
# 回归 ：有序Logit（TRRI为因变量）
#         + 多项Logit（stype四类为因变量）
#         全国 + 分气候区；仅站点尺度
# 方差 ：varpart（投资 vs 地理+控制）
#
# 输出：
#   G1. 中国地图：Köppen气候区底色 + 站点CCM响应类型
#   G2. 热图：站点 × time lag（S-Map系数），气候区分面板
#   G3. 热图：城市 × stype站点占比
#   I.  有序Logit：全部变量系数表（全国+分气候区）
#   J.  多项Logit：全部变量系数表（全国+分气候区）
#   K.  方差分解：各分量独立+共享解释力（全国+分气候区）
# =============================================================================

pacman::p_load(
  dplyr, tidyr, purrr, ggplot2, stringr, readr,
  vegan, tibble, scales, showtext, sysfonts,
  targets, MASS, nnet, ggnewscale, patchwork,
  terra, sf, rnaturalearth, rnaturalearthdata
)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT     <- "data_proc/output_10y_built_up_05_01"
CCM_RDS <- "data_proc/output_10y_built_up_05_01/ccm_05_results.rds"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
BS <- 14
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

stype_levels <- c("always_inhibit", "inhibit_promote",
                  "promote_inhibit", "always_promote")
stype_labels <- c("全程抑制", "抑制→促进", "促进→抑制", "全程促进")
stype_colors <- c(
  always_inhibit  = "#4575B4",
  inhibit_promote = "#74ADD1",
  promote_inhibit = "#FDAE61",
  always_promote  = "#D73027"
)
koppen_colors <- c(A = "#2CA25F", B = "#D95F0E", C = "#3182BD", D = "#756BB1")
koppen_labels <- c(A = "A 热带", B = "B 干旱", C = "C 温带", D = "D 大陆")

winsorize <- function(x, p = 0.99) {
  q <- quantile(x, p, na.rm = TRUE); ifelse(x > q, q, x)
}
sig_star <- function(p) case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "ns"
)


# =============================================================================
# B. 读取 CCM 结果
# =============================================================================

cat("\n=== B. 读取CCM结果 ===\n")
ccm_results <- readRDS(CCM_RDS)
cat(sprintf("站点数: %d | tp: %d..%d | 行数: %d\n",
            length(unique(ccm_results$meteo_stat_id)),
            min(ccm_results$tp), max(ccm_results$tp), nrow(ccm_results)))


# =============================================================================
# C. TRRI / stype 分类
# =============================================================================

cat("\n=== C. TRRI / stype 分类 ===\n")

N_TP     <- length(unique(ccm_results$tp))
MAX_TP   <- max(ccm_results$tp)
TRRI_MAX <- 2L * N_TP
cat(sprintf("N_TP=%d, TRRI范围: 1..%d\n", N_TP, TRRI_MAX))

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

cat("站点类型分布:\n"); print(count(trri_df, stype))
write_csv(trri_df %>% dplyr::select(-coef_seq),
          file.path(OUT, "trri_station.csv"))


# =============================================================================
# D. 投资变量（pa_built_10y）
# =============================================================================

cat("\n=== D. 投资变量构建 ===\n")

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
year_cols <- paste0("inv_", 2002:2030)
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

cat(sprintf("pa_built_10y 非NA城市: %d\n", sum(!is.na(invest_tbl$pa_built_10y))))


# =============================================================================
# E. 控制变量（人均GDP；CGI政府治理指数2021）
# =============================================================================

cat("\n=== E. 控制变量 ===\n")

city_db_raw <- readxl::read_excel(
  "data_raw/china_city_db.xlsx",
  sheet = 1, col_names = TRUE
) %>%
  rename(year = 1, city_name = 3, pgdp = 15, urban_rate = 26) %>%
  dplyr::select(year, city_name, pgdp, urban_rate) %>%
  filter(year >= 2011, year <= 2020) %>%
  mutate(across(c(pgdp, urban_rate), as.numeric))

city_ctrl <- city_db_raw %>%
  group_by(city_name) %>%
  summarise(pgdp_10y       = mean(pgdp,       na.rm = TRUE),
            urban_rate_10y = mean(urban_rate, na.rm = TRUE),
            .groups = "drop")

# CGI 政府治理与公共服务指数（2021）
cgi_raw <- readxl::read_excel(
  "data_raw/green_invest/cgi_2021.xlsx",
  col_names = TRUE
) %>%
  rename(cgi_rank = 1, city_name = 2, province = 3, cgi_score = 4) %>%
  dplyr::select(city_name, cgi_score) %>%
  mutate(
    # 统一城市名格式：末尾加"市"（部分可能已有）
    city_name = ifelse(str_detect(city_name, "市$"), city_name,
                       paste0(city_name, "市"))
  )

cat(sprintf("CGI数据：%d 城市，得分范围 %.1f ~ %.1f\n",
            nrow(cgi_raw), min(cgi_raw$cgi_score), max(cgi_raw$cgi_score)))

city_vars <- invest_tbl %>%
  left_join(city_ctrl, by = "city_name") %>%
  left_join(cgi_raw,   by = "city_name")

invest_station <- station_city_map %>% left_join(city_vars, by = "city_name")


# =============================================================================
# F. 合并分析数据框（仅站点级）
# =============================================================================

cat("\n=== F. 合并分析数据框 ===\n")

anal_df <- trri_df %>%
  left_join(invest_station, by = "meteo_stat_id") %>%
  filter(!is.na(TRRI), !is.na(longitude), !is.na(latitude),
         !is.na(koppen_group)) %>%
  mutate(
    koppen_B         = as.integer(koppen_group == "B"),
    koppen_C         = as.integer(koppen_group == "C"),
    koppen_D         = as.integer(koppen_group == "D"),
    pa_built_10y_w   = winsorize(pa_built_10y),
    pgdp_10y_w       = winsorize(pgdp_10y),
    urban_rate_10y_w = winsorize(urban_rate_10y),
    cgi_score_w      = winsorize(cgi_score),
    TRRI_ord         = factor(TRRI, levels = 1:TRRI_MAX, ordered = TRUE),
    stype_fct        = factor(stype, levels = stype_levels, labels = stype_labels)
  )

cat(sprintf("站点总数: %d\n", nrow(anal_df)))
cat(sprintf("pa_built_10y 非NA: %d\n",    sum(!is.na(anal_df$pa_built_10y))))
cat(sprintf("pgdp_10y 非NA: %d\n",        sum(!is.na(anal_df$pgdp_10y))))
cat(sprintf("cgi_score 非NA: %d\n",       sum(!is.na(anal_df$cgi_score))))
cat("气候组分布:\n"); print(count(anal_df, koppen_group))

# 城市级（仅用于G3热图，不做回归）
city_agg <- anal_df %>%
  group_by(city_name) %>%
  summarise(
    n_stations        = n(),
    koppen_group      = first(koppen_group),
    n_always_inhibit  = sum(stype == "always_inhibit"),
    n_inhibit_promote = sum(stype == "inhibit_promote"),
    n_promote_inhibit = sum(stype == "promote_inhibit"),
    n_always_promote  = sum(stype == "always_promote"),
    .groups = "drop"
  )


# =============================================================================
# G. 可视化
# =============================================================================

cat("\n=== G. 可视化 ===\n")

# ---------- G1. 中国地图：Köppen底色（仅中国境内）+ 站点CCM类型 + 边际密度 ----------

koppen_tif  <- "data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif"
china_bbox  <- c(xmin = 72, xmax = 136, ymin = 17, ymax = 54)

koppen_to_group <- function(x) {
  case_when(
    x >= 1  & x <= 4  ~ "A",
    x >= 5  & x <= 9  ~ "B",
    x >= 10 & x <= 17 ~ "C",
    x >= 18 & x <= 28 ~ "D",
    TRUE ~ NA_character_
  )
}

# --- 底图矢量 ---
world_land  <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"),
                        error = function(e) NULL)
china_sf    <- tryCatch(
  ne_countries(country = c("China","Hong Kong S.A.R.","Macao S.A.R.",
                            "Taiwan"), scale = "medium", returnclass = "sf") %>%
    st_union(),
  error = function(e) NULL)
china_prov  <- tryCatch(ne_states(country = "China", returnclass = "sf"),
                        error = function(e) NULL)

# --- Köppen 栅格：仅保留中国境内像元 ---
koppen_rast <- tryCatch({
  r <- terra::rast(koppen_tif)
  terra::crop(r, terra::ext(china_bbox))
}, error = function(e) { cat("  [Köppen raster读取失败]\n"); NULL })

koppen_df_china <- NULL
if (!is.null(koppen_rast) && !is.null(china_sf)) {
  # 用中国边界 mask 栅格
  china_vect <- tryCatch(terra::vect(china_sf), error = function(e) NULL)
  if (!is.null(china_vect)) {
    kr_masked <- tryCatch(terra::mask(koppen_rast, china_vect),
                          error = function(e) koppen_rast)  # 退回无mask版
  } else {
    kr_masked <- koppen_rast
  }
  koppen_df_china <- terra::as.data.frame(kr_masked, xy = TRUE) %>%
    rename(koppen_code = 3) %>%
    mutate(koppen_grp = koppen_to_group(koppen_code)) %>%
    filter(!is.na(koppen_grp))
}

map_df <- trri_df %>%
  mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))

stype_colors_cn <- setNames(stype_colors, stype_labels)

# Köppen色系：绿/橙/蓝/紫（暖色系，区别于CCM的红蓝系）
koppen_colors_map <- c(A = "#1B9E77", B = "#D95F02", C = "#7570B3", D = "#E7298A")
koppen_labels_map <- c(A = "A 热带/亚热带", B = "B 干旱", C = "C 温带", D = "D 大陆")

# --- 主地图 ---
theme_map <- theme_minimal(base_size = 13) +
  theme(
    text             = element_text(family = "heiti"),
    panel.background = element_rect(fill = "#D6EAF8", color = NA),  # 海洋浅蓝
    panel.grid.major = element_line(color = "white", linewidth = 0.3),
    axis.text        = element_text(size = 10),
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11, face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                    family = "heiti"),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11,
                                    family = "heiti")
  )

p_map_main <- ggplot()

# 世界陆地（浅灰）
if (!is.null(world_land))
  p_map_main <- p_map_main +
    geom_sf(data = world_land, fill = "grey88", color = "grey75",
            linewidth = 0.15, inherit.aes = FALSE)

# 中国境内 Köppen 底色
if (!is.null(koppen_df_china))
  p_map_main <- p_map_main +
    geom_raster(data = koppen_df_china,
                aes(x = x, y = y, fill = koppen_grp), alpha = 0.65) +
    scale_fill_manual(values = koppen_colors_map, labels = koppen_labels_map,
                      name = "Köppen气候区", na.value = "grey88")

# 省界
if (!is.null(china_prov))
  p_map_main <- p_map_main +
    geom_sf(data = china_prov, fill = NA, color = "grey55",
            linewidth = 0.2, inherit.aes = FALSE)

# 站点（CCM类型：红橙黄绿蓝系与Köppen色系区分）
stype_colors_point <- c(
  "全程抑制"  = "#1A6FBF",   # 深蓝
  "抑制→促进" = "#74C6E8",   # 浅蓝
  "促进→抑制" = "#F4A261",   # 橙
  "全程促进"  = "#C1121F"    # 深红
)

p_map_main <- p_map_main +
  ggnewscale::new_scale_color() +
  geom_point(data = map_df,
             aes(x = longitude, y = latitude, color = stype_label),
             size = 1.8, alpha = 0.88, shape = 16) +
  scale_color_manual(values = stype_colors_point, name = "CCM响应类型") +
  coord_sf(xlim = c(72, 136), ylim = c(17, 54), expand = FALSE) +
  labs(title    = "城市植被热响应类型（CCM）与Köppen气候区",
       subtitle = sprintf("X=VPD均值（去趋势）；tp=0..%d；共%d站", MAX_TP, nrow(map_df)),
       x = "经度", y = "纬度") +
  theme_map

# --- 保存主地图 ---
ggsave(file.path(OUT, "map_stype_koppen.png"),
       p_map_main, width = 14, height = 9, dpi = 300)
cat("-> map_stype_koppen.png\n")

# --- 经纬度 × CCM类型占比图（独立输出）---
make_prop_df <- function(data, coord_col, bin_size, min_n = 3) {
  data %>%
    filter(!is.na(.data[[coord_col]]), !is.na(stype_label)) %>%
    mutate(bin = floor(.data[[coord_col]] / bin_size) * bin_size + bin_size / 2) %>%
    group_by(bin, stype_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(bin) %>%
    mutate(pct = n / sum(n), total = sum(n)) %>%
    ungroup() %>%
    filter(total >= min_n)
}

lon_prop <- make_prop_df(map_df, "longitude", 3)
lat_prop <- make_prop_df(map_df, "latitude",  2)

theme_prop <- theme_minimal(base_size = 13) +
  theme(text            = element_text(family = "heiti"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 11),
        axis.title       = element_text(size = 12),
        legend.text      = element_text(size = 11, family = "heiti"),
        legend.title     = element_text(size = 12, face = "bold", family = "heiti"),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                        family = "heiti"),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11,
                                        family = "heiti"))

# 图2：经度区间 × CCM占比（横轴=经度，纵轴=占比，堆叠条形）
n_lon_total <- n_distinct(map_df$meteo_stat_id[!is.na(map_df$longitude)])
p_lon_prop <- ggplot(lon_prop,
                     aes(x = factor(bin), y = pct, fill = stype_label)) +
  geom_col(position = "stack", alpha = 0.88, width = 0.85) +
  geom_text(data = lon_prop %>% group_by(bin) %>% slice(1),
            aes(x = factor(bin), y = 1.03, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3, color = "grey40", family = "heiti") +
  scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  scale_x_discrete(labels = function(x)
    paste0(as.numeric(x) - 1.5, "°–", as.numeric(x) + 1.5, "°E")) +
  labs(title    = "不同经度区间各类CCM响应类型的站点占比",
       subtitle = sprintf("每3°经度一组；共%d站；n=各组站点数", n_lon_total),
       x = "经度区间", y = "占比") +
  theme_prop +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(file.path(OUT, "ccm_prop_by_lon.png"),
       p_lon_prop, width = 13, height = 6, dpi = 300)
cat("-> ccm_prop_by_lon.png\n")

# 图3：纬度区间 × CCM占比（纵轴=纬度，横轴=占比，水平堆叠条形）
n_lat_total <- n_distinct(map_df$meteo_stat_id[!is.na(map_df$latitude)])
p_lat_prop <- ggplot(lat_prop,
                     aes(y = factor(bin), x = pct, fill = stype_label)) +
  geom_col(position = "stack", alpha = 0.88, width = 0.85) +
  geom_text(data = lat_prop %>% group_by(bin) %>% slice(1),
            aes(y = factor(bin), x = 1.03, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3, color = "grey40", family = "heiti",
            hjust = 0) +
  scale_fill_manual(values = stype_colors_point, name = "CCM响应类型") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.12), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x)
    paste0(as.numeric(x) - 1, "°–", as.numeric(x) + 1, "°N")) +
  labs(title    = "不同纬度区间各类CCM响应类型的站点占比",
       subtitle = sprintf("每2°纬度一组；共%d站；n=各组站点数", n_lat_total),
       y = "纬度区间", x = "占比（累加至100%）") +
  theme_prop +
  theme(legend.position = "right")

ggsave(file.path(OUT, "ccm_prop_by_lat.png"),
       p_lat_prop, width = 9, height = 10, dpi = 300)
cat("-> ccm_prop_by_lat.png\n")

# ---------- G2. 热图：站点 × time lag，按stype+气候区 ----------

coef_long <- ccm_results %>%
  dplyr::select(meteo_stat_id, tp, mean_coef) %>%
  left_join(trri_df %>% dplyr::select(meteo_stat_id, stype, koppen_group, TRRI),
            by = "meteo_stat_id") %>%
  filter(!is.na(stype)) %>%
  mutate(coef_clip = pmax(pmin(mean_coef, 0.5), -0.5))

MAX_PER_ZONE <- 120
koppen_groups_plot <- sort(unique(coef_long$koppen_group))

p_heatmap_list <- map(koppen_groups_plot, function(grp) {
  d <- coef_long %>% filter(koppen_group == grp)
  stations <- d %>% dplyr::select(meteo_stat_id, stype, TRRI) %>%
    distinct() %>% arrange(factor(stype, levels = stype_levels), TRRI)
  if (nrow(stations) > MAX_PER_ZONE) {
    stations <- stations %>%
      group_by(stype) %>%
      slice_sample(prop = MAX_PER_ZONE / nrow(stations)) %>%
      ungroup() %>% arrange(factor(stype, levels = stype_levels), TRRI)
  }
  stations <- stations %>% mutate(y_rank = row_number())
  d <- d %>%
    filter(meteo_stat_id %in% stations$meteo_stat_id) %>%
    left_join(dplyr::select(stations, meteo_stat_id, y_rank, stype),
              by = c("meteo_stat_id", "stype"))

  # stype分隔线
  sep_lines <- stations %>% group_by(stype) %>%
    summarise(y_max = max(y_rank), .groups = "drop") %>%
    filter(y_max < max(stations$y_rank))

  # stype标注位置
  stype_anno <- stations %>% group_by(stype) %>%
    summarise(y_mid = mean(y_rank), .groups = "drop") %>%
    mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))

  koppen_name <- c(A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")[grp]

  ggplot(d, aes(x = factor(tp), y = factor(y_rank), fill = coef_clip)) +
    geom_tile() +
    geom_hline(data = sep_lines, aes(yintercept = y_max + 0.5),
               color = "white", linewidth = 1.5, inherit.aes = FALSE) +
    annotate("text",
             x    = rep(0.3, nrow(stype_anno)),
             y    = stype_anno$y_mid,
             label = as.character(stype_anno$stype_label),
             hjust = 0, size = 3.2, family = "heiti", color = "grey20") +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish,
      name = "S-Map系数"
    ) +
    scale_x_discrete(labels = paste0("tp", 0:MAX_TP)) +
    labs(
      title    = sprintf("%s气候区：站点热响应热图（n=%d）",
                         koppen_name, length(unique(d$meteo_stat_id))),
      subtitle = "行=站点（全程抑制→…→全程促进）；列=时间滞后；蓝=抑制 红=促进",
      x = "Time lag", y = NULL
    ) +
    theme_cn(bs = 12) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank())
})

iwalk(p_heatmap_list, function(p, i) {
  grp   <- koppen_groups_plot[i]
  fname <- file.path(OUT, sprintf("heatmap_coef_%s.png", grp))
  ggsave(fname, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("-> heatmap_coef_%s.png\n", grp))
})

# ---------- G3. 热图：城市 × stype站点占比 ----------

city_stype_long <- city_agg %>%
  pivot_longer(
    cols      = c(n_always_inhibit, n_inhibit_promote, n_promote_inhibit, n_always_promote),
    names_to  = "stype",
    values_to = "count"
  ) %>%
  mutate(
    stype       = str_remove(stype, "^n_"),
    stype_label = factor(stype, levels = stype_levels, labels = stype_labels),
    pct         = count / n_stations * 100
  )

city_order <- city_agg %>%
  arrange(koppen_group, desc(n_stations)) %>%
  mutate(city_rank = row_number())

city_stype_long <- city_stype_long %>%
  left_join(dplyr::select(city_order, city_name, city_rank), by = "city_name")

koppen_sep <- city_order %>%
  group_by(koppen_group) %>% summarise(y_max = max(city_rank), .groups = "drop") %>%
  filter(y_max < max(city_order$city_rank))

p_city_stype <- ggplot(city_stype_long,
                       aes(x = stype_label, y = factor(city_rank), fill = pct)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_hline(data = koppen_sep, aes(yintercept = y_max + 0.5),
             color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_fill_gradient(low = "white", high = "#D73027",
                      name = "占比(%)", limits = c(0, 100)) +
  labs(
    title    = "城市各热响应类型站点占比",
    subtitle = sprintf("共%d城市，按气候区（A/B/C/D）和站点数降序排列",
                       nrow(city_agg)),
    x = "响应类型", y = "城市（↑占更多站点）"
  ) +
  theme_cn() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid  = element_blank())

ggsave(file.path(OUT, "heatmap_city_stype.png"),
       p_city_stype, width = 9, height = 14, dpi = 300)
cat("-> heatmap_city_stype.png\n")

write_csv(
  city_stype_long %>%
    dplyr::select(city_name, koppen_group, stype_label, count, pct) %>%
    arrange(koppen_group, city_name),
  file.path(OUT, "city_stype_counts.csv")
)
cat("-> city_stype_counts.csv\n")


# =============================================================================
# I. 有序Logit（TRRI ~ 所有变量）
#    仅站点级；全国 + 分气候区；单一控制变量集
#    全国：inv_w + cgi_score_w + lon + lat + Köppen哑变量
#    分区：inv_w + cgi_score_w + lon + lat
#    输出：所有回归变量的系数表
# =============================================================================

cat("\n=== I. 有序Logit（TRRI）===\n")

# 全国：地理+Köppen哑变量；分区内：仅地理（已按区分层，不需要Köppen哑变量）
ctrl_all   <- c("cgi_score_w", "longitude", "latitude",
                "koppen_B", "koppen_C", "koppen_D")
ctrl_inner <- c("cgi_score_w", "longitude", "latitude")

run_olr_full <- function(df, inv_var, ctrl_vars, grp_label) {
  d <- df %>%
    filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0,
           !is.na(TRRI),
           if_all(all_of(ctrl_vars[ctrl_vars %in% names(df)]), ~!is.na(.x)))
  if (nrow(d) < 15) {
    cat(sprintf("  [跳过 %s] n=%d\n", grp_label, nrow(d)))
    return(NULL)
  }
  d$inv_w  <- d[[inv_var]]   # 仅截尾，不标准化
  d$Y_ord  <- factor(d$TRRI, levels = sort(unique(d$TRRI)), ordered = TRUE)
  ctrl_use <- ctrl_vars[ctrl_vars %in% names(d)]
  ctrl_use <- ctrl_use[sapply(ctrl_use, function(v) var(d[[v]], na.rm = TRUE) > 0)]
  all_pred <- c("inv_w", ctrl_use)
  fml <- as.formula(paste("Y_ord ~", paste(all_pred, collapse = "+")))
  m <- tryCatch(
    MASS::polr(fml, data = d, Hess = TRUE, method = "logistic"),
    error = function(e) { cat(sprintf("  [polr error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s <- tryCatch(summary(m),
                error = function(e) { cat(sprintf("  [summary error: %s]\n", e$message)); NULL })
  if (is.null(s)) return(NULL)

  coef_mat <- s$coefficients[all_pred, , drop = FALSE]
  z_vec    <- coef_mat[, "t value"]
  pval_vec <- 2 * pnorm(abs(z_vec), lower.tail = FALSE)

  m0  <- tryCatch(MASS::polr(Y_ord ~ 1, data = d, Hess = FALSE),
                  error = function(e) NULL)
  mcf <- if (!is.null(m0))
    round(1 - as.numeric(logLik(m)) / as.numeric(logLik(m0)), 4)
  else NA_real_

  cat(sprintf("  [%s] n=%d  inv_w: β=%.4f p=%.4f%s  cgi: β=%.4f p=%.4f%s  McF=%.4f\n",
              grp_label, nrow(d),
              coef_mat["inv_w", "Value"], pval_vec["inv_w"],
              sig_star(pval_vec["inv_w"]),
              ifelse("cgi_score_w" %in% all_pred, coef_mat["cgi_score_w","Value"], NA),
              ifelse("cgi_score_w" %in% all_pred, pval_vec["cgi_score_w"], NA),
              ifelse("cgi_score_w" %in% all_pred, sig_star(pval_vec["cgi_score_w"]), ""),
              mcf))

  var_labels <- c(
    inv_w       = "投资强度（pa_built_10y）",
    cgi_score_w = "CGI政府治理指数（2021，截尾）",
    longitude   = "经度",
    latitude    = "纬度",
    koppen_B    = "Köppen B（干旱区）",
    koppen_C    = "Köppen C（温带）",
    koppen_D    = "Köppen D（大陆）"
  )

  tibble(
    koppen    = grp_label,
    n         = nrow(d),
    mcfadden  = mcf,
    variable  = rownames(coef_mat),
    var_label = var_labels[rownames(coef_mat)],
    coef      = round(coef_mat[, "Value"],      6),
    se        = round(coef_mat[, "Std. Error"], 6),
    z_val     = round(z_vec,    3),
    p_val     = round(pval_vec, 4),
    sig       = sig_star(pval_vec),
    direction = ifelse(coef_mat[, "Value"] > 0, "+", "-")
  )
}

koppen_all <- c("ALL", sort(unique(anal_df$koppen_group)))

olr_tbl <- map_dfr(koppen_all, function(grp) {
  df_g   <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
  run_olr_full(df_g, "pa_built_10y_w", ctrl_v, grp)
})

cat("\n=== I. 有序Logit系数汇总 ===\n")
print(olr_tbl %>%
        filter(variable %in% c("inv_w", "cgi_score_w")) %>%
        dplyr::select(koppen, variable, n, mcfadden, coef, se, p_val, sig, direction))

write_csv(olr_tbl, file.path(OUT, "olr_trri_full_coef.csv"))
cat("-> olr_trri_full_coef.csv（含所有变量系数）\n")

# 投资+CGI p值热图（气候区为y轴，变量为x轴）
if (nrow(olr_tbl) > 0) {
  koppen_nm <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  var_nm    <- c(inv_w = "投资（pa_built_10y）", cgi_score_w = "CGI治理指数")

  p_olr_pval <- olr_tbl %>%
    filter(variable %in% c("inv_w", "cgi_score_w")) %>%
    mutate(
      grp_label  = factor(koppen_nm[koppen],
                          levels = rev(koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)])),
      var_nm_lbl = factor(var_nm[variable], levels = var_nm),
      cell_label = sprintf("β=%.3f\np=%.3f%s\nn=%d", coef, p_val, sig, n),
      p_fill     = pmin(p_val, 0.5)
    ) %>%
    filter(!is.na(grp_label), !is.na(var_nm_lbl)) %>%
    ggplot(aes(x = var_nm_lbl, y = grp_label, fill = p_fill)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.8, family = "heiti") +
    scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#F0F0F0",
                         midpoint = 0.05, limits = c(0, 0.5),
                         name = "p值（红=显著）") +
    labs(title    = "有序Logit：管理变量对TRRI的效应",
         subtitle = "因变量=TRRI（有序因子，1-18）；管理=投资强度+CGI治理指数",
         x = NULL, y = NULL) +
    theme_cn() + theme(panel.grid = element_blank())

  ggsave(file.path(OUT, "olr_trri_invest_summary.png"),
         p_olr_pval, width = 8, height = 6, dpi = 300)
  cat("-> olr_trri_invest_summary.png\n")
}


# =============================================================================
# I2. 有序Logit 可视化
#     (a) 森林图：所有预测变量系数 ± 95%CI，分气候区
#     (b) 散点+箱线：投资/CGI 分位 × TRRI（原始数据关系）
#     (c) 预测概率曲线：固定其他变量，投资/CGI 沿范围变化时 P(TRRI=k) 变化
# =============================================================================

cat("\n=== I2. 有序Logit 可视化 ===\n")

if (nrow(olr_tbl) > 0) {
  koppen_nm  <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  koppen_col <- c("全部"="#555555","A(热带)"="#1B9E77","B(干旱)"="#D95F02",
                  "C(温带)"="#7570B3","D(大陆)"="#E7298A")

  # (a) 森林图：系数 ± 1.96*SE -------------------------------------------
  forest_df <- olr_tbl %>%
    mutate(
      grp_label = factor(koppen_nm[koppen],
                         levels = koppen_nm[c("ALL","B","C","D")]),
      ci_lo     = coef - 1.96 * se,
      ci_hi     = coef + 1.96 * se,
      var_short = case_when(
        variable == "inv_w"       ~ "单位面积投资",
        variable == "cgi_score_w" ~ "CGI治理指数",
        variable == "longitude"   ~ "经度",
        variable == "latitude"    ~ "纬度",
        variable == "koppen_B"    ~ "Köppen B",
        variable == "koppen_C"    ~ "Köppen C",
        variable == "koppen_D"    ~ "Köppen D",
        TRUE ~ variable
      ),
      var_short = factor(var_short,
                         levels = rev(c("单位面积投资","CGI治理指数","经度","纬度",
                                        "Köppen B","Köppen C","Köppen D")))
    ) %>%
    filter(!is.na(grp_label))

  p_forest <- ggplot(forest_df,
                     aes(x = coef, y = var_short, color = grp_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.6) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(0.6), linewidth = 0.8) +
    geom_point(aes(shape = sig %in% c("*","**","***",".")),
               position = position_dodge(0.6), size = 2.8) +
    scale_color_manual(values = koppen_col, name = "气候区") +
    scale_shape_manual(values = c("TRUE"=16, "FALSE"=1),
                       labels = c("TRUE"="显著(p<0.1)", "FALSE"="不显著"),
                       name = "") +
    labs(title    = "有序Logit各预测变量系数（森林图）",
         subtitle = "水平线=95%CI；实心圆=p<0.1显著；纵轴=预测变量；横轴=log-odds系数",
         x = "系数（Log-odds）", y = NULL) +
    theme_cn() +
    theme(legend.position = "right",
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  ggsave(file.path(OUT, "olr_forest.png"),
         p_forest, width = 10, height = 6, dpi = 300)
  cat("-> olr_forest.png\n")

  # (b) 散点+箱线：将连续投资/CGI 分成5分位组，展示各组 TRRI 分布 --------
  make_quantile_box <- function(df, xvar, xlabel) {
    df %>%
      filter(!is.na(.data[[xvar]]), .data[[xvar]] > 0,
             !is.na(TRRI), !is.na(koppen_group)) %>%
      mutate(
        q5     = cut(.data[[xvar]],
                     breaks = quantile(.data[[xvar]], probs = 0:5/5, na.rm = TRUE),
                     include.lowest = TRUE, labels = paste0("Q", 1:5)),
        kgroup = factor(koppen_nm[koppen_group],
                        levels = koppen_nm[c("B","C","D")])
      ) %>%
      filter(!is.na(q5), !is.na(kgroup)) %>%
      ggplot(aes(x = q5, y = TRRI, fill = kgroup)) +
      geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5,
                   linewidth = 0.5, alpha = 0.75) +
      geom_smooth(aes(x = as.integer(q5), group = kgroup, color = kgroup),
                  method = "loess", se = FALSE, linewidth = 1.2,
                  show.legend = FALSE) +
      scale_fill_manual(values  = c("B(干旱)"="#D95F02","C(温带)"="#7570B3",
                                    "D(大陆)"="#E7298A"), name = "气候区") +
      scale_color_manual(values = c("B(干旱)"="#D95F02","C(温带)"="#7570B3",
                                    "D(大陆)"="#E7298A")) +
      scale_y_continuous(breaks = seq(1, 18, 2)) +
      labs(title    = sprintf("%s 五分位组 × TRRI 分布", xlabel),
           subtitle = "箱线=各分位组TRRI分布；折线=loess趋势",
           x = sprintf("%s（Q1=最低，Q5=最高）", xlabel),
           y = "TRRI（1=全程抑制，18=全程促进）") +
      theme_cn() +
      theme(legend.position = "right")
  }

  p_box_inv <- make_quantile_box(anal_df, "pa_built_10y_w", "单位面积投资")
  p_box_cgi <- make_quantile_box(anal_df, "cgi_score_w",    "CGI治理指数")

  ggsave(file.path(OUT, "olr_box_invest.png"),
         p_box_inv, width = 10, height = 5, dpi = 300)
  cat("-> olr_box_invest.png\n")

  ggsave(file.path(OUT, "olr_box_cgi.png"),
         p_box_cgi, width = 10, height = 5, dpi = 300)
  cat("-> olr_box_cgi.png\n")

  # (c) 预测概率曲线：固定其他变量于中位数，沿自变量范围预测 TRRI --------
  # 用全国模型（ALL）预测；展示 predicted mean TRRI（加权期望值）
  make_pred_curve <- function(xvar, xlabel, n_pts = 50) {
    d_all <- anal_df %>%
      filter(!is.na(.data[[xvar]]), .data[[xvar]] > 0,
             !is.na(TRRI), !is.na(cgi_score_w),
             !is.na(longitude), !is.na(latitude),
             !is.na(koppen_B))
    if (nrow(d_all) < 30) return(NULL)

    # 重建全国模型（与 run_olr_full 一致，但用于预测）
    d_all$inv_w <- d_all[[xvar]]
    d_all$Y_ord <- factor(d_all$TRRI, levels = sort(unique(d_all$TRRI)), ordered = TRUE)
    ctrl_use <- ctrl_all[ctrl_all %in% names(d_all)]
    ctrl_use <- ctrl_use[sapply(ctrl_use, function(v) var(d_all[[v]], na.rm=TRUE) > 0)]
    # 投资变量名替换
    pred_vars <- c("inv_w", ctrl_use[ctrl_use != xvar])
    fml <- as.formula(paste("Y_ord ~", paste(pred_vars, collapse = "+")))
    m <- tryCatch(
      MASS::polr(fml, data = d_all, Hess = FALSE, method = "logistic"),
      error = function(e) NULL)
    if (is.null(m)) return(NULL)

    # 预测数据框：xvar 沿范围变化，其余固定于中位数
    xseq   <- seq(min(d_all[[xvar]], na.rm=TRUE),
                  max(d_all[[xvar]], na.rm=TRUE), length.out = n_pts)
    newdat <- d_all[rep(1, n_pts), pred_vars, drop = FALSE]
    for (v in pred_vars) {
      if (v != "inv_w") newdat[[v]] <- median(d_all[[v]], na.rm = TRUE)
    }
    newdat$inv_w <- xseq
    # 强制转换 Y_ord levels 与训练集一致
    levels_ord <- levels(d_all$Y_ord)
    probs <- tryCatch(predict(m, newdat, type = "probs"), error = function(e) NULL)
    if (is.null(probs)) return(NULL)
    if (!is.matrix(probs)) probs <- matrix(probs, nrow = 1)

    # 计算预测 E[TRRI]
    trri_vals <- as.integer(levels_ord)
    pred_mean <- as.numeric(probs %*% trri_vals)

    # 95% CI via delta method（简化：用 ±1 TRRI 单位作为示意）
    tibble(x = xseq, pred_trri = pred_mean,
           xlabel = xlabel, xvar = xvar)
  }

  pred_inv <- make_pred_curve("pa_built_10y_w", "单位面积投资（万元/公顷）")
  pred_cgi <- make_pred_curve("cgi_score_w",    "CGI治理指数")

  pred_all <- bind_rows(pred_inv, pred_cgi) %>% filter(!is.na(pred_trri))

  if (nrow(pred_all) > 0) {
    # 同时展示两条预测曲线（分面）
    p_pred <- ggplot(pred_all, aes(x = x, y = pred_trri)) +
      geom_line(color = "#D73027", linewidth = 1.2) +
      geom_rug(data = anal_df %>%
                 pivot_longer(c(pa_built_10y_w, cgi_score_w),
                              names_to = "xvar", values_to = "x") %>%
                 filter(!is.na(x), x > 0) %>%
                 mutate(xlabel = ifelse(xvar == "pa_built_10y_w",
                                        "单位面积投资（万元/公顷）", "CGI治理指数")),
               aes(x = x), sides = "b", alpha = 0.15, color = "grey40",
               inherit.aes = FALSE) +
      facet_wrap(~xlabel, scales = "free_x") +
      labs(title    = "有序Logit：预测 TRRI 随管理变量的变化（全国模型）",
           subtitle = "其余变量固定于中位数；竖线密度=实际数据分布",
           x = "自变量取值", y = "预测 TRRI（期望值）") +
      theme_cn()

    ggsave(file.path(OUT, "olr_pred_curve.png"),
           p_pred, width = 11, height = 5, dpi = 300)
    cat("-> olr_pred_curve.png\n")
  }
}


# =============================================================================
# J. 多项Logit：stype（四类名义变量）~ pa_built_10y + cgi + 控制
#    仅站点级；全国 + 分气候区；输出所有变量系数
# =============================================================================

cat("\n=== J. 多项Logit（stype四类）===\n")

run_mnl_full <- function(df, inv_var, ctrl_vars, grp_label,
                         ref = "always_inhibit") {
  d <- df %>%
    filter(!is.na(.data[[inv_var]]), .data[[inv_var]] > 0, !is.na(stype),
           if_all(all_of(ctrl_vars[ctrl_vars %in% names(df)]), ~!is.na(.x)))
  if (nrow(d) < 20) {
    cat(sprintf("  [跳过 %s] n=%d\n", grp_label, nrow(d)))
    return(NULL)
  }
  d$inv_w   <- d[[inv_var]]  # 仅截尾，不标准化
  d$stype_f <- relevel(factor(d$stype, levels = stype_levels), ref = ref)
  ctrl_use  <- ctrl_vars[ctrl_vars %in% names(d)]
  ctrl_use  <- ctrl_use[sapply(ctrl_use, function(v) var(d[[v]], na.rm = TRUE) > 0)]
  all_pred  <- c("inv_w", ctrl_use)
  fml <- as.formula(paste("stype_f ~", paste(all_pred, collapse = "+")))
  m <- tryCatch(
    suppressWarnings(nnet::multinom(fml, data = d, trace = FALSE, MaxNWts = 10000)),
    error = function(e) { cat(sprintf("  [multinom error: %s]\n", e$message)); NULL }
  )
  if (is.null(m)) return(NULL)
  s      <- summary(m)
  coef_m <- s$coefficients
  se_m   <- s$standard.errors
  z_m    <- coef_m / se_m
  pval_m <- 2 * pnorm(abs(z_m), lower.tail = FALSE)

  cat(sprintf("  MNL [%s] n=%d  inv_w p: %s\n",
              grp_label, nrow(d),
              paste(sprintf("%s:%.3f%s", rownames(coef_m),
                            pval_m[, "inv_w"], sig_star(pval_m[, "inv_w"])),
                    collapse = " ")))

  var_labels <- c(
    inv_w       = "投资强度（pa_built_10y）",
    cgi_score_w = "CGI政府治理指数（2021，截尾）",
    longitude   = "经度",
    latitude    = "纬度",
    koppen_B    = "Köppen B（干旱区）",
    koppen_C    = "Köppen C（温带）",
    koppen_D    = "Köppen D（大陆）"
  )

  map_dfr(rownames(coef_m), function(cat) {
    vars <- colnames(coef_m)
    tibble(
      koppen      = grp_label,
      reference   = ref,
      outcome_cat = cat,
      n           = nrow(d),
      variable    = vars,
      var_label   = var_labels[vars],
      coef        = round(coef_m[cat, vars], 6),
      se          = round(se_m[cat,   vars], 6),
      z_val       = round(z_m[cat,    vars], 3),
      p_val       = round(pval_m[cat,  vars], 4),
      sig         = sig_star(pval_m[cat, vars]),
      direction   = ifelse(coef_m[cat, vars] > 0,
                           "↑更可能落入此类", "↓更不可能落入此类")
    )
  })
}

mnl_tbl <- map_dfr(koppen_all, function(grp) {
  df_g   <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  ctrl_v <- if (grp == "ALL") ctrl_all else ctrl_inner
  run_mnl_full(df_g, "pa_built_10y_w", ctrl_v, grp)
})

cat("\n=== J. 多项Logit：管理变量系数汇总 ===\n")
print(mnl_tbl %>%
        filter(variable %in% c("inv_w", "cgi_score_w")) %>%
        dplyr::select(koppen, outcome_cat, variable, n, coef, se, p_val, sig))

write_csv(mnl_tbl, file.path(OUT, "mnl_stype_full_coef.csv"))
cat("-> mnl_stype_full_coef.csv（含所有变量系数）\n")

# 投资+CGI系数热图（气候区 × 响应类型，两个变量并列）
if (nrow(mnl_tbl) > 0) {
  koppen_nm  <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")
  stype_lmap <- setNames(stype_labels, stype_levels)
  var_nm     <- c(inv_w = "投资", cgi_score_w = "CGI治理")

  p_mnl <- mnl_tbl %>%
    filter(variable %in% c("inv_w", "cgi_score_w")) %>%
    mutate(
      grp_label  = factor(koppen_nm[koppen],
                          levels = koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)]),
      cat_label  = factor(stype_lmap[outcome_cat], levels = stype_labels),
      var_lbl    = factor(var_nm[variable], levels = var_nm),
      cell_label = sprintf("β=%.2f\np=%.3f%s", coef, p_val, sig),
      coef_clip  = pmax(pmin(coef, 1), -1)
    ) %>%
    filter(!is.na(grp_label), !is.na(cat_label)) %>%
    ggplot(aes(x = cat_label, y = grp_label, fill = coef_clip)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.2, family = "heiti") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                         midpoint = 0, name = "β（截断±1）") +
    facet_wrap(~var_lbl, nrow = 1) +
    labs(title    = "多项Logit：管理变量对响应类型归属的效应",
         subtitle = "参照类别=全程抑制；β>0表示该变量增加→更可能落入该类",
         x = "对比类别（vs 全程抑制）", y = NULL) +
    theme_cn() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          panel.grid  = element_blank())

  ggsave(file.path(OUT, "mnl_stype_invest_beta.png"),
         p_mnl, width = 14, height = 6, dpi = 300)
  cat("-> mnl_stype_invest_beta.png\n")
}


# =============================================================================
# K. 方差分解（varpart）
#    分组：[管理] = pa_built_10y_w + cgi_score_w
#          [地理] = lon/lat/Köppen/pgdp（基础版不含pgdp）
#    输出四分量：管理独立[a]、地理独立[b]、共享[c]、未解释[d]
#    仅站点级；全国 + 分气候区；基础控制 + 扩展控制（含pgdp）
# =============================================================================

cat("\n=== K. 方差分解 ===\n")

# 管理变量组（需要两者均非NA且>0）
mgmt_vars <- c("pa_built_10y_w", "cgi_score_w")

run_vp_full <- function(df, outcome_var, mgmt_set, geo_set, grp_label, ctrl_label) {
  # 过滤：管理变量均有值，且投资 > 0
  d <- df %>%
    filter(!is.na(pa_built_10y), pa_built_10y > 0,
           !is.na(.data[[outcome_var]])) %>%
    filter(if_all(all_of(mgmt_set[mgmt_set %in% names(.)]), ~!is.na(.x))) %>%
    filter(if_all(all_of(geo_set[geo_set %in% names(.)]),  ~!is.na(.x)))
  if (nrow(d) < 15) {
    cat(sprintf("  [跳过 %s|%s] n=%d\n", grp_label, ctrl_label, nrow(d)))
    return(NULL)
  }
  mgmt_use <- mgmt_set[mgmt_set %in% names(d)]
  mgmt_use <- mgmt_use[sapply(mgmt_use, function(v) var(d[[v]], na.rm=TRUE) > 0)]
  geo_use  <- geo_set[geo_set %in% names(d)]
  geo_use  <- geo_use[sapply(geo_use,  function(v) var(d[[v]], na.rm=TRUE) > 0)]
  if (length(mgmt_use) == 0 || length(geo_use) == 0) return(NULL)

  vp <- tryCatch(
    vegan::varpart(d[[outcome_var]],
                   dplyr::select(d, all_of(mgmt_use)),
                   dplyr::select(d, all_of(geo_use))),
    error = function(e) { cat(sprintf("  [varpart error: %s]\n", e$message)); NULL })
  if (is.null(vp)) return(NULL)

  fr <- vp$part$indfract
  cat(sprintf(
    "  [%s|%s] n=%d  管理=%.2f%%  地理=%.2f%%  共享=%.2f%%  未解释=%.2f%%\n",
    grp_label, ctrl_label, nrow(d),
    fr$Adj.R.square[1]*100, fr$Adj.R.square[2]*100,
    fr$Adj.R.square[3]*100, fr$Adj.R.square[4]*100))

  tibble(
    koppen        = grp_label,
    ctrl_type     = ctrl_label,
    n             = nrow(d),
    mgmt_vars_used = paste(mgmt_use, collapse = "+"),
    geo_vars_used  = paste(geo_use,  collapse = "+"),
    mgmt_only     = round(fr$Adj.R.square[1]*100, 2),  # [a] 管理独立
    geo_only      = round(fr$Adj.R.square[2]*100, 2),  # [b] 地理独立
    shared        = round(fr$Adj.R.square[3]*100, 2),  # [c] 共享
    unexplained   = round(fr$Adj.R.square[4]*100, 2),  # [d] 未解释
    mgmt_total    = round((fr$Adj.R.square[1]+fr$Adj.R.square[3])*100, 2),
    geo_total     = round((fr$Adj.R.square[2]+fr$Adj.R.square[3])*100, 2)
  )
}

# 地理变量组：全国含Köppen哑变量；分区内仅经纬度
geo_vp_all   <- c("longitude", "latitude", "koppen_B", "koppen_C", "koppen_D")
geo_vp_inner <- c("longitude", "latitude")

vp_tbl <- map_dfr(koppen_all, function(grp) {
  df_g  <- if (grp == "ALL") anal_df else filter(anal_df, koppen_group == grp)
  geo_v <- if (grp == "ALL") geo_vp_all else geo_vp_inner
  run_vp_full(df_g, "TRRI", mgmt_vars, geo_v, grp, "")
})

cat("\n=== K. 方差分解汇总 ===\n")
print(vp_tbl %>%
        dplyr::select(koppen, n, mgmt_only, geo_only, shared, unexplained))

write_csv(vp_tbl, file.path(OUT, "varpart_results.csv"))
cat("-> varpart_results.csv\n")

koppen_nm <- c(ALL="全部", A="A(热带)", B="B(干旱)", C="C(温带)", D="D(大陆)")

if (nrow(vp_tbl) > 0) {
  vp_plot_df <- vp_tbl %>%
    mutate(grp_label = factor(koppen_nm[koppen],
                              levels = koppen_nm[intersect(c("ALL","A","B","C","D"), koppen)])) %>%
    filter(!is.na(grp_label))

  # 图1：各分量堆叠条形图
  p_vp_stack <- vp_plot_df %>%
    pivot_longer(c(mgmt_only, shared, geo_only, unexplained),
                 names_to = "component", values_to = "r2") %>%
    mutate(
      r2_show    = pmax(r2, 0),
      comp_label = factor(component,
                          levels = c("mgmt_only","shared","geo_only","unexplained"),
                          labels = c("[a] 管理（独立）","[c] 共享",
                                     "[b] 地理（独立）","[d] 未解释"))
    ) %>%
    ggplot(aes(x = grp_label, y = r2_show, fill = comp_label)) +
    geom_col(position = "stack", alpha = 0.88, width = 0.65) +
    geom_text(aes(label = ifelse(r2_show >= 0.5, sprintf("%.1f%%", r2_show), "")),
              position = position_stack(vjust = 0.5),
              size = 3.5, family = "heiti", color = "white") +
    scale_fill_manual(
      values = c("[a] 管理（独立）"="#D73027", "[c] 共享"="#FDAE61",
                 "[b] 地理（独立）"="#4575B4", "[d] 未解释"="#CCCCCC"),
      name = "方差分量") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title    = "方差分解：管理（投资+CGI）vs 地理",
         subtitle = "[a]=管理独立；[b]=地理独立；[c]=共享；[d]=未解释\n（调整R²，负值显示为0；管理=pa_built_10y_w + cgi_score_w）",
         x = NULL, y = "调整R²（%）") +
    theme_cn() +
    theme(legend.position = "right")

  ggsave(file.path(OUT, "varpart_stack.png"),
         p_vp_stack, width = 11, height = 5, dpi = 300)
  cat("-> varpart_stack.png\n")

  # 图2：管理独立解释力[a]各气候区
  p_vp_mgmt <- vp_plot_df %>%
    mutate(
      r2_show  = pmax(mgmt_only, 0),
      r2_label = sprintf("%.2f%%", mgmt_only)
    ) %>%
    ggplot(aes(x = grp_label, y = r2_show, fill = grp_label)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = r2_label), vjust = -0.4, size = 4, family = "heiti") +
    scale_fill_manual(
      values = c("全部"="#888888", "A(热带)"="#2CA25F", "B(干旱)"="#D95F0E",
                 "C(温带)"="#3182BD", "D(大陆)"="#756BB1"),
      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title    = "方差分解：管理独立解释力（[a]分量）",
         subtitle = "管理 = pa_built_10y_w + cgi_score_w；负值代表独立贡献低于随机水平（显示为0）",
         x = NULL, y = "管理独立调整R²（%）") +
    theme_cn()

  ggsave(file.path(OUT, "varpart_mgmt_only.png"),
         p_vp_mgmt, width = 9, height = 5, dpi = 300)
  cat("-> varpart_mgmt_only.png\n")
}


# =============================================================================
# 完成
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("05_1_analysis.R 运行完成\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("输出目录  : %s\n", OUT))
cat(sprintf("CCM站点数 : %d\n", length(unique(ccm_results$meteo_stat_id))))
cat(sprintf("分析站点数: %d\n", nrow(anal_df)))
cat(sprintf("分析城市数: %d（仅用于G3热图）\n", nrow(city_agg)))
cat("\n输出文件:\n")
for (f in sort(list.files(OUT, pattern = "\\.(csv|png)$")))
  cat(sprintf("  %s\n", f))
