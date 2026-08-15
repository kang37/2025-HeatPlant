#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 16_hcsif_figures.R
#   基于全量 HCSIF VPD->SIF(buf750) CCM 结果，出两张图：
#   (1) 站点×tp 符号热力图：行=站点，列=滞后tp，颜色=促进/抑制；
#       站点自下往上排列：全促进 → 促进转抑制 → 抑制转促进 → 全抑制
#   (2) 中国地图：Köppen 气候区底色 + 4 种 VPD-SIF 模式站点分布
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(stringr)
})
# 全英文标签，不依赖中文字体，避免非交互 Rscript 下 showtext 缺字形显示为 "..."

PROJ    <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
CCM_RDS <- file.path(PROJ, "data_proc/ccm_hcsif_500m/ccm_hcsif_vpd_20260726_2018.rds")
OUT     <- file.path(PROJ, "data_proc/output_hcsif_new")
SIF_Y   <- "SIF_buf750_dt"; VPD_X <- "vpd_mean_dt"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# 1. 读取 + 模式分类（与 15 脚本一致）
# ===========================================================================
ccm <- readRDS(CCM_RDS) %>% as.data.frame() %>%
  filter(y_var == SIF_Y, x_var == VPD_X)
N_TP <- length(unique(ccm$tp))

# 严格 5 类：按符号变号次数分类，中间掺杂即归为 other
classify5 <- function(coefs) {
  s <- sign(coefs)
  if (any(is.na(s)) || any(s == 0)) return("other")
  nch <- sum(diff(s) != 0)                 # 符号变化次数
  if (nch == 0) return(if (s[1] > 0) "always_promote" else "always_inhibit")
  if (nch == 1) return(if (s[1] < 0) "inhibit_promote" else "promote_inhibit")
  "other"
}
# 单次干净翻转的位置（变号后第一个 tp 的索引, 1-based）；非单次翻转则 NA
trans_pos <- function(coefs) {
  s <- sign(coefs)
  if (any(is.na(s)) || any(s == 0)) return(NA_integer_)
  ch <- which(diff(s) != 0)
  if (length(ch) != 1) return(NA_integer_)
  as.integer(ch + 1L)
}

pat <- ccm %>%
  arrange(meteo_stat, tp) %>%
  group_by(meteo_stat) %>%
  summarise(coef_seq = list(mean_coef),
            longitude = first(longitude), latitude = first(latitude),
            .groups = "drop") %>%
  filter(map_lgl(coef_seq, ~length(.x) == N_TP)) %>%
  mutate(stype = map_chr(coef_seq, classify5),
         tpos  = map_int(coef_seq, trans_pos))

cat("严格5类分布:\n"); print(table(pat$stype))

# ===========================================================================
# 2. 热力图：站点×tp 符号
# ===========================================================================
# 站点排序：自下往上 全促进→促进转抑制→抑制转促进→全抑制
# 自下往上：全促进 → 先促进再抑制 → 先抑制再促进 → 全抑制 → 其他(顶部单独一带)
grp_rank <- c(always_promote = 1L, promote_inhibit = 2L,
              inhibit_promote = 3L, always_inhibit = 4L, other = 5L)
order_tbl <- pat %>%
  mutate(
    n_neg = map_int(coef_seq, ~sum(.x < 0, na.rm = TRUE)),
    g = grp_rank[stype],
    # 组内键 w：小=更靠下(更红/促进)，大=更靠上(更蓝/抑制)
    w = case_when(
      stype == "always_promote"  ~ 0,
      stype == "promote_inhibit" ~ -as.numeric(tpos),    # 翻转越晚(红得越久)越靠下
      stype == "inhibit_promote" ~  as.numeric(tpos),    # 翻转越晚(越晚转促进)越靠上(接全抑制)
      stype == "always_inhibit"  ~ 0,
      stype == "other"           ~ as.numeric(n_neg))    # other内按蓝格多少排
  ) %>%
  arrange(g, w) %>%              # g=1(全促进)最前 → station_order 最小 → y 最底
  mutate(station_order = row_number())

# 长表：station × tp × 符号
heat_long <- ccm %>%
  select(meteo_stat, tp, mean_coef) %>%
  inner_join(order_tbl %>% select(meteo_stat, station_order, stype),
             by = "meteo_stat") %>%
  mutate(sign = case_when(mean_coef > 0 ~ "Promote",
                          mean_coef < 0 ~ "Inhibit",
                          TRUE ~ NA_character_))

# 组分界线（用于在图上画白线分隔4类）
grp_bounds <- order_tbl %>% group_by(stype) %>%
  summarise(ymax = max(station_order), .groups = "drop") %>%
  arrange(ymax) %>% pull(ymax)
grp_bounds <- grp_bounds[-length(grp_bounds)] + 0.5

grp_label_pos <- order_tbl %>% group_by(stype) %>%
  summarise(y = mean(station_order), n = n(), .groups = "drop") %>%
  mutate(lab = c(always_promote = "Always promote", promote_inhibit = "Promote->Inhibit",
                 inhibit_promote = "Inhibit->Promote", always_inhibit = "Always inhibit",
                 other = "Other (mixed)")[stype])

p_heat <- ggplot(heat_long, aes(x = factor(tp), y = station_order, fill = sign)) +
  geom_tile() +
  geom_hline(yintercept = grp_bounds, color = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("Promote" = "#C1121F", "Inhibit" = "#1A6FBF"),
                    na.value = "grey85", name = "VPD->SIF effect") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = grp_label_pos$y,
                     labels = sprintf("%s\n(n=%d)", grp_label_pos$lab, grp_label_pos$n)) +
  labs(title = "Sign of VPD->SIF causal effect across CCM lag, per station",
       subtitle = sprintf("SIF = HCSIF buf750 (750m); each row = 1 station (n=%d), each column = 1 lag step (8 days);\nbottom->top: always promote -> promote-then-inhibit -> inhibit-then-promote -> always inhibit; within block sorted red->blue (transition timing)",
                          nrow(order_tbl)),
       x = "CCM lag tp (each step = 8 days)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 10, lineheight = 0.9),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT, "heatmap_station_tp_sign.png"), p_heat,
       width = 9, height = 11, dpi = 300)
cat("-> heatmap_station_tp_sign.png\n")

# ===========================================================================
# 3. 地图：Köppen 气候区 + 4 模式站点分布
# ===========================================================================
ok <- requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("terra", quietly = TRUE) &&
      requireNamespace("sf", quietly = TRUE) &&
      requireNamespace("ggnewscale", quietly = TRUE)

if (!ok) {
  cat("[跳过地图] 缺少 rnaturalearth/terra/sf/ggnewscale 之一\n")
} else {
  suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
  koppen_tif <- file.path(PROJ, "data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif")
  china_bbox <- c(xmin = 72, xmax = 136, ymin = 17, ymax = 54)
  koppen_to_group <- function(x) case_when(
    x >= 1 & x <= 4 ~ "A", x >= 5 & x <= 9 ~ "B",
    x >= 10 & x <= 17 ~ "C", x >= 18 & x <= 28 ~ "D", TRUE ~ NA_character_)

  world_land <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"), error = function(e) NULL)
  china_sf   <- tryCatch(ne_countries(country = c("China","Hong Kong S.A.R.","Macao S.A.R.","Taiwan"),
                                      scale = "medium", returnclass = "sf") %>% st_union(),
                         error = function(e) NULL)
  china_prov <- tryCatch(ne_states(country = "China", returnclass = "sf"), error = function(e) NULL)

  koppen_df <- NULL
  koppen_rast <- tryCatch(terra::crop(terra::rast(koppen_tif), terra::ext(china_bbox)),
                          error = function(e) NULL)
  if (!is.null(koppen_rast) && !is.null(china_sf)) {
    cv <- tryCatch(terra::vect(china_sf), error = function(e) NULL)
    kr <- if (!is.null(cv)) tryCatch(terra::mask(koppen_rast, cv), error = function(e) koppen_rast) else koppen_rast
    koppen_df <- terra::as.data.frame(kr, xy = TRUE) %>%
      rename(koppen_code = 3) %>% mutate(koppen_grp = koppen_to_group(koppen_code)) %>%
      filter(!is.na(koppen_grp))
  }

  stype_levels <- c("always_inhibit","inhibit_promote","promote_inhibit","always_promote","other")
  stype_labels <- c("Always inhibit","Inhibit->Promote","Promote->Inhibit","Always promote","Other (mixed)")
  map_df <- pat %>% mutate(stype_label = factor(stype, levels = stype_levels, labels = stype_labels))

  koppen_colors_map <- c(A = "#5A8F76", B = "#EAD5A0", C = "#A3B86C", D = "#C2DFCD")
  koppen_labels_map <- c(A = "A Tropical", B = "B Arid", C = "C Temperate", D = "D Continental")
  stype_colors_point <- c("Always inhibit" = "#1A6FBF","Inhibit->Promote" = "#74C6E8",
                          "Promote->Inhibit" = "#F4A261","Always promote" = "#C1121F",
                          "Other (mixed)" = "#9E9E9E")

  theme_map <- theme_minimal(base_size = 13) +
    theme(panel.background = element_rect(fill = "#D6EAF8", color = NA),
          panel.grid.major = element_line(color = "white", linewidth = 0.3),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  p_map <- ggplot()
  if (!is.null(world_land))
    p_map <- p_map + geom_sf(data = world_land, fill = "grey88", color = "grey75",
                             linewidth = 0.15, inherit.aes = FALSE)
  if (!is.null(koppen_df))
    p_map <- p_map + geom_raster(data = koppen_df, aes(x = x, y = y, fill = koppen_grp), alpha = 0.65) +
      scale_fill_manual(values = koppen_colors_map, labels = koppen_labels_map,
                        name = "Köppen气候区", na.value = "grey88")
  if (!is.null(china_prov))
    p_map <- p_map + geom_sf(data = china_prov, fill = NA, color = "grey55",
                             linewidth = 0.2, inherit.aes = FALSE)
  p_map <- p_map + ggnewscale::new_scale_color() +
    geom_point(data = map_df, aes(x = longitude, y = latitude, color = stype_label),
               size = 1.8, alpha = 0.88, shape = 16) +
    scale_color_manual(values = stype_colors_point, name = "VPD-SIF pattern") +
    coord_sf(xlim = c(72, 136), ylim = c(17, 54), expand = FALSE) +
    labs(title = "VPD-SIF causal response patterns and Koppen climate zones",
         subtitle = sprintf("SIF = HCSIF buf750; tp = 0..%d; n = %d stations", max(ccm$tp), nrow(map_df)),
         x = "Longitude", y = "Latitude") + theme_map

  ggsave(file.path(OUT, "map_stype_koppen_hcsif.png"), p_map, width = 14, height = 9, dpi = 300)
  cat("-> map_stype_koppen_hcsif.png\n")
}

cat("完成\n")
