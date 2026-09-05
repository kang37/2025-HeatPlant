#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 09_fig_classify_134.R — 134站"方向 x tp速度分箱"分类的两张图
#
#   fig1_bubble: x=optimal tp, y=S-map系数(median_coef_2v), 颜色=方向,
#                点大小=nsig_0to8(信号广度)。一张图装下三个分类因子。
#   fig2_map:    Köppen气候区底图(仿05_fig_maps_causal_criteria.R) + 站点,
#                颜色=方向, facet=tp速度分箱三档。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(sf); library(terra)
  library(rnaturalearth); library(ggnewscale); library(showtext); library(sysfonts)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
font_add("cjk", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)

OUT_DIR <- file.path(PROJ, "data_proc/smap_bivar_134")
st <- readRDS(file.path(OUT_DIR, "classify_134.rds"))

dir_colors <- c(Promote = "#C2703A", Inhibit = "#2E7D74", Ambiguous = "grey65")

# ===========================================================================
# 1. 气泡图: optimal tp x S-map系数, 颜色=方向, 大小=nsig_0to8
# ===========================================================================
set.seed(1)
st[, tp_jit := optimal_tp + runif(.N, -0.15, 0.15)]

p_bubble <- ggplot(st, aes(x = tp_jit, y = median_coef_2v)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_point(aes(color = direction_2v, size = nsig_0to8), alpha = 0.75) +
  scale_color_manual(values = dir_colors, name = "方向(二变量S-map)") +
  scale_size_continuous(range = c(1.5, 6), name = "tp=0..8内\n显著个数") +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "134 站分类概览: optimal tp x S-map 系数 x 显著广度",
       subtitle = "颜色=方向(主导符号占比>=75%判定, 否则Ambiguous); 点大小=显著tp个数(信号广度)",
       x = "optimal tp (每步8天; 已加抖动便于区分重叠点)",
       y = expression(paste("S-map 系数中位数  ", partialdiff, "SIF/", partialdiff, "VPD"))) +
  theme_minimal(base_size = 12, base_family = "cjk") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig1_bubble_tp_coef.png"), p_bubble, width = 9, height = 6, dpi = 300)
cat("已保存 fig1_bubble_tp_coef.png\n")

# ===========================================================================
# 2. 地图: Köppen底图 + 站点(颜色=方向) facet=tp速度分箱
# ===========================================================================
coords <- fread("data_raw/hcsif/stations_924.csv")
st <- merge(st, coords, by = "meteo_stat")

koppen_tif <- "data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif"
china_bbox <- c(xmin = 72, xmax = 136, ymin = 17, ymax = 54)
koppen_to_group <- function(x) {
  fifelse(x >= 1  & x <= 4,  "A",
  fifelse(x >= 5  & x <= 9,  "B",
  fifelse(x >= 10 & x <= 17, "C",
  fifelse(x >= 18 & x <= 28, "D", NA_character_))))
}

world_land <- tryCatch(ne_countries(scale = "medium", returnclass = "sf"), error = function(e) NULL)
china_sf   <- tryCatch(
  st_union(ne_countries(country = c("China","Hong Kong S.A.R.","Macao S.A.R.","Taiwan"),
                        scale = "medium", returnclass = "sf")),
  error = function(e) NULL)
china_prov <- tryCatch(ne_states(country = "China", returnclass = "sf"), error = function(e) NULL)

koppen_rast <- tryCatch({
  rr <- terra::rast(koppen_tif); terra::crop(rr, terra::ext(china_bbox))
}, error = function(e) NULL)

koppen_df_china <- NULL
if (!is.null(koppen_rast) && !is.null(china_sf)) {
  china_vect <- tryCatch(terra::vect(china_sf), error = function(e) NULL)
  kr_masked <- if (!is.null(china_vect)) tryCatch(terra::mask(koppen_rast, china_vect),
                                                  error = function(e) koppen_rast) else koppen_rast
  koppen_df_china <- as.data.table(terra::as.data.frame(kr_masked, xy = TRUE))
  setnames(koppen_df_china, 3, "koppen_code")
  koppen_df_china[, koppen_grp := koppen_to_group(koppen_code)]
  koppen_df_china <- koppen_df_china[!is.na(koppen_grp)]
}
koppen_colors_map <- c(A = "#5A8F76", B = "#EAD5A0", C = "#A3B86C", D = "#C2DFCD")
koppen_labels_map <- c(A = "A 热带/亚热带", B = "B 干旱", C = "C 温带", D = "D 大陆")

theme_map <- theme_minimal(base_size = 12, base_family = "cjk") +
  theme(panel.background = element_rect(fill = "#D6EAF8", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.3),
        legend.position  = "right",
        legend.key.size  = unit(0.4, "cm"),
        strip.text       = element_text(face = "bold", size = 11),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10))

base_map <- function() {
  p <- ggplot()
  if (!is.null(world_land))
    p <- p + geom_sf(data = world_land, fill = "grey88", color = "grey75",
                      linewidth = 0.15, inherit.aes = FALSE)
  if (!is.null(koppen_df_china))
    p <- p + geom_raster(data = koppen_df_china, aes(x = x, y = y, fill = koppen_grp), alpha = 0.55) +
      scale_fill_manual(values = koppen_colors_map, labels = koppen_labels_map,
                        name = "Köppen气候区", na.value = "grey88")
  if (!is.null(china_prov))
    p <- p + geom_sf(data = china_prov, fill = NA, color = "grey55", linewidth = 0.2, inherit.aes = FALSE)
  p + coord_sf(xlim = c(72, 136), ylim = c(17, 54), expand = FALSE) +
    labs(x = "经度", y = "纬度") + theme_map
}

p_map <- base_map() +
  ggnewscale::new_scale_color() +
  geom_point(data = st, aes(x = longitude, y = latitude, color = direction_2v,
                             size = nsig_0to8), alpha = 0.85, shape = 16) +
  scale_color_manual(values = dir_colors, name = "方向") +
  scale_size_continuous(range = c(1.2, 4), name = "tp=0..8内\n显著个数") +
  facet_wrap(~tp_bin) +
  labs(title = "134 站分类地图: 方向 x tp 响应速度分箱",
       subtitle = "二变量S-map方向(2026-09-05主线); 面板=optimal tp所在速度分箱")

ggsave(file.path(OUT_DIR, "fig2_map_direction_tpbin.png"), p_map, width = 13, height = 6, dpi = 300)
cat("已保存 fig2_map_direction_tpbin.png\n")
