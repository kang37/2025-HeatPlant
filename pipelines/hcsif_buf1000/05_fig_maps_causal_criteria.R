#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 05_fig_maps_causal_criteria.R — "确认 VPD→SIF 方向"两步判据的站点地图
#
# 数据源：data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv
# (904站 x tp=-8..8 全量 fwd CCM，见 docs/02_ccm.md 2026-09-04 条目)。
# 判据：p_surr<0.05 且 Δρ=rho-rho_min>0(单tp显著)；第二步 optimal tp(显著tp里
# rho最高者)>=0。底图仿 05_1_analysis.R 的 G1 节(Köppen栅格 + rnaturalearth
# 省界，terra+sf+ggnewscale)。
#
# 4 张图：
#   第一步(≥1个tp显著)：map1_first_sig_tp   点色=显著tp里最小的那个(-8..8范围)
#                        map2_nsig_0to8     点色=在tp=0..8子范围内有几个tp显著
#                                           (0/1/2/3/>3，不满足第一步的点灰色)
#   第二步(optimal tp>=0)：map3_optimal_tp  点色=optimal tp(仅通过第二步的站)
#                          map4_optimal_rho 点色=optimal tp处的rho(仅通过第二步的站)
#                                           (未通过第一步或第二步的点灰色)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(sf); library(terra)
  library(rnaturalearth); library(ggnewscale); library(showtext); library(sysfonts)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
font_add("cjk", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)

PROJ <- "."
OUT  <- "data_proc/ccm_hcsif_buf1000_fwd_negtp_full"
IN_CSV <- file.path(OUT, "fwd_negtp_full_20260904_1126.csv")

# ===========================================================================
# 1. 读数据、算两步判据
# ===========================================================================
r <- fread(IN_CSV)
r[, sig := !is.na(p_surr) & p_surr < 0.05 & (rho - rho_min) > 0]

coords <- fread("data_raw/hcsif/stations_924.csv")
r <- merge(r, coords, by = "meteo_stat")

# --- 第一步汇总：至少1个tp(-8..8全范围)显著 ---
step1 <- r[sig == TRUE, .(
  first_sig_tp = min(tp),                                  # map1: 显著里最小的tp
  n_sig_0to8   = sum(tp >= 0 & tp <= 8)                     # map2 的计数子范围: 0..8
), by = .(meteo_stat, latitude, longitude)]
step1[, n_sig_0to8_cat := fifelse(n_sig_0to8 == 0, "0(仅负tp显著)",
                            fifelse(n_sig_0to8 == 1, "1",
                            fifelse(n_sig_0to8 == 2, "2",
                            fifelse(n_sig_0to8 == 3, "3", ">3"))))]
step1[, n_sig_0to8_cat := factor(n_sig_0to8_cat,
        levels = c("0(仅负tp显著)", "1", "2", "3", ">3"))]

# --- 第二步汇总：在第一步基础上, optimal tp(argmax rho)>=0 ---
opt <- r[sig == TRUE, .SD[which.max(rho)], by = .(meteo_stat, latitude, longitude)][
  , .(meteo_stat, latitude, longitude, opt_tp = tp, opt_rho = rho)]
step2 <- opt[opt_tp >= 0]

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
log_msg("第一步(全量-8..8至少1个tp显著): ", nrow(step1), " 站")
log_msg("第二步(+optimal tp>=0): ", nrow(step2), " 站")

# 全部站点坐标(含未通过判据的，用来画灰点)
all_stations <- unique(r[, .(meteo_stat, latitude, longitude)])

# ===========================================================================
# 2. Köppen 底图(仿 05_1_analysis.R G1 节)
# ===========================================================================
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
  rr <- terra::rast(koppen_tif)
  terra::crop(rr, terra::ext(china_bbox))
}, error = function(e) { log_msg("Köppen栅格读取失败"); NULL })

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
        legend.key.size  = unit(0.45, "cm"),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10))

base_map <- function() {
  p <- ggplot()
  if (!is.null(world_land))
    p <- p + geom_sf(data = world_land, fill = "grey88", color = "grey75",
                      linewidth = 0.15, inherit.aes = FALSE)
  if (!is.null(koppen_df_china))
    p <- p + geom_raster(data = koppen_df_china, aes(x = x, y = y, fill = koppen_grp), alpha = 0.6) +
      scale_fill_manual(values = koppen_colors_map, labels = koppen_labels_map,
                        name = "Köppen气候区", na.value = "grey88")
  if (!is.null(china_prov))
    p <- p + geom_sf(data = china_prov, fill = NA, color = "grey55", linewidth = 0.2, inherit.aes = FALSE)
  p + coord_sf(xlim = c(72, 136), ylim = c(17, 54), expand = FALSE) +
    labs(x = "经度", y = "纬度") + theme_map
}

save_map <- function(p, file, extra_h = 0) {
  ggsave(file.path(OUT, file), p, width = 10, height = 7 + extra_h, dpi = 300)
  log_msg("已保存: ", file)
}

# ===========================================================================
# 3. 第一步：map1 首个显著tp、map2 tp=0..8内显著个数
# ===========================================================================
grey_pts <- function(p, data_grey) {
  p + ggnewscale::new_scale_color() +
    geom_point(data = data_grey, aes(x = longitude, y = latitude),
               color = "grey60", size = 1.1, alpha = 0.5, shape = 16)
}

not_step1 <- all_stations[!step1, on = "meteo_stat"]

p1 <- grey_pts(base_map(), not_step1) +
  ggnewscale::new_scale_color() +
  geom_point(data = step1, aes(x = longitude, y = latitude, color = first_sig_tp),
             size = 1.8, alpha = 0.9, shape = 16) +
  scale_color_gradient2(low = "#2E7D74", mid = "grey85", high = "#C2703A", midpoint = 0,
                        name = "首个显著 tp") +
  labs(title = "第一步：首个显著 tp 的取值",
       subtitle = sprintf("p_surr<0.05 且 Δρ>0，tp=-8..8；%d/%d 站通过(灰点=未通过)",
                          nrow(step1), nrow(all_stations)))
save_map(p1, "map1_first_sig_tp.png")

p2 <- grey_pts(base_map(), not_step1) +
  ggnewscale::new_scale_color() +
  geom_point(data = step1, aes(x = longitude, y = latitude, color = n_sig_0to8_cat),
             size = 1.8, alpha = 0.9, shape = 16) +
  scale_color_manual(values = c("0(仅负tp显著)" = "#B0B0D8", "1" = "#FDD49E",
                                "2" = "#FC8D59", "3" = "#D7301F", ">3" = "#7F0000"),
                     name = "tp=0..8内\n显著个数") +
  labs(title = "第一步：tp=0..8 范围内显著 tp 的个数",
       subtitle = sprintf("同一批第一步通过站(%d站)，仅按 0..8 子范围计数；灰点=未通过第一步",
                          nrow(step1)))
save_map(p2, "map2_nsig_0to8.png")

# ===========================================================================
# 4. 第二步：map3 optimal tp、map4 optimal tp 处的 rho
# ===========================================================================
not_step2 <- all_stations[!step2, on = "meteo_stat"]

p3 <- grey_pts(base_map(), not_step2) +
  ggnewscale::new_scale_color() +
  geom_point(data = step2, aes(x = longitude, y = latitude, color = opt_tp),
             size = 1.8, alpha = 0.9, shape = 16) +
  scale_color_gradient(low = "#FFF7BC", high = "#8C2D04", name = "optimal tp") +
  labs(title = "第二步：optimal tp 的取值(仅通过站)",
       subtitle = sprintf("optimal tp = 显著tp中rho最高者，且要求>=0；%d/%d 站通过(灰点=未通过)",
                          nrow(step2), nrow(all_stations)))
save_map(p3, "map3_optimal_tp.png")

p4 <- grey_pts(base_map(), not_step2) +
  ggnewscale::new_scale_color() +
  geom_point(data = step2, aes(x = longitude, y = latitude, color = opt_rho),
             size = 1.8, alpha = 0.9, shape = 16) +
  scale_color_viridis_c(name = expression(rho), option = "plasma") +
  labs(title = "第二步：optimal tp 处的 rho(仅通过站)",
       subtitle = sprintf("%d/%d 站通过(灰点=未通过)", nrow(step2), nrow(all_stations)))
save_map(p4, "map4_optimal_rho.png")

log_msg("四张图全部完成")
