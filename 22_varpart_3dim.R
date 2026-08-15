#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 22_varpart_3dim.R — 按"土地利用与空间结构 / 气候与自然环境 / 社会经济与人类活动"
#                     三个维度，对 VPD->SIF 因果强度做方差分解与影响因素分析
#
# 站点到城市用 city.shp 做空间连接(旧协变量表只覆盖 683 站, 空间连接可覆盖全部)。
# 城市名在三套表里格式不一("齐齐哈尔市" / "七台河"), 统一去掉行政后缀再匹配。
#
# 每个变量都报覆盖率。方差分解对样本量敏感, 所以跑两套:
#   CORE  只用覆盖率 >=90% 的变量, 样本接近全量
#   FULL  纳入全部变量, 样本按完整记录缩减
# 两套并报, 避免把"样本变了"误读成"解释力变了"。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(terra); library(ranger); library(readxl)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)

# ===== 1. 响应变量 ========================================================
cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
resp <- cm[, .(rho_tp0 = rho[tp == 0][1], rho_max = max(rho, na.rm = TRUE),
               coef_tp0 = mean_coef[tp == 0][1]), by = .(stat_id = meteo_stat)]
resp <- resp[is.finite(rho_max) & !is.na(rho_tp0)]

st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")

# ===== 2. 维度一: 土地利用与空间结构 ======================================
g <- fread("data_raw/covariates_1km/glc_station_2020.csv")
L <- function(...) rowSums(as.matrix(g[, sprintf("LC%02d", c(...)), with = FALSE]))
d1 <- data.table(
  stat_id = g$stat_id,
  imperv     = g$LC23,                      # 不透水面(GLC_FCS30D)
  forest     = L(4:13), grass = L(17), crop = L(0:3), water = L(27),
  # "树种": 森林按叶型/物候拆分, 均以缓冲区面积为分母(避免无林站点除零)
  needleleaf = L(8:11), broadleaf = L(4:7), mixedleaf = L(12:13),
  evergreen  = L(4, 5, 8, 9), deciduous = L(6, 7, 10, 11))

ws <- fread("data_raw/WS土地利用比例_1000m.csv", encoding = "UTF-8")
setnames(ws, names(ws), sub("^﻿", "", names(ws)))
setnames(ws, names(ws), sub("^(LC[0-9]{2})_.*$", "\\1", names(ws)))
d1 <- merge(d1, ws[, .(stat_id, imperv_ws = LC23)], by = "stat_id", all.x = TRUE)

# ===== 3. 维度二: 气候与自然环境 ==========================================
mt <- fread("data_raw/covariates_1km/meteo_station_agg.csv")[
  , .(stat_id, tavg, rh, cloud, precip)]
rs <- rbindlist(lapply(2000:2022, function(y)
  fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rs <- rs[, .(rsds = mean(rsds_mean)), by = .(stat_id, year)][
  , .(rsds_mean = mean(rsds), rsds_sd = sd(rsds)), by = stat_id]
dem <- fread("data_raw/covariates_1km/dem_cell_1km.csv")[
  , .(elev = mean(elev, na.rm = TRUE)), by = stat_id]
d2 <- Reduce(function(a, b) merge(a, b, by = "stat_id", all = TRUE), list(mt, rs, dem))

# ===== 4. 维度三: 社会经济与人类活动 ======================================
nt <- fread("data_raw/mvnl_monthly_buffer1000_2000_2024.csv", encoding = "UTF-8")
setnames(nt, names(nt), sub("^﻿", "", names(nt)))
nt <- nt[year == 2020 & quality_flag == "OK",
         .(ntl = log1p(mean(area_weighted_mean))), by = stat_id]

# 站点 -> 城市: 空间连接
cv <- vect("data_raw/china_cities/city.shp")
pt <- vect(st, geom = c("longitude", "latitude"), crs = crs(cv))
# city.shp 只含地级市城区多边形, 大量站点落在城区之外; 先取落入的,
# 落不进去的按最近城区赋值——社会经济属性本来就是按所属城市赋给站点的。
hit <- terra::relate(pt, cv, "intersects")
idx <- apply(hit, 1, function(r) { w <- which(r); if (length(w)) w[1] else NA_integer_ })
n_in <- sum(!is.na(idx))
if (anyNA(idx)) idx[is.na(idx)] <- terra::nearest(pt[is.na(idx)], cv)$to_id
smap <- data.table(stat_id = st$stat_id, city = strip_city(cv$ct_name[idx]))
smap <- unique(smap, by = "stat_id")
log_msg("落入城区多边形的站: ", n_in, " / ", nrow(smap), "; 其余按最近城区赋值")

iv <- fread("data_raw/green_invest/city_invest_data.csv", encoding = "UTF-8")
setnames(iv, 1, "city_raw"); iv[, city := strip_city(city_raw)]
ycol <- grep("^园林绿化_20(0[3-9]|1[0-9]|2[0-2])$", names(iv), value = TRUE)
iv[, invest := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = ycol]
iv <- iv[is.finite(invest), .(invest = mean(invest)), by = city]

# 城镇化率列前几千行(1990 年代)全空, readxl 会把整列猜成 logical, 需放大猜测行数
cd <- as.data.table(read_excel("data_raw/中国城市数据库1990-2023.xlsx",
                               sheet = "原始数据", guess_max = 1e6))
# 现成的"常住人口城镇化率"在 2000-2022 内全空(只有 2023 一年 297 条),
# 改用户籍口径的非农业人口比重作代理: 覆盖 288 城 x 2488 城年。
# 这是与常住口径不同的定义(偏低), 但作横截面排序用是一致的。
num <- function(x) suppressWarnings(as.numeric(x))
cd <- cd[as.integer(num(年份)) %in% 2000:2022,
         .(city = strip_city(城市),
           nong = num(`非农业人口数(万人)`), hu = num(`户籍人口(万人)`))]
cd <- cd[is.finite(nong) & is.finite(hu) & hu > 0]
cd <- cd[, .(urban_rate = mean(100 * nong / hu)), by = city]

d3 <- Reduce(function(a, b) merge(a, b, by = "city", all.x = TRUE),
             list(smap, iv, cd))[, .(stat_id, invest = log1p(invest), urban_rate)]
d3 <- unique(merge(nt, d3, by = "stat_id", all = TRUE), by = "stat_id")

# ===== 5. 旧协变量(建筑/路网/气候区划) ====================================
old <- as.data.table(readRDS("data_proc/output_10y_built_up_05_01/station_covariates.rds"))
setnames(old, "meteo_stat_id", "stat_id"); old[, stat_id := as.integer(stat_id)]
old <- unique(old, by = "stat_id")     # 旧表有 6 个重复站号
old <- old[, .(stat_id, koppen_group, road_density, building_footprint,
               mean_height, building_vol_density)]

dt <- Reduce(function(a, b) merge(a, b, by = "stat_id", all.x = TRUE),
             list(resp, d1, d2, d3, old))
log_msg("总表: ", nrow(dt), " 站\n")

# ===== 6. 覆盖率 ==========================================================
VARS <- c("imperv","imperv_ws","forest","grass","crop","water","needleleaf","broadleaf",
          "mixedleaf","evergreen","deciduous","road_density","building_footprint",
          "mean_height","building_vol_density",
          "tavg","rh","cloud","precip","rsds_mean","rsds_sd","elev",
          "ntl","invest","urban_rate")
cov_tab <- data.table(var = VARS,
                      cover = sapply(VARS, function(v) mean(!is.na(dt[[v]]))))
setorder(cov_tab, cover)
log_msg("-- 各变量覆盖率 --")
print(cov_tab[, .(var, cover = round(cover, 3))])
fwrite(cov_tab, file.path(OUT, "varpart3_coverage.csv"))

stopifnot(!any(duplicated(dt$stat_id)))

DIM <- list(
  landuse = c("imperv","forest","grass","crop","water","needleleaf","broadleaf",
              "evergreen","road_density","building_footprint","mean_height"),
  climate = c("tavg","rh","cloud","precip","rsds_mean","rsds_sd","elev"),
  socioec = c("ntl","invest","urban_rate"))
CORE <- lapply(DIM, function(v) v[sapply(v, function(x) mean(!is.na(dt[[x]])) >= .90)])

winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

run <- function(blocks, yv, tag) {
  vars <- unlist(blocks, use.names = FALSE)
  d <- dt[, c(yv, vars), with = FALSE][complete.cases(dt[, c(yv, vars), with = FALSE])]
  for (v in vars) set(d, j = v, value = winz(as.numeric(d[[v]])))
  f <- as.formula(paste(yv, "~", paste(vars, collapse = "+")))
  R2 <- summary(lm(f, d))$r.squared
  vp <- rbindlist(lapply(names(blocks), function(b) {
    rest <- setdiff(vars, blocks[[b]])
    r2w <- if (!length(rest)) 0 else summary(lm(as.formula(paste(yv,"~",paste(rest,collapse="+"))), d))$r.squared
    r2a <- summary(lm(as.formula(paste(yv,"~",paste(blocks[[b]],collapse="+"))), d))$r.squared
    data.table(dim = b, uniq = R2 - r2w, alone = r2a)
  }))
  rf <- ranger(f, d, num.trees = 1000, importance = "permutation", seed = 42)
  imp <- data.table(var = names(rf$variable.importance),
                    rel = as.numeric(rf$variable.importance))
  imp[, rel := rel / sum(pmax(rel, 0))]; imp <- imp[order(-rel)]
  z <- copy(d); for (v in vars) set(z, j = v, value = as.numeric(scale(z[[v]])))
  co <- coef(summary(lm(f, z)))
  std <- data.table(var = rownames(co)[-1], beta = co[-1,1], p = co[-1,4])[order(-abs(beta))]

  log_msg("\n============ ", tag, " | ", yv, " | n=", nrow(d), " ============")
  log_msg("线性 R2 = ", round(R2,4), " | 随机森林 R2 = ", round(rf$r.squared,4))
  log_msg("-- 三维度方差分解 --")
  print(vp[order(-uniq)][, .(dim, uniq = round(uniq,4), alone = round(alone,4))])
  log_msg("独占之和 ", round(sum(vp$uniq),4), " | 总 R2 ", round(R2,4),
          " | 维度间共享 ", round(R2 - sum(vp$uniq),4))
  log_msg("-- 标准化系数(前10) --")
  print(std[1:min(10,.N), .(var, beta = round(beta,4), p = signif(p,3),
                            sig = ifelse(p<.05,"*",""))])
  log_msg("-- 随机森林重要性(前10) --")
  print(imp[1:min(10,.N), .(var, rel = round(rel,4))])
  fwrite(vp,  file.path(OUT, sprintf("varpart3_%s_%s_dim.csv",  tag, yv)))
  fwrite(std, file.path(OUT, sprintf("varpart3_%s_%s_coef.csv", tag, yv)))
  fwrite(imp, file.path(OUT, sprintf("varpart3_%s_%s_rfimp.csv",tag, yv)))
}

for (yv in c("rho_tp0","rho_max","coef_tp0")) run(CORE, yv, "CORE")
for (yv in c("rho_tp0","rho_max","coef_tp0")) run(DIM,  yv, "FULL")

# ===== 7. 不透水面来源敏感性: GLC vs 王琳 =================================
log_msg("\n============ 不透水面来源敏感性 ============")
for (yv in c("rho_tp0","rho_max")) {
  d <- dt[!is.na(imperv) & !is.na(imperv_ws) & is.finite(get(yv))]
  log_msg(sprintf("%s: 与 GLC 不透水面 r=%.4f | 与王琳不透水面 r=%.4f  (n=%d)",
                  yv, cor(d[[yv]], d$imperv), cor(d[[yv]], d$imperv_ws), nrow(d)))
}
log_msg("\n完成")
