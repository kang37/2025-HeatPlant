#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20_drivers_all.R — 把本轮新采集的协变量并入 VPD->SIF 因果强度的影响因素分析
#
# 响应变量(站点级, 来自 1000 m 缓冲区的 CCM 结果):
#   rho_tp0   同期因果强度
#   rho_max   各滞后中的最大因果强度
#   coef_tp0  同期因果系数(带符号: VPD 对 SIF 是抑制还是促进)
#
# 解释变量按来源分块, 便于做方差分解:
#   地形      DEM
#   气候      ERA5 辐射(23 年均值与年际标准差)、降水、柯本气候带
#   地表覆盖  GLC_FCS30D 2020 (不透水面/林/草/耕地/水体)
#   城市结构  夜间灯光、路网密度、建筑基底与体积、10 年建成区扩张
#   社会经济  人口、人均 GDP、绿地投资
#   土壤      粘粒、砂粒
#
# 两套样本:
#   A 全量  只用本轮新变量(地形+气候+地表覆盖+灯光), 904 站
#   B 完整  并入旧协变量(社会经济/建筑/土壤), 样本降到旧表覆盖的站点
# 两套都报, 因为 B 的解释力提升可能只是样本变了而非变量变好。
#
# 方差分解用分层分割: 每个块的"独占贡献" = 全模型 R2 - 去掉该块的 R2。
# 块间相关会使独占贡献之和小于总 R2, 差额即共享部分, 一并报出。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ranger)
})

PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
setwd(PROJ)
OUT <- "data_proc/output_hcsif_buf1000"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CCM <- "data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv"
GLC_YEAR <- 2020L; NTL_YEAR <- 2020L
set.seed(42)
log_msg <- function(...) cat(...,"\n", sep="")

# ===== 1. 响应变量 =========================================================
cm <- fread(CCM)
resp <- cm[, .(rho_tp0  = rho[tp == 0][1],
               rho_max  = max(rho, na.rm = TRUE),
               coef_tp0 = mean_coef[tp == 0][1],
               n_tp     = .N), by = .(stat_id = meteo_stat)]
resp <- resp[is.finite(rho_max) & !is.na(rho_tp0)]
log_msg("响应变量: ", nrow(resp), " 站")

# ===== 2. 新协变量 =========================================================
dem <- fread("data_raw/covariates_1km/dem_cell_1km.csv")[, .(elev = mean(elev, na.rm = TRUE)),
                                                         by = stat_id]

rs <- rbindlist(lapply(2000:2022, function(y)
  fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rs <- rs[, .(rsds = mean(rsds_mean, na.rm = TRUE)), by = .(stat_id, year)]
rs <- rs[, .(rsds_mean = mean(rsds), rsds_sd = sd(rsds)), by = stat_id]

g <- fread(sprintf("data_raw/covariates_1km/glc_station_%d.csv", GLC_YEAR))
LCn <- function(...) rowSums(as.matrix(g[, sprintf("LC%02d", c(...)), with = FALSE]))
glc <- data.table(stat_id = g$stat_id,
                  imperv = g$LC23,
                  forest = LCn(4:13),
                  grass  = LCn(17),
                  crop   = LCn(0:3),
                  water  = LCn(27),
                  shrub  = LCn(14:16))
glc[, veg := forest + grass + crop + shrub]

nt <- fread(sprintf("data_raw/covariates_1km/ntl_station_%d.csv", NTL_YEAR))
nt <- nt[, .(stat_id, ntl = log1p(ntl_mean))]      # 灯光重尾, 取对数

new <- Reduce(function(a, b) merge(a, b, by = "stat_id", all = TRUE),
              list(resp, dem, rs, glc, nt))
log_msg("并入新变量后: ", nrow(new), " 站; 完整记录 ", sum(complete.cases(new)), " 站")

# ===== 3. 旧协变量 =========================================================
old <- as.data.table(readRDS("data_proc/output_10y_built_up_05_01/station_covariates.rds"))
setnames(old, "meteo_stat_id", "stat_id")
old[, stat_id := as.integer(stat_id)]              # 旧表里站号存成了字符
old <- old[, .(stat_id, koppen_group, precip_mean, soil_clay, soil_sand,
               pa_built_10y, cgi_score, pop_10y, pgdp_10y, road_density,
               building_footprint, mean_height, building_vol_density)]
full <- merge(new, old, by = "stat_id")
log_msg("并入旧变量后: ", nrow(full), " 站")

# ===== 4. 变量分块 =========================================================
BLK_A <- list("terrain" = "elev", "climate" = c("rsds_mean", "rsds_sd"),
              "landcover" = c("imperv", "forest", "grass", "crop", "water"),
              "urban" = "ntl")
BLK_B <- list("terrain" = "elev",
              "climate" = c("rsds_mean", "rsds_sd", "precip_mean"),
              "landcover" = c("imperv", "forest", "grass", "crop", "water"),
              "urban" = c("ntl", "road_density", "building_footprint",
                           "mean_height", "building_vol_density", "pa_built_10y"),
              "socioecon" = c("pop_10y", "pgdp_10y", "cgi_score"),
              "soil" = c("soil_clay", "soil_sand"))

winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

# ===== 5. 方差分解 + 回归 + 随机森林 ======================================
run_one <- function(dt, blocks, yv, tag) {
  vars <- unlist(blocks, use.names = FALSE)
  d <- dt[, c(yv, vars), with = FALSE]
  d <- d[complete.cases(d)]
  for (v in vars) set(d, j = v, value = winz(as.numeric(d[[v]])))
  n <- nrow(d)
  f  <- as.formula(paste(yv, "~", paste(vars, collapse = "+")))
  m  <- lm(f, d); R2 <- summary(m)$r.squared

  vp <- rbindlist(lapply(names(blocks), function(b) {
    rest <- setdiff(vars, blocks[[b]])
    r2_wo <- if (!length(rest)) 0 else
      summary(lm(as.formula(paste(yv, "~", paste(rest, collapse = "+"))), d))$r.squared
    r2_al <- summary(lm(as.formula(paste(yv, "~", paste(blocks[[b]], collapse = "+"))), d))$r.squared
    data.table(blk = b, uniq = R2 - r2_wo, alone = r2_al)
  }))

  rf <- ranger(f, d, num.trees = 1000, importance = "permutation", seed = 42)
  imp <- data.table(var = names(rf$variable.importance),
                    impt = as.numeric(rf$variable.importance))
  imp[, rel := impt / sum(pmax(impt, 0))]
  setorder(imp, -impt)

  z <- copy(d); for (v in vars) set(z, j = v, value = as.numeric(scale(z[[v]])))
  mz <- lm(f, z); co <- coef(summary(mz))
  std <- data.table(var = rownames(co)[-1], beta = co[-1, 1], p = co[-1, 4])
  std <- std[order(-abs(beta))]

  log_msg("\n================ ", tag, " | 响应=", yv, " | n=", n, " ================")
  log_msg("线性模型 R2 = ", round(R2, 4), "  |  随机森林 R2 = ", round(rf$r.squared, 4))
  log_msg("\n-- 方差分解(块) --")
  vp2 <- copy(vp); vp2[, `:=`(uniq = round(uniq, 4), alone = round(alone, 4))]
  print(vp2[order(-uniq)])
  log_msg("独占之和 ", round(sum(vp$uniq), 4), " vs 总 R2 ", round(R2, 4),
          " -> 块间共享 ", round(R2 - sum(vp$uniq), 4))
  log_msg("\n-- 标准化回归系数(前8) --")
  print(std[1:min(8, .N), .(var, beta = round(beta, 4),
                            p = signif(p, 3), sig = ifelse(p < .05, "*", ""))])
  log_msg("\n-- 随机森林置换重要性(前8) --")
  print(imp[1:min(8, .N), .(var, rel = round(rel, 4))])

  fwrite(vp,  file.path(OUT, sprintf("drivers_all_%s_%s_varpart.csv", tag, yv)))
  fwrite(std, file.path(OUT, sprintf("drivers_all_%s_%s_coef.csv",    tag, yv)))
  fwrite(imp, file.path(OUT, sprintf("drivers_all_%s_%s_rfimp.csv",   tag, yv)))
  invisible(list(n = n, R2 = R2, rf = rf$r.squared))
}

for (yv in c("rho_tp0", "rho_max", "coef_tp0")) run_one(new,  BLK_A, yv, "A_new_only")
for (yv in c("rho_tp0", "rho_max", "coef_tp0")) run_one(full, BLK_B, yv, "B_with_old")

# ===== 6. 按柯本气候带的单因素方差分析 ====================================
log_msg("\n================ 柯本气候带 单因素方差分析 ================")
for (yv in c("rho_tp0", "rho_max", "coef_tp0")) {
  d <- full[!is.na(koppen_group) & is.finite(get(yv))]
  d[, kg := factor(koppen_group)]
  keep <- d[, .N, by = kg][N >= 10]$kg          # 样本过少的气候带并不稳健, 剔除
  d <- d[kg %in% keep]; d[, kg := droplevels(kg)]
  a <- anova(lm(as.formula(paste(yv, "~ kg")), d))
  log_msg(sprintf("%s: F(%d,%d)=%.2f, p=%.3g, eta2=%.4f  (%d 组, n=%d)",
                  yv, a$Df[1], a$Df[2], a$`F value`[1], a$`Pr(>F)`[1],
                  a$`Sum Sq`[1] / sum(a$`Sum Sq`), nlevels(d$kg), nrow(d)))
  print(d[, .(n = .N, mu = round(mean(get(yv)), 4)), by = kg][order(-mu)])
}
log_msg("\n完成")
