#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 15_covariates_1km.R
#
# 基于站点 1000 m 缓冲区提取自变量，输出两级结果：
#   (1) 网格级 —— 每个 HCSIF 500 m 像元一行(与 SIF 的取值网格一一对应)
#   (2) 站点级 —— 每个气象站一行(网格级按面积汇总)
#
# 当前覆盖:
#   DEM       本地 data_raw/DEM1km.tif (EPSG:4326, ~1 km)
#   CLCD      Zenodo COG 直读(HTTP range)，30 m，9 类土地利用，含不透水面
#
# CLCD 说明: 文件是 Cloud Optimized GeoTIFF，用 /vsicurl/ 按需读取窗口，
# 无需下载 18.4 GB 全量(实测 Zenodo 仅 221 KB/s，全下要 23 小时)。
# Zenodo 的 HEAD 不声明 Accept-Ranges，必须设 CPL_VSIL_CURL_USE_HEAD=NO，
# 否则 GDAL 判定不可分块读、退化成下载整个文件。
#
# 用法:
#   Rscript 15_covariates_1km.R dem
#   Rscript 15_covariates_1km.R clcd <年份> [站点数(测试用)]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
GRID_SHP <- file.path(PROJ, "data_raw/hcsif/grid/hcsif_grid_buf1000.shp")
OUT_DIR  <- file.path(PROJ, "data_raw/covariates_1km")
ZEN      <- "https://zenodo.org/api/records/15853565/files/CLCD_v01_%d_albert.tif/content"

CLCD_CLASS <- c("1" = "cropland", "2" = "forest",  "3" = "shrub",
                "4" = "grassland","5" = "water",   "6" = "snowice",
                "7" = "barren",   "8" = "impervious", "9" = "wetland")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# 1000 m 缓冲区覆盖的 HCSIF 像元多边形(每站约 15 个)
load_cells <- function() {
  v <- terra::vect(GRID_SHP)
  log_msg("读入网格: ", nrow(v), " 个像元多边形, ", length(unique(v$stat_id)), " 个站点")
  v
}

# ===========================================================================
# DEM
# ===========================================================================
run_dem <- function() {
  v <- load_cells()
  r <- terra::rast(file.path(PROJ, "data_raw/DEM1km.tif"))
  log_msg("DEM: ", paste(dim(r)[1:2], collapse = " x "), " | 分辨率 ", signif(terra::res(r)[1], 6))

  # DEM 约 1 km，而单个 HCSIF 像元仅 500 m —— 多数像元只覆盖 1 个 DEM 像元。
  # 因此网格级取该像元中心的高程，站点级再按缓冲区内所有像元汇总。
  ctr <- terra::centroids(v)
  if (!terra::same.crs(ctr, terra::crs(r))) ctr <- terra::project(ctr, terra::crs(r))
  elev <- terra::extract(r, ctr, ID = FALSE)[[1]]

  cell <- data.table(stat_id = v$stat_id, cell = v$cell, elev = elev)
  fwrite(cell, file.path(OUT_DIR, "dem_cell_1km.csv"))
  log_msg("网格级写出: dem_cell_1km.csv (", nrow(cell), " 行)")

  st <- cell[, .(
    elev_mean  = mean(elev, na.rm = TRUE),
    elev_min   = min(elev,  na.rm = TRUE),
    elev_max   = max(elev,  na.rm = TRUE),
    elev_sd    = sd(elev,   na.rm = TRUE),      # 缓冲区内起伏度(1 km DEM 下偏保守)
    n_cell     = .N,
    n_valid    = sum(!is.na(elev))
  ), by = stat_id]
  st[is.nan(elev_mean), elev_mean := NA_real_]
  fwrite(st, file.path(OUT_DIR, "dem_station_1km.csv"))
  log_msg("站点级写出: dem_station_1km.csv (", nrow(st), " 站)")
  print(summary(st$elev_mean))
}

# ===========================================================================
# CLCD
# ===========================================================================
run_clcd <- function(year, n_test = NA_integer_) {
  terra::setGDALconfig("CPL_VSIL_CURL_USE_HEAD", "NO")
  terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
  terra::setGDALconfig("GDAL_HTTP_MAX_RETRY", "5")
  terra::setGDALconfig("GDAL_HTTP_RETRY_DELAY", "3")

  v <- load_cells()
  if (!is.na(n_test)) {
    keep <- sort(unique(v$stat_id))[seq_len(n_test)]
    v <- v[v$stat_id %in% keep]
    log_msg("测试模式: 仅 ", n_test, " 个站点, ", nrow(v), " 个像元")
  }

  url <- sprintf(ZEN, year)
  r   <- terra::rast(paste0("/vsicurl/", url))
  log_msg("CLCD ", year, " 打开成功 | 分辨率 ", terra::res(r)[1], " m")

  # 把像元多边形投到 CLCD 的 Albers 上，保持 CLCD 原生网格不重采样
  # (类别栅格重采样会制造不存在的类别)
  vp <- terra::project(v, terra::crs(r))

  sids  <- sort(unique(vp$stat_id))
  t0    <- Sys.time()
  parts <- vector("list", length(sids))
  names(parts) <- as.character(sids)

  # 单站提取；失败(多为 Zenodo 429 限流)返回 NULL 交由外层重试
  grab <- function(sid) {
    sub <- vp[vp$stat_id == sid]
    win <- try(terra::crop(r, terra::ext(sub) * 1.02), silent = TRUE)  # 略放大避免边界丢失
    if (inherits(win, "try-error")) return(NULL)
    ex <- try(terra::extract(win, sub), silent = TRUE)
    if (inherits(ex, "try-error")) return(NULL)
    ex <- as.data.table(ex)
    setnames(ex, 1:2, c("ID", "cls"))
    ex <- ex[!is.na(cls)]
    if (!nrow(ex)) return(NULL)

    tab <- ex[, .N, by = .(ID, cls)]
    tot <- tab[, .(n_px = sum(N)), by = ID]
    tab <- merge(tab, tot, by = "ID")
    tab[, frac := N / n_px]
    tab[, cls_name := CLCD_CLASS[as.character(cls)]]

    w <- dcast(tab, ID + n_px ~ cls_name, value.var = "frac", fill = 0)
    w[, `:=`(stat_id = sub$stat_id[ID], cell = sub$cell[ID])]
    w
  }

  # 多轮重试: Zenodo 对密集 range 请求会返回 429，单轮必然漏站。
  # 每轮只补上一轮失败的站，轮间退避。
  todo <- as.character(sids)
  for (pass in 1:4) {
    for (j in seq_along(todo)) {
      sid <- todo[j]
      parts[[sid]] <- grab(as.integer(sid))
      if (j %% 50 == 0)
        log_msg(sprintf("  [第%d轮 %3d/%3d] 累计 %.1f 分钟", pass, j, length(todo),
                        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      Sys.sleep(0.05)                      # 轻微节流，降低 429 概率
    }
    todo <- names(parts)[vapply(parts, is.null, logical(1))]
    log_msg(sprintf("第 %d 轮后: 成功 %d/%d 站", pass, length(sids) - length(todo), length(sids)))
    if (!length(todo)) break
    if (pass < 4) { log_msg("  退避 90 秒后重试 ", length(todo), " 个站"); Sys.sleep(90) }
  }

  if (length(todo)) {
    log_msg("!! ", year, " 仍有 ", length(todo), " 个站失败，不写出残缺结果")
    stop("年份 ", year, " 提取不完整，请重跑")
  }

  cell <- rbindlist(parts[!vapply(parts, is.null, logical(1))], fill = TRUE)

  # 某个类别在部分网格里不出现时，dcast 不会生成该列，rbindlist(fill=TRUE)
  # 于是填入 NA。但"该类别占比为 0"才是正确语义 —— 留着 NA 会让后续
  # 加权汇总整列变成 NA，也会让各类占比之和不等于 1。
  for (cn in unname(CLCD_CLASS)) {
    if (!cn %in% names(cell)) data.table::set(cell, j = cn, value = 0)
    else data.table::set(cell, i = which(is.na(cell[[cn]])), j = cn, value = 0)
  }
  data.table::set(cell, j = "ID", value = NULL)
  data.table::set(cell, j = "year", value = year)
  setcolorder(cell, c("stat_id", "cell", "n_px", unname(CLCD_CLASS)))

  rs <- rowSums(cell[, unname(CLCD_CLASS), with = FALSE])
  if (max(abs(rs - 1)) > 1e-6)
    log_msg("警告: ", sum(abs(rs - 1) > 1e-6), " 个网格的类别占比之和不为 1")

  suffix <- if (is.na(n_test)) "" else "_test"
  fp1 <- file.path(OUT_DIR, sprintf("clcd_cell_%d%s.csv", year, suffix))
  fwrite(cell, fp1)
  log_msg("网格级写出: ", basename(fp1), " (", nrow(cell), " 行)")

  # 站点级: 按像元内 CLCD 像元数加权，等价于整个缓冲区的类别占比
  st <- cell[, c(list(n_px = sum(n_px)),
                 lapply(.SD, function(x) sum(x * n_px) / sum(n_px))),
             by = stat_id, .SDcols = unname(CLCD_CLASS)]
  data.table::set(st, j = "year", value = year)
  st_rs <- rowSums(st[, unname(CLCD_CLASS), with = FALSE])
  log_msg(sprintf("站点级占比之和: 均值 %.6f  范围 [%.6f, %.6f]",
                  mean(st_rs), min(st_rs), max(st_rs)))
  fp2 <- file.path(OUT_DIR, sprintf("clcd_station_%d%s.csv", year, suffix))
  fwrite(st, fp2)
  log_msg("站点级写出: ", basename(fp2), " (", nrow(st), " 站)")

  cat("\n各类占比(站点级均值):\n")
  print(round(sapply(st[, unname(CLCD_CLASS), with = FALSE], mean, na.rm = TRUE), 4))
  invisible(st)
}

# ===========================================================================
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("用法: Rscript 15_covariates_1km.R dem | clcd <年份> [站点数]")
if (args[1] == "dem") {
  run_dem()
} else if (args[1] == "clcd") {
  yr <- as.integer(args[2])
  nt <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
  run_clcd(yr, nt)
} else stop("未知参数: ", args[1])
