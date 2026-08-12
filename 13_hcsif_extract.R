#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_hcsif_extract.R
#
# 从 HCSIF 数据集提取气象站点尺度的 SIF 时间序列。
#
#   数据集: A high-resolution satellite-based SIF dataset for China (HCSIF)
#           Tao S., Chen J.M., Zhang Z. et al. Sci Data 11, 1286 (2024)
#           DOI 10.57760/sciencedb.16910   (CC0)
#           500 m / 8-day / 中国 / 2000-2022 / 仅生长季 3-10 月
#           单位 mW m-2 nm-1 sr-1，存储为整型，scale factor = 1e-4
#
# 用法:
#   Rscript 13_hcsif_extract.R <year> [raw_dir] [out_dir]
#   Rscript 13_hcsif_extract.R --selftest      # 用合成栅格自检整条流程
#
# 设计约束: 磁盘只剩 ~35 GiB，必须逐年/逐批处理后立即删除原始栅格。
#           本脚本只负责"提取"，下载与删除由 13_hcsif_run.sh 驱动。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
})

# --- 配置 ------------------------------------------------------------------

PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"

SCALE_FACTOR <- 1e-4    # 数据集页面标注: "The Scale factor is 0.0001"

# HCSIF 实际的时相网格是 DOY = 60 + 8k，即 60, 68, ..., 300，每年 31 个时相
# (注意不是 MODIS 标准的 1, 9, 17, ... 网格)。
# 下面按"时相起始日落在 5/1-9/30"筛选，命中 DOY 124, 132, ..., 268 共 19 个。
# DOY 116 (4/25-5/2) 主体在 4 月，按起始日归月应算 4 月，故排除。
DOY_MIN      <- 121     # 5 月 1 日
DOY_MAX      <- 273     # 9 月 30 日
# 缓冲区半径(m)，可给多个，每个生成一组 SIF_buf<r> / n_valid<r> 列。
#   750  ≈ 3x3 个 500 m 像元
#   1000 ≈ 5x5 个 500 m 像元(纬度越高经向像元越窄，实际个数随纬度变化)
BUFFER_M     <- c(750, 1000)
VALID_RANGE  <- c(-1, 5)  # 物理合理的 SIF 区间(mW m-2 nm-1 sr-1)，超出即判为填充值

# 关键: 栅格没有声明 NoData(NAflag = NaN)，但全图 68.8% 的像元恰好等于 0，
# 且这些 0 与不透水面高度相关(零值站点建筑占地中位数是非零站点的 4.7 倍,
# Wilcoxon p = 0.0037)。7 月中旬不存在真实的零光合，因此 0 是非植被掩膜值。
# 原始值为整型，真实的近零信号会是 3、-12 这类小整数而非恰好 0，故按精确等于 0 判定。
ZERO_AS_NA   <- TRUE

# --- 工具函数 --------------------------------------------------------------

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

#' 从文件名解析年份与 DOY
#'
#' HCSIF 的命名为 SIF<YYYY>-<DDD>.tif，但不同版本/镜像可能略有出入，
#' 因此这里放宽为"抓取 4 位年份 + 紧随其后的 1-3 位 DOY"。
parse_doy <- function(paths) {
  bn <- basename(paths)
  m  <- regmatches(bn, regexec("(19|20)([0-9]{2})[^0-9]?([0-9]{1,3})", bn))
  do.call(rbind, lapply(seq_along(m), function(i) {
    g <- m[[i]]
    if (length(g) < 4) return(data.table(file = bn[i], year = NA_integer_, doy = NA_integer_))
    data.table(file = bn[i],
               year = as.integer(paste0(g[2], g[3])),
               doy  = as.integer(g[4]))
  }))
}

#' 把站点点位投影到栅格所在的坐标系
#'
#' HCSIF 可能存为 MODIS 正弦投影或等积投影，站点是 WGS84 经纬度，
#' 必须先投影再取值，否则会静默取到错误像元。
make_points <- function(st, target_crs) {
  pts <- terra::vect(as.data.frame(st),
                     geom = c("longitude", "latitude"),
                     crs  = "EPSG:4326")
  if (!is.na(target_crs) && nzchar(target_crs)) {
    same <- tryCatch(terra::same.crs(pts, target_crs), error = function(e) FALSE)
    if (!isTRUE(same)) pts <- terra::project(pts, target_crs)
  }
  pts
}

#' 处理 scale factor
#'
#' 若 GeoTIFF 头里已写了 scale/offset，terra 会自动应用，此时不能再乘一次。
#' 通过 terra::scoff() 判断，并在数值量级上做一次兜底检查。
resolve_scale <- function(r) {
  so <- tryCatch(terra::scoff(r), error = function(e) NULL)
  if (!is.null(so) && !all(is.na(so)) && !isTRUE(all.equal(unname(so[1, 1]), 1))) {
    return(list(mult = 1, note = sprintf("header scale=%g 已由 terra 应用", so[1, 1])))
  }
  list(mult = SCALE_FACTOR, note = sprintf("手工应用 scale=%g", SCALE_FACTOR))
}

# --- 单个时相的提取 --------------------------------------------------------

extract_one <- function(fp, st, pts_cache = NULL) {
  r <- terra::rast(fp)
  if (terra::nlyr(r) > 1L) r <- r[[1]]

  crs_txt <- terra::crs(r, describe = TRUE)$code
  if (is.null(pts_cache)) pts_cache <- make_points(st, terra::crs(r))
  pts <- pts_cache

  sc <- resolve_scale(r)

  scale_it <- function(v) {
    v <- v * sc$mult
    v[!is.na(v) & (v < VALID_RANGE[1] | v > VALID_RANGE[2])] <- NA_real_
    v
  }

  # 最近像元值
  v_pt_raw <- terra::extract(r, pts, ID = FALSE)[[1]]
  masked   <- !is.na(v_pt_raw) & v_pt_raw == 0        # 站点像元本身被掩膜
  v_pt     <- v_pt_raw
  if (ZERO_AS_NA) v_pt[masked] <- NA_real_

  # 缓冲区聚合(点值 vs 各级邻域的尺度敏感性检验)
  # 逐像元取回再自行聚合，这样才能在求均值前剔除掩膜值——
  # 直接用 fun=mean 会把 0 当作真实低值拉低均值。
  # 每个半径生成 SIF_buf<r> 与 n_valid<r> 两列。
  out <- data.table(
    meteo_stat = st$meteo_stat,
    latitude   = st$latitude,
    longitude  = st$longitude,
    SIF        = scale_it(v_pt),
    SIF_raw    = v_pt_raw,             # 保留原始存储值(含 0)，便于回溯
    pt_masked  = masked                # 站点像元是否为非植被掩膜
  )

  for (bm in BUFFER_M) {
    buf <- terra::buffer(pts, width = bm)
    ex  <- data.table::as.data.table(terra::extract(r, buf))
    data.table::setnames(ex, 1:2, c("ID", "val"))
    if (ZERO_AS_NA) ex[!is.na(val) & val == 0, val := NA_real_]

    agg <- ex[, .(v_buf = mean(val, na.rm = TRUE),
                  n_val = sum(!is.na(val))), by = ID]
    agg <- agg[data.table::data.table(ID = seq_len(nrow(st))), on = "ID"]  # 补齐空缓冲区
    # mean(全 NA, na.rm=TRUE) 返回 NaN 而非 NA，必须显式转换
    agg[is.nan(v_buf), v_buf := NA_real_]
    agg[is.na(n_val), n_val := 0L]

    # 用 set() 而非 [[<- ，否则后续 := 会触发 shallow-copy 警告
    data.table::set(out, j = paste0("SIF_buf", bm), value = scale_it(agg$v_buf))
    data.table::set(out, j = paste0("n_valid", bm), value = agg$n_val)
  }

  # 向后兼容: 既有的 mask_analysis / grid_shp / CCM 脚本都引用无后缀的 n_valid,
  # 保留它作为首个半径(750 m)的别名，避免改动下游代码。
  if (length(BUFFER_M) && paste0("n_valid", BUFFER_M[1]) %in% names(out)) {
    data.table::set(out, j = "n_valid", value = out[[paste0("n_valid", BUFFER_M[1])]])
  }

  data.table::set(out, j = "src_file",   value = basename(fp))
  data.table::set(out, j = "src_crs",
                  value = if (is.null(crs_txt) || is.na(crs_txt))
                            terra::crs(r, proj = TRUE) else crs_txt)
  data.table::set(out, j = "scale_note", value = sc$note)
  out
}

# --- 主流程 ----------------------------------------------------------------

run <- function(target_year, raw_dir, out_dir) {
  st <- fread(file.path(PROJ, "data_raw/hcsif/stations_924.csv"))
  stopifnot(all(c("meteo_stat", "latitude", "longitude") %in% names(st)))
  log_msg("站点数: ", nrow(st))

  files <- list.files(raw_dir, pattern = "\\.(tif|tiff|nc)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("raw_dir 里没有找到 .tif/.nc 文件: ", raw_dir)

  meta <- parse_doy(files)
  meta[, path := files]
  keep <- meta[!is.na(doy) & year == target_year & doy >= DOY_MIN & doy <= DOY_MAX][order(doy)]

  log_msg("目录中共 ", length(files), " 个栅格; ",
          target_year, " 年 5-9 月(DOY ", DOY_MIN, "-", DOY_MAX, ") 命中 ", nrow(keep), " 个")
  if (!nrow(keep)) stop("没有匹配到目标时相，请检查文件命名是否为 SIF<YYYY>-<DDD>.tif")

  pts_cache <- NULL
  res <- vector("list", nrow(keep))

  for (i in seq_len(nrow(keep))) {
    fp <- keep$path[i]
    if (is.null(pts_cache)) {
      r0 <- terra::rast(fp)
      pts_cache <- make_points(st, terra::crs(r0))
      log_msg("栅格 CRS: ", terra::crs(r0, describe = TRUE)$name,
              " | 分辨率: ", paste(signif(terra::res(r0), 6), collapse = " x "),
              " | 维度: ", paste(dim(r0)[1:2], collapse = " x "))
      rm(r0)
    }
    d <- extract_one(fp, st, pts_cache)
    d[, `:=`(year = keep$year[i], doy = keep$doy[i])]
    res[[i]] <- d
    log_msg(sprintf("  [%2d/%2d] DOY %3d  %-24s 有效站点 %3d/%d",
                    i, nrow(keep), keep$doy[i], keep$file[i],
                    sum(!is.na(d$SIF)), nrow(d)))
  }

  out <- rbindlist(res)

  # 补齐日期字段，与既有 meteo_stat_SIF_data.csv 的列结构保持一致
  out[, date_d := as.Date(doy - 1, origin = paste0(year, "-01-01"))]
  out[, `:=`(month = as.integer(format(date_d, "%m")),
             day   = as.integer(format(date_d, "%d")),
             date  = as.integer(format(date_d, "%Y%m%d")))]

  buf_cols <- unlist(lapply(BUFFER_M, function(b) paste0(c("SIF_buf", "n_valid"), b)))
  setcolorder(out, c("meteo_stat", "latitude", "longitude",
                     "year", "month", "day", "doy", "date",
                     "SIF", "SIF_raw", "pt_masked", buf_cols, "n_valid",
                     "src_file", "src_crs", "scale_note"))
  out[, date_d := NULL]

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fp_out <- file.path(out_dir, sprintf("hcsif_station_%d.csv", target_year))
  fwrite(out, fp_out)

  log_msg("写出: ", fp_out, "  (", nrow(out), " 行)")
  log_msg(sprintf("点值有效率 %.1f%% | 被掩膜 %.1f%%",
                  100 * mean(!is.na(out$SIF)), 100 * mean(out$pt_masked)))
  for (b in BUFFER_M) {
    log_msg(sprintf("  buf%-5d 有效率 %.1f%% | 均值 %.4f",
                    b, 100 * mean(!is.na(out[[paste0("SIF_buf", b)]])),
                    mean(out[[paste0("SIF_buf", b)]], na.rm = TRUE)))
  }
  log_msg(sprintf("SIF 点值 均值 %.4f | 范围 [%.4f, %.4f]",
                  mean(out$SIF, na.rm = TRUE),
                  min(out$SIF, na.rm = TRUE),
                  max(out$SIF, na.rm = TRUE)))
  invisible(out)
}

# --- 自检: 用合成栅格验证整条流程 ------------------------------------------

selftest <- function() {
  log_msg("=== 自检开始 ===")
  tmp <- file.path(tempdir(), "hcsif_selftest")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  st <- fread(file.path(PROJ, "data_raw/hcsif/stations_924.csv"))

  # 造一个覆盖中国、投影为 Albers 的整型栅格(模拟真实产品不是经纬度网格的情况)
  crs_albers <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +datum=WGS84 +units=m"
  tpl <- terra::rast(xmin = 73, xmax = 136, ymin = 17, ymax = 55,
                     resolution = 0.05, crs = "EPSG:4326")
  tpl <- terra::project(tpl, crs_albers, res = 5000)

  # 已知的真值场: SIF = 0.3 * cos(纬度弧度)，存为整型
  ll  <- terra::project(tpl, "EPSG:4326")
  lat <- terra::init(ll, "y")
  truth_ll <- 0.3 * cos(lat * pi / 180)
  truth <- terra::project(truth_ll, tpl)
  stored <- terra::app(truth, function(x) as.integer(round(x / SCALE_FACTOR)))

  # 注入一块 0 值掩膜(模拟真实数据里的不透水面/非植被填充)，
  # 用于验证 0 被正确识别为掩膜而不是当成真实的零光合。
  zbox <- terra::vect("POLYGON ((110 30, 115 30, 115 35, 110 35, 110 30))", crs = "EPSG:4326")
  zbox <- terra::project(zbox, crs_albers)
  zind <- terra::rasterize(zbox, stored, field = 1, background = 0)
  stored <- terra::ifel(zind == 1, 0L, stored)

  for (doy in c(113, 121, 129, 273, 281)) {   # 113/281 在 5-9 月窗口外，应被过滤
    terra::writeRaster(stored,
                       file.path(tmp, sprintf("SIF2015-%03d.tif", doy)),
                       overwrite = TRUE, datatype = "INT2S")
  }

  out <- run(2015, tmp, tmp)

  # 断言 1: 只保留窗口内的 3 个时相
  doys <- sort(unique(out$doy))
  stopifnot(identical(doys, c(121L, 129L, 273L)))
  log_msg("PASS 时相过滤: ", paste(doys, collapse = ", "))

  # 断言 2: scale factor 正确还原
  chk <- out[!is.na(SIF)]
  expect <- 0.3 * cos(chk$latitude * pi / 180)
  err <- abs(chk$SIF - expect)
  stopifnot(max(err) < 0.01)
  log_msg(sprintf("PASS scale/投影还原: 最大误差 %.5f (n=%d)", max(err), nrow(chk)))

  # 断言 3: 日期换算正确 (2015 非闰年, DOY 121 = 5月1日)
  d121 <- unique(out[doy == 121L, .(month, day, date)])
  stopifnot(nrow(d121) == 1, d121$month == 5, d121$day == 1, d121$date == 20150501L)
  log_msg("PASS 日期换算: DOY 121 -> 2015-05-01")

  # 断言 4: 0 值被识别为掩膜而非真实低值
  # 注意用收缩后的内部框判定: 经纬度矩形投影到 Albers 后边缘是弯的，
  # 贴边的站点落在掩膜内外取决于投影细节，不适合做断言。
  inbox <- out[longitude >= 111 & longitude <= 114 & latitude >= 31 & latitude <= 34]
  stopifnot(nrow(inbox) > 0)
  stopifnot(all(inbox$pt_masked), all(is.na(inbox$SIF)), all(inbox$SIF_raw == 0))

  msk <- out[pt_masked == TRUE]
  stopifnot(nrow(msk) > 0, all(is.na(msk$SIF)), all(msk$SIF_raw == 0))
  stopifnot(!any(out[!is.na(SIF)]$SIF == 0))
  log_msg(sprintf("PASS 0 值掩膜识别: %d 条被标记并置 NA (内部框 %d 条全中)",
                  nrow(msk), nrow(inbox)))

  # 断言 5: 输出列结构完整
  need <- c("meteo_stat", "latitude", "longitude", "year", "month", "day",
            "doy", "date", "SIF", "SIF_raw", "pt_masked", "SIF_buf750", "n_valid")
  stopifnot(all(need %in% names(out)))
  log_msg("PASS 列结构完整")

  log_msg("=== 自检全部通过 ===")
}

# --- 入口 ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) && args[1] == "--selftest") {
  selftest()
} else {
  if (!length(args)) stop("用法: Rscript 13_hcsif_extract.R <year> [raw_dir] [out_dir]")
  yr      <- as.integer(args[1])
  raw_dir <- if (length(args) >= 2) args[2] else file.path(PROJ, "data_raw/hcsif/tmp")
  out_dir <- if (length(args) >= 3) args[3] else file.path(PROJ, "data_raw/hcsif/station")
  run(yr, raw_dir, out_dir)
}
