#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 17_glc_merge.R — 合并 46 个瓦片分片，算 1000 m 缓冲区的 GLC_FCS30D 类别占比
#
# 分片是 (stat_id, cell, year, cls, n) 的像元计数。跨瓦片边界的缓冲区网格会在
# 两个瓦片各得一部分像元，故按 (stat_id, cell, year, cls) 求和即自动按像元数加权。
#
# 输出(与 15_covariates_1km.R / 16_ntl_1km.R 同构):
#   glc_cell_<年>.csv     每个 HCSIF 500 m 像元一行, 各类占比
#   glc_station_<年>.csv  每个站点一行, 缓冲区内各类占比
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_raw/covariates_1km"
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# GLC_FCS30D 精细分类体系, 29 类。顺序即 WS 文件里 LC00-LC28 的顺序。
LC <- data.table(
  cls  = c(10L, 11L, 12L, 20L, 51L, 52L, 61L, 62L, 71L, 72L, 81L, 82L, 91L, 92L,
           120L, 121L, 122L, 130L, 140L, 150L, 152L, 153L, 180L, 190L,
           200L, 201L, 202L, 210L, 220L),
  code = c("LC00","LC01","LC02","LC03","LC04","LC05","LC06","LC07","LC08","LC09",
           "LC10","LC11","LC12","LC13","LC14","LC15","LC16","LC17","LC18","LC19",
           "LC20","LC21","LC22","LC23","LC24","LC25","LC26","LC27","LC28"))

fs <- list.files("data_raw/glc/parts", "^glc_long_.*[.]csv$", full.names = TRUE)
log_msg("分片 ", length(fs), " 个")
d <- rbindlist(lapply(fs, fread))
log_msg("原始 ", nrow(d), " 行, ", uniqueN(d$stat_id), " 站, ",
        uniqueN(d$year), " 年, ", uniqueN(d$cls), " 个类别码")

# GLC_FCS30D 把湿地(180)细分为 181 沼泽 / 182 草本沼泽 / 183 泛洪地 / 185 红树林 /
# 186 盐沼 / 187 潮滩 等亚类。WS 文件只有一个 LC22_湿地，故先并回 180 再统计。
WETLAND <- 181:187
if (d[cls %in% WETLAND, .N]) {
  log_msg("湿地亚类 ", paste(sort(intersect(unique(d$cls), WETLAND)), collapse = ","),
          " 合并为 180 (", d[cls %in% WETLAND, sum(n)], " 像元)")
  d[cls %in% WETLAND, cls := 180L]
}

unk <- setdiff(unique(d$cls), LC$cls)
if (length(unk)) log_msg("注意: 类表外的编码 ", paste(sort(unk), collapse = ","),
                         " (合计 ", d[cls %in% unk, sum(n)], " 像元), 已计入分母但不单列")

# 跨瓦片边界的同一网格会重复出现，求和合并
d <- d[, .(n = sum(n)), by = .(stat_id, cell, year, cls)]

for (lvl in c("cell", "station")) {
  key <- if (lvl == "cell") c("stat_id", "cell", "year") else c("stat_id", "year")
  a   <- d[, .(n = sum(n)), by = c(key, "cls")]
  a[, frac := n / sum(n), by = key]                 # 分母含类表外编码，占比之和仍为 1
  a <- merge(a, LC, by = "cls")
  w <- dcast(a, as.formula(paste(paste(key, collapse = "+"), "~ code")),
             value.var = "frac", fill = 0)
  for (cc in setdiff(LC$code, names(w))) w[, (cc) := 0]   # 全域不出现的类别补 0
  setcolorder(w, c(key, LC$code))
  w[, proportion_sum := rowSums(.SD), .SDcols = LC$code]
  for (y in sort(unique(w$year)))
    fwrite(w[year == y], file.path(OUT, sprintf("glc_%s_%d.csv", lvl, y)))
  log_msg(lvl, " 级写出 ", uniqueN(w$year), " 年, 占比和 [",
          paste(round(range(w$proportion_sum), 6), collapse = ", "), "]")
}
log_msg("完成")
