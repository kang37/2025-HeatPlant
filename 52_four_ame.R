#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 52_four_ame.R — 把两层分解的系数翻译成"每个变量把站推向哪一类"
#
# P(始终抑制)   = P(抑制) * (1 - P(翻|抑制))
# P(抑制转促进) = P(抑制) *      P(翻|抑制)
# P(始终促进)   = P(促进) * (1 - P(翻|促进))
# P(促进转抑制) = P(促进) *      P(翻|促进)
#
# 第 2 层 a(抑制组内翻转)已判定为准完全分离且折外 AUC 0.516, 其系数不可用。
# 故给两个版本:
#   完整版  三层全用       —— 仅供对照, 抑制支的拆分不可信
#   保守版  第2层a 固定为基础发生率 90.7% —— 可发表版本
# 在保守版里, 变量对"始终抑制 vs 抑制转促进"的相对推力恒为 0,
# 也就是说: 协变量只能告诉你站落在抑制支还是促进支, 不能告诉你抑制支内部落在哪一类。
#
# 报的是平均边际效应(AME): 每个自变量 +1 个标准差, 四类概率各变多少个百分点。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"; set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

obj <- readRDS(file.path(OUT, "class_model_data.rds"))
PRED <- obj$PRED
D <- obj$D[strict == TRUE]
D[, flip := as.integer(stype %in% c("inhibit_promote", "promote_inhibit"))]
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }
z <- copy(D); for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
z[, wt := winz(w, .01)][, wt := wt/mean(wt)][, ww := wt]
CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔", ntl="夜间灯光",
        invest="绿地投资", urban_rate="城镇化率")
f <- function(y) as.formula(paste(y, "~", paste(PRED, collapse = "+")))
m1  <- suppressWarnings(glm(f("grp2"), binomial(), z, weights = ww))
z[, ww := 1]
m2i <- suppressWarnings(glm(f("flip"), binomial(), z[grp2 == 1]))
m2p <- suppressWarnings(glm(f("flip"), binomial(), z[grp2 == 0]))
BASE_I <- mean(z[grp2 == 1]$flip)          # 抑制支的翻转基础发生率

probs <- function(d, conservative = TRUE) {
  pi_ <- predict(m1, d, type = "response")
  qi  <- if (conservative) rep(BASE_I, nrow(d)) else predict(m2i, d, type = "response")
  qp  <- predict(m2p, d, type = "response")
  cbind(始终抑制 = pi_*(1-qi), 抑制转促进 = pi_*qi,
        始终促进 = (1-pi_)*(1-qp), 促进转抑制 = (1-pi_)*qp)
}
ame <- function(conservative) {
  P0 <- probs(z, conservative)
  rbindlist(lapply(PRED, function(v) {
    d <- copy(z); set(d, j = v, value = d[[v]] + 1)      # +1 个标准差
    dd <- colMeans(probs(d, conservative) - P0) * 100
    data.table(var = v, 变量 = CN[v], t(dd))
  }))
}
co <- fread(file.path(OUT, "class4_layer_coef.csv"))[var != "(Intercept)"]
p1 <- co[layer %like% "第1层", .(var, p_dir = p)]
p2 <- co[layer %like% "第2层b", .(var, p_flip_p = p)]

A <- ame(TRUE)
A <- merge(merge(A, p1, by = "var"), p2, by = "var")
A[, 总推力 := abs(始终抑制) + abs(抑制转促进) + abs(始终促进) + abs(促进转抑制)]
setorder(A, -总推力)
log_msg("=== 保守版: 每个自变量 +1 SD, 四类概率各变几个百分点 ===")
log_msg("(第2层a 固定为基础发生率 ", round(100*BASE_I, 1), "%, 因其系数不可用)\n")
sig <- function(p) fifelse(p < .001, "***", fifelse(p < .01, "**", fifelse(p < .05, "*", "")))
print(A[, .(变量, 始终抑制 = sprintf("%+.1f", 始终抑制),
            抑制转促进 = sprintf("%+.1f", 抑制转促进),
            始终促进 = sprintf("%+.1f", 始终促进),
            促进转抑制 = sprintf("%+.1f", 促进转抑制),
            方向层p = paste0(signif(p_dir, 2), sig(p_dir)),
            翻转层p = paste0(signif(p_flip_p, 2), sig(p_flip_p)))])
fwrite(A, file.path(OUT, "class4_ame_conservative.csv"))

B <- ame(FALSE); B[, 总推力 := abs(始终抑制)+abs(抑制转促进)+abs(始终促进)+abs(促进转抑制)]
setorder(B, -总推力)
fwrite(B, file.path(OUT, "class4_ame_full.csv"))
log_msg("\n=== 完整版(三层全用)与保守版的差别: 仅供对照 ===")
cmp <- merge(A[, .(变量, 保守_始终抑制 = 始终抑制)], B[, .(变量, 完整_始终抑制 = 始终抑制)], by = "变量")
cmp <- cmp[order(-abs(完整_始终抑制))]
print(head(cmp[, .(变量, 保守_始终抑制 = sprintf("%+.1f", 保守_始终抑制),
                   完整_始终抑制 = sprintf("%+.1f", 完整_始终抑制))], 6))
log_msg("\n四类基础占比: ", paste(sprintf("%s %.1f%%", c("始终抑制","抑制转促进","始终促进","促进转抑制"),
        100*colMeans(probs(z, TRUE))), collapse = " | "))
