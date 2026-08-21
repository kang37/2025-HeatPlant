#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 35_survival_deep.R — VPD->SIF 因果符号"翻转时间"的生存分析(深化版)
#
# 23_hazard_all.R 只做了离散时间风险模型, 且把 437 个"多次变号"站直接删掉。
# 本脚本在同一份数据上做四件事:
#
#   1. 前提检验     —— 生存分析该查什么、不该查什么(见文件末尾的说明块)
#   2. KM + log-rank —— 非参数地看两组翻转得快不快, 并给空间稳健版本
#   3. Cox 比例风险  —— 半参数, 带 18 协变量; PH 假设检验; Conley 空间 HAC
#   4. 竞争风险      —— 把"持续翻转"与"瞬时翻转"设为互斥事件, 437 站全部收回
#
# 关键的重新设定(相对 23 号脚本):
#   23 号: 事件 = 持续翻转; 多次变号站被"删除"。删除等价于假定这些站与保留站
#          同分布, 但它们恰恰是因果信号弱、S-map 系数在零附近抖动的站, 不随机。
#   本脚本: 时间 = 首次变号所在 tp(1..8), 全部 898 站都进模型;
#          status 0 = 到 tp8 从未变号(右删失, 即 always_inhibit/always_promote)
#          status 1 = 首次变号且此后不再变号(持续翻转, 与 23 号的事件完全一致)
#          status 2 = 首次变号但此后又变(瞬时翻转) —— 竞争事件
#   核对: status 分布 = 26/217/214 与 83/135/223, 与严格判据的分类逐个吻合。
#   "事件类型在事件发生后才知道"不是问题 —— 竞争风险里的死因也是事后判定的。
#
# 输出: data_proc/output_hcsif_buf1000/surv_*.csv 与 surv_base.rds
# 用法: Rscript 35_survival_deep.R [--cut 200]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(survival); library(splines)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
aa  <- commandArgs(trailingOnly = TRUE)
CUT <- if (any(aa == "--cut")) as.numeric(aa[which(aa == "--cut") + 1]) else 200
MAXTP <- 8L
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) log_msg("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74))

# ===========================================================================
# 0. 数据构建
# ===========================================================================
hdr("0. 数据构建")

env <- new.env(); lines <- readLines("22_varpart_3dim.R")
kcut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(
  eval(parse(text = paste(lines[1:kcut], collapse = "\n")), envir = env)))
dtc  <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)
log_msg("协变量 ", length(PRED), " 个: ", paste(PRED, collapse = ", "))

cm <- fread("data_proc/ccm_hcsif_buf1000/ccm_hcsif_buf1000_vpd_20260727_1031.csv")
cm <- cm[tp <= MAXTP]; setorder(cm, meteo_stat, tp)
N_TP <- length(unique(cm$tp))

pat <- cm[, .(co = list(mean_coef), rho_tp0 = rho[tp == 0][1],
              rho_max = max(rho, na.rm = TRUE)), by = .(stat_id = meteo_stat)]
n_all <- nrow(pat)
pat <- pat[sapply(co, length) == N_TP]
log_msg("站点 ", n_all, " -> 完整 ", N_TP, " 个滞后者 ", nrow(pat))

# 首次变号时间 + 事件类型
build_surv <- function(co) {
  s <- sign(co)
  if (anyNA(s) || any(s == 0)) return(list(t = NA_real_, st = NA_integer_, d0 = NA_character_))
  k  <- which(diff(s) != 0)                       # 变号发生在 tp=k 与 tp=k-1 之间
  d0 <- if (s[1] < 0) "inhibit_first" else "promote_first"
  if (!length(k)) return(list(t = as.numeric(MAXTP), st = 0L, d0 = d0))
  list(t = as.numeric(k[1]), st = if (length(k) == 1L) 1L else 2L, d0 = d0)
}
b <- lapply(pat$co, build_surv)
pat[, `:=`(time      = sapply(b, `[[`, "t"),
           status    = as.integer(sapply(b, `[[`, "st")),
           start_dir = sapply(b, `[[`, "d0"),
           n_flip    = sapply(co, function(x) sum(diff(sign(x)) != 0)))]
pat[, co := NULL]
pat <- pat[!is.na(status)]

# 与 23 号严格判据的对照: status==1 应恰为 inhibit_promote/promote_inhibit
pat[, stype := fcase(status == 0L & start_dir == "inhibit_first", "always_inhibit",
                     status == 0L & start_dir == "promote_first", "always_promote",
                     status == 1L & start_dir == "inhibit_first", "inhibit_promote",
                     status == 1L & start_dir == "promote_first", "promote_inhibit",
                     default = "multi_flip")]
log_msg("\n-- status(0=从未变号 1=持续翻转 2=瞬时翻转) x 起始方向 --")
print(dcast(pat[, .N, by = .(start_dir, status)], start_dir ~ status,
            value.var = "N", fill = 0L))
log_msg("\n-- 与 23 号严格判据分类的对照 --")
print(pat[, .N, by = stype][order(-N)])
log_msg("23 号可用样本 ", pat[status != 2L, .N], " 站; 本脚本 ", nrow(pat),
        " 站, 收回 ", pat[status == 2L, .N], " 站(", 
        round(100 * pat[status == 2L, .N] / nrow(pat), 1), "%)")

# 被剔除者是否随机? 用 CCM 技能比一比
log_msg("\n-- 被 23 号剔除者 vs 保留者的 CCM 技能(检验'剔除非随机') --")
cmp <- pat[, .(n = .N, rho_tp0 = round(mean(rho_tp0, na.rm = TRUE), 4),
               rho_max = round(mean(rho_max, na.rm = TRUE), 4)),
           by = .(dropped = status == 2L)]
print(cmp)
wt <- wilcox.test(rho_max ~ (status == 2L), data = pat)
log_msg("rho_max 两组 Wilcoxon 检验: W=", signif(wt$statistic, 5),
        ", p=", signif(wt$p.value, 3))

st <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
winz <- function(x, p = .01) { q <- quantile(x, c(p, 1-p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

D <- merge(merge(pat, dtc[, c("stat_id", "koppen_group", PRED), with = FALSE], by = "stat_id"),
           st[, .(stat_id, longitude, latitude)], by = "stat_id")
D <- D[complete.cases(D[, c("time", "status", PRED), with = FALSE])]
# 标准化: 系数读作"该变量增加一个标准差"
for (v in PRED) set(D, j = v, value = as.numeric(scale(winz(as.numeric(D[[v]])))))
D[, ev_any := as.integer(status > 0L)]
D[, status_f := factor(status, 0:2, c("censor", "sustained", "transient"))]
log_msg("\n建模样本 ", nrow(D), " 站(协变量完整)")
saveRDS(D, file.path(OUT, "surv_base.rds"))

# ---- 空间工具: 距离矩阵 / Conley 三明治 / Moran's I -----------------------
dist_km <- function(lon, lat) {
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  sqrt(dx^2 + (outer(lat, lat, "-") * 111.32)^2)
}
# coxph 的稳健方差 = t(dfbeta) %*% dfbeta。把中间换成空间核就是 Conley 版:
#   V = t(Dfb) %*% K %*% Dfb, K 为 Bartlett 核, 截断距离外为 0。
# 截断 -> 0 时 K -> I, 退化为 coxph(robust=TRUE), 故它是稳健方差的推广。
conley_cox <- function(fit, lon, lat, cut, id = NULL, key = NULL) {
  if (is.null(id)) {
    Dfb <- as.matrix(residuals(fit, type = "dfbeta"))
  } else {
    # collapse 后行序为 sort(unique(id)) 而非原始行序, 必须按 rownames 回对坐标,
    # 否则 Conley 用错了距离——这类错位不会报错, 只会悄悄给出错的标准误。
    Dfb <- as.matrix(residuals(fit, type = "dfbeta", collapse = id))
    ord <- match(rownames(Dfb), as.character(key))
    stopifnot(!anyNA(ord))
    lon <- lon[ord]; lat <- lat[ord]
  }
  DM <- dist_km(lon, lat)
  stopifnot(nrow(Dfb) == nrow(DM))
  K <- DM; K[] <- pmax(0, 1 - DM / cut)   # pmax 会丢掉 dim, 必须写回矩阵
  V <- t(Dfb) %*% K %*% Dfb
  list(se = sqrt(pmax(diag(V), 0)),
       eff_n = nrow(DM) / max(1, mean(rowSums(DM <= cut))))
}
moran_I <- function(x, lon, lat, cut, nperm = 999) {
  DM <- dist_km(lon, lat); W <- (DM > 0 & DM <= cut) * 1
  n <- length(x); z <- x - mean(x); S0 <- sum(W)
  I <- (n / S0) * sum(W * outer(z, z)) / sum(z^2)
  Ip <- replicate(nperm, { zz <- sample(z); (n / S0) * sum(W * outer(zz, zz)) / sum(zz^2) })
  list(I = I, p = (1 + sum(abs(Ip) >= abs(I))) / (nperm + 1))
}

# ===========================================================================
# 1. 前提检验
# ===========================================================================
hdr("1. 前提检验")

log_msg("-- 1.1 事件数与每自变量事件数(EPV) --")
epv <- D[, .(n = .N, ev_any = sum(status > 0), ev1 = sum(status == 1),
             ev2 = sum(status == 2), cens = sum(status == 0)), by = start_dir]
epv[, `:=`(epv_any = round(ev_any / length(PRED), 1), epv_1 = round(ev1 / length(PRED), 1))]
print(epv)
log_msg("经验下限 EPV>=10。低于 10 时系数偏倚、区间过窄, 须折扣解读。")

log_msg("\n-- 1.2 共线性(VIF, 站点层面) --")
vif <- sapply(PRED, function(v) {
  r2 <- summary(lm(as.formula(paste(v, "~", paste(setdiff(PRED, v), collapse = "+"))),
                   D))$r.squared
  1 / (1 - r2) })
print(round(sort(vif, decreasing = TRUE)[1:6], 2))
log_msg("最大 VIF = ", round(max(vif), 2), if (max(vif) < 10) "  [<10, 可解释]" else "  <<< 共线警告")

log_msg("\n-- 1.3 事件时间的并结(ties) --")
tt <- D[status > 0, .N, by = time][order(time)]
print(tt)
log_msg("时间只有 ", nrow(tt), " 个离散取值, 并结极重(每个时点数十至上百站)。")
log_msg("=> Cox 的精确离散偏似然不可行(组合爆炸), 用 Efron 近似;")
log_msg("   并以离散时间 cloglog 模型作为分组时间下的精确对应模型交叉核对。")

log_msg("\n-- 1.4 函数形式: 线性 vs 自然样条(df=3), 逐变量 LRT --")
f_lin <- as.formula(paste("Surv(time, ev_any) ~", paste(PRED, collapse = "+")))
m_lin <- coxph(f_lin, data = D, ties = "efron")
# 五类森林占比在多数站为 0(零膨胀), 内节点与边界重合, 样条无法构造。
# 这类变量本身近似"有/无"的二值信息, 线性设定即已足够, 单独标出而非硬拟合。
spl_ok <- sapply(PRED, function(v) {
  q <- quantile(D[[v]], c(.25, .5, .75)); length(unique(q)) == 3L })
log_msg("零膨胀、无法构造样条而跳过的变量: ",
        if (any(!spl_ok)) paste(PRED[!spl_ok], collapse = ", ") else "无")
ff <- rbindlist(lapply(PRED[spl_ok], function(v) {
  f2 <- as.formula(paste("Surv(time, ev_any) ~",
        paste(c(setdiff(PRED, v), sprintf("ns(%s,3)", v)), collapse = "+")))
  m2 <- try(coxph(f2, data = D, ties = "efron"), silent = TRUE)
  if (inherits(m2, "try-error")) return(NULL)
  data.table(var = v, lrt = 2 * (m2$loglik[2] - m_lin$loglik[2]), df = 2)
}))
ff[, p := pchisq(lrt, df, lower.tail = FALSE)]
ff[, p_bonf := pmin(1, p * .N)]
print(ff[order(p)][1:6, .(var, lrt = round(lrt, 2), p = signif(p, 3), p_bonf = signif(p_bonf, 3))])
log_msg("Bonferroni 后仍显著者需改用样条: ", 
        if (ff[p_bonf < .05, .N]) paste(ff[p_bonf < .05, var], collapse = ", ") else "无")

log_msg("\n-- 1.5 影响点(标准化 dfbeta) --")
dfb <- residuals(m_lin, type = "dfbeta")
sdf <- sweep(abs(dfb), 2, sqrt(diag(vcov(m_lin))), "/")
log_msg("最大标准化 dfbeta = ", round(max(sdf), 3),
        " (经验阈值 0.2); 超阈站点数 ", sum(apply(sdf, 1, max) > 0.2))

log_msg("\n-- 1.6 残差空间自相关(Moran's I, ", CUT, " km 邻域) --")
mr <- residuals(m_lin, type = "martingale")
mi <- moran_I(mr, D$longitude, D$latitude, CUT)
log_msg("Moran's I = ", round(mi$I, 4), ", 置换 p = ", signif(mi$p, 3),
        if (mi$p < .05) "  <<< 存在空间自相关, 标准误必须用 Conley" else "  [不显著]")

log_msg("\n-- 1.7 不需要检验的项 --")
log_msg("正态性: KM/log-rank/Cox 都不对结局或残差作正态假定(Cox 的基线风险完全非参数),")
log_msg("        故不做 Shapiro-Wilk 之类的检验; 正态性只在推断依赖大样本正态近似时",
        "以'样本够不够大'的形式出现。")
log_msg("等方差: 同理不适用。异方差在生存模型里表现为 PH 假设不成立(见 3.2)。")
log_msg("独立性: 需要, 且本数据违反(见 1.6) —— 用 Conley 空间 HAC 修正方差。")
log_msg("非信息删失: 无法从数据检验。此处删失 = 到 tp8 仍未变号, 由观测窗口决定,")
log_msg("        与站点特征无关, 属设计性删失, 假定成立较安全。")

# ===========================================================================
# 2. Kaplan-Meier 与 log-rank
# ===========================================================================
hdr("2. Kaplan-Meier 与 log-rank(事件 = 首次变号, 不分持续/瞬时)")

km <- survfit(Surv(time, ev_any) ~ start_dir, data = D)
log_msg("-- 生存函数 S(t) = 到第 t 步仍未变号的比例 --")
ktab <- data.table(time = km$time, n_risk = km$n.risk, n_event = km$n.event,
                   surv = round(km$surv, 4), lo = round(km$lower, 4), hi = round(km$upper, 4),
                   grp = rep(names(km$strata), km$strata))
print(ktab)
fwrite(ktab, file.path(OUT, "surv_km_bystartdir.csv"))
mq <- quantile(km, probs = c(.25, .5))
log_msg("\n中位翻转时间(步, 每步8天):")
print(round(mq$quantile, 2))

sd0 <- survdiff(Surv(time, ev_any) ~ start_dir, data = D, rho = 0)
sd1 <- survdiff(Surv(time, ev_any) ~ start_dir, data = D, rho = 1)
log_msg("\nlog-rank (rho=0, 各时点等权):  chisq=", round(sd0$chisq, 3),
        ", df=", length(sd0$n) - 1, ", p=", signif(pchisq(sd0$chisq, length(sd0$n)-1, lower.tail=FALSE), 3))
log_msg("Peto-Wilcoxon (rho=1, 早期加权): chisq=", round(sd1$chisq, 3),
        ", p=", signif(pchisq(sd1$chisq, length(sd1$n)-1, lower.tail=FALSE), 3))

# log-rank = Cox 中二值协变量的得分检验。故把方差换成 Conley,
# 就得到"空间稳健的 log-rank"——观测独立性被违反时这才是可信的 p 值。
m_g <- coxph(Surv(time, ev_any) ~ start_dir, data = D, ties = "efron")
cg  <- conley_cox(m_g, D$longitude, D$latitude, CUT)
z_naive <- coef(m_g) / sqrt(diag(vcov(m_g)))
z_conley <- coef(m_g) / cg$se
log_msg("\n空间稳健版: beta=", round(coef(m_g), 4),
        " | 普通 SE ", round(sqrt(diag(vcov(m_g))), 4), " p=", signif(2*pnorm(-abs(z_naive)), 3),
        " | Conley", CUT, "km SE ", round(cg$se, 4), " p=", signif(2*pnorm(-abs(z_conley)), 3))
log_msg("有效独立样本量 ~", round(cg$eff_n, 1), " (名义 ", nrow(D), ")")

log_msg("\n-- 分气候区的 log-rank(组内) --")
kres <- rbindlist(lapply(c("inhibit_first", "promote_first"), function(g) {
  d <- D[start_dir == g & !is.na(koppen_group)]
  s <- survdiff(Surv(time, ev_any) ~ koppen_group, data = d)
  data.table(grp = g, n = nrow(d), k = length(s$n), chisq = s$chisq,
             p = pchisq(s$chisq, length(s$n) - 1, lower.tail = FALSE))
}))
print(kres[, .(grp, n, k, chisq = round(chisq, 3), p = signif(p, 3))])
fwrite(kres, file.path(OUT, "surv_logrank_koppen.csv"))

# ===========================================================================
# 3. Cox 比例风险模型
# ===========================================================================
hdr("3. Cox 比例风险模型(18 协变量, 事件 = 首次变号)")

cox_table <- function(fit, d, tag, id = NULL) {
  b <- coef(fit); vn <- names(b)
  se_n <- sqrt(diag(vcov(fit)))
  cs <- conley_cox(fit, d$longitude, d$latitude, CUT, id = id, key = d$stat_id)
  se_c <- cs$se; names(se_c) <- vn
  out <- data.table(tag = tag, var = vn, beta = b, HR = exp(b),
                    se_naive = se_n, se_conley = se_c)
  out[, `:=`(zstat = beta / se_conley)]
  out[, `:=`(p_conley = 2 * pnorm(-abs(zstat)), p_naive = 2 * pnorm(-abs(beta / se_naive)),
             lo = exp(beta - 1.96 * se_conley), hi = exp(beta + 1.96 * se_conley),
             eff_n = cs$eff_n)]
  out[]
}

# strata(start_dir): 两组允许有各自的基线风险, 只共享协变量效应。
# 这比放一个 start_dir 主效应更弱的假定 —— 无需假定两组风险成比例。
f_full <- as.formula(paste("Surv(time, ev_any) ~ strata(start_dir) +",
                          paste(PRED, collapse = "+")))
m_full <- coxph(f_full, data = D, ties = "efron")
t_full <- cox_table(m_full, D, "cox_any_strata")
log_msg("n = ", m_full$n, ", 事件 ", m_full$nevent, ", 有效独立样本 ~",
        round(t_full$eff_n[1], 1))
log_msg("HR>1 = 该变量越大, 每步变号风险越高 -> 变号越早")
print(t_full[order(p_conley)][1:10, .(var, HR = round(HR, 3),
      CI = sprintf("[%.2f,%.2f]", lo, hi), p_conley = signif(p_conley, 3),
      p_naive = signif(p_naive, 3), sig = fifelse(p_conley < .05, "*", ""))])
log_msg("Conley 下显著 ", t_full[p_conley < .05, .N], " 个 | 不修正会有 ",
        t_full[p_naive < .05, .N], " 个")

log_msg("\n-- 3.2 比例风险假设检验(Schoenfeld 残差) --")
zp <- cox.zph(m_full)
zt <- data.table(var = rownames(zp$table), chisq = zp$table[, "chisq"],
                 df = zp$table[, "df"], p = zp$table[, "p"])
print(zt[order(p)][1:8, .(var, chisq = round(chisq, 2), df, p = signif(p, 3),
                          sig = fifelse(p < .05, "违反", ""))])
log_msg("全局检验 p = ", signif(zt[var == "GLOBAL", p], 3),
        if (zt[var == "GLOBAL", p] < .05) "  <<< PH 假设整体被拒" else "  [不拒绝]")
fwrite(zt, file.path(OUT, "surv_cox_phtest.csv"))

viol <- zt[var != "GLOBAL" & p < .05, var]
if (length(viol)) {
  log_msg("违反 PH 的变量: ", paste(viol, collapse = ", "))
  log_msg("处理: 对全部违反变量同时加时间交互项 tt(x) = x * log(t)。")
  log_msg("      这不只是'补丁'——交互项本身就是结论: 系数为正 = 该变量的作用")
  log_msg("      随滞后增大而增强, 为负 = 只在短滞后起作用、之后衰减。")
  f_tt <- as.formula(paste("Surv(time, ev_any) ~ strata(start_dir) +",
            paste(PRED, collapse = "+"), "+",
            paste(sprintf("tt(%s)", viol), collapse = "+")))
  m_tt <- coxph(f_tt, data = D, ties = "efron",
                tt = function(x, t, ...) x * log(t))
  ctt <- summary(m_tt)$coefficients
  k <- grep("^tt\\(", rownames(ctt))
  tv <- data.table(var = sub("^tt\\((.*)\\)$", "\\1", rownames(ctt)[k]),
                   beta_tt = ctt[k, "coef"], p_tt = ctt[k, "Pr(>|z|)"])
  tv[, beta_main := ctt[match(var, rownames(ctt)), "coef"]]
  tv[, 走向 := fifelse(beta_tt > 0, "随滞后增强", "随滞后衰减")]
  print(tv[order(p_tt), .(var, beta_main = round(beta_main, 3),
        beta_tt = round(beta_tt, 3), p_tt = signif(p_tt, 3),
        走向 = fifelse(p_tt < .05, 走向, ""))])
  log_msg("LRT(加时间交互 vs 不加): chisq=",
          round(2 * (m_tt$loglik[2] - m_full$loglik[2]), 2), ", df=", length(k),
          ", p=", signif(pchisq(2 * (m_tt$loglik[2] - m_full$loglik[2]),
                                length(k), lower.tail = FALSE), 3))
  fwrite(tv, file.path(OUT, "surv_cox_timevarying.csv"))
}

# 1.4 查出 elev / rsds_mean 线性设定不足, 这里回代样条作稳健性检查:
# 若主结论(系数符号与显著性)在样条设定下不变, 说明非线性不影响解读。
nl <- ff[p_bonf < .05, var]
if (length(nl)) {
  f_ns <- as.formula(paste("Surv(time, ev_any) ~ strata(start_dir) +",
      paste(c(setdiff(PRED, nl), sprintf("ns(%s,3)", nl)), collapse = "+")))
  m_ns <- coxph(f_ns, data = D, ties = "efron")
  keep <- setdiff(PRED, nl)
  cmpn <- data.table(var = keep, lin = coef(m_full)[keep], spl = coef(m_ns)[keep])
  log_msg("\n-- 3.2b 非线性稳健性: ", paste(nl, collapse = ", "), " 改用自然样条 --")
  log_msg("其余 ", length(keep), " 个系数与线性设定的相关 r = ",
          round(cor(cmpn$lin, cmpn$spl), 4), "; 最大绝对变动 ",
          round(max(abs(cmpn$lin - cmpn$spl)), 4))
  # 符号翻转只在系数接近零处才值得担心, 故与显著性一并看
  cmpn <- merge(cmpn, t_full[, .(var, p_lin = p_conley)], by = "var")
  cs2 <- conley_cox(m_ns, D$longitude, D$latitude, CUT)
  cmpn[, p_spl := 2 * pnorm(-abs(spl / cs2$se[match(var, names(coef(m_ns)))]))]
  cmpn[, 翻转 := fifelse(sign(lin) != sign(spl), "是", "")]
  print(cmpn[order(p_lin)][1:8, .(var, lin = round(lin, 3), spl = round(spl, 3),
        p_lin = signif(p_lin, 3), p_spl = signif(p_spl, 3), 翻转)])
  fl <- cmpn[翻转 == "是"]
  log_msg("符号翻转 ", nrow(fl), " 个: ",
          if (nrow(fl)) paste(sprintf("%s(|beta|<=%.3f, 两侧 p 均>%.2f)", fl$var,
              pmax(abs(fl$lin), abs(fl$spl)), pmin(fl$p_lin, fl$p_spl)),
              collapse = ", ") else "无")
  log_msg("在线性或样条设定下任一显著(p<.05)且符号翻转的变量: ",
          if (cmpn[翻转 == "是" & (p_lin < .05 | p_spl < .05), .N])
            paste(cmpn[翻转 == "是" & (p_lin < .05 | p_spl < .05), var], collapse = ", ")
          else "无  [主结论对非线性设定稳健]")
}

# 分组时间下的精确对应: 离散时间 cloglog 模型(与 Cox 同一 PH 参数)
log_msg("\n-- 3.3 与离散时间模型的交叉核对 --")
pp <- D[, .(tp = seq_len(time)), by = .(stat_id, start_dir, time, ev_any)]
pp[, ev := as.integer(tp == time & ev_any == 1)]
pp <- merge(pp, D[, c("stat_id", PRED), with = FALSE], by = "stat_id")
f_d <- as.formula(paste("ev ~ factor(tp) + start_dir +", paste(PRED, collapse = "+")))
g_cll <- glm(f_d, data = pp, family = binomial(link = "cloglog"))
g_lgt <- glm(f_d, data = pp, family = binomial(link = "logit"))
cmp3 <- data.table(var = PRED, cox = coef(m_full)[PRED],
                   cloglog = coef(g_cll)[PRED], logit = coef(g_lgt)[PRED])
log_msg("Cox(Efron) 与离散 cloglog 系数相关 r = ",
        round(cor(cmp3$cox, cmp3$cloglog), 4),
        "; 与 logit r = ", round(cor(cmp3$cox, cmp3$logit), 4))
print(cmp3[order(-abs(cox))][1:6, .(var, cox = round(cox, 3),
      cloglog = round(cloglog, 3), logit = round(logit, 3))])
fwrite(cmp3, file.path(OUT, "surv_cox_vs_discrete.csv"))

# 分组各自拟合(便于与 23 号逐组结果对照)
res_grp <- rbindlist(lapply(c("inhibit_first", "promote_first"), function(g) {
  d <- D[start_dir == g]
  m <- coxph(as.formula(paste("Surv(time, ev_any) ~", paste(PRED, collapse = "+"))),
             data = d, ties = "efron")
  zph <- cox.zph(m)
  o <- cox_table(m, d, paste0("cox_any_", g))
  set(o, j = "ph_p", value = zph$table[match(o$var, rownames(zph$table)), "p"])
  log_msg("\n### ", g, " | n=", m$n, " 事件=", m$nevent,
          " | PH 全局 p=", signif(zph$table["GLOBAL", "p"], 3))
  print(o[order(p_conley)][1:6, .(var, HR = round(HR, 3),
        CI = sprintf("[%.2f,%.2f]", lo, hi), p = signif(p_conley, 3),
        sig = fifelse(p_conley < .05, "*", ""))])
  o
}))
fwrite(rbind(t_full, res_grp, fill = TRUE), file.path(OUT, "surv_cox_coef.csv"))

# ===========================================================================
# 4. 竞争风险
# ===========================================================================
hdr("4. 竞争风险: 持续翻转(事件1) vs 瞬时翻转(事件2)")

log_msg("为什么需要它: 23 号把事件2的 437 站删掉。删除隐含'它们与保留站同分布',")
log_msg("而 1.0 节已显示两者 CCM 技能不同 -> 删除有偏。竞争风险把事件2 保留在")
log_msg("风险集里, 只是记为另一种结局, 不作任何同分布假定。")

aj <- survfit(Surv(time, status_f) ~ start_dir, data = D)
cif <- data.table(time = rep(aj$time, 2),
                  state = rep(colnames(aj$pstate)[2:3], each = length(aj$time)),
                  cif = c(aj$pstate[, 2], aj$pstate[, 3]),
                  grp = rep(rep(names(aj$strata), aj$strata), 2))
fwrite(cif, file.path(OUT, "surv_cif.csv"))
log_msg("\n-- Aalen-Johansen 累积发生率 CIF(tp=8 处) --")
print(cif[time == MAXTP][order(grp, state), .(grp, state, cif = round(cif, 4))])
log_msg("注意: 用 1-KM 会高估 —— 把竞争事件当删失, 等于假定它们'以后还会持续翻转'。")

kmc <- survfit(Surv(time, status == 1L) ~ start_dir, data = D)
log_msg("对照: 1-KM(把事件2 当删失) 在 tp=8 处 = ",
        paste(round(1 - summary(kmc, times = MAXTP)$surv, 4), collapse = " / "),
        "  (高估量见与上表之差)")

# Fine-Gray 子分布风险模型。Gray 检验就是子分布风险的组间检验,
# 这里用 Fine-Gray 回归中组指示变量的 Wald 检验实现同一件事。
fg1 <- finegray(Surv(time, status_f) ~ ., data = D, etype = "sustained")
m_gray <- coxph(Surv(fgstart, fgstop, fgstatus) ~ start_dir,
                weights = fgwt, data = fg1, robust = TRUE, id = stat_id)
log_msg("\n-- Gray 型检验(事件1 的子分布风险, 两组) --")
log_msg("beta=", round(coef(m_gray), 4), ", 稳健 p=",
        signif(summary(m_gray)$coefficients[, "Pr(>|z|)"], 3))

# 同一批协变量, 两种风险各估一遍:
#   cause-specific  = 问"在还没变号的站里, 谁更快发生持续翻转"(病因机制)
#   subdistribution = 问"这个变量如何改变最终持续翻转的累积概率"(预测/占比)
# 二者方向不一致处, 说明该变量的作用有相当一部分是通过竞争事件实现的。
fg_res <- rbindlist(lapply(c("inhibit_first", "promote_first"), function(g) {
  d <- D[start_dir == g]
  m_cs1 <- coxph(as.formula(paste("Surv(time, status == 1L) ~", paste(PRED, collapse = "+"))),
                 data = d, ties = "efron")
  m_cs2 <- coxph(as.formula(paste("Surv(time, status == 2L) ~", paste(PRED, collapse = "+"))),
                 data = d, ties = "efron")
  fg <- finegray(Surv(time, status_f) ~ ., data = d, etype = "sustained")
  m_fg <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                                 paste(PRED, collapse = "+"))),
                weights = fgwt, data = fg, robust = TRUE, id = stat_id)
  # Conley: dfbeta 按站折叠后再加空间核
  cs_fg <- conley_cox(m_fg, d$longitude, d$latitude, CUT,
                      id = fg$stat_id, key = d$stat_id)
  t1 <- cox_table(m_cs1, d, paste0("cause_specific_ev1_", g))
  t2 <- cox_table(m_cs2, d, paste0("cause_specific_ev2_", g))
  tf <- data.table(tag = paste0("fine_gray_ev1_", g), var = names(coef(m_fg)),
                   beta = coef(m_fg), HR = exp(coef(m_fg)),
                   se_naive = sqrt(diag(m_fg$naive.var)), se_conley = cs_fg$se)
  tf[, zstat := beta / se_conley][, `:=`(p_conley = 2 * pnorm(-abs(zstat)),
      p_naive = 2 * pnorm(-abs(beta / se_naive)),
      lo = exp(beta - 1.96 * se_conley), hi = exp(beta + 1.96 * se_conley))]

  log_msg("\n### ", g, " | 事件1 ", m_cs1$nevent, " 例, 事件2 ", m_cs2$nevent, " 例")
  mg <- merge(t1[, .(var, cs1 = beta, p1 = p_conley)],
              tf[, .(var, fg = beta, pf = p_conley)], by = "var")
  mg <- merge(mg, t2[, .(var, cs2 = beta, p2 = p_conley)], by = "var")
  mg[, flip := sign(cs1) != sign(fg)]
  print(mg[order(p1)][1:8, .(var, cs_ev1 = round(cs1, 3), p1 = signif(p1, 3),
        FG_ev1 = round(fg, 3), pFG = signif(pf, 3), cs_ev2 = round(cs2, 3),
        symbol_flip = fifelse(flip, "方向相反", ""))])
  log_msg("cause-specific 与 Fine-Gray 方向相反的变量: ",
          if (mg[flip == TRUE, .N]) paste(mg[flip == TRUE, var], collapse = ", ") else "无")
  rbind(t1, t2, tf, fill = TRUE)
}))
fwrite(fg_res, file.path(OUT, "surv_competing_coef.csv"))

# ===========================================================================
# 5. 分气候区
# ===========================================================================
hdr("5. 分气候区")

# A 区只有 2 站, 进模型只会制造不收敛的交互项, 直接剔除
kn <- D[!is.na(koppen_group), .N, by = koppen_group][N >= 20, koppen_group]
Dk <- D[koppen_group %in% kn]
log_msg("保留站数 >= 20 的气候区: ", paste(sort(kn), collapse = ", "),
        " (剔除 ", D[!is.na(koppen_group) & !koppen_group %in% kn, .N], " 站)")
log_msg("有柯本分区的站 ", nrow(Dk), " / ", nrow(D),
        " (", round(100 * nrow(Dk) / nrow(D), 1), "%)")
print(dcast(Dk[, .N, by = .(koppen_group, start_dir)], koppen_group ~ start_dir,
            value.var = "N", fill = 0L))

# 5.1 分层基线: 各区自有基线风险, 协变量效应仍共享(自由度代价小)
m_ks <- coxph(as.formula(paste("Surv(time, ev_any) ~ strata(start_dir, koppen_group) +",
                               paste(PRED, collapse = "+"))), data = Dk, ties = "efron")
t_ks <- cox_table(m_ks, Dk, "cox_strata_koppen")
log_msg("\n-- 5.1 按气候区分层的基线风险(协变量效应共享) --")
print(t_ks[order(p_conley)][1:8, .(var, HR = round(HR, 3),
      CI = sprintf("[%.2f,%.2f]", lo, hi), p = signif(p_conley, 3),
      sig = fifelse(p_conley < .05, "*", ""))])

# 5.2 效应是否随气候区变化: 对每个变量做 koppen 交互的 LRT
log_msg("\n-- 5.2 协变量效应的气候区异质性(交互 LRT) --")
base_ll <- m_ks$loglik[2]
het <- rbindlist(lapply(PRED, function(v) {
  f2 <- as.formula(paste("Surv(time, ev_any) ~ strata(start_dir, koppen_group) +",
        paste(PRED, collapse = "+"), "+", v, ":koppen_group"))
  m2 <- try(coxph(f2, data = Dk, ties = "efron"), silent = TRUE)
  if (inherits(m2, "try-error")) return(NULL)
  dfd <- length(coef(m2)) - length(coef(m_ks))
  data.table(var = v, lrt = 2 * (m2$loglik[2] - base_ll), df = dfd)
}))
het[, p := pchisq(lrt, df, lower.tail = FALSE)][, p_bonf := pmin(1, p * .N)]
print(het[order(p)][1:6, .(var, lrt = round(lrt, 2), df, p = signif(p, 3),
                           p_bonf = signif(p_bonf, 3))])
log_msg("Bonferroni 后有异质性的变量: ",
        if (het[p_bonf < .05, .N]) paste(het[p_bonf < .05, var], collapse = ", ") else "无")
fwrite(het, file.path(OUT, "surv_koppen_heterogeneity.csv"))

# 5.3 逐格单独拟合(精简变量集), 并报每格的有效独立样本量
PRED_S <- c("imperv", "grass", "elev", "rsds_sd", "ntl", "urban_rate")
log_msg("\n-- 5.3 逐格单独拟合 | 精简集: ", paste(PRED_S, collapse = ", "), " --")
zres <- list()
for (kg in sort(unique(Dk$koppen_group))) for (g in c("inhibit_first", "promote_first")) {
  d <- Dk[koppen_group == kg & start_dir == g]
  n_ev <- sum(d$ev_any); lab <- sprintf("%s | %s", kg, g)
  if (n_ev < 5 * length(PRED_S)) {
    log_msg(sprintf("[跳过] %s: 站 %d, 事件 %d (每自变量 %.1f < 5)",
                    lab, nrow(d), n_ev, n_ev / length(PRED_S))); next }
  m <- coxph(as.formula(paste("Surv(time, ev_any) ~", paste(PRED_S, collapse = "+"))),
             data = d, ties = "efron")
  o <- cox_table(m, d, lab); en <- o$eff_n[1]
  o[, `:=`(kg = kg, grp = g, n = nrow(d), n_ev = n_ev)]
  log_msg(sprintf("\n=== %s | 站 %d | 事件 %d | 每自变量 %.1f | 有效独立样本 ~%.1f%s ===",
                  lab, nrow(d), n_ev, n_ev / length(PRED_S), en,
                  if (en < 10) "  <<< Conley 在该格已失效, 结果不可用" else ""))
  print(o[order(p_conley), .(var, HR = round(HR, 3),
        CI = sprintf("[%.2f,%.2f]", lo, hi), p = signif(p_conley, 3),
        sig = fifelse(p_conley < .05, "*", ""))])
  zres[[lab]] <- o
}
if (length(zres)) fwrite(rbindlist(zres), file.path(OUT, "surv_koppen_cox.csv"))

hdr("完成")
log_msg("输出目录: ", OUT)
