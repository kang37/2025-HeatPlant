#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 50_tables.R — 把 45/46/47/49 的结果汇成可直接贴进论文的表
# 输出: data_proc/output_hcsif_buf1000/class_results_tables.md
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
OUT <- "data_proc/output_hcsif_buf1000"
CN <- c(imperv="不透水面", grass="草地", water="水体", crop="耕地",
        ever_needle="常绿针叶林", deci_needle="落叶针叶林",
        ever_broad="常绿阔叶林", deci_broad="落叶阔叶林", mixedleaf="混交林",
        tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率",
        `(Intercept)`="截距")
DIMOF <- c(rep("土地利用",9), rep("气候环境",7), rep("社会经济",3))
names(DIMOF) <- c("imperv","grass","water","crop","ever_needle","deci_needle",
  "ever_broad","deci_broad","mixedleaf","tavg","rh","cloud","precip",
  "rsds_mean","rsds_sd","elev","ntl","invest","urban_rate")
ORD <- c("imperv","grass","water","ever_needle","deci_needle","ever_broad",
         "deci_broad","mixedleaf","tavg","rh","cloud","precip","rsds_mean",
         "rsds_sd","elev","ntl","invest","urban_rate")
stars <- function(p) fifelse(p < .001, "***", fifelse(p < .01, "**",
                     fifelse(p < .05, "*", fifelse(p < .1, "†", ""))))
fp <- function(p) fifelse(p < .001, format(p, digits = 2, scientific = TRUE),
                          sprintf("%.3f", p))
md <- c()
add <- function(...) md <<- c(md, paste0(...))
tbl <- function(dt) {
  add("| ", paste(names(dt), collapse = " | "), " |")
  add("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  for (i in seq_len(nrow(dt)))
    add("| ", paste(sapply(dt[i], as.character), collapse = " | "), " |")
  add("")
}

add("# VPD→SIF 因果方向的站点分类：模型结果表")
add("")
add("因变量 = 站点在 tp=0 的 S-map 同期偏导数符号（1 = 负 = 抑制组，0 = 正 = 促进组）。")
add("自变量全部标准化（先 1% 缩尾再 z 标准化），故系数可横向比较大小。")
add("地表覆盖为成分数据，以耕地为参照类剔除，各地类系数读作「替代耕地的效应」。")
add("显著性: *** p<.001, ** p<.01, * p<.05, † p<.1。")
add("")

# ===== 表 1 / 2: 二分类 logit ==============================================
co <- fread(file.path(OUT, "class_logit_coef.csv"))
co[, `:=`(cn = CN[var], dim = DIMOF[var])]
mk <- function(smp, wm) {
  d <- co[sample == smp & wmode == wm & var != "(Intercept)"]
  setorderv(d, "var"); d <- d[match(ORD, var)]
  data.table(维度 = d$dim, 变量 = d$cn,
             `β` = sprintf("%.3f", d$beta),
             `SE(Conley)` = sprintf("%.3f", d$se_conley),
             z = sprintf("%.2f", d$z_c),
             p = paste0(fp(d$p_conley), stars(d$p_conley)),
             `OR [95% CI]` = sprintf("%.2f [%.2f, %.2f]", d$or, d$or_lo, d$or_hi))
}
add("## 表 1　主模型：加权二元 logistic 回归（权重 = |mean_coef(tp=0)|）")
add("")
add("**样本 A（严格四分类站，n = 385；抑制 214 / 促进 171）**")
add(""); tbl(mk("A", "wgt"))
add("**样本 B（全部完整站，n = 762；抑制 406 / 促进 356）**")
add(""); tbl(mk("B", "wgt"))
add("> 标准误为 Conley 空间 HAC（Bartlett 核，截断 200 km）。最大 VIF：样本 A 5.38（平均气温），样本 B 5.73。")
add("> 权重有效样本量 ESS：样本 A 225.4（58.5% of n），样本 B 403.6（53.0% of n）。")
add("")
add("## 表 2　敏感性：不加权 logistic 回归")
add(""); add("**样本 A**"); add(""); tbl(mk("A", "unw"))
add("**样本 B**"); add(""); tbl(mk("B", "unw"))

# 稳健性小结
w <- dcast(co[wmode == "wgt" & var != "(Intercept)"], var ~ sample,
           value.var = c("beta", "p_conley"))
w[, `:=`(cn = CN[var], 同号 = sign(beta_A) == sign(beta_B),
         sigA = p_conley_A < .05, sigB = p_conley_B < .05)]
add("## 表 3　跨样本一致性")
add("")
tbl(data.table(
  类别 = c("两样本同号且都显著", "仅样本 A 显著", "仅样本 B 显著", "两样本都不显著", "两样本变号"),
  个数 = c(w[sigA & sigB, .N], w[sigA & !sigB, .N], w[!sigA & sigB, .N],
           w[!sigA & !sigB, .N], w[同号 == FALSE, .N]),
  变量 = c(paste(w[sigA & sigB]$cn, collapse = "、"),
           paste(w[sigA & !sigB]$cn, collapse = "、"),
           paste(w[!sigA & sigB]$cn, collapse = "、"),
           paste(w[!sigA & !sigB]$cn, collapse = "、"),
           paste(w[同号 == FALSE]$cn, collapse = "、"))))
add("> 最后一行与前四行不互斥：落叶针叶林两样本都不显著，同时两样本符号相反（β = +0.099 / −0.059），")
add("> 这正说明它没有稳定效应，而不是存在口径冲突。")

# ===== 表 3b: 地表覆盖稀疏度诊断 ==========================================
obj <- readRDS(file.path(OUT, "class_model_data.rds"))
DA <- obj$D[strict == TRUE]
LC <- c("imperv","crop","water","grass","ever_needle","deci_broad",
        "ever_broad","deci_needle","mixedleaf")
sp <- rbindlist(lapply(LC, function(v) { x <- DA[[v]]
  data.table(var = v, 非零站 = sum(x > 0), 非零占比 = sprintf("%.1f%%", 100*mean(x > 0)),
             P90 = sprintf("%.4f", quantile(x, .9)), 最大 = sprintf("%.4f", max(x))) }))
sp[, cn := CN[var]]
sp <- merge(sp, co[sample == "A" & wmode == "wgt", .(var, p_conley)], by = "var", all.x = TRUE)
setorder(sp, 非零站)
add("## 表 3b　地表覆盖变量的稀疏度诊断（样本 A）")
add("")
add("成分数据里多数地类在多数站为 0。占比越稀疏，标准化后少数站对系数的杠杆越大，")
add("「显著」越可能是十几个站撑出来的假象。这张表决定哪些地类系数可以引用。")
add("")
tbl(data.table(变量 = sp$cn, 非零站 = sp$非零站, 非零占比 = sp$非零占比,
               `P90 占比` = sp$P90, `最大占比` = sp$最大,
               `主模型 p` = fifelse(is.na(sp$p_conley), "参照类",
                                    paste0(fp(sp$p_conley), stars(sp$p_conley))),
               判断 = fifelse(is.na(sp$p_conley), "参照类",
                        fifelse(sp$非零站 < 30, "不可引用",
                        fifelse(as.numeric(sp$P90) < 0.01, "谨慎", "可引用")))))
add("> **混交林（13 站、最大占比 0.14%）与落叶针叶林（25 站）不可引用**：")
add("> 缩尾后混交林只剩 10 个不同取值、标准差 8e−5，其 p = 0.002 是数值假象而非生态信号。")
add("> 落叶阔叶林非零站 143 个但 P90 仅 0.44%，同属需谨慎之列。")
add("> 因此表 3 里「两样本同号且都显著」的 6 个变量，真正稳健的是 5 个：")
add("> 海拔、辐射均值、城镇化率、常绿阔叶林、水体。")
add("")

# ===== 表 4: 模型比较 ======================================================
auc_w <- function(y, s, wt = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(wt[i]*wt[n]*(d>0)) +
                   .5*sum(wt[i]*wt[n]*(d==0)); den <- den + sum(wt[i]*wt[n]) }
  num/den }
lcv <- fread(file.path(OUT, "class_logit_cv.csv"))
xcv <- fread(file.path(OUT, "class_xgb_cv.csv"))[, .(sample, model = "XGBoost", pred, y, w, absco)]
cv <- rbind(lcv[, .(sample, model = ifelse(model == "logit_wgt", "logit（加权）", "logit（不加权）"),
                    pred, y, w, absco)], xcv)
cv[, tert := cut(absco, quantile(absco, 0:3/3), include.lowest = TRUE,
                 labels = c("弱","中","强")), by = .(sample, model)]
a <- cv[, .(整体 = auc_w(y, pred)), by = .(sample, model)]
b <- dcast(cv[, .(A = auc_w(y, pred)), by = .(sample, model, tert)],
           sample + model ~ tert, value.var = "A")
ab <- merge(a, b, by = c("sample", "model"))
add("## 表 4　空间分块交叉验证的 AUC（折 = 经纬度 k-means 10 类）")
add("")
tbl(data.table(样本 = ab$sample, 模型 = ab$model,
               整体 = sprintf("%.3f", ab$整体),
               `弱信号层` = sprintf("%.3f", ab$弱),
               `中信号层` = sprintf("%.3f", ab$中),
               `强信号层` = sprintf("%.3f", ab$强)))
add("> 三分位按 |mean_coef(tp=0)| 划分。AUC = 0.5 为随机猜测。")
add("")

# ===== 表 5: 非线性回检 ====================================================
nl <- fread(file.path(OUT, "class_logit_nonlin_cmp.csv"))
add("## 表 5　非线性设定的回检")
add("")
tbl(data.table(样本 = nl$sample, 模型 = nl$model, 参数数 = nl$k,
               AIC = round(nl$aic), `空间分块 AUC` = sprintf("%.3f", nl$auc),
               `ΔAUC(比 M0)` = sprintf("%+.3f", nl$dAUC),
               `LR χ²` = ifelse(nl$ddf > 0, sprintf("%.1f", nl$LRchi), "—"),
               `LR p` = ifelse(nl$ddf > 0, paste0(fp(nl$LRp), stars(nl$LRp)), "—")))
add("> M1 加海拔自然样条(df=3)；M2 再加辐射均值样条；M3 再加五个林型的存在虚拟变量；M4 再加海拔×辐射交互。")
add("")

# ===== 表 6: SHAP 重要性 ===================================================
imp <- fread(file.path(OUT, "class_xgb_imp.csv"))
imp[, cn := CN[var]]
sh <- dcast(imp, var + cn ~ sample, value.var = c("mean_abs", "rel", "cor_xs"))
setorder(sh, -mean_abs_A)
add("## 表 6　XGBoost 折外 TreeSHAP 重要性（单位：对数几率）")
add("")
tbl(data.table(变量 = sh$cn,
               `平均绝对SHAP(A)` = sprintf("%.3f", sh$mean_abs_A),
               `占比 A` = sprintf("%.1f%%", 100*sh$rel_A),
               `平均绝对SHAP(B)` = sprintf("%.3f", sh$mean_abs_B),
               `占比 B` = sprintf("%.1f%%", 100*sh$rel_B),
               方向 = fifelse(is.na(sh$cor_xs_A), "树从未用它分裂",
                        fifelse(sh$cor_xs_A > 0, "越大越偏抑制组", "越大越偏促进组"))))
add("> 树模型不受成分共线之困，故耕地放回，共 19 个特征。方向由 SHAP 值与变量取值的 Spearman 相关定号。")
add("> 落叶针叶林与混交林的 SHAP 恒为 0：两者分别只有 6.5% / 3.4% 的站非零，")
add("> 在 min_child_weight = 8 的正则下树从未在它们上分裂。混交林在 logit 里显著而在树里为零，")
add("> 正是「线性模型被少数非零站撑起」与「树模型拒绝在稀疏变量上分裂」两种行为的直接对照。")
add("")
ia <- fread(file.path(OUT, "class_xgb_inter.csv"))[sample == "A"][1:6]
ia[, pair := sapply(strsplit(key, " "), function(v) paste(CN[v], collapse = " × "))]
add("**样本 A 最强的 6 个成对交互（|SHAP 交互值| 均值之和）**")
add(""); tbl(data.table(变量对 = ia$pair, 强度 = sprintf("%.3f", ia$strength)))

# ===== 表 7: 四分类两层 ====================================================
L <- fread(file.path(OUT, "class4_layer_coef.csv"))[var != "(Intercept)"]
L[, cn := CN[var]]
add("## 表 7　四分类的两层分解（样本 A）")
add("")
add("四类 = 起始方向（第 1 层）× 后续是否翻转（第 2 层）。两层都是二元 logit，标准误为 Conley 200 km。")
add("")
for (ly in unique(L$layer)) {
  d <- L[layer == ly][match(ORD, var)]
  add("**", ly, "**"); add("")
  tbl(data.table(变量 = d$cn, `β` = sprintf("%.3f", d$beta),
                 `SE(Conley)` = sprintf("%.3f", d$se_conley),
                 p = paste0(fp(d$p), stars(d$p)),
                 OR = sprintf("%.2f", d$or)))
}
add("> 第 2 层 a 的样本为 214 站、其中 194 站翻转（少数类仅 9.3%），拟合概率有 14.5% 贴边，")
add("> 属准完全分离，系数与 p 值不可解读；混交林在全样本仅 3.4% 的站非零。")
add("")

# ===== 表 8: 四分类 XGBoost ================================================
cm <- fread(file.path(OUT, "class4_xgb_cv.csv"))
LV <- c("始终抑制","抑制转促进","始终促进","促进转抑制")
KEY <- c(always_inhibit="始终抑制", inhibit_promote="抑制转促进",
         always_promote="始终促进", promote_inhibit="促进转抑制")
cm[, 真实 := KEY[stype]]
cm[, 预测 := LV[max.col(as.matrix(.SD))], .SDcols = c("p1","p2","p3","p4")]
ova <- sapply(1:4, function(k) auc_w(as.integer(cm$真实 == LV[k]),
                                     as.matrix(cm[, .(p1,p2,p3,p4)])[, k]))
add("## 表 8　四分类 XGBoost 的折外表现（样本 A）")
add("")
tbl(data.table(类别 = LV, n = as.integer(table(factor(cm$真实, LV))),
               `占比` = sprintf("%.1f%%", 100*as.numeric(table(factor(cm$真实, LV)))/nrow(cm)),
               `one-vs-rest AUC` = sprintf("%.3f", ova),
               `召回率` = sprintf("%.0f%%", 100*sapply(LV, function(v)
                 mean(cm[真实 == v]$预测 == v)))))
add("> 总体准确率 0.681（全猜最大类的基准 0.504），宏平均 AUC 0.735。")
add("> 「始终抑制」只有 20 站，折外没有任何一站被预测为该类。")
add("")

# ===== 表 9: 多项 logit ====================================================
mo <- fread(file.path(OUT, "class4_multinom_coef.csv"))[var != "(Intercept)"]
mo[, `:=`(cn = CN[var], 对比 = KEY[vs])]
mo <- mo[order(p)][1:15]
add("## 表 9　多项 logistic 回归（参照类 = 始终抑制，仅作定性参考）")
add("")
tbl(data.table(`对比类别` = mo$对比, 变量 = mo$cn,
               `β` = sprintf("%.2f", mo$beta), SE = sprintf("%.2f", mo$se),
               p = paste0(fp(mo$p), stars(mo$p))))
add("> 只列 p 值最小的 15 行。该模型的标准误未做空间修正（多项 logit 上 Conley 无现成实现），")
add("> 且参照类只有 20 站，所有系数的绝对水平都不可靠，仅用于核对方向是否与两层分解一致。")

# ===== 表 10: 二分类 vs 四分类的取舍 ======================================
fp10 <- file.path(OUT, "class_flip_predictability.csv")
if (file.exists(fp10)) {
  r <- fread(fp10)
  add("## 表 10　该用二分类还是四分类：逐层可预测性")
  add("")
  add("四分类 = 起始方向 × 是否翻转。二分类只用第一维。下表把两维各自的折外可预测性分开量：")
  add("")
  tbl(data.table(层 = r$层, n = r$n, 阳性率 = sprintf("%.3f", r$阳性率),
                 `AUC logit` = sprintf("%.3f", r$AUC_logit),
                 `AUC XGBoost` = sprintf("%.3f", r$AUC_xgb),
                 `伪R² logit` = sprintf("%.3f", r$`伪R2_logit`),
                 `伪R² XGBoost` = sprintf("%.3f", r$`伪R2_xgb`)))
  add("> 伪 R² = 1 − 折外对数似然 / 常数模型对数似然；≤ 0 表示还不如直接用基础发生率。")
  add("> 折的划分与表 4 完全相同（经纬度 k-means 10 类，seed 42），故数值可直接对比。")
  add("> 第 2 层 a 的 XGBoost AUC 0.166 远低于 0.5，不是「反向可预测」，是 20 个少数类站散在 10 个空间块里")
  add("> 造成的退化现象——该层没有可恢复的信号。")
  add("> 第 1 层 logit 的 AUC 0.864 与伪 R² 0.053 并存，说明 logit 跨空间块<b>排序好但校准差</b>：")
  add("> 谁更可能属抑制组排得准，概率的绝对水平偏极端。若下游要用概率值而非排序，用 XGBoost 的输出。")
  add("")
  add("**样本量代价**：二分类可用 762 站，四分类只有严格分类的 385 站，丢掉 377 站（49.5%）；")
  add("四类中最小的「始终抑制」仅 20 站，按 10 折分平均每折 2 站。")
  add("")
}

writeLines(md, file.path(OUT, "class_results_tables.md"))
cat("已写出 class_results_tables.md（", length(md), " 行）\n", sep = "")
