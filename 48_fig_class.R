#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 48_fig_class.R — 分类分析的全套图
#   fig_class_coef      加权 logit 的 Conley 系数森林图(A / B 两样本)
#   fig_class_auc       空间分块 CV 的 AUC 对比 + 按 |coef0| 分层
#   fig_class_shap      SHAP 蜂群图(样本 A)
#   fig_class_shap_dep  前 6 个变量的 SHAP 依赖图(看非线性与阈值)
#   fig_class_compare   logit 标准化系数 vs SHAP 重要性
#   fig_class_map       站点地图: 实测组别 vs 折外预测概率
#   fig_class4          四分类: 两层分解系数 / 逐类 SHAP / 混淆矩阵
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(sysfonts)
  library(patchwork) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
font_add("heiti", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"
FIG <- function(n) file.path(OUT, paste0(n, ".png"))
th <- function(b = 11) theme_minimal(base_size = b, base_family = "heiti") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = .35),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold", size = b + 2),
        plot.subtitle = element_text(colour = "grey35", size = b - 1.5, lineheight = 1.15),
        plot.caption = element_text(colour = "grey45", size = b - 2.5, hjust = 0))

CN <- c(imperv="不透水面", grass="草地", water="水体", crop="耕地",
        ever_needle="常绿针叶林", deci_needle="落叶针叶林",
        ever_broad="常绿阔叶林", deci_broad="落叶阔叶林", mixedleaf="混交林",
        tavg="平均气温", rh="相对湿度", cloud="云量", precip="降水量",
        rsds_mean="辐射均值", rsds_sd="辐射年际变率", elev="海拔",
        ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")
DIMOF <- c(rep("土地利用", 9), rep("气候环境", 7), rep("社会经济", 3))
names(DIMOF) <- c("imperv","grass","water","crop","ever_needle","deci_needle",
                  "ever_broad","deci_broad","mixedleaf",
                  "tavg","rh","cloud","precip","rsds_mean","rsds_sd","elev",
                  "ntl","invest","urban_rate")
PALD <- c(土地利用 = "#5B8C5A", 气候环境 = "#2E6E8E", 社会经济 = "#C2703A")
SMP <- c(A = "样本 A：严格四分类站 n=385", B = "样本 B：全部完整站 n=762")

# ===== 1. 系数森林图 =======================================================
co <- fread(file.path(OUT, "class_logit_coef.csv"))[var != "(Intercept)" & wmode == "wgt"]
co[, `:=`(cn = CN[var], dim = DIMOF[var],
          lo = beta - 1.96*se_conley, hi = beta + 1.96*se_conley,
          sig = p_conley < .05, smp = SMP[sample])]
ord <- co[sample == "A"][order(beta)]$cn
co[, cn := factor(cn, levels = ord)]
p1 <- ggplot(co, aes(beta, cn, colour = dim)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .65) +
  geom_point(aes(shape = sig), size = 2.5, stroke = .9) +
  facet_wrap(~ smp, nrow = 1) +
  scale_colour_manual(values = PALD, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                     labels = c("Conley p ≥ 0.05", "Conley p < 0.05"), name = NULL) +
  labs(title = "什么决定一个站属于抑制组还是促进组",
       subtitle = paste0("加权 logistic 回归（权重 = |mean_coef(tp=0)|），",
                         "自变量已标准化，误差棒为 Conley 空间 HAC 95% 区间（截断 200 km）\n",
                         "系数在零线右边 = 该变量越大越可能属于抑制组，左边 = 越可能属于促进组；",
                         "实心点 = p < 0.05，空心点 = 不显著"),
       x = "标准化 logit 系数（对数几率）", y = NULL,
       caption = "地表覆盖以耕地为参照类剔除，故各地类系数读作「替代耕地的效应」") +
  th() + theme(legend.position = "top")
ggsave(FIG("fig_class_coef"), p1, width = 10.5, height = 6.4, dpi = 300)

# ===== 2. AUC 对比 =========================================================
auc_w <- function(y, s, w = rep(1, length(y))) {
  p <- which(y == 1); n <- which(y == 0); if (!length(p) || !length(n)) return(NA_real_)
  num <- 0; den <- 0
  for (i in p) { d <- sign(s[i] - s[n]); num <- num + sum(w[i]*w[n]*(d>0)) +
                   .5*sum(w[i]*w[n]*(d==0)); den <- den + sum(w[i]*w[n]) }
  num/den }
lcv <- fread(file.path(OUT, "class_logit_cv.csv"))
xcv <- fread(file.path(OUT, "class_xgb_cv.csv"))[, .(sample, model = "XGBoost", pred, y, w, absco)]
cv <- rbind(lcv[, .(sample, model = ifelse(model == "logit_wgt", "logit（加权）", "logit（不加权）"),
                    pred, y, w, absco)], xcv)
a1 <- cv[, .(AUC = auc_w(y, pred)), by = .(sample, model)]
a1[, smp := SMP[sample]]
pa <- ggplot(a1, aes(AUC, model, fill = model)) +
  geom_col(width = .6) + geom_vline(xintercept = .5, colour = "grey40", linetype = 2) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -.15, size = 3.2, family = "heiti") +
  facet_wrap(~ smp, nrow = 1) + coord_cartesian(xlim = c(.5, 1)) +
  scale_fill_manual(values = c("logit（不加权）" = "#9BA8A5", "logit（加权）" = "#1C6B66",
                               XGBoost = "#C2703A"), guide = "none") +
  labs(subtitle = "整体 AUC（空间分块 10 折）", x = NULL, y = NULL) + th()
cv[, tert := cut(absco, quantile(absco, 0:3/3), include.lowest = TRUE,
                 labels = c("弱（|coef| 小）","中","强（|coef| 大）")), by = .(sample, model)]
a2 <- cv[, .(AUC = auc_w(y, pred)), by = .(sample, model, tert)]
a2[, smp := SMP[sample]]
pb <- ggplot(a2, aes(tert, AUC, colour = model, group = model)) +
  geom_hline(yintercept = .5, colour = "grey40", linetype = 2) +
  geom_line(linewidth = .8) + geom_point(size = 2.4) +
  facet_wrap(~ smp, nrow = 1) +
  scale_colour_manual(values = c("logit（不加权）" = "#9BA8A5", "logit（加权）" = "#1C6B66",
                                 XGBoost = "#C2703A"), name = NULL) +
  labs(subtitle = "按 |mean_coef(tp=0)| 三分位分层的 AUC：标签越可信越好预测",
       x = "标签信噪比分层", y = "AUC") + th() + theme(legend.position = "top")
p2 <- pa / pb + plot_annotation(
  title = "预测能力：空间分块交叉验证",
  subtitle = paste0("折 = 对站点经纬度做 k-means 聚成 10 类。随机分折会把空间相邻站拆到不同折，",
                    "AUC 会偏乐观约 0.05\nAUC = 0.5 是抛硬币"),
  theme = th())
ggsave(FIG("fig_class_auc"), p2, width = 10.5, height = 8, dpi = 300)

# ===== 3. SHAP 蜂群图 ======================================================
sh <- fread(file.path(OUT, "class_xgb_shap.csv"))
imp <- fread(file.path(OUT, "class_xgb_imp.csv"))
beeswarm <- function(smp, topn = 14) {
  s <- sh[sample == smp]
  iv <- imp[sample == smp][order(-mean_abs)]$var[1:topn]
  s <- s[var %in% iv]
  s[, xr := (frank(x, ties.method = "average") - 1)/(.N - 1), by = var]
  s[, cn := factor(CN[var], levels = rev(CN[iv]))]
  ggplot(s, aes(shap, cn, colour = xr)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
    geom_jitter(height = .22, size = .75, alpha = .65) +
    scale_colour_gradient(low = "#2E6E8E", high = "#C2703A",
                          breaks = c(0, 1), labels = c("低", "高"),
                          name = "该变量的取值") +
    labs(subtitle = SMP[smp], x = "SHAP 值（对数几率，>0 推向「抑制组」）", y = NULL) +
    th() + theme(legend.position = "right", legend.key.height = unit(1.6, "cm"))
}
p3 <- (beeswarm("A") | beeswarm("B")) + plot_layout(guides = "collect") +
  plot_annotation(title = "XGBoost 的 SHAP 归因：每个点是一个站",
    subtitle = paste0("SHAP 由折外模型算出（每个站的归因来自没见过它的那个模型）。",
                      "横轴 = 该变量把这个站推离基准的幅度\n",
                      "颜色 = 该站在这个变量上的取值高低；一列点若「左蓝右橙」说明变量越大越偏抑制组"),
    theme = th())
ggsave(FIG("fig_class_shap"), p3, width = 13, height = 6.8, dpi = 300)

# ===== 4. SHAP 依赖图 ======================================================
iv6 <- imp[sample == "A"][order(-mean_abs)]$var[1:6]
sd6 <- sh[sample == "A" & var %in% iv6]
sd6[, cn := factor(CN[var], levels = CN[iv6])]
p4 <- ggplot(sd6, aes(x, shap)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_point(aes(colour = factor(grp2)), size = .8, alpha = .55) +
  geom_smooth(method = "loess", se = FALSE, colour = "#1C1C1C", linewidth = .8, span = .8) +
  facet_wrap(~ cn, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = c(`0` = "#C2703A", `1` = "#2E6E8E"),
                      labels = c("促进组", "抑制组"), name = "实测组别") +
  labs(title = "SHAP 依赖图：变量与效应之间是不是直线",
       subtitle = paste0("样本 A。黑线是 loess 平滑。若明显弯折或有平台，说明 logit 里应该加样条；",
                         "若基本是直线，logit 的线性设定就够用\n",
                         "纵轴 > 0 表示该变量在这个站把预测推向抑制组"),
       x = "变量原始取值", y = "SHAP 值（对数几率）") +
  th() + theme(legend.position = "top")
ggsave(FIG("fig_class_shap_dep"), p4, width = 11, height = 6.6, dpi = 300)

# ===== 5. logit 系数 vs SHAP 重要性 ========================================
cmp <- merge(co[sample == "A", .(var, beta, p_conley, sig)],
             imp[sample == "A", .(var, rel, cor_xs)], by = "var")
cmp[, `:=`(cn = CN[var], dim = DIMOF[var])]
p5 <- ggplot(cmp, aes(abs(beta), rel, colour = dim)) +
  geom_point(aes(shape = sig), size = 3) +
  ggrepel::geom_text_repel(aes(label = cn), size = 3, family = "heiti",
                           max.overlaps = 30, seed = 1) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("Conley p ≥ .05", "Conley p < .05"), name = NULL) +
  scale_colour_manual(values = PALD, name = NULL) +
  labs(title = "两条路线给出的重要性排序对不对得上",
       subtitle = paste0("横轴 = logit 标准化系数的绝对值（线性、可推断）；",
                         "纵轴 = SHAP 相对重要性（含非线性与交互，不可推断）\n",
                         "落在左上 = 效应主要是非线性的，线性模型会低估；",
                         "落在右下 = 线性模型里显著但树模型不太用它"),
       x = "|logit 标准化系数|", y = "SHAP 相对重要性") +
  th() + theme(legend.position = "top")
if (requireNamespace("ggrepel", quietly = TRUE))
  ggsave(FIG("fig_class_compare"), p5, width = 8.5, height = 6.6, dpi = 300)

# ===== 6. 地图 =============================================================
obj <- readRDS(file.path(OUT, "class_model_data.rds"))
mp <- merge(fread(file.path(OUT, "class_xgb_cv.csv"))[sample == "A"],
            obj$D[, .(stat_id, longitude, latitude)], by = "stat_id")
pm1 <- ggplot(mp, aes(longitude, latitude)) +
  geom_point(aes(colour = factor(y), size = absco), alpha = .8) +
  scale_colour_manual(values = c(`0` = "#C2703A", `1` = "#2E6E8E"),
                      labels = c("促进组", "抑制组"), name = "实测组别") +
  scale_size_continuous(range = c(.6, 3.2), name = "|coef(tp=0)|") +
  coord_quickmap() + labs(subtitle = "实测：点越大标签越可信") + th()
pm2 <- ggplot(mp, aes(longitude, latitude)) +
  geom_point(aes(colour = pred), size = 1.7) +
  scale_colour_gradient2(low = "#C2703A", mid = "grey88", high = "#2E6E8E",
                         midpoint = .5, name = "P(抑制组)") +
  coord_quickmap() + labs(subtitle = "折外预测概率（该站所在空间块未参与训练）") + th()
p6 <- (pm1 | pm2) + plot_annotation(
  title = "空间格局：抑制组集中在华北平原与东北，促进组在南方与西部高地",
  subtitle = "样本 A。预测图与实测图的相似度就是空间分块 AUC 想量化的东西",
  theme = th())
ggsave(FIG("fig_class_map"), p6, width = 12.5, height = 5.6, dpi = 300)

# ===== 7. 四分类 ===========================================================
L <- fread(file.path(OUT, "class4_layer_coef.csv"))[var != "(Intercept)"]
L[, `:=`(cn = CN[var], lo = beta - 1.96*se_conley, hi = beta + 1.96*se_conley,
         sig = p < .05)]
L[, layer := factor(layer, levels = unique(layer))]
L[, cn := factor(cn, levels = rev(CN[c("imperv","grass","water","ever_needle",
  "deci_needle","ever_broad","deci_broad","mixedleaf","tavg","rh","cloud","precip",
  "rsds_mean","rsds_sd","elev","ntl","invest","urban_rate")]))]
q1 <- ggplot(L, aes(beta, cn, colour = DIMOF[as.character(var)])) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .6) +
  geom_point(aes(shape = sig), size = 2.3, stroke = .85) +
  facet_wrap(~ layer, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = PALD, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                     labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  labs(subtitle = paste0("两层分解：四类 = 起始方向 × 后续是否翻转（各层都能上 Conley 标准误）。",
                         "实心 = p < 0.05，空心 = 不显著"),
       x = "标准化 logit 系数", y = NULL,
       caption = paste0("第2层a 里「混交林」的巨大系数是伪信号：全样本只有 3.4% 的站混交林占比非零，",
                        "且该层 194/214 都翻转，属准完全分离，不要解读")) +
  th(10) + theme(legend.position = "top")
i4 <- fread(file.path(OUT, "class4_xgb_imp.csv"))
i4[, cor_xs := fifelse(is.na(cor_xs), 0, cor_xs)]   # 该类内某地类恒为 0 时相关系数为 NA
i4[, `:=`(cn = CN[var], signed = rel * sign(cor_xs))]
ordv <- i4[, .(s = sum(rel)), by = cn][order(s)]$cn
i4[, cn := factor(cn, levels = ordv)]
i4[, cls := factor(cls, levels = c("始终抑制","抑制转促进","始终促进","促进转抑制"))]
q2 <- ggplot(i4, aes(cls, cn, fill = signed)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = ifelse(rel > .06, sprintf("%.0f%%", 100*rel), "")),
            size = 2.6, family = "heiti", colour = "grey15") +
  scale_fill_gradient2(low = "#2E6E8E", mid = "white", high = "#C2703A", midpoint = 0,
                       name = "带符号的\nSHAP 相对重要性") +
  labs(subtitle = "逐类 SHAP 重要性（橙 = 该变量越大越像这一类，蓝 = 越不像）",
       x = NULL, y = NULL) + th(10)
cm <- fread(file.path(OUT, "class4_xgb_cv.csv"))
LV <- c("始终抑制","抑制转促进","始终促进","促进转抑制")
KEY <- c(always_inhibit="始终抑制", inhibit_promote="抑制转促进",
         always_promote="始终促进", promote_inhibit="促进转抑制")
cm[, 真实 := KEY[stype]]
cm[, 预测 := LV[max.col(as.matrix(.SD))], .SDcols = c("p1","p2","p3","p4")]
tab <- cm[, .N, by = .(真实, 预测)]
tab[, share := N/sum(N), by = 真实]
tab[, `:=`(真实 = factor(真实, levels = rev(LV)), 预测 = factor(预测, levels = LV))]
q3 <- ggplot(tab, aes(预测, 真实, fill = share)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%d\n%.0f%%", N, 100*share)), size = 2.9, family = "heiti") +
  scale_fill_gradient(low = "white", high = "#2E6E8E", guide = "none") +
  labs(subtitle = "折外混淆矩阵：行内百分比", x = "预测类别", y = "实际类别") + th(10)
p7 <- q1 / (q2 | q3) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(title = "四分类：始终抑制 / 抑制转促进 / 始终促进 / 促进转抑制",
    subtitle = paste0("只用样本 A（385 站），因为四类只在「至多变号一次」的站上有定义。",
                      "四类样本量 20 / 194 / 69 / 102，「始终抑制」只有 20 站"),
    theme = th())
ggsave(FIG("fig_class4"), p7, width = 13, height = 11, dpi = 300)

cat("已输出 7 张图到 ", OUT, "\n", sep = "")

# ===== 8. 非线性检验 =======================================================
nl <- fread(file.path(OUT, "class_logit_nonlin_cmp.csv"))
nl[, `:=`(smp = SMP[sample], model = factor(model, levels = unique(model)))]
n1 <- ggplot(nl, aes(model, auc, group = smp, colour = smp)) +
  geom_line(linewidth = .8) + geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.1, size = 2.8, family = "heiti",
            show.legend = FALSE) +
  scale_colour_manual(values = c("#1C6B66", "#C2703A"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(.08, .16))) +
  labs(subtitle = "空间分块 CV 的 AUC：加了样条、存在虚拟、交互项之后并没有变好",
       x = NULL, y = "AUC") + th(10) +
  theme(legend.position = "top", axis.text.x = element_text(size = 8))
n2 <- ggplot(nl, aes(model, aic, group = smp, colour = smp)) +
  geom_line(linewidth = .8) + geom_point(size = 2.6) +
  scale_colour_manual(values = c("#1C6B66", "#C2703A"), guide = "none") +
  facet_wrap(~ smp, scales = "free_y", nrow = 1) +
  labs(subtitle = "样本内 AIC：样本 B 一路下降；样本 A 只有辐射样条真正改善了拟合",
       x = NULL, y = "AIC") + th(10) + theme(axis.text.x = element_text(size = 7))
cu <- fread(file.path(OUT, "class_logit_nonlin_curve.csv"))
cu[, `:=`(smp = SMP[sample], cn = CN[var])]
n3 <- ggplot(cu, aes(x, fit)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = smp), alpha = .18, colour = NA) +
  geom_line(aes(colour = smp), linewidth = .85) +
  facet_wrap(~ cn, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("#1C6B66", "#C2703A"), name = NULL) +
  scale_fill_manual(values = c("#1C6B66", "#C2703A"), guide = "none") +
  labs(subtitle = "自然样条的偏效应曲线（以左端为 0 基准），阴影为 Conley 95% 区间",
       x = "变量原始取值", y = "对数几率") + th(10) + theme(legend.position = "top")
p8 <- n1 / n3 / n2 + plot_layout(heights = c(1, 1, .9)) +
  plot_annotation(title = "把 SHAP 看到的非线性搬回 logit：样本内变好，样本外没变好",
    subtitle = paste0("M1 加海拔自然样条，M2 再加辐射样条，M3 再加五个林型的「有/无」虚拟变量，",
                      "M4 再加海拔×辐射交互\n",
                      "似然比检验全部显著、样本 B 的 AIC 一路下降，但空间分块 AUC 不升反降\n",
                      "= 这些弯曲是本地的、跨不过空间块，推断仍应以线性 logit 为准"),
    theme = th())
ggsave(FIG("fig_class_nonlin"), p8, width = 11.5, height = 11, dpi = 300)
cat("已输出 fig_class_nonlin\n")
