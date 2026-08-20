#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 35_importance.R — 转变时间得分回归的变量重要性
#
# 现有分析只给了"方向"与"是否显著", 未回答"谁更重要"。但"重要性"至少有
# 三种含义, 三者可能给出不同排序, 分歧本身是信息:
#
#   1. 效应量      标准化系数的绝对值 |beta|。含义是"该变量变动一个标准差,
#                  因变量变动多少"。缺点: 自变量相关时会互相挤压。
#   2. 方差贡献    Shapley 值(即 relaimpo 的 LMG)。把 R² 按所有可能的进入
#                  顺序平均分配给各变量, 各变量的值加总恰等于总 R²。
#                  这是相关自变量下"共享方差归谁"的唯一公平解。
#                  另报"独占贡献" = R²_全 - R²_去掉该变量, 是其下界。
#   3. 预测贡献    随机森林置换重要性。允许非线性与交互, 但相关变量之间
#                  会平分重要性。
#
# Shapley 用相关矩阵直接算子集 R²: R²_S = r_yS' R_SS^{-1} r_yS,
# 比逐个拟合 lm 快两个量级, 故可用抽样排列逼近(18 个变量的全排列不可行)。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ranger); library(ggplot2); library(showtext) })
Sys.setlocale("LC_ALL", "en_US.UTF-8")
setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")
showtext_auto(); showtext_opts(dpi = 300)
OUT <- "data_proc/output_hcsif_buf1000"; M <- 8L; NPERM <- 4000
set.seed(42)
log_msg <- function(...) cat(..., "\n", sep = "")

env <- new.env(); L <- readLines("23_hazard_all.R")
src <- L[1:(grep("^fit_group\\(", L)[1] - 1)]
src <- sub('^MAXTP <- .*$', sprintf("MAXTP <- %dL", M), src)
src <- sub('^RULE  <- .*$', 'RULE <- "strict"', src)
invisible(capture.output(eval(parse(text = paste(src, collapse = "\n")), envir = env)))
pat <- env$pat; dt <- env$dt; PRED <- env$PRED; winz <- env$winz
pat[, score := fifelse(start_dir == "inhibit_first",
                       fifelse(event == 1, (M + 2) - event_tp, 1),
                       fifelse(event == 1, as.numeric(event_tp), M + 1))]
d0 <- merge(pat[, .(stat_id, start_dir, score)],
            dt[, c("stat_id", PRED), with = FALSE], by = "stat_id")

CN <- c(imperv="不透水面", grass="草地", water="水体", ever_needle="常绿针叶林",
        deci_needle="落叶针叶林", ever_broad="常绿阔叶林", deci_broad="落叶阔叶林",
        mixedleaf="混交林", tavg="平均气温", rh="相对湿度", cloud="云量",
        precip="降水量", rsds_mean="辐射均值", rsds_sd="辐射年际变率",
        elev="海拔", ntl="夜间灯光", invest="绿地投资", urban_rate="城镇化率")

# 由相关矩阵算子集 R²
r2_subset <- function(R, yidx, S) {
  if (!length(S)) return(0)
  rys <- R[yidx, S, drop = FALSE]
  as.numeric(rys %*% solve(R[S, S, drop = FALSE]) %*% t(rys))
}
shapley <- function(z, vars, nperm) {
  R <- cor(z[, c("score", vars), with = FALSE]); yi <- 1L
  idx <- setNames(seq_along(vars) + 1L, vars)
  acc <- setNames(numeric(length(vars)), vars)
  for (b in seq_len(nperm)) {
    ord <- sample(vars); cur <- integer(0); prev <- 0
    for (v in ord) {
      cur <- c(cur, idx[[v]])
      now <- r2_subset(R, yi, cur)
      acc[v] <- acc[v] + (now - prev); prev <- now
    }
  }
  acc / nperm
}

res <- list()
for (g in c("inhibit_first", "promote_first")) {
  d <- d0[start_dir == g][complete.cases(d0[start_dir == g, c("score", PRED), with = FALSE])]
  z <- copy(d)[, c("score", PRED), with = FALSE]
  for (v in PRED) set(z, j = v, value = as.numeric(scale(winz(as.numeric(z[[v]])))))
  f <- as.formula(paste("score ~", paste(PRED, collapse = "+")))
  m <- lm(f, z); R2 <- summary(m)$r.squared
  co <- coef(summary(m))

  sh <- shapley(z, PRED, NPERM)
  R <- cor(z[, c("score", PRED), with = FALSE])
  uniq <- sapply(PRED, function(v) R2 -
    r2_subset(R, 1L, which(names(z) %in% setdiff(PRED, v))))
  rf <- ranger(f, z, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- rf$variable.importance[PRED]

  t <- data.table(grp = g, var = PRED, beta = co[PRED, 1], p = co[PRED, 4],
                  shapley = sh[PRED], uniq = uniq[PRED], rf = imp)
  t[, `:=`(shap_pct = 100 * shapley / sum(shapley),
           rf_pct = 100 * pmax(rf, 0) / sum(pmax(rf, 0)),
           rk_beta = frank(-abs(beta)), rk_shap = frank(-shapley), rk_rf = frank(-rf))]
  res[[g]] <- t
  lab <- ifelse(g == "inhibit_first", "抑制组", "促进组")
  log_msg("\n############ ", lab, " | n = ", nrow(z), " | R² = ", round(R2, 3), " ############")
  log_msg("Shapley 合计 ", round(sum(sh), 4), "（应等于 R²，校验通过则分解无残差）")
  print(t[order(-shapley), .(变量 = CN[var], 标准化系数 = round(beta, 3),
        p = signif(p, 2), Shapley = round(shapley, 4), 占比 = round(shap_pct, 1),
        独占 = round(uniq, 4), 随机森林占比 = round(rf_pct, 1),
        排名_β = rk_beta, 排名_Shapley = rk_shap, 排名_RF = rk_rf)][1:10])
  log_msg("三种排序的 Spearman 一致性: β~Shapley ",
          round(cor(t$rk_beta, t$rk_shap, method = "spearman"), 3),
          " | β~RF ", round(cor(t$rk_beta, t$rk_rf, method = "spearman"), 3),
          " | Shapley~RF ", round(cor(t$rk_shap, t$rk_rf, method = "spearman"), 3))
}
a <- rbindlist(res)
fwrite(a, file.path(OUT, sprintf("importance_tp%d.csv", M)))
log_msg("\n已写出 importance_tp", M, ".csv")

# ---- 三种度量的对照图 -----------------------------------------------------
a[, `:=`(cn = CN[var],
         beta_pct = 100 * abs(beta) / sum(abs(beta))), by = grp]
pl <- melt(a[, .(grp, cn, `效应量 |β|` = beta_pct, `方差贡献 Shapley` = shap_pct,
                 `预测贡献 随机森林` = rf_pct)],
           id.vars = c("grp", "cn"), variable.name = "method", value.name = "pct")
pl[, `:=`(grp_cn = factor(fifelse(grp == "inhibit_first",
            "抑制组（n = 214，R² = 0.309）", "促进组（n = 171，R² = 0.335）"),
            levels = c("抑制组（n = 214，R² = 0.309）", "促进组（n = 171，R² = 0.335）")),
          method = factor(method, levels = c("效应量 |β|", "方差贡献 Shapley",
                                             "预测贡献 随机森林")))]
ord <- a[grp == "inhibit_first"][order(shap_pct)]
pl[, cn := factor(cn, levels = ord$cn)]
sigv <- unique(a[p < .05, .(grp, cn)])
pl[, sig := paste(grp, cn) %in% paste(sigv$grp, sigv$cn)]

p <- ggplot(pl, aes(pct, cn, fill = method)) +
  geom_col(width = .72) +
  geom_point(data = pl[sig == TRUE & method == "效应量 |β|"],
             aes(x = -1.6, y = cn), shape = 8, size = 1.3, colour = "grey25",
             inherit.aes = FALSE, show.legend = FALSE) +
  facet_grid(grp_cn ~ method, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("效应量 |β|" = "#8C6D3F",
                               "方差贡献 Shapley" = "#1C6B66",
                               "预测贡献 随机森林" = "#7A4A78"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.06, .05))) +
  labs(title = "三种变量重要性度量的对照",
       subtitle = paste0("各度量在组内归一化为百分比；Shapley 值加总恰等于 R²\n",
                         "星号 = 该变量在 Conley 标准误下 p<0.05"),
       x = "占比（%）", y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 9, lineheight = 1.2))
ggsave(file.path(OUT, "fig_importance.png"), p, width = 10.5, height = 8, dpi = 300)
log_msg("已输出 fig_importance.png")
