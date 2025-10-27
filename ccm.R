# ============================================================================
# CCM分析准备（修正版 - 解决库长度问题）
# ============================================================================

library(rEDM)

# ============================================================================
# 1. 单个站点的CCM分析示例
# ============================================================================

# ============================================================================
# 完整的滞后因果分析：t月温度异常 → t+1月、t+2月SIF
# ============================================================================

# ============================================================================
# 1. 数据准备
# ============================================================================

# 选择一个数据点充足的站点
example_station <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  slice(1) %>%
  pull(meteo_stat_id)

# 提取单个站点数据（去趋势后的数据）
single_station_data <- data_for_ccm %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  select(temp_anomaly, sif) %>%
  as.data.frame()

# 添加时间索引列
single_station_data <- data.frame(
  time = 1:nrow(single_station_data),
  temp_anomaly = single_station_data$temp_anomaly,
  sif = single_station_data$sif
)

cat("=== 站点信息 ===\n")
cat("站点ID:", example_station, "\n")
cat("数据点数:", nrow(single_station_data), "\n")
cat("时间跨度: 6-8月的数据\n\n")

n_data <- nrow(single_station_data)

# ============================================================================
# 2. 创建滞后变量（用于相关性分析）
# ============================================================================

lag_data <- single_station_data %>%
  mutate(
    # t+1月的SIF（lead函数：向前看1个位置）
    sif_lag1 = lead(sif, 1),
    # t+2月的SIF
    sif_lag2 = lead(sif, 2),
    # 同时也保留当月SIF（t+0）
    sif_lag0 = sif
  )

cat("=== 滞后变量创建完成 ===\n")
cat("可用于分析的数据点:\n")
cat("- Tp=0 (同期):", sum(!is.na(lag_data$sif_lag0)), "个\n")
cat("- Tp=1 (滞后1月):", sum(!is.na(lag_data$sif_lag1)), "个\n")
cat("- Tp=2 (滞后2月):", sum(!is.na(lag_data$sif_lag2)), "个\n\n")

# ============================================================================
# 3. 确定最优嵌入维度
# ============================================================================

lib_pred_size <- min(100, n_data)

embed_temp <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  maxE = 10,
  showPlot = FALSE
)

embed_sif <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "sif",
  target = "sif",
  maxE = 10,
  showPlot = FALSE
)

best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
E_for_ccm <- max(best_E_temp, best_E_sif)

cat("=== 嵌入维度 ===\n")
cat("Temperature Anomaly 最优 E:", best_E_temp, "\n")
cat("SIF 最优 E:", best_E_sif, "\n")
cat("CCM使用 E:", E_for_ccm, "\n\n")

# ============================================================================
# 4. CCM分析：测试不同滞后期的因果关系
# ============================================================================

# 计算库长度参数
tau <- 1
embedding_loss <- (E_for_ccm - 1) * tau
max_available_lib <- n_data - embedding_loss
lib_start <- max(E_for_ccm + 2, 10)
lib_end <- max_available_lib
lib_step <- max(2, floor((lib_end - lib_start) / 15))
libSizes_str <- paste(lib_start, lib_end, lib_step)

cat("=== 开始CCM分析 ===\n")
cat("测试滞后期: Tp = 0, 1, 2\n")
cat("嵌入维度 E:", E_for_ccm, "\n")
cat("库长度范围:", libSizes_str, "\n\n")

# 存储CCM结果
ccm_results_list <- list()
ccm_summary_list <- list()

for (tp in 0:2) {
  cat("正在分析 Tp =", tp, "...\n")
  
  # 调整最大库长度（因为滞后会损失数据点）
  adjusted_max_lib <- max_available_lib - tp
  adjusted_lib_end <- min(lib_end, adjusted_max_lib)
  adjusted_libSizes_str <- paste(lib_start, adjusted_lib_end, lib_step)
  
  tryCatch({
    ccm_result <- CCM(
      dataFrame = single_station_data,
      E = E_for_ccm,
      Tp = tp,
      columns = "temp_anomaly",
      target = "sif",
      libSizes = adjusted_libSizes_str,
      sample = 100,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 保存原始结果
    ccm_results_list[[paste0("Tp", tp)]] <- ccm_result
    
    # 提取汇总
    ccm_summary <- ccm_result %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Tp = tp)
    
    ccm_summary_list[[paste0("Tp", tp)]] <- ccm_summary
    
  }, error = function(e) {
    cat("  ✗ Tp =", tp, "分析失败:", e$message, "\n")
  })
}

# 合并所有Tp的汇总数据
ccm_all_summary <- bind_rows(ccm_summary_list)

cat("\n✓ CCM分析完成\n\n")

# ============================================================================
# 5. 相关性分析：测试不同滞后期的相关性
# ============================================================================

cat("=== 相关性分析 ===\n\n")

correlation_results <- list()

# Tp = 0 (同期)
cor_tp0 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag0)
correlation_results[["Tp0"]] <- tibble(
  Tp = 0,
  r = cor_tp0$estimate,
  p_value = cor_tp0$p.value,
  ci_lower = cor_tp0$conf.int[1],
  ci_upper = cor_tp0$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag0))
)

cat("Tp = 0 (同期: t月温度 vs t月SIF)\n")
cat("  Pearson r =", round(cor_tp0$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp0$p.value), "\n")
cat("  95% CI: [", round(cor_tp0$conf.int[1], 3), ",", 
    round(cor_tp0$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp0$estimate < 0, 
                 "✓ 负相关：温度↑ SIF↓", 
                 "✗ 正相关：温度↑ SIF↑"), "\n\n")

# Tp = 1 (滞后1月)
cor_tp1 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag1)
correlation_results[["Tp1"]] <- tibble(
  Tp = 1,
  r = cor_tp1$estimate,
  p_value = cor_tp1$p.value,
  ci_lower = cor_tp1$conf.int[1],
  ci_upper = cor_tp1$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag1))
)

cat("Tp = 1 (滞后1月: t月温度 vs t+1月SIF)\n")
cat("  Pearson r =", round(cor_tp1$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp1$p.value), "\n")
cat("  95% CI: [", round(cor_tp1$conf.int[1], 3), ",", 
    round(cor_tp1$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp1$estimate < 0, 
                 "✓ 负相关：t月温度↑ → t+1月SIF↓", 
                 "✗ 正相关：t月温度↑ → t+1月SIF↑"), "\n\n")

# Tp = 2 (滞后2月)
cor_tp2 <- cor.test(lag_data$temp_anomaly, lag_data$sif_lag2)
correlation_results[["Tp2"]] <- tibble(
  Tp = 2,
  r = cor_tp2$estimate,
  p_value = cor_tp2$p.value,
  ci_lower = cor_tp2$conf.int[1],
  ci_upper = cor_tp2$conf.int[2],
  n = sum(!is.na(lag_data$temp_anomaly) & !is.na(lag_data$sif_lag2))
)

cat("Tp = 2 (滞后2月: t月温度 vs t+2月SIF)\n")
cat("  Pearson r =", round(cor_tp2$estimate, 3), "\n")
cat("  p-value =", format.pval(cor_tp2$p.value), "\n")
cat("  95% CI: [", round(cor_tp2$conf.int[1], 3), ",", 
    round(cor_tp2$conf.int[2], 3), "]\n")
cat("  ", ifelse(cor_tp2$estimate < 0, 
                 "✓ 负相关：t月温度↑ → t+2月SIF↓", 
                 "✗ 正相关：t月温度↑ → t+2月SIF↑"), "\n\n")

# 合并相关性结果
correlation_summary <- bind_rows(correlation_results)

# ============================================================================
# 6. 提取CCM的最终ρ值（用于综合报告）
# ============================================================================

ccm_final_rho <- ccm_all_summary %>%
  group_by(Tp) %>%
  filter(lib_size == max(lib_size)) %>%
  summarise(
    rho_final = mean(rho_mean),
    .groups = "drop"
  )

# 计算CCM趋势（收敛性）
ccm_trends <- ccm_all_summary %>%
  group_by(Tp) %>%
  summarise(
    trend = cor(lib_size, rho_mean),
    .groups = "drop"
  )

ccm_final_summary <- ccm_final_rho %>%
  left_join(ccm_trends, by = "Tp")

# ============================================================================
# 7. 综合结果表格
# ============================================================================

comprehensive_results <- correlation_summary %>%
  left_join(ccm_final_summary, by = "Tp") %>%
  mutate(
    causality = case_when(
      rho_final > 0.15 & trend > 0 ~ "有因果",
      TRUE ~ "无因果"
    ),
    effect_direction = case_when(
      r < 0 ~ "负向(温度↑ SIF↓)",
      r > 0 ~ "正向(温度↑ SIF↑)",
      TRUE ~ "无关"
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("=== 综合结果汇总 ===\n\n")
print(comprehensive_results %>%
        select(Tp, r, p_value, significance, rho_final, trend, 
               causality, effect_direction))

# ============================================================================
# 8. 可视化：相关性散点图
# ============================================================================

# 准备绘图数据
plot_data_tp0 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag0) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=0 (同期)")

plot_data_tp1 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag1) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=1 (滞后1月)")

plot_data_tp2 <- lag_data %>%
  select(temp_anomaly, sif = sif_lag2) %>%
  filter(!is.na(sif)) %>%
  mutate(Tp = "Tp=2 (滞后2月)")

plot_data_all <- bind_rows(plot_data_tp0, plot_data_tp1, plot_data_tp2) %>%
  mutate(Tp = factor(Tp, levels = c("Tp=0 (同期)", 
                                    "Tp=1 (滞后1月)", 
                                    "Tp=2 (滞后2月)")))

# 创建散点图
p_scatter <- ggplot(plot_data_all, aes(x = temp_anomaly, y = sif)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  facet_wrap(~ Tp, ncol = 3) +
  labs(
    title = paste("站点", example_station, "- 相关性分析"),
    subtitle = "温度异常 vs SIF（不同滞后期）",
    x = "温度异常 (°C)",
    y = "SIF",
    caption = "蓝线：线性拟合；阴影：95%置信区间"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11)
  )

# 添加相关系数标注
correlation_labels <- comprehensive_results %>%
  mutate(
    label = paste0("r = ", round(r, 3), significance, 
                   "\np = ", format.pval(p_value, digits = 2)),
    Tp = factor(paste0("Tp=", Tp, " (", 
                       c("同期", "滞后1月", "滞后2月")[Tp + 1], ")"),
                levels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)"))
  )

p_scatter_final <- p_scatter +
  geom_text(data = correlation_labels,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, 
            color = "red", fontface = "bold")

print(p_scatter_final)

# ============================================================================
# 9. 可视化：CCM收敛图
# ============================================================================

# CCM图：显示收敛性
p_ccm <- ggplot(ccm_all_summary, aes(x = lib_size, y = rho_mean, 
                                     color = factor(Tp))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = rho_mean - rho_sd, 
                  ymax = rho_mean + rho_sd,
                  fill = factor(Tp)),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.15, linetype = "dashed", 
             color = "gray", alpha = 0.7) +
  scale_color_manual(
    values = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A"),
    labels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)")
  ) +
  scale_fill_manual(
    values = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A"),
    labels = c("Tp=0 (同期)", "Tp=1 (滞后1月)", "Tp=2 (滞后2月)")
  ) +
  labs(
    title = paste("站点", example_station, "- CCM因果分析"),
    subtitle = "温度异常 → SIF（不同滞后期的因果强度）",
    x = "Library Size",
    y = "Cross Map Skill (ρ)",
    color = "滞后期",
    fill = "滞后期",
    caption = "ρ > 0.15 且上升趋势 → 有因果关系"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_ccm)


# ============================================================================
# 10. 可视化：相关系数和CCM ρ随滞后期变化
# ============================================================================

# 准备综合对比数据
comparison_data <- comprehensive_results %>%
  select(Tp, r, rho_final) %>%
  pivot_longer(cols = c(r, rho_final),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    "r" = "相关系数 (Pearson r)",
                    "rho_final" = "CCM因果强度 (ρ)")
  )

p_comparison <- ggplot(comparison_data, aes(x = Tp, y = value, 
                                            color = metric, 
                                            group = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0.15, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = 0:2, 
                     labels = c("0\n(同期)", "1\n(滞后1月)", "2\n(滞后2月)")) +
  scale_color_manual(values = c("相关系数 (Pearson r)" = "#E41A1C",
                                "CCM因果强度 (ρ)" = "#377EB8")) +
  labs(
    title = paste("站点", example_station, "- 滞后效应综合对比"),
    subtitle = "相关性 vs 因果性随滞后期的变化",
    x = "滞后期 (月)",
    y = "强度",
    color = "指标",
    caption = "红色虚线(0.15)：CCM因果显著性阈值"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_comparison)


# ============================================================================
# 11. 创建组合图
# ============================================================================

combined_plot <- (p_scatter_final) / (p_ccm + p_comparison) +
  plot_annotation(
    title = paste("站点", example_station, 
                  "- 温度异常对SIF的滞后因果效应完整分析"),
    subtitle = "上：相关性分析（因果性质） | 左下：CCM收敛性（因果方向） | 右下：综合对比",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)


# ============================================================================
# 12. 最终结论报告
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("                   最终因果关系结论\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("站点:", example_station, "\n")
cat("数据: 去趋势后的6-8月温度异常和SIF\n\n")

for (tp in 0:2) {
  result <- comprehensive_results %>% filter(Tp == tp)
  
  cat("【Tp =", tp, "】", 
      c("同期效应", "滞后1月效应", "滞后2月效应")[tp + 1], "\n")
  
  # 因果方向（CCM）
  if (result$causality == "有因果") {
    cat("  ✓ 因果关系: 温度异常 → SIF (ρ =", 
        round(result$rho_final, 3), ")\n")
  } else {
    cat("  ✗ 无显著因果 (ρ =", round(result$rho_final, 3), ")\n")
  }
  
  # 因果性质（相关性）
  if (result$p_value < 0.05) {
    if (result$r < 0) {
      cat("  ✓ 效应方向: 负向影响 - 温度异常↑ → SIF↓\n")
      cat("    (r =", round(result$r, 3), result$significance, ")\n")
      cat("  → 符合研究预期！\n")
    } else {
      cat("  ✗ 效应方向: 正向影响 - 温度异常↑ → SIF↑\n")
      cat("    (r =", round(result$r, 3), result$significance, ")\n")
      cat("  → 不符合研究预期\n")
    }
  } else {
    cat("  - 相关性不显著 (r =", round(result$r, 3), ", p > 0.05)\n")
  }
  
  cat("\n")
}

# 找出最强的滞后效应
best_lag <- comprehensive_results %>%
  filter(causality == "有因果", r < 0, p_value < 0.05) %>%
  arrange(desc(abs(r))) %>%
  slice(1)

if (nrow(best_lag) > 0) {
  cat("*** 最佳滞后期 ***\n")
  cat("Tp =", best_lag$Tp, "月滞后效应最强\n")
  cat("- CCM因果强度: ρ =", round(best_lag$rho_final, 3), "\n")
  cat("- 负向影响强度: r =", round(best_lag$r, 3), 
      best_lag$significance, "\n")
  cat("→ t月温度异常每升高1°C，t+", best_lag$Tp, 
      "月SIF预期降低约", round(abs(best_lag$r), 3), "个标准差\n")
} else {
  cat("*** 警告 ***\n")
  cat("未发现显著的负向滞后因果效应！\n")
  cat("建议：\n")
  cat("1. 检查是否只分析了真正的热月数据\n")
  cat("2. 尝试更长的滞后期（Tp = 3, 4, 5...）\n")
  cat("3. 检查数据质量和站点选择\n")
}

cat("\n", rep("=", 70), "\n", sep = "")


# ============================================================================
# 2. 批量分析所有站点（修正版）
# ============================================================================
# 函数：对单个站点进行CCM分析。
# ============================================================================
# 2. 批量分析所有站点（增强版 - 保存原始结果和图表）
# ============================================================================

perform_ccm_analysis_enhanced <- function(station_id, data, min_points = 30) {
  
  # 提取该站点数据
  station_data <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    select(temp_anomaly, sif) %>%
    as.data.frame()
  
  n_data <- nrow(station_data)
  
  # 检查数据点数
  if (n_data < min_points) {
    return(NULL)
  }
  
  tryCatch({
    # 添加时间列
    station_data <- data.frame(
      time = 1:n_data,
      temp_anomaly = station_data$temp_anomaly,
      sif = station_data$sif
    )
    
    lib_pred_size <- min(100, n_data)
    
    # 确定最优E
    embed_temp <- EmbedDimension(
      dataFrame = station_data,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "temp_anomaly",
      target = "temp_anomaly",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    embed_sif <- EmbedDimension(
      dataFrame = station_data,
      lib = paste("1", lib_pred_size),
      pred = paste("1", lib_pred_size),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    E_ccm <- max(best_E_temp, best_E_sif)
    
    # 计算实际可用的最大库长度
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    max_available_lib <- n_data - embedding_loss
    
    # 设置库长度
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    # 检查是否有足够的数据进行分析
    if (lib_end <= lib_start || lib_end < 15) {
      return(NULL)
    }
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # CCM分析
    ccm_result <- CCM(
      dataFrame = station_data,
      E = E_ccm,
      Tp = 0,
      columns = "temp_anomaly",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 提取汇总结果
    ccm_summary <- ccm_result %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_temp_xmap_sif_mean = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_temp_xmap_sif_sd = sd(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sif_xmap_temp_mean = mean(`sif:temp_anomaly`, na.rm = TRUE),
        rho_sif_xmap_temp_sd = sd(`sif:temp_anomaly`, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 获取最终rho值
    final_rho_temp_xmap_sif <- ccm_summary$rho_temp_xmap_sif_mean[nrow(ccm_summary)]
    final_rho_sif_xmap_temp <- ccm_summary$rho_sif_xmap_temp_mean[nrow(ccm_summary)]
    
    # 计算趋势
    trend_temp_xmap_sif <- cor(ccm_summary$lib_size, 
                               ccm_summary$rho_temp_xmap_sif_mean)
    trend_sif_xmap_temp <- cor(ccm_summary$lib_size, 
                               ccm_summary$rho_sif_xmap_temp_mean)
    
    # 创建CCM图表
    ccm_plot <- ggplot(ccm_summary) +
      geom_line(aes(x = lib_size, y = rho_temp_xmap_sif_mean, 
                    color = "Temp xmap SIF"), 
                linewidth = 1) +
      geom_point(aes(x = lib_size, y = rho_temp_xmap_sif_mean, 
                     color = "Temp xmap SIF"), 
                 size = 2) +
      geom_ribbon(aes(x = lib_size, 
                      ymin = rho_temp_xmap_sif_mean - rho_temp_xmap_sif_sd,
                      ymax = rho_temp_xmap_sif_mean + rho_temp_xmap_sif_sd,
                      fill = "Temp xmap SIF"),
                  alpha = 0.2) +
      geom_line(aes(x = lib_size, y = rho_sif_xmap_temp_mean, 
                    color = "SIF xmap Temp"), 
                linewidth = 1) +
      geom_point(aes(x = lib_size, y = rho_sif_xmap_temp_mean, 
                     color = "SIF xmap Temp"), 
                 size = 2) +
      geom_ribbon(aes(x = lib_size, 
                      ymin = rho_sif_xmap_temp_mean - rho_sif_xmap_temp_sd,
                      ymax = rho_sif_xmap_temp_mean + rho_sif_xmap_temp_sd,
                      fill = "SIF xmap Temp"),
                  alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", 
                 alpha = 0.3) +
      scale_color_manual(values = c("Temp xmap SIF" = "red", 
                                    "SIF xmap Temp" = "blue")) +
      scale_fill_manual(values = c("Temp xmap SIF" = "red", 
                                   "SIF xmap Temp" = "blue")) +
      labs(
        title = paste("CCM分析 - 站点", station_id),
        subtitle = paste0("E = ", E_ccm, " | n = ", n_data, 
                          " | ρ(T→S) = ", round(final_rho_sif_xmap_temp, 2),
                          " | ρ(S→T) = ", round(final_rho_temp_xmap_sif, 2)),
        x = "Library Size",
        y = "Cross Map Skill (ρ)",
        color = "Direction",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # 返回增强结果（包含原始数据和图表）
    tibble(
      meteo_stat_id = station_id,
      n_points = n_data,
      max_lib_size = max_available_lib,
      E = E_ccm,
      best_E_temp = best_E_temp,
      best_E_sif = best_E_sif,
      rho_temp_xmap_sif = final_rho_temp_xmap_sif,
      trend_temp_xmap_sif = trend_temp_xmap_sif,
      rho_sif_xmap_temp = final_rho_sif_xmap_temp,
      trend_sif_xmap_temp = trend_sif_xmap_temp,
      # 嵌套列：保存原始CCM结果
      ccm_raw_data = list(ccm_result),
      # 嵌套列：保存汇总数据
      ccm_summary_data = list(ccm_summary),
      # 嵌套列：保存图表对象
      ccm_plot = list(ccm_plot)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# 获取所有站点
all_stations <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  pull(meteo_stat_id)

cat("\n准备分析", length(all_stations), "个站点...\n")

# 批量分析
library(progress)
pb <- progress_bar$new(
  format = "  CCM分析 [:bar] :percent eta: :eta",
  total = length(all_stations), 
  clear = FALSE
)

ccm_results_all <- map_dfr(all_stations, function(sid) {
  pb$tick()
  perform_ccm_analysis_enhanced(sid, data_for_ccm, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_all), "/", length(all_stations), "\n")

# 保存结果（注意：列表列需要特殊处理）
# 保存基本结果（不含列表列）
ccm_results_basic <- ccm_results_all %>%
  select(-ccm_raw_data, -ccm_summary_data, -ccm_plot)

write_csv(ccm_results_basic, "data_proc/ccm_results.csv")

# 保存完整结果为RDS格式（可以保存列表列）
saveRDS(ccm_results_all, "data_proc/ccm_results_full.rds")

cat("✓ 基本结果已保存到: data_proc/ccm_results.csv\n")
cat("✓ 完整结果已保存到: data_proc/ccm_results_full.rds\n")

# ============================================================================
# 3. 结果分析和可视化
# ============================================================================

# 因果关系分类
ccm_summary_final <- ccm_results_all %>%
  mutate(
    sif_causes_temp = rho_temp_xmap_sif > 0.1 & trend_temp_xmap_sif > 0,
    temp_causes_sif = rho_sif_xmap_temp > 0.1 & trend_sif_xmap_temp > 0,
    causality_type = case_when(
      temp_causes_sif & !sif_causes_temp ~ "Temp → SIF",
      !temp_causes_sif & sif_causes_temp ~ "SIF → Temp",
      temp_causes_sif & sif_causes_temp ~ "双向因果",
      TRUE ~ "无显著因果"
    )
  )

# 统计
cat("\n=== CCM分析总结 ===\n")
cat("成功分析的站点数:", nrow(ccm_summary_final), "\n\n")
print(table(ccm_summary_final$causality_type))

# 汇总图
ggplot(ccm_summary_final, aes(x = rho_temp_xmap_sif, y = rho_sif_xmap_temp)) +
  geom_point(aes(color = causality_type), size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(
    "Temp → SIF" = "red",
    "SIF → Temp" = "blue",
    "双向因果" = "purple",
    "无显著因果" = "gray"
  )) +
  labs(
    title = "CCM分析结果：热月温度异常与SIF的因果关系",
    subtitle = "虚线表示 ρ = 0.1 的阈值",
    x = "ρ (Temp xmap SIF)\n[SIF影响Temp的证据]",
    y = "ρ (SIF xmap Temp)\n[Temp影响SIF的证据]",
    color = "因果关系类型"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# ============================================================================
# 4. 查看和导出特定站点的CCM结果和图表
# ============================================================================

# 4.1 查看双向因果的站点
bidirectional_stations <- ccm_summary_final %>%
  filter(causality_type == "双向因果") %>%
  arrange(desc(rho_temp_xmap_sif + rho_sif_xmap_temp))

cat("\n双向因果站点数:", nrow(bidirectional_stations), "\n")
print(head(bidirectional_stations %>% 
             select(meteo_stat_id, n_points, E, 
                    rho_temp_xmap_sif, rho_sif_xmap_temp,
                    trend_temp_xmap_sif, trend_sif_xmap_temp), 
           10))

# 4.2 查看单个站点的原始CCM数据
view_station_ccm <- function(station_id, results_df) {
  station_row <- results_df %>% filter(meteo_stat_id == station_id)
  
  if (nrow(station_row) == 0) {
    cat("站点", station_id, "未找到\n")
    return(NULL)
  }
  
  cat("\n=== 站点", station_id, "的CCM详细结果 ===\n")
  cat("数据点数:", station_row$n_points, "\n")
  cat("嵌入维度 E:", station_row$E, "\n")
  cat("最大库长度:", station_row$max_lib_size, "\n\n")
  
  # 显示汇总数据
  cat("CCM汇总数据:\n")
  print(station_row$ccm_summary_data[[1]])
  
  # 显示图表
  cat("\n显示CCM图表...\n")
  print(station_row$ccm_plot[[1]])
  
  # 返回原始数据（如果需要进一步分析）
  return(list(
    summary = station_row$ccm_summary_data[[1]],
    raw_data = station_row$ccm_raw_data[[1]],
    plot = station_row$ccm_plot[[1]]
  ))
}

# 示例：查看第一个双向因果站点
if (nrow(bidirectional_stations) > 0) {
  example_bidir_id <- bidirectional_stations$meteo_stat_id[1]
  station_details <- view_station_ccm(example_bidir_id, ccm_results_all)
}

# 4.3 批量导出双向因果站点的图表
export_bidirectional_plots <- function(results_df, output_dir = "data_proc/ccm_plots") {
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 筛选双向因果站点
  bidir_results <- results_df %>%
    mutate(
      sif_causes_temp = rho_temp_xmap_sif > 0.1 & trend_temp_xmap_sif > 0,
      temp_causes_sif = rho_sif_xmap_temp > 0.1 & trend_sif_xmap_temp > 0
    ) %>%
    filter(sif_causes_temp & temp_causes_sif)
  
  cat("\n导出", nrow(bidir_results), "个双向因果站点的图表...\n")
  
  pb <- progress_bar$new(
    format = "  导出图表 [:bar] :percent",
    total = nrow(bidir_results)
  )
  
  for (i in 1:nrow(bidir_results)) {
    pb$tick()
    
    station_id <- bidir_results$meteo_stat_id[i]
    plot_obj <- bidir_results$ccm_plot[[i]]
    
    # 保存为PNG
    filename <- file.path(output_dir, paste0("ccm_", station_id, ".png"))
    ggsave(filename, plot_obj, width = 8, height = 6, dpi = 300)
  }
  
  cat("✓ 图表已保存到:", output_dir, "\n")
}

# 导出所有双向因果站点的图表
export_bidirectional_plots(ccm_results_all)

# 4.4 创建多站点对比图（显示前9个双向因果站点）
create_multipanel_plot <- function(results_df, n_plots = 9) {
  
  # 筛选双向因果站点
  bidir_results <- results_df %>%
    mutate(
      sif_causes_temp = rho_temp_xmap_sif > 0.1 & trend_temp_xmap_sif > 0,
      temp_causes_sif = rho_sif_xmap_temp > 0.1 & trend_sif_xmap_temp > 0
    ) %>%
    filter(sif_causes_temp & temp_causes_sif) %>%
    arrange(desc(rho_temp_xmap_sif + rho_sif_xmap_temp)) %>%
    head(n_plots)
  
  if (nrow(bidir_results) == 0) {
    cat("没有双向因果的站点\n")
    return(NULL)
  }
  
  # 合并所有站点的汇总数据
  combined_data <- map_dfr(1:nrow(bidir_results), function(i) {
    bidir_results$ccm_summary_data[[i]] %>%
      mutate(meteo_stat_id = bidir_results$meteo_stat_id[i])
  })
  
  # 转换为长格式
  plot_data <- combined_data %>%
    pivot_longer(
      cols = c(rho_temp_xmap_sif_mean, rho_sif_xmap_temp_mean),
      names_to = "direction",
      values_to = "rho"
    ) %>%
    mutate(
      direction = recode(direction,
                         "rho_temp_xmap_sif_mean" = "Temp xmap SIF (S→T)",
                         "rho_sif_xmap_temp_mean" = "SIF xmap Temp (T→S)")
    )
  
  # 创建分面图
  p <- ggplot(plot_data, aes(x = lib_size, y = rho, color = direction)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0.1, linetype = "dashed", 
               color = "red", alpha = 0.5) +
    scale_color_manual(values = c(
      "Temp xmap SIF (S→T)" = "red",
      "SIF xmap Temp (T→S)" = "blue"
    )) +
    facet_wrap(~ meteo_stat_id, scales = "free", ncol = 3) +
    labs(
      title = "双向因果站点CCM分析对比",
      subtitle = paste("前", nrow(bidir_results), "个ρ总和最高的站点"),
      x = "Library Size",
      y = "Cross Map Skill (ρ)",
      color = "方向"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 9)
    )
  
  return(p)
}

# 创建多站点对比图
multipanel_plot <- create_multipanel_plot(ccm_results_all, n_plots = 9)
if (!is.null(multipanel_plot)) {
  print(multipanel_plot)
  ggsave("data_proc/ccm_bidirectional_comparison.png", 
         multipanel_plot, width = 12, height = 10, dpi = 300)
}

# 4.5 创建汇总表格
create_summary_table <- function(results_df) {
  results_df %>%
    mutate(
      sif_causes_temp = rho_temp_xmap_sif > 0.1 & trend_temp_xmap_sif > 0,
      temp_causes_sif = rho_sif_xmap_temp > 0.1 & trend_sif_xmap_temp > 0,
      causality_type = case_when(
        temp_causes_sif & !sif_causes_temp ~ "Temp → SIF",
        !temp_causes_sif & sif_causes_temp ~ "SIF → Temp",
        temp_causes_sif & sif_causes_temp ~ "双向因果",
        TRUE ~ "无显著因果"
      ),
      rho_diff = rho_sif_xmap_temp - rho_temp_xmap_sif,
      dominant_direction = case_when(
        causality_type == "双向因果" & abs(rho_diff) < 0.05 ~ "均衡",
        causality_type == "双向因果" & rho_diff > 0.05 ~ "T→S主导",
        causality_type == "双向因果" & rho_diff < -0.05 ~ "S→T主导",
        TRUE ~ as.character(causality_type)
      )
    ) %>%
    select(meteo_stat_id, n_points, E, causality_type, dominant_direction,
           rho_temp_xmap_sif, rho_sif_xmap_temp, rho_diff,
           trend_temp_xmap_sif, trend_sif_xmap_temp) %>%
    arrange(desc(causality_type == "双向因果"), desc(abs(rho_diff)))
}

summary_table <- create_summary_table(ccm_results_all)
write_csv(summary_table, "data_proc/ccm_summary_table.csv")

cat("\n✓ 汇总表格已保存到: data_proc/ccm_summary_table.csv\n")

# 显示双向因果的汇总
cat("\n双向因果站点详细信息:\n")
print(summary_table %>% filter(causality_type == "双向因果") %>% head(20))
