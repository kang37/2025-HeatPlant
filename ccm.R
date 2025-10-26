# ============================================================================
# CCM分析准备（修正版 - 解决库长度问题）
# ============================================================================

library(rEDM)

# ============================================================================
# 1. 单个站点的CCM分析示例
# ============================================================================

# 选择一个数据点充足的站点
example_station <- data_for_ccm %>%
  group_by(meteo_stat_id) %>%
  summarise(n = n()) %>%
  filter(n >= 30) %>%
  slice(1) %>%
  pull(meteo_stat_id)

# 提取单个站点数据，必须是data.frame格式
single_station_data <- data_for_ccm %>%
  filter(meteo_stat_id == example_station) %>%
  arrange(year, month) %>%
  select(temp_anomaly, sif) %>%
  as.data.frame()

# 添加时间索引列（rEDM要求第一列是时间）
single_station_data <- data.frame(
  time = 1:nrow(single_station_data),
  temp_anomaly = single_station_data$temp_anomaly,
  sif = single_station_data$sif
)

cat("示例站点:", example_station, "\n")
cat("数据点数:", nrow(single_station_data), "\n")

n_data <- nrow(single_station_data)
lib_pred_size <- min(100, n_data)

# 1.1 确定最优嵌入维度
embed_temp <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  maxE = 10,
  showPlot = TRUE
)

embed_sif <- EmbedDimension(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "sif",
  target = "sif",
  maxE = 10,
  showPlot = TRUE
)

# 选择最优E
best_E_temp <- embed_temp$E[which.max(embed_temp$rho)]
best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]

cat("\nTemperature Anomaly 最优 E:", best_E_temp, 
    "(ρ =", round(max(embed_temp$rho), 3), ")\n")
cat("SIF 最优 E:", best_E_sif, 
    "(ρ =", round(max(embed_sif$rho), 3), ")\n")

# 1.2 检测非线性
nonlinear_temp <- PredictNonlinear(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "temp_anomaly",
  target = "temp_anomaly",
  E = best_E_temp,
  showPlot = TRUE
)

nonlinear_sif <- PredictNonlinear(
  dataFrame = single_station_data,
  lib = paste("1", lib_pred_size),
  pred = paste("1", lib_pred_size),
  columns = "sif",
  target = "sif",
  E = best_E_sif,
  showPlot = TRUE
)

# 1.3 进行CCM分析
E_for_ccm <- max(best_E_temp, best_E_sif)

# 关键修改：计算实际可用的最大库长度
# 嵌入后会损失 (E-1)*tau 个数据点
tau <- 1  # 默认时间延迟
embedding_loss <- (E_for_ccm - 1) * tau
max_available_lib <- n_data - embedding_loss

cat("\n数据点信息:\n")
cat("原始数据点数:", n_data, "\n")
cat("嵌入维度 E:", E_for_ccm, "\n")
cat("嵌入损失:", embedding_loss, "\n")
cat("最大可用库长度:", max_available_lib, "\n")

# 设置库长度范围（确保不超过最大可用长度）
lib_start <- max(E_for_ccm + 2, 10)
lib_end <- max_available_lib  # 使用实际可用的最大长度
lib_step <- max(2, floor((lib_end - lib_start) / 15))

# 确保至少有足够的库长度进行分析
if (lib_end <= lib_start) {
  stop("数据点数不足以进行CCM分析")
}

libSizes_str <- paste(lib_start, lib_end, lib_step)

cat("\n开始CCM分析...\n")
cat("使用 E =", E_for_ccm, "\n")
cat("Library sizes:", libSizes_str, "\n")

# CCM分析
ccm_result <- CCM(
  dataFrame = single_station_data,
  E = E_for_ccm,
  Tp = 0,
  columns = "temp_anomaly",
  target = "sif",
  libSizes = libSizes_str,
  sample = 100,
  random = TRUE,
  showPlot = FALSE
)

# 查看结果
cat("\nCCM结果列名:", colnames(ccm_result), "\n")
head(ccm_result)

# 提取结果
ccm_summary <- ccm_result %>%
  rename(lib_size = LibSize) %>%
  group_by(lib_size) %>%
  summarise(
    rho_temp_xmap_sif = mean(`temp_anomaly:sif`, na.rm = TRUE),
    rho_sif_xmap_temp = mean(`sif:temp_anomaly`, na.rm = TRUE),
    .groups = "drop"
  )

# 可视化CCM结果
ggplot(ccm_summary) +
  geom_line(aes(x = lib_size, y = rho_temp_xmap_sif, color = "Temp xmap SIF"), 
            linewidth = 1) +
  geom_point(aes(x = lib_size, y = rho_temp_xmap_sif, color = "Temp xmap SIF"), 
             size = 2) +
  geom_line(aes(x = lib_size, y = rho_sif_xmap_temp, color = "SIF xmap Temp"), 
            linewidth = 1) +
  geom_point(aes(x = lib_size, y = rho_sif_xmap_temp, color = "SIF xmap Temp"), 
             size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", alpha = 0.3) +
  scale_color_manual(values = c("Temp xmap SIF" = "red", 
                                "SIF xmap Temp" = "blue")) +
  labs(
    title = paste("CCM分析 - 站点", example_station),
    subtitle = paste("E =", E_for_ccm, "| 数据点数 =", n_data),
    x = "Library Size",
    y = "Cross Map Skill (ρ)",
    color = "Direction"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 因果关系判断
final_rho_temp_xmap_sif <- ccm_summary$rho_temp_xmap_sif[nrow(ccm_summary)]
final_rho_sif_xmap_temp <- ccm_summary$rho_sif_xmap_temp[nrow(ccm_summary)]

trend_temp_xmap_sif <- cor(ccm_summary$lib_size, ccm_summary$rho_temp_xmap_sif)
trend_sif_xmap_temp <- cor(ccm_summary$lib_size, ccm_summary$rho_sif_xmap_temp)

cat("\n=== CCM结果解读 ===\n")
cat("站点:", example_station, "\n")
cat("数据点数:", n_data, "\n")
cat("嵌入维度 E:", E_for_ccm, "\n\n")

cat("Temp xmap SIF:\n")
cat("  最终 ρ:", round(final_rho_temp_xmap_sif, 3), "\n")
cat("  趋势（收敛性）:", round(trend_temp_xmap_sif, 3), "\n")
cat("  → ", ifelse(final_rho_temp_xmap_sif > 0.2 & trend_temp_xmap_sif > 0,
                   "✓ SIF 影响 Temp", "✗ 无强证据"), "\n\n")

cat("SIF xmap Temp:\n")
cat("  最终 ρ:", round(final_rho_sif_xmap_temp, 3), "\n")
cat("  趋势（收敛性）:", round(trend_sif_xmap_temp, 3), "\n")
cat("  → ", ifelse(final_rho_sif_xmap_temp > 0.2 & trend_sif_xmap_temp > 0,
                   "✓ Temp 影响 SIF", "✗ 无强证据"), "\n")

# ============================================================================
# 2. 批量分析所有站点（修正版）
# ============================================================================

perform_ccm_analysis_batch <- function(station_id, data, min_points = 30) {
  
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
    
    # 关键修改：计算实际可用的最大库长度
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    max_available_lib <- n_data - embedding_loss
    
    # 设置库长度（确保不超过最大可用长度）
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib  # 使用实际可用的最大长度
    
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
    
    # 提取结果
    ccm_summary <- ccm_result %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_temp_xmap_sif = mean(`temp_anomaly:sif`, na.rm = TRUE),
        rho_sif_xmap_temp = mean(`sif:temp_anomaly`, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 获取最终rho值
    final_rho_temp_xmap_sif <- ccm_summary$rho_temp_xmap_sif[nrow(ccm_summary)]
    final_rho_sif_xmap_temp <- ccm_summary$rho_sif_xmap_temp[nrow(ccm_summary)]
    
    # 计算趋势
    trend_temp_xmap_sif <- cor(ccm_summary$lib_size, 
                               ccm_summary$rho_temp_xmap_sif)
    trend_sif_xmap_temp <- cor(ccm_summary$lib_size, 
                               ccm_summary$rho_sif_xmap_temp)
    
    # 返回结果
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
      trend_sif_xmap_temp = trend_sif_xmap_temp
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
  perform_ccm_analysis_batch(sid, data_for_ccm, min_points = 30)
})

cat("\n成功分析的站点数:", nrow(ccm_results_all), "/", length(all_stations), "\n")

# 保存结果
write_csv(ccm_results_all, "data_proc/ccm_results.csv")

# ============================================================================
# 3. 结果分析和可视化
# ============================================================================

# 因果关系分类
ccm_summary_final <- ccm_results_all %>%
  mutate(
    # 判断因果关系（ρ > 0.1 且 趋势 > 0）
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

# 可视化
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

# 查看结果
head(ccm_summary_final, 20)
