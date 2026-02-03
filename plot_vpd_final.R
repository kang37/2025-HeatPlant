# ============================================================================
# Script: plot_vpd_final.R
# Description: 从预处理好的RDS文件中读取数据，并生成一个字体更大、正确显示中文的图表。
# ============================================================================

# --- 0. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  ggplot2, dplyr, showtext, scales
)

# --- 1. 配置中文字体 ---
# 使用showtext包来解决ggplot中文乱码问题
# Noto Sans SC 是一个开源且高质量的中文字体
# 第一次运行时，font_add_google可能需要几秒钟来下载字体
cat("--> 正在配置中文字体...\n")
font_add_google("Noto Sans SC", "notosans-sc")
showtext_auto()

# --- 2. 读取预处理好的数据 ---
cat("--> 正在从 vpd_timeseries_agg.rds 读取数据...\n")
vpd_timeseries_agg <- readRDS("vpd_timeseries_agg.rds")

# --- 3. 可视化 (大字体版本) ---
cat("--> 正在生成最终图表...\n")

plot_vpd_final <- ggplot(vpd_timeseries_agg, aes(x = date, y = vpd_mean, group = 1)) +
  # 阴影区域
  geom_ribbon(aes(ymin = vpd_mean - vpd_sd, ymax = vpd_mean + vpd_sd, fill = region), alpha = 0.3) +
  # 均值折线 (加粗)
  geom_line(aes(color = region), linewidth = 0.8) +
  # 按区域分面
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  # 格式化X轴
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  # 添加标题和标签
  labs(
    title = "各区域VPD时间序列 (2000-2023)",
    subtitle = "实线为区域月平均VPD，阴影为区域内站点的标准差范围",
    x = "年份",
    y = "区域月平均 VPD (kPa)"
  ) +
  # 设置主题和字体大小
  theme_bw() +
  theme(
    # 指定所有文本元素使用刚刚加载的字体
    text = element_text(family = "notosans-sc"),
    plot.title = element_text(face = "bold", size = 100, hjust = 0.5),
    plot.subtitle = element_text(size = 100, hjust = 0.5, margin = margin(b = 15)),
    # 加大分面标题的字体
    strip.text = element_text(face = "bold", size = 18),
    # 加大坐标轴标签和刻度的字体
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100),
    legend.position = "none"
  )

# --- 4. 保存图表 ---
ggsave("vpd_timeseries_by_region_final.png", plot_vpd_final, width = 14, height = 22, dpi = 300, bg = "white")

