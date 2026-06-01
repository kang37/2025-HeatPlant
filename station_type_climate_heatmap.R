pacman::p_load(dplyr, ggplot2, tidyr, purrr, targets)

setwd("/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant")

results <- readRDS("data_proc/results_weekly_0_5.rds")

# =============================================================================
# 1. 站点分类：4 类
# =============================================================================
classify_station <- function(coefs) {
  # coefs: tp=0..5 排序的 mean_coef 向量
  tp0_pos <- coefs[1] > 0   # tp=0 促进
  any_switch_to_neg <- any(coefs[-1] < 0) & mean(coefs[-1] < 0) >= 0.5
  any_switch_to_pos <- any(coefs[-1] > 0) & mean(coefs[-1] > 0) >= 0.5

  if ( tp0_pos & !any_switch_to_neg) return("全程促进")
  if (!tp0_pos & !any_switch_to_pos) return("全程抑制")
  if ( tp0_pos &  any_switch_to_neg) return("促进→抑制")
  if (!tp0_pos &  any_switch_to_pos) return("抑制→促进")
  return("全程抑制")  # 兜底
}

station_type <- results %>%
  arrange(meteo_stat_id, tp) %>%
  group_by(meteo_stat_id) %>%
  summarise(
    coef_seq     = list(mean_coef),
    koppen_class = first(koppen_class),
    koppen_group = first(koppen_group),
    .groups = "drop"
  ) %>%
  filter(map_lgl(coef_seq, ~length(.x) == 6)) %>%
  mutate(type = map_chr(coef_seq, classify_station))

cat("站点类型分布:\n")
print(count(station_type, type))

# =============================================================================
# 2. 气候类型聚合标签（保留二级，加描述）
# =============================================================================
climate_labels <- c(
  "Am"  = "Am\n热带季风",
  "As"  = "As\n热带草原",
  "BSh" = "BSh\n热半干旱",
  "BWh" = "BWh\n热荒漠",
  "BWk" = "BWk\n冷荒漠",
  "Csc" = "Csc\n温带地中海",
  "Cwa" = "Cwa\n温带季风",
  "Cwc" = "Cwc\n温带高原",
  "Dsd" = "Dsd\n大陆性干旱",
  "Dwa" = "Dwa\n大陆性季风",
  "Dwd" = "Dwd\n大陆严寒"
)

# 气候带分组颜色（用于行标签着色）
group_colors <- c(A = "#E41A1C", B = "#FF7F00", C = "#4DAF4A", D = "#377EB8")

# 类型顺序（逻辑排列）
type_order <- c("全程促进", "促进→抑制", "抑制→促进", "全程抑制")
type_colors <- c(
  "全程促进"   = "#1B7837",
  "促进→抑制"  = "#74C476",
  "抑制→促进"  = "#FC8D59",
  "全程抑制"   = "#B2182B"
)

# =============================================================================
# 3. 构建热力图数据
# =============================================================================
hm_data <- station_type %>%
  mutate(
    climate_label = climate_labels[koppen_class],
    type = factor(type, levels = type_order)
  ) %>%
  count(climate_label, koppen_group, type, .drop = FALSE) %>%
  group_by(climate_label) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) %>%
  ungroup() %>%
  filter(!is.na(climate_label))

# 气候带行排序：A → B → C → D，组内按总站点数降序
climate_order <- hm_data %>%
  distinct(climate_label, koppen_group, total) %>%
  arrange(koppen_group, desc(total)) %>%
  pull(climate_label)

hm_data <- hm_data %>%
  mutate(
    climate_label = factor(climate_label, levels = rev(climate_order)),
    type          = factor(type, levels = type_order)
  )

# =============================================================================
# 4. 热力图：填色 = 百分比，数字标注 = 数量
# =============================================================================
p_hm <- ggplot(hm_data, aes(x = type, y = climate_label, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(n > 0, sprintf("%d\n(%.0f%%)", n, pct), "")),
    size = 20, color = "white", fontface = "bold", lineheight = 1.1
  ) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#2171B5", "#08306B"),
    name = "占比 (%)",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "站点因果响应类型 × 柯本气候分区",
    subtitle = "格内数字：站点数（括号内为该气候带内占比）",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 30) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 20),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11),
    axis.text.x      = element_text(face = "bold", size = 11,
                                    color = type_colors[type_order]),
    axis.text.y      = element_text(size = 20),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.key.height = unit(1.2, "cm")
  )

# 在 y 轴标签旁加气候带色块（用 geom_tile 在负 x 位置）
# 提取行对应气候组颜色，用边栏注释
group_strip <- hm_data %>%
  distinct(climate_label, koppen_group) %>%
  mutate(x_strip = 0.35)  # 占位用

p_hm_final <- p_hm +
  # 右侧加气候大组标注
  geom_text(
    data = hm_data %>%
      group_by(koppen_group) %>%
      summarise(climate_label = last(levels(climate_label)[
        levels(climate_label) %in% unique(as.character(climate_label))
      ]), .groups = "drop") %>%
      mutate(type = factor(type_order[4], levels = type_order),
             label = paste0("（", koppen_group, "类）")),
    aes(x = type, y = climate_label, label = label),
    hjust = -0.05, size = 20, color = "grey50",
    inherit.aes = FALSE, nudge_x = 0.52
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 60, 5, 5))

ggsave("data_proc/station_type_climate_heatmap.png",
       p_hm_final, width = 10, height = 7, dpi = 300)
cat("-> 图已保存: data_proc/station_type_climate_heatmap.png\n")

# =============================================================================
# 5. 补充：百分比堆叠条形图（更直观比较各气候带的类型构成）
# =============================================================================
p_bar <- hm_data %>%
  mutate(climate_label = factor(climate_label, levels = climate_order)) %>%
  ggplot(aes(x = climate_label, y = pct, fill = type)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(pct >= 8, sprintf("%.0f%%", pct), "")),
    position = position_stack(vjust = 0.5),
    size = 20, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = type_colors, name = "响应类型") +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "各气候分区的站点响应类型构成",
    subtitle = "百分比 = 该气候带内各类型占比",
    x = "柯本气候类型", y = "站点占比 (%)"
  ) +
  theme_minimal(base_size = 30) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.text.x   = element_text(angle = 30, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave("data_proc/station_type_climate_barplot.png",
       p_bar, width = 11, height = 6, dpi = 300)
cat("-> 图已保存: data_proc/station_type_climate_barplot.png\n")

# =============================================================================
# 6. 数值汇总
# =============================================================================
cat("\n各气候带站点类型数量矩阵:\n")
station_type %>%
  mutate(type = factor(type, levels = type_order)) %>%
  count(koppen_class, type) %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  mutate(合计 = rowSums(across(where(is.numeric)))) %>%
  print()
