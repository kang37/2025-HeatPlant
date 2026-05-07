# 加载包。
library(targets)

# Set target options:
tar_option_set(
  packages = c(
    "tibble", "dplyr", "ggplot2", "lubridate", "purrr", "data.table", "stringr",
    "readr", "tidyr", "showtext", "rEDM", "sf", "rnaturalearth", 
    "rnaturalearthdata", "knitr", "readxl"
  )
)

# 自定义函数。
# 函数：读取每日气象数据。
read_one_meteo_daily <- function(path) {
  station_id <- str_extract(basename(path), "\\d+")
  read_csv(path, skip = 1, show_col_types = FALSE, na = c("", "NA")) %>%
    mutate(across(where(is.numeric), ~ ifelse(.x >= 999990, NA_real_, .x))) %>%
    mutate(meteo_stat_id = station_id, .before = 1) %>%
    return()
}

# 函数：安全去趋势
safe_detrend <- function(x, time_idx) {
  # 检查是否有足够的非NA值
  valid_data <- !is.na(x)
  n_valid <- sum(valid_data)
  
  if (n_valid < 3) {
    # 数据点太少，返回原值或NA
    return(rep(NA_real_, length(x)))
  }
  
  tryCatch({
    # 拟合线性模型
    model <- lm(x ~ time_idx)
    # 返回残差
    residuals(model)
  }, error = function(e) {
    # 如果出错，返回NA
    rep(NA_real_, length(x))
  })
}

# 函数：进行CCM分析。
perform_ccm_heat_sif <- function(
    station_id, data, min_points = 30, 
    sif_col = "sif_detrended", heat_col = "heat_index_composite_detrended", 
    tp_x = 0) {
  # 提取特定站点的数据。
  station_data_ccm <- data %>%
    filter(meteo_stat_id == station_id) %>%
    arrange(year, month) %>%
    group_by(year) %>% 
    mutate(heat_index = dplyr::lag(!!sym(heat_col), n = tp_x)) %>% 
    ungroup() %>% 
    select(sif = all_of(sif_col), heat_index) %>%
    filter(!is.na(sif), !is.na(heat_index)) %>%
    mutate(time = row_number(), .before = 1) %>%
    as.data.frame()
  
  # 如果样本不足，返回空值。
  n_data <- nrow(station_data_ccm)
  if (n_data < min_points) return(NULL)
  
  tryCatch({
    # 确定最优E
    embed_sif <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      columns = "sif",
      target = "sif",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    embed_heat <- EmbedDimension(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      columns = "heat_index",
      target = "heat_index",
      maxE = min(8, floor(n_data/10)),
      showPlot = FALSE
    )
    
    best_E_sif <- embed_sif$E[which.max(embed_sif$rho)]
    best_E_heat <- embed_heat$E[which.max(embed_heat$rho)]
    E_ccm <- max(best_E_sif, best_E_heat)
    
    # CCM参数
    tau <- 1
    embedding_loss <- (E_ccm - 1) * tau
    # Tp设置为0为同步，否则为滞后关系。
    tp <- tp_x
    max_available_lib <- n_data - embedding_loss - tp
    
    lib_start <- max(E_ccm + 2, 10)
    lib_end <- max_available_lib
    
    if (lib_end <= lib_start || lib_end < 15) return(NULL)
    
    lib_step <- max(2, floor((lib_end - lib_start) / 10))
    libSizes_str <- paste(lib_start, lib_end, lib_step)
    
    # CCM: 热事件指数 → SIF
    ccm_heat_to_sif <- CCM(
      dataFrame = station_data_ccm,
      E = E_ccm,
      Tp = tp,
      columns = "heat_index",
      target = "sif",
      libSizes = libSizes_str,
      sample = 50,
      random = TRUE,
      showPlot = FALSE
    )
    
    # 汇总
    ccm_summary_heat <- ccm_heat_to_sif %>%
      rename(lib_size = LibSize) %>%
      group_by(lib_size) %>%
      summarise(
        rho_mean = mean(`heat_index:sif`, na.rm = TRUE),
        rho_sd = sd(`heat_index:sif`, na.rm = TRUE),
        .groups = "drop"
      )
    
    final_rho_heat <- ccm_summary_heat$rho_mean[nrow(ccm_summary_heat)]
    trend_heat <- cor(ccm_summary_heat$lib_size, ccm_summary_heat$rho_mean)
    
    # S-map分析
    smap_heat <- SMap(
      dataFrame = station_data_ccm,
      lib = paste("1", n_data),
      pred = paste("1", n_data),
      E = E_ccm,
      theta = 2,
      columns = "heat_index",
      target = "sif",
      embedded = FALSE
    )
    
    smap_coeffs <- smap_heat$coefficients
    
    if (!is.null(smap_coeffs) && ncol(smap_coeffs) >= 2) {
      heat_coef_col <- 
        which(grepl("heat", colnames(smap_coeffs), ignore.case = TRUE))
      
      if (length(heat_coef_col) > 0) {
        smap_coef_heat <- smap_coeffs[, heat_coef_col[1]]
      } else {
        smap_coef_heat <- smap_coeffs[, 2]
      }
      
      smap_coef_heat <- smap_coef_heat[!is.na(smap_coef_heat)]
      
      mean_smap_coef <- mean(smap_coef_heat, na.rm = TRUE)
      sd_smap_coef <- sd(smap_coef_heat, na.rm = TRUE)
      
    } else {
      mean_smap_coef <- NA
      sd_smap_coef <- NA
      smap_coef_heat <- NA
    }
    
    # 判断因果
    ccm_threshold_rho <- 0.1
    ccm_threshold_trend <- 0
    
    heat_causes_sif <- (
      final_rho_heat > ccm_threshold_rho & 
        trend_heat > ccm_threshold_trend
    )
    
    effect_type <- case_when(
      !heat_causes_sif ~ "无因果",
      is.na(mean_smap_coef) ~ "S-map失败",
      mean_smap_coef > 0 ~ "促进",
      mean_smap_coef < 0 ~ "抑制"
    )
    
    effect_stability <- case_when(
      is.na(mean_smap_coef) | is.na(sd_smap_coef) ~ "未知",
      abs(mean_smap_coef) > sd_smap_coef ~ "稳定",
      TRUE ~ "波动大"
    )
    
    # 返回结果
    tibble(
      meteo_stat_id = station_id,
      tp = tp_x, 
      n_points = n_data,
      E = E_ccm,
      
      # 热事件指数 → SIF
      rho_heat_to_sif = final_rho_heat,
      trend_heat_to_sif = trend_heat,
      heat_causes_sif = heat_causes_sif,
      effect_index_heat = mean_smap_coef,
      effect_index_sd_heat = sd_smap_coef,
      effect_type_heat = effect_type,
      effect_stability_heat = effect_stability,
      
      # 保存数据
      smap_coefficients_heat = list(smap_coef_heat),
      ccm_summary_heat = list(ccm_summary_heat)
    )
    
  }, error = function(e) {
    message("Error in station ", station_id, ": ", e$message)
    return(NULL)
  })
}

# Replace the target list below with your own:
list(
  # 气象数据文件列表。
  tar_target(
    meteo_file_list, 
    list.files(
      path = "data_raw/meteo_data_1961-2023",
      pattern = "\\.txt$",
      full.names = TRUE
    ) %>%
      .[!grepl("sta_lonlat_china.txt", .)]
  ),
  # 读取SIF站点列表。
  # Bug: 为何SIF数据只有2000-2023年？LHSIF是否覆盖更多年份？
  tar_target(
    meteo_sif_data, 
    read.csv("data_raw/meteo_stat_SIF_data.csv") %>% 
      tibble() %>% 
      rename_with(~tolower(.x)) %>%
      rename(meteo_stat_id = meteo_stat) %>% 
      filter(!is.na(sif)) %>% 
      mutate(meteo_stat_id = as.character(meteo_stat_id))
  ), 
  tar_target(
    target_stations, 
    unique(meteo_sif_data$meteo_stat_id)
  ),
  # 读取每日气象数据并初步筛选。
  tar_target(
    meteo_data_daily,
    map(
      meteo_file_list[str_extract(basename(meteo_file_list), "\\d+") %in% target_stations],
      read_one_meteo_daily
    ) %>%
      list_rbind() %>%
      rename_with(~ tolower(.x)) %>%
      mutate(
        date = as.Date(date),
        year = year(date),
        month = month(date),
        day = day(date)
      ) %>%
      # Bug：仅分析每年的5-9月是否合理？
      filter(year %in% c(unique(meteo_sif_data$year)), month %in% 5:9)
  ),
  # 计算每日VPD。
  tar_target(
    meteo_data_daily_vpd,
    meteo_data_daily %>%
      mutate(
        # 饱和水汽压 (kPa)
        svp = 0.6112 * exp((17.67 * tavg) / (tavg + 243.5)),
        # 实际水汽压 (kPa)
        avp = (rh / 100) * svp,
        # VPD (kPa)
        vpd = svp - avp
      ) %>%
      mutate(
        vpd = case_when(
          vpd < 0 ~ NA_real_,
          vpd > 10 ~ NA_real_,
          is.na(tavg) | is.na(rh) ~ NA_real_,
          TRUE ~ vpd
        )
      )
  ),
  # 月度热胁迫指标计算。
  tar_target(
    monthly_heat_metrics,
    meteo_data_daily_vpd %>%
      mutate(
        # Bug：是否需要基于不同阈值做敏感性分析？2.0 kPa是常用的高VPD阈值。
        is_heat_event = vpd > 2.0,
        heat_intensity = pmax(vpd - 2.0, 0)
      ) %>%
      group_by(meteo_stat_id, year, month) %>%
      summarise(
        n_days = n(),
        vpd_mean = mean(vpd, na.rm = TRUE),
        heat_event_days = sum(is_heat_event, na.rm = TRUE),
        heat_event_freq = heat_event_days / n_days, 
        heat_over_sum = sum(heat_intensity, na.rm = TRUE),  
        .groups = "drop"
      ) %>%
      group_by(meteo_stat_id) %>%
      mutate(
        heat_freq_z = scale(heat_event_freq)[,1],
        heat_intensity_z = scale(heat_over_sum)[,1],
        heat_index_composite = 0.5 * heat_freq_z + 0.5 * heat_intensity_z
      ) %>%
      filter(n() > 30) %>% 
      ungroup()
  ),
  # 汇总月度SIF。
  tar_target(
    sif_monthly,
    meteo_sif_data %>%
      group_by(meteo_stat_id, year, month) %>%
      summarise(sif = mean(sif, na.rm = TRUE), .groups = "drop")
  ),
  # 合并数据并进行去趋势。
  tar_target(
    data_heat_sif,
    monthly_heat_metrics %>%
      inner_join(sif_monthly, by = c("meteo_stat_id", "year", "month")) %>%
      filter(!is.na(sif)) %>% 
      group_by(meteo_stat_id) %>%
      arrange(year, month) %>%
      mutate(
        time_idx = row_number(),
        sif_detrended = safe_detrend(sif, time_idx),
        heat_index_composite_detrended = safe_detrend(heat_index_composite, time_idx),
        heat_event_freq_detrended = safe_detrend(heat_event_freq, time_idx),
        heat_intensity_sum_detrended = safe_detrend(heat_over_sum, time_idx)
      ) %>%
      ungroup() %>%
      group_by(meteo_stat_id) %>%
      filter(sum(!is.na(sif_detrended)) >= 30) %>%  
      ungroup()
  ), 
  # 分析维度 ----
  tar_target(
    # 2020年GDP (万元)
    gdp_data_2020, 
    read_excel("data_raw/中国城市数据库1990-2023.xlsx") %>%
      filter(年份 == 2020) %>%
      dplyr::select(city_name = 城市, gdp = "地区生产总值(万元)") %>%
      mutate(gdp = as.numeric(gdp))
  ), 
  tar_target(
    # A. 2020年绿化投资 (万元)
    green_invest_2020, 
    read_csv("data_raw/green_invest/中国城市数据.csv", show_col_types = FALSE) %>%
      dplyr::select(city_name = 城市名称, invest = 园林绿化_2020) %>%
      mutate(invest = as.numeric(invest),
             city_name = ifelse(str_detect(city_name, "市$"), city_name, paste0(city_name, "市")))
  ), 
  # 绿地面积。
  tar_target(
    green_area_2020, 
    read_excel(
      "data_raw/城市绿地面积数据_2003-2023.xlsx", sheet = "绿地面积_明细数据"
    ) %>%
      filter(年份 == 2020) %>%
      dplyr::select(
        city_name = 城市, 
        area_green_tot = "绿地面积(公顷)", 
        area_green_built = "建成区绿地面积(公顷)", 
        area_green_park = "公园绿地面积(公顷)"
      ) %>%
      mutate(
        area_green_tot = as.numeric(area_green_tot),
        area_green_built = as.numeric(area_green_built),
        area_green_park = as.numeric(area_green_park),
        city_name = ifelse(
          str_detect(city_name, "市$"), city_name, paste0(city_name, "市")
        )
      )
  ), 
  # 空间站点匹配到城市。
  tar_target(
    china_cities_shp, 
    st_read("data_raw/china_cities/city.shp", quiet = TRUE) %>% 
      st_transform(crs = 4326)
  ), 
  tar_target(
    pts_sf, 
    ccm_results_heat %>% 
      distinct(meteo_stat_id, longitude, latitude) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  ), 
  tar_target(
    station_city_map, 
    pts_sf %>% 
      st_join(china_cities_shp, join = st_within) %>% 
      as.data.frame() %>%
      dplyr::select(meteo_stat_id, city_name = ct_name)
  ), 
  tar_target(
    invest_metrics, 
    station_city_map %>%
      left_join(green_invest_2020, by = "city_name") %>%
      left_join(gdp_data_2020, by = "city_name") %>%
      left_join(green_area_2020, by = "city_name") %>%
      mutate(
        invest_ratio = invest / gdp * 100,
        invest_pa_tot = invest / area_green_tot, 
        invest_pa_built = invest / area_green_built, 
        invest_pa_park = invest / area_green_park
      ) %>%
      dplyr::select(
        meteo_stat_id, invest_ratio, 
        invest_pa_park, invest_pa_built, invest_pa_park
      )
  ), 
  # 站点坐标。
  # Bug：可以合并它和target_stations？
  tar_target(
    station_coords, 
    read_csv("data_raw/meteo_stat_SIF_data.csv") %>%
      rename_with(~tolower(.x)) %>%
      select(meteo_stat_id = meteo_stat, longitude, latitude) %>%
      distinct(meteo_stat_id, .keep_all = TRUE) %>%
      mutate(meteo_stat_id = as.character(meteo_stat_id))
  ), 
  # CCM ----
  tar_target(
    all_stations_heat, 
    data_heat_sif %>%
      group_by(meteo_stat_id) %>%
      summarise(n = n()) %>%
      filter(n >= 30) %>%
      pull(meteo_stat_id)
  ), 
  # 不同Tp下各站点CCM结果。
  tar_target(
    ccm_results_heat, 
    map_dfr(c(0:4), function(current_tp) {
      # 内部循环遍历所有站点
      map_dfr(seq_along(all_stations_heat), function(i) {
        perform_ccm_heat_sif(
          all_stations_heat[i], data_heat_sif, tp_x = current_tp
        )
      })
    }) %>% 
      # 加入站点的经纬度数据。
      left_join(station_coords, by = "meteo_stat_id") %>%
      filter(!is.na(longitude), !is.na(latitude)) %>%
      mutate(effect_strength = abs(effect_index_heat))
  )
)
