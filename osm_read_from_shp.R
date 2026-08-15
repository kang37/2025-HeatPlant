# 加载必要的库
library(tidyverse)
library(sf)

# 读取站点元数据
stationMetadata <- read_csv('E:/GEE/airquality_stations_locations.csv')

# 读取福建省边界的shapefile
fujian_boundary <- st_read('F:/fujian/fjshape/fjshape/fujian.shp')

# 将站点数据转换为sf对象，并设置坐标参考系为WGS84（EPSG:4326）
stations_sf <- st_as_sf(stationMetadata, coords = c('Longitude', 'Latitude'), crs = 4326)

# 确保福建边界的坐标系与站点一致
fujian_boundary <- st_transform(fujian_boundary, st_crs(stations_sf))

# 提取位于福建省边界内的站点
fujian_stations <- stations_sf %>%
  st_intersection(fujian_boundary)

# 输出福建站点的结果以检查
print(fujian_stations)

# 初始化一个空的tibble数据框用于存储输出
osmOutput <- tibble()

# 从本地读取福建的建筑物shapefile
buildings_sf <- st_read('E:/WeChat Files/wl479415435/FileStorage/File/2024-11/fujian-latest-free.shp')

# 确保建筑物数据和站点数据的坐标系一致
buildings_sf <- st_transform(buildings_sf, st_crs(fujian_stations))

# 从福建站点开始逐个计算建筑物相关信息
for (i in seq(1, nrow(fujian_stations), 1)) {
  
  print(paste0('Iteration number...', i, '  -  ', fujian_stations$AirQualityStation[i]))
  
  # 选择当前迭代的站点数据
  select <- fujian_stations[i, ]
  
  # 计算站点所在的UTM区号
  zone <- (floor((st_coordinates(select)[1] + 180) / 6) %% 60) + 1
  utmCrs <- st_crs(paste0("+proj=utm +zone=", zone, " +ellps=WGS84"))
  
  # 将站点数据转换为UTM坐标系，并设置缓冲区为30米
  selectUTM <- st_transform(select, utmCrs)
  selectUTM_buff <- selectUTM %>% st_buffer(30)
  
  # 筛选建筑物数据，并计算与站点缓冲区的交集
  buildingFootprint <- buildings_sf %>%
    st_transform(utmCrs) %>%
    st_intersection(selectUTM_buff)
  
  if (nrow(buildingFootprint) > 0) {
    # 如果缓冲区内存在建筑物，计算建筑物的总面积
    buildingFootprint <- buildingFootprint %>%
      mutate(area = as.numeric(st_area(geometry))) %>%
      summarise(area = sum(area))
    buildingFootprint <- buildingFootprint$area
  } else {
    print('no buildings....')
    buildingFootprint <- 0
  }
  
  # 提取当前站点的经纬度
  coords <- st_coordinates(select)
  
  # 将结果存储在osmInner tibble中
  osmInner <- tibble(
    AirQualityStation = fujian_stations$AirQualityStation[i],
    Longitude = coords[1],
    Latitude = coords[2],
    buildingFootprint = buildingFootprint
  )
  
  # 将osmInner绑定到osmOutput
  osmOutput <- osmOutput %>% bind_rows(osmInner)
}

# 将福建站点及其建筑物信息保存为CSV文件
osmOutput %>%
  write_csv('E:/GEE/fujian_stations_buildings_with_coords30.csv')

# 将福建站点元数据（含坐标）另存为CSV文件
fujian_stations %>%
  st_drop_geometry() %>%
  write_csv('E:/GEE/fujian_stations_with_coords.csv')
