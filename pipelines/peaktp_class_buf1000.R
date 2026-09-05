#!/usr/bin/env Rscript
# buf1000: 按 rho(交叉映射技巧)最大的 tp 处 VPD->SIF 的方向(促进/抑制) 二分类。
# 出图：中国气候区分布地图 + 气候区构成条形 + 峰值tp直方。全英文。
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(sf); library(rnaturalearth)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
OUT  <- file.path(PROJ,"data_proc/output_hcsif_buf1000"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
latest<-function(dir,pat){f<-list.files(dir,pattern=pat,full.names=TRUE);f<-f[!grepl("_raw",f)];tail(sort(f),1)}

ccm <- as.data.frame(readRDS(latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_v3"),"buf1000_vpd_.*rds$")))
ccm <- ccm[ccm$y_var=="SIF_buf1000_dt" & ccm$x_var=="vpd_mean_dt",]

# 每站：rho 最大的 tp，及该 tp 处方向
cls <- ccm %>% arrange(meteo_stat,tp) %>% group_by(meteo_stat) %>%
  filter(n()==9) %>%
  summarise(peak_tp   = tp[which.max(rho)],
            rho_peak  = max(rho),
            mc_peak   = mean_coef[which.max(rho)],
            longitude = first(longitude), latitude = first(latitude),
            .groups="drop") %>%
  mutate(class = ifelse(mc_peak>0,"Promote","Inhibit"))
cat("总站点(完整9tp):",nrow(cls),"\n")
print(cls %>% count(class) %>% mutate(pct=round(100*n/sum(n),1)))
cat("峰值tp分布:\n"); print(table(cls$peak_tp))

# 气候区
cov <- readRDS(file.path(PROJ,"data_proc/output_10y_built_up_05_01/station_covariates.rds")) %>%
  transmute(meteo_stat=as.integer(meteo_stat_id), koppen_group)
cls <- cls %>% mutate(meteo_stat=as.integer(meteo_stat)) %>% left_join(cov,by="meteo_stat")
fwrite(cls, file.path(OUT,"peaktp_class.csv"))

col2 <- c("Promote"="#C1121F","Inhibit"="#1A6FBF")

# ---- 图1：中国气候区 + 两类站点分布地图 ----
koppen_tif<-file.path(PROJ,"data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif")
bb<-c(72,136,17,54)
k2g<-function(x)case_when(x>=1&x<=4~"A",x>=5&x<=9~"B",x>=10&x<=17~"C",x>=18&x<=28~"D",TRUE~NA_character_)
world<-tryCatch(ne_countries(scale="medium",returnclass="sf"),error=function(e)NULL)
chn<-tryCatch(ne_countries(country=c("China","Taiwan","Hong Kong S.A.R.","Macao S.A.R."),scale="medium",returnclass="sf")%>%st_union(),error=function(e)NULL)
prov<-tryCatch(ne_states(country="China",returnclass="sf"),error=function(e)NULL)
kdf<-NULL
kr<-tryCatch(terra::crop(terra::rast(koppen_tif),terra::ext(bb[1],bb[2],bb[3],bb[4])),error=function(e)NULL)
if(!is.null(kr)&&!is.null(chn)){cv<-tryCatch(terra::vect(chn),error=function(e)NULL)
  krm<-if(!is.null(cv))tryCatch(terra::mask(kr,cv),error=function(e)kr) else kr
  kdf<-terra::as.data.frame(krm,xy=TRUE)%>%rename(code=3)%>%mutate(g=k2g(code))%>%filter(!is.na(g))}
kcol<-c(A="#5A8F76",B="#EAD5A0",C="#A3B86C",D="#C2DFCD")
klab<-c(A="A Tropical",B="B Arid",C="C Temperate",D="D Continental")

pm<-ggplot()
if(!is.null(world))pm<-pm+geom_sf(data=world,fill="grey88",color="grey75",linewidth=.15,inherit.aes=FALSE)
if(!is.null(kdf))pm<-pm+geom_raster(data=kdf,aes(x,y,fill=g),alpha=.6)+
  scale_fill_manual(values=kcol,labels=klab,name="Koppen zone",na.value="grey88")
if(!is.null(prov))pm<-pm+geom_sf(data=prov,fill=NA,color="grey55",linewidth=.2,inherit.aes=FALSE)
pm<-pm+ggnewscale::new_scale_color()+
  geom_point(data=cls,aes(longitude,latitude,color=class),size=1.8,alpha=.85)+
  scale_color_manual(values=col2,name="VPD->SIF at peak-rho tp")+
  coord_sf(xlim=bb[1:2],ylim=bb[3:4],expand=FALSE)+
  labs(title="VPD->SIF effect at the lag of strongest cross-map skill (buffer 1 km)",
       subtitle=sprintf("Per station: direction of mean_coef at tp with max rho; n=%d (Promote=%d, Inhibit=%d)",
                        nrow(cls),sum(cls$class=="Promote"),sum(cls$class=="Inhibit")),
       x="Longitude",y="Latitude")+
  theme_minimal(base_size=12)+theme(panel.background=element_rect(fill="#D6EAF8",color=NA),
    plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,color="grey40"))
ggsave(file.path(OUT,"peaktp_class_map.png"),pm,width=14,height=9,dpi=300)
cat("-> peaktp_class_map.png\n")

# ---- 图2：各气候区两类构成 ----
zb <- cls %>% filter(!is.na(koppen_group), koppen_group %in% c("A","B","C","D")) %>%
  count(koppen_group,class)
p2<-ggplot(zb,aes(koppen_group,n,fill=class))+geom_col(position="fill",width=.65)+
  geom_text(aes(label=n),position=position_fill(vjust=.5),color="white",size=3.2)+
  scale_fill_manual(values=col2,name="Direction")+scale_y_continuous(labels=scales::percent)+
  labs(title="Promote vs Inhibit (at peak-rho tp) by Koppen zone — buffer 1 km",
       subtitle="Bar label = station count",x="Koppen zone",y="Share")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"peaktp_class_by_zone.png"),p2,width=8,height=5,dpi=300)
cat("-> peaktp_class_by_zone.png\n")

# ---- 图3：峰值tp分布(按类) ----
p3<-ggplot(cls,aes(factor(peak_tp),fill=class))+geom_bar(position="stack")+
  scale_fill_manual(values=col2,name="Direction")+
  labs(title="Distribution of the peak-rho lag (buffer 1 km)",
       subtitle="tp at which VPD->SIF cross-map skill (rho) is maximal; each step = 8 days",
       x="Peak-rho tp (each step = 8 days)",y="Number of stations")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"peaktp_distribution.png"),p3,width=8,height=5,dpi=300)
cat("-> peaktp_distribution.png\n")
cat("\n各气候区 x 类别:\n"); print(as.data.frame(tidyr::pivot_wider(zb,names_from=class,values_from=n,values_fill=0)))
