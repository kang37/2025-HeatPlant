#!/usr/bin/env Rscript
# 1km: Ushio式归一化 vs 线性去趋势 —— 严格5类对比 + 转移 + 气候区地图。全英文。
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggalluvial); library(sf); library(rnaturalearth)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
OUT<-file.path(PROJ,"data_proc/output_hcsif_buf1000_norm"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
latest<-function(dir,pat){f<-list.files(dir,pattern=pat,full.names=TRUE);f<-f[!grepl("_raw",f)];tail(sort(f),1)}
classify5<-function(c){s<-sign(c);if(any(is.na(s))||any(s==0))return("other");nch<-sum(diff(s)!=0)
 if(nch==0)return(if(s[1]>0)"always_promote" else "always_inhibit")
 if(nch==1)return(if(s[1]<0)"inhibit_promote" else "promote_inhibit");"other"}
lev<-c("always_promote","promote_inhibit","inhibit_promote","always_inhibit","other")
lab<-c("Always promote","Promote->Inhibit","Inhibit->Promote","Always inhibit","Other (mixed)")
col<-setNames(c("#C1121F","#F4A261","#74C6E8","#1A6FBF","#9E9E9E"),lab)

styp<-function(rds,keepll=FALSE){SIF<-"SIF_buf1000_dt"
 d<-as.data.frame(readRDS(rds));d<-d[d$y_var==SIF&d$x_var=="vpd_mean_dt",]
 g<-d|>arrange(meteo_stat,tp)|>group_by(meteo_stat)|>filter(n()==9)
 if(keepll) g<-g|>summarise(stype=classify5(mean_coef),longitude=first(longitude),latitude=first(latitude),.groups="drop")
 else g<-g|>summarise(stype=classify5(mean_coef),.groups="drop")
 g|>mutate(meteo_stat=as.integer(meteo_stat))}

nrm<-styp(latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_norm"),"buf1000_vpd_.*rds$"),keepll=TRUE)
det<-styp(latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_v3"),"buf1000_vpd_.*rds$"))
cat("=== 5类分布对比 (1km) ===\n")
cmp<-full_join(count(nrm,stype)|>rename(normalize=n), count(det,stype)|>rename(detrend=n), by="stype")
cmp$stype<-factor(cmp$stype,levels=lev); cmp<-cmp|>arrange(stype); print(as.data.frame(cmp))

# 转移(去趋势 -> 归一化)
m<-inner_join(det|>rename(detrend=stype), nrm|>select(meteo_stat,normalize=stype), by="meteo_stat")
same<-mean(m$detrend==m$normalize)
cat(sprintf("\n共同站点 %d；去趋势与归一化 类别一致: %d (%.1f%%)\n", nrow(m), sum(m$detrend==m$normalize), 100*same))
tm<-m|>count(detrend,normalize)|>mutate(detrend=factor(lab[match(detrend,lev)],levels=lab),
                                        normalize=factor(lab[match(normalize,lev)],levels=lab))
fwrite(as.data.table(tidyr::pivot_wider(tm,names_from=normalize,values_from=n,values_fill=0)),
       file.path(OUT,"transition_detrend_to_norm.csv"))

# 冲积图
al<-m|>mutate(id=row_number())|>pivot_longer(c(detrend,normalize),names_to="method",values_to="stype")|>
  mutate(method=factor(method,levels=c("detrend","normalize"),labels=c("Detrend","Normalize")),
         stype=factor(lab[match(stype,lev)],levels=lab))
pa<-ggplot(al,aes(x=method,stratum=stype,alluvium=id,fill=stype))+
  geom_flow(stat="alluvium",alpha=.55,color=NA,decreasing=FALSE)+
  geom_stratum(width=.4,color="white",decreasing=FALSE)+
  geom_text(stat="stratum",aes(label=after_stat(count)),size=3,decreasing=FALSE)+
  scale_fill_manual(values=col,name="Pattern")+
  labs(title="1 km strict 5-class: Detrend vs Ushio-style Normalize",
       subtitle=sprintf("n=%d stations classifiable both ways; %.0f%% keep the same class",nrow(m),100*same),
       x="Preprocessing",y="Number of stations")+
  theme_minimal(base_size=12)+theme(plot.title=element_text(face="bold"),panel.grid.major.x=element_blank())
ggsave(file.path(OUT,"sankey_detrend_vs_norm.png"),pa,width=10,height=8,dpi=300)
cat("-> sankey_detrend_vs_norm.png\n")

# 归一化版 气候区地图
cov<-readRDS(file.path(PROJ,"data_proc/output_10y_built_up_05_01/station_covariates.rds"))|>
  transmute(meteo_stat=as.integer(meteo_stat_id),koppen_group)
nrm2<-nrm|>left_join(cov,by="meteo_stat")|>mutate(sl=factor(lab[match(stype,lev)],levels=lab))
fwrite(as.data.table(nrm2),file.path(OUT,"norm_5class.csv"))
bb<-c(72,136,17,54);k2g<-function(x)case_when(x>=1&x<=4~"A",x>=5&x<=9~"B",x>=10&x<=17~"C",x>=18&x<=28~"D",TRUE~NA_character_)
world<-tryCatch(ne_countries(scale="medium",returnclass="sf"),error=function(e)NULL)
chn<-tryCatch(ne_countries(country=c("China","Taiwan","Hong Kong S.A.R.","Macao S.A.R."),scale="medium",returnclass="sf")|>st_union(),error=function(e)NULL)
prov<-tryCatch(ne_states(country="China",returnclass="sf"),error=function(e)NULL)
kdf<-NULL;kr<-tryCatch(terra::crop(terra::rast(file.path(PROJ,"data_raw/koppen_geiger_tif/1991_2020/koppen_geiger_0p5.tif")),terra::ext(bb[1],bb[2],bb[3],bb[4])),error=function(e)NULL)
if(!is.null(kr)&&!is.null(chn)){cv<-tryCatch(terra::vect(chn),error=function(e)NULL);krm<-if(!is.null(cv))tryCatch(terra::mask(kr,cv),error=function(e)kr) else kr
 kdf<-terra::as.data.frame(krm,xy=TRUE)|>rename(code=3)|>mutate(g=k2g(code))|>filter(!is.na(g))}
kcol<-c(A="#5A8F76",B="#EAD5A0",C="#A3B86C",D="#C2DFCD");klab<-c(A="A Tropical",B="B Arid",C="C Temperate",D="D Continental")
pm<-ggplot()
if(!is.null(world))pm<-pm+geom_sf(data=world,fill="grey88",color="grey75",linewidth=.15,inherit.aes=FALSE)
if(!is.null(kdf))pm<-pm+geom_raster(data=kdf,aes(x,y,fill=g),alpha=.6)+scale_fill_manual(values=kcol,labels=klab,name="Koppen zone",na.value="grey88")
if(!is.null(prov))pm<-pm+geom_sf(data=prov,fill=NA,color="grey55",linewidth=.2,inherit.aes=FALSE)
pm<-pm+ggnewscale::new_scale_color()+
  geom_point(data=nrm2,aes(longitude,latitude,color=sl),size=1.7,alpha=.85)+
  scale_color_manual(values=col,name="VPD->SIF pattern (5-class)")+
  coord_sf(xlim=bb[1:2],ylim=bb[3:4],expand=FALSE)+
  labs(title="VPD->SIF strict 5-class pattern (1 km, Ushio-style normalization)",
       subtitle=sprintf("n=%d; normalize to unit variance (no detrend/deseasonalize); S-map sign sequence over tp=0..8",nrow(nrm2)),
       x="Longitude",y="Latitude")+
  theme_minimal(base_size=12)+theme(panel.background=element_rect(fill="#D6EAF8",color=NA),
    plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,color="grey40"))
ggsave(file.path(OUT,"norm_5class_map.png"),pm,width=14,height=9,dpi=300)
cat("-> norm_5class_map.png\n")
cat("\n=== 转移矩阵(行=去趋势, 列=归一化) ===\n")
print(as.data.frame(tidyr::pivot_wider(tm,names_from=normalize,values_from=n,values_fill=0)))
