#!/usr/bin/env Rscript
# =============================================================================
# province_ccm_fig2.R
#   仿 Nat Commun 2024 (s41467-024-48199-z) Fig.2a: 州级 CCM 替代检验小多图。
#   本研究: 每个省份一个子图, 子图内每个站点=一个圆(观测CCM技巧rho),
#   竖线=该站季节替代零模型的95%包络(0 -> rho_null95)。
#   观测rho 超过包络且收敛 → 显著(实心, 红=促进/蓝=抑制); 否则空心灰。
#   数据: 归一化(Ushio式)+年块替代 CCM。每站取峰值rho的tp。全英文。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(sf); library(rnaturalearth)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
SIF_COL<-"SIF_buf1000"

# ---- CCM: 每站峰值tp ----
d<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
d<-d[y_var==paste0(SIF_COL,"_dt") & x_var=="vpd_mean_dt"]
d[, drho:=rho-rho_min]
full9<-d[,.N,by=meteo_stat][N==9,meteo_stat]
d<-d[meteo_stat %in% full9]
# 每站取"最强证据"的tp: 在收敛(drho>0)的滞后中 p_surr 最小者; 若无收敛则取峰值rho
d[, cand:=drho>0]
pick_row<-function(sd){
  s<-sd[cand==TRUE]
  if(nrow(s)) s[order(p_surr,-rho)][1] else sd[order(-rho)][1]
}
pk<-d[, pick_row(.SD), by=meteo_stat]
pk[, sig:=(p_surr<0.1) & (drho>0)]             # 与全文一致的宽松判据
pk[, dir:=fifelse(mean_coef>0,"Promote","Inhibit")]
pk[, status:=fifelse(!sig,"Not significant", paste("Sig:",dir))]

# ---- 站点 -> 省份 (空间连接) ----
prov<-ne_states(country=c("China","Taiwan"),returnclass="sf")
pts<-st_as_sf(pk[!is.na(longitude)&!is.na(latitude)],coords=c("longitude","latitude"),crs=4326,remove=FALSE)
ji<-st_join(pts, prov[,c("name_en","region")], join=st_intersects)
pk2<-as.data.table(st_drop_geometry(ji))
# 落海/边界外的最近邻补齐
miss<-which(is.na(pk2$name_en))
if(length(miss)){
  nr<-st_nearest_feature(pts[miss,],prov)
  pk2$name_en[miss]<-prov$name_en[nr]; pk2$region[miss]<-prov$region[nr]
}
pk2<-pk2[!is.na(name_en)]
pk2[name_en %in% c("Paracel Islands","Spratly Islands"), name_en:=NA]; pk2<-pk2[!is.na(name_en)]

# 省内按rho排序给x位置; 省按显著站数排序(多在前)
prov_ord<-pk2[, .(nsig=sum(sig), n=.N), by=name_en][order(-nsig,-n)]
pk2[, prov_f:=factor(name_en, levels=prov_ord$name_en)]
pk2[, prov_lab:=sprintf("%s  (%d/%d sig)", name_en,
        prov_ord$nsig[match(name_en,prov_ord$name_en)],
        prov_ord$n[match(name_en,prov_ord$name_en)])]
lev_lab<-sprintf("%s  (%d/%d sig)",prov_ord$name_en,prov_ord$nsig,prov_ord$n)
pk2[, prov_lab:=factor(prov_lab, levels=lev_lab)]
setorder(pk2, prov_f, rho)
pk2[, xpos:=seq_len(.N), by=prov_f]

cat("省份数:",uniqueN(pk2$name_en)," 站点数:",nrow(pk2)," 显著:",sum(pk2$sig),"\n")
cat("\n各省 显著/总:\n"); print(prov_ord)
fwrite(pk2[,.(meteo_stat,name_en,region,tp,rho,rho_min,rho_null95,drho,p_surr,mean_coef,dir,sig)],
       file.path(OUT,"province_ccm_stations.csv"))

# ---- 图: 省级小多图 ----
cols<-c("Sig: Promote"="#C1121F","Sig: Inhibit"="#1A6FBF","Not significant"="grey70")
p<-ggplot(pk2,aes(xpos,rho))+
  geom_linerange(aes(ymin=0,ymax=rho_null95),color="grey75",linewidth=.5)+  # 季节替代95%包络
  geom_point(aes(color=status,fill=status,shape=sig),size=1.7,stroke=.5)+
  scale_color_manual(values=cols,name=NULL)+
  scale_fill_manual(values=cols,name=NULL)+
  scale_shape_manual(values=c(`TRUE`=21,`FALSE`=1),guide="none")+
  facet_wrap(~prov_lab,scales="free_x",ncol=5)+
  labs(title="Province-level VPD->SIF CCM surrogate test (per station)",
       subtitle=paste0("Each circle = one station at its strongest-evidence lag; grey bar = 0 to seasonal-surrogate 95% null (rho_null95).\n",
                       "Filled = passes seasonal-surrogate test (p<0.1) & converges; red=promote, blue=inhibit. Ushio-normalized, 1km. Panel title: sig/total stations."),
       x="Stations within province (ordered by cross-map skill)",y="Cross-map skill rho")+
  theme_bw(base_size=10)+
  theme(plot.title=element_text(face="bold"),legend.position="top",
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(size=7.5,face="bold"),panel.grid.minor=element_blank())
ng<-uniqueN(pk2$prov_lab); nr<-ceiling(ng/5)
ggsave(file.path(OUT,"province_ccm_fig2.png"),p,width=15,height=max(8,1.6*nr),dpi=300,limitsize=FALSE)
cat("\n-> province_ccm_fig2.png  (", ng, "provinces,", nr, "rows )\n")
