#!/usr/bin/env Rscript
# =============================================================================
# sig1to3_smap_boxplots_grouped.R
#   版本二: 无 jitter; 横轴按站点类别分组(促进/抑制/混合各一子图), 子图共用Y轴范围。
#   同站的多个 tp 相邻; 各子图内站点按系数中位数排序。用缓存的逐点系数。全英文。
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
cf<-readRDS(file.path(OUT,"sig1to3_coef_full.rds"))   # meteo_stat,tp,nsig,coef (逐点)

# 记录方向(主导>=75%) + 站点类别
recdir<-cf[,.(frac_pos=mean(coef>0)),by=.(meteo_stat,tp,nsig)]
recdir[, dir:=fifelse(frac_pos>=0.75,"Promote",fifelse(frac_pos<=0.25,"Inhibit","Ambiguous"))]
stcat<-recdir[,.(cat={if(all(dir=="Promote"))"Always promote" else if(all(dir=="Inhibit"))"Always inhibit" else "Mixed"}),by=meteo_stat]
stmed<-cf[,.(st_med=median(coef)),by=meteo_stat]
cf<-merge(merge(cf,stcat,by="meteo_stat"),stmed,by="meteo_stat")
cf[, cat:=factor(cat,levels=c("Always promote","Always inhibit","Mixed"))]
catcol<-c(`Always promote`="#C1121F",`Always inhibit`="#1A6FBF",`Mixed`="#B0A160")

plot_n<-function(NT,wd){
  x<-cf[nsig==NT]
  # 站点排序: 类别分组内按站点中位数; 站点为一个x刻度
  sord<-unique(x[,.(meteo_stat,cat,st_med)])[order(cat,st_med)]
  x[, stf:=factor(meteo_stat,levels=sord$meteo_stat)]
  nc<-stcat[meteo_stat%in%x$meteo_stat][,.N,by=cat]
  base<-ggplot(x,aes(stf,coef))+
    geom_hline(yintercept=0,color="grey40",linewidth=.4)+
    facet_grid(~cat,scales="free_x",space="free_x")+
    labs(x="Station ID (grouped by type; ordered by median within group)",
         y="S-map coefficient (local dSIF/dVPD)")+
    theme_bw(base_size=11)+
    theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"),
          axis.text.x=element_text(angle=90,vjust=.5,hjust=1,
                     size=if(NT==1)3.0 else if(NT==2)4.2 else 6.5))
  if(NT==1){
    # 1 tp: 每站一个盒, 按类别着色
    p<-base+geom_boxplot(aes(fill=cat),outlier.shape=NA,color="grey25",linewidth=.25)+
      scale_fill_manual(values=catcol,guide="none")+
      labs(title="S-map partial derivative dSIF/dVPD per station — 1 significant tp (grouped by station type, no jitter)",
           subtitle=sprintf("Dominant >=75%% station type; %d stations; same Y range; ordered by median coefficient",uniqueN(x$meteo_stat)))
  }else{
    # 2/3 tp: 同站的多个tp在同一x刻度上并排, 按tp着色
    x[, tpf:=factor(tp)]
    p<-base%+%x+aes(fill=tpf)+
      geom_boxplot(outlier.shape=NA,color="grey25",linewidth=.25,
                   position=position_dodge2(preserve="single"))+
      scale_fill_brewer(palette="Set2",name="Lag tp")+
      labs(title=sprintf("S-map partial derivative dSIF/dVPD per station — %d significant tp (grouped by station type; colored by tp; no jitter)",NT),
           subtitle=sprintf("Dominant >=75%% station type; %d stations x %d tp side-by-side per station; same Y range; ordered by median",uniqueN(x$meteo_stat),NT))+
      theme(legend.position="top")
  }
  fn<-sprintf("sig%dtp_smap_boxplots_grouped.png",NT)
  ggsave(file.path(OUT,fn),p,width=wd,height=6.5,dpi=300,limitsize=FALSE)
  cat("->",fn,"| ",paste(sprintf("%s=%d",nc$cat,nc$N),collapse=", "),"\n")
}
plot_n(1,26); plot_n(2,22); plot_n(3,12)
