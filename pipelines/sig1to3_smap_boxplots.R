#!/usr/bin/env Rscript
# =============================================================================
# sig1to3_smap_boxplots.R
#   对 1/2/3 tp 显著的站, 每个 站×tp 一个箱形图(逐时间点 S-map 偏导数)+jitter。
#   方向按"主导符号(>=75%)"定站点类别: 总是促进/总是抑制/混合, 用于箱体着色。
#   站点按其系数中位数排序; 2/3 tp 时同站相邻。三张图分开(1/2/3 tp)。全英文。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(rEDM); library(ggplot2)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
BUF_R<-1000L; SIF_COL<-sprintf("SIF_buf%d",BUF_R); SIF_DT<-paste0(SIF_COL,"_dt")
STAT_DIR<-file.path(PROJ,"data_raw/hcsif/station_v3"); CACHE<-file.path(PROJ,"data_raw/hcsif/vpd_8day_cache.rds")
MIN_PTS<-40L; MIN_DAYS_WIN<-5L; SEED_BASE<-20260824L
CACHE_FULL<-file.path(OUT,"sig1to3_coef_full.rds")

ccm<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
ccm<-ccm[y_var==SIF_DT & x_var=="vpd_mean_dt"]; ccm[, sig:=p_surr<0.1 & (rho-rho_min)>0]
full9<-ccm[,.N,by=meteo_stat][N==9,meteo_stat]
st<-ccm[meteo_stat%in%full9,.(nsig=sum(sig)),by=meteo_stat]
sel_st<-st[nsig>=1&nsig<=3,meteo_stat]
recs<-merge(ccm[meteo_stat%in%sel_st & sig==TRUE,.(meteo_stat,tp,E)],st,by="meteo_stat")

if(file.exists(CACHE_FULL)){
  cf<-readRDS(CACHE_FULL); cat("读取全序列缓存\n")
}else{
  ids<-unique(recs$meteo_stat)
  fs<-list.files(STAT_DIR,pattern="^hcsif_station_[0-9]{4}\\.csv$",full.names=TRUE)
  sif<-rbindlist(lapply(fs,fread),fill=TRUE)[meteo_stat%in%ids,c("meteo_stat","year","doy","date",SIF_COL),with=FALSE]
  vpd<-readRDS(CACHE)[n_days>=MIN_DAYS_WIN & meteo_stat%in%ids]
  d<-merge(sif,vpd,by=c("meteo_stat","year","doy")); setorder(d,meteo_stat,year,doy)
  d[, idx8:=as.integer(round(as.numeric(as.Date(as.character(date),format="%Y%m%d")-as.Date("2000-01-01"))/8))]
  znorm<-function(x){s<-sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)return(rep(NA_real_,length(x)));(x-mean(x,na.rm=TRUE))/s}
  for(v in c(SIF_COL,"vpd_mean")) d[, paste0(v,"_dt"):=znorm(get(v)), by=meteo_stat]
  get_coefs<-function(dd,tp_x,best_E){
    dd<-copy(dd); dd[, x_lag:=shift(vpd_mean_dt,n=tp_x,type="lag"), by=year]; dd<-dd[!is.na(idx8)]
    tmp<-dd[,.(idx8,y=get(SIF_DT),x=x_lag)]
    full<-tmp[data.table(idx8=seq.int(min(dd$idx8),max(dd$idx8))),on="idx8"]; setorder(full,idx8)
    if(full[!is.na(y)&!is.na(x),.N]<MIN_PTS) return(NULL)
    df<-data.frame(time=seq_len(nrow(full)),sif=full$y,heat=full$x)
    set.seed(SEED_BASE+as.integer(dd$meteo_stat[1])+tp_x*1000L)
    sm<-tryCatch(rEDM::SMap(dataFrame=df,E=best_E,theta=2,lib=paste("1",nrow(df)),
                  pred=paste("1",nrow(df)),columns="heat",target="sif",embedded=FALSE),error=function(e)NULL)
    if(is.null(sm))return(NULL)
    co<-as.data.table(sm$coefficients); cc<-which(grepl("heat",colnames(co),ignore.case=TRUE))[1]; if(is.na(cc))cc<-3L
    tt<-co[["time"]]; cv<-co[[cc]]; ok<-!is.na(cv)&!is.nan(cv)&tt>=1&tt<=nrow(full); cv[ok]
  }
  cf<-rbindlist(lapply(split(recs,seq_len(nrow(recs))),function(r){
    v<-get_coefs(d[meteo_stat==r$meteo_stat],r$tp,r$E); if(is.null(v)||!length(v))return(NULL)
    data.table(meteo_stat=r$meteo_stat,tp=r$tp,nsig=r$nsig,coef=v)}),fill=TRUE)
  saveRDS(cf,CACHE_FULL)
}
cat("系数点总数:",nrow(cf)," 记录数:",cf[,uniqueN(paste(meteo_stat,tp))],"\n")

# ---- 记录方向(主导>=75%) + 站点类别 ----
recdir<-cf[,.(med=median(coef),frac_pos=mean(coef>0),.N),by=.(meteo_stat,tp,nsig)]
recdir[, dir:=fifelse(frac_pos>=0.75,"Promote",fifelse(frac_pos<=0.25,"Inhibit","Ambiguous"))]
stcat<-recdir[,.(cat={u<-unique(dir)
  if(all(dir=="Promote"))"Always promote" else if(all(dir=="Inhibit"))"Always inhibit" else "Mixed"},
  st_med=median(unlist(cf[meteo_stat==.BY$meteo_stat,coef]))),by=meteo_stat]
cf<-merge(cf,recdir[,.(meteo_stat,tp,rec_med=med,dir)],by=c("meteo_stat","tp"))
cf<-merge(cf,stcat,by="meteo_stat")

catcol<-c(`Always promote`="#C1121F",`Always inhibit`="#1A6FBF",`Mixed`="#B0A160")

plot_n<-function(NT,wd){
  x<-cf[nsig==NT]
  # 排序: 站点按 st_med, 站内按 tp 相邻
  ordkey<-unique(x[,.(meteo_stat,tp,st_med)])[order(st_med,tp)]
  ordkey[, pos:=.I]; ordkey[, xlab:=as.character(meteo_stat)]
  x<-merge(x,ordkey[,.(meteo_stat,tp,pos,xlab)],by=c("meteo_stat","tp"))
  x[, pos:=factor(pos,levels=ordkey$pos,labels=ordkey$xlab)]
  ncat<-stcat[meteo_stat%in%x$meteo_stat,.N,by=cat]
  sub<-sprintf("%d stations x %d tp = %d boxes; box color = station type (dominant >=75%%): %s",
               uniqueN(x$meteo_stat),NT,nrow(ordkey),
               paste(sprintf("%s=%d",ncat$cat,ncat$N),collapse=", "))
  p<-ggplot(x,aes(pos,coef,fill=cat))+
    geom_hline(yintercept=0,color="grey40",linewidth=.4)+
    geom_boxplot(outlier.shape=NA,color="grey25",linewidth=.25)+
    geom_jitter(width=.18,size=.35,alpha=.25,color="grey15")+
    scale_fill_manual(values=catcol,name="Station type")+
    labs(title=sprintf("Time-varying S-map partial derivative dSIF/dVPD per station x tp — %d significant tp",NT),
         subtitle=sub, x="Station ID (ordered by median coefficient; same station's tp adjacent)",
         y="S-map coefficient (local dSIF/dVPD)")+
    theme_bw(base_size=11)+
    theme(plot.title=element_text(face="bold"),legend.position="top",
          axis.text.x=element_text(angle=90,vjust=.5,hjust=1,
                     size=if(NT==1)3.2 else if(NT==2)4.5 else 7))
  fn<-sprintf("sig%dtp_smap_boxplots.png",NT)
  ggsave(file.path(OUT,fn),p,width=wd,height=6.5,dpi=300,limitsize=FALSE)
  cat("->",fn,"|",sub,"\n")
}
plot_n(1, 26)
plot_n(2, 22)
plot_n(3, 12)
