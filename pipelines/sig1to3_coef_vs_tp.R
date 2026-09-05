#!/usr/bin/env Rscript
# =============================================================================
# sig1to3_coef_vs_tp.R
#   对 1/2/3 tp 显著的站(分3子图), 展示 S-map 偏导数 ∂SIF/∂VPD 随滞后 tp(0..8)
#   的变化: 每个 tp 一个 boxplot + jitter。用全 tp 扫描(非仅显著 tp)。归一化。全英文。
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(rEDM); library(ggplot2) })
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
BUF_R<-1000L; SIF_COL<-sprintf("SIF_buf%d",BUF_R); SIF_DT<-paste0(SIF_COL,"_dt")
STAT_DIR<-file.path(PROJ,"data_raw/hcsif/station_v3"); CACHE<-file.path(PROJ,"data_raw/hcsif/vpd_8day_cache.rds")
MIN_PTS<-40L; MIN_DAYS_WIN<-5L; SEED_BASE<-20260824L; TP_SEQ<-0:8
CACHE_ALL<-file.path(OUT,"sig1to3_coef_alltp.rds")

ccm<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
ccm<-ccm[y_var==SIF_DT & x_var=="vpd_mean_dt"]; ccm[, sig:=p_surr<0.1 & (rho-rho_min)>0]
full9<-ccm[,.N,by=meteo_stat][N==9,meteo_stat]
st<-ccm[meteo_stat%in%full9,.(nsig=sum(sig)),by=meteo_stat]
sel_st<-st[nsig>=1&nsig<=3,meteo_stat]
Etab<-ccm[meteo_stat%in%sel_st,.(meteo_stat,tp,E)]     # 每站每tp的E

if(file.exists(CACHE_ALL)){
  cf<-readRDS(CACHE_ALL); cat("读取全tp缓存\n")
}else{
  fs<-list.files(STAT_DIR,pattern="^hcsif_station_[0-9]{4}\\.csv$",full.names=TRUE)
  sif<-rbindlist(lapply(fs,fread),fill=TRUE)[meteo_stat%in%sel_st,c("meteo_stat","year","doy","date",SIF_COL),with=FALSE]
  vpd<-readRDS(CACHE)[n_days>=MIN_DAYS_WIN & meteo_stat%in%sel_st]
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
  grid<-CJ(meteo_stat=sel_st,tp=TP_SEQ)
  grid<-merge(grid,Etab,by=c("meteo_stat","tp"))
  grid<-merge(grid,st,by="meteo_stat")
  cf<-rbindlist(lapply(split(grid,seq_len(nrow(grid))),function(r){
    v<-get_coefs(d[meteo_stat==r$meteo_stat],r$tp,r$E); if(is.null(v)||!length(v))return(NULL)
    data.table(meteo_stat=r$meteo_stat,tp=r$tp,nsig=r$nsig,coef=v)}),fill=TRUE)
  saveRDS(cf,CACHE_ALL)
}
cat("系数点:",nrow(cf)," 记录:",cf[,uniqueN(paste(meteo_stat,tp))]," (站",cf[,uniqueN(meteo_stat)],"x tp)\n")

# 每 tp 汇总(便于看趋势线)
cf[, grp:=factor(nsig,levels=1:3,labels=c("1 significant tp (n=260)","2 significant tp (n=90)","3 significant tp (n=25)"))]
med<-cf[,.(med=median(coef)),by=.(grp,tp)]

# jitter 抽样(点太多), 每(grp,tp)最多2000点
set.seed(7)
jit<-cf[, .SD[sample(.N, min(.N,2000))], by=.(grp,tp)]

p<-ggplot(cf,aes(factor(tp),coef))+
  geom_hline(yintercept=0,color="grey40",linewidth=.4)+
  geom_jitter(data=jit,width=.28,size=.2,alpha=.08,color="grey20")+
  geom_boxplot(aes(fill=factor(tp)),outlier.shape=NA,color="grey20",linewidth=.3,alpha=.85)+
  geom_line(data=med,aes(factor(tp),med,group=1),color="black",linewidth=.6)+
  geom_point(data=med,aes(factor(tp),med),color="black",size=1.4)+
  facet_wrap(~grp,nrow=1)+
  scale_fill_brewer(palette="Spectral",guide="none")+
  labs(title="How the S-map partial derivative dSIF/dVPD changes with lag tp",
       subtitle="Full tp sweep (0..8) for stations with 1/2/3 significant tp; box+jitter per tp; black line = median trend",
       x="Lag tp (each step = 8 days)",y="S-map coefficient (local dSIF/dVPD)")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"))
ggsave(file.path(OUT,"sig1to3_coef_vs_tp.png"),p,width=15,height=6,dpi=300)
cat("-> sig1to3_coef_vs_tp.png\n")
cat("\n各组各tp的系数中位数:\n"); print(dcast(med,tp~grp,value.var="med"))
