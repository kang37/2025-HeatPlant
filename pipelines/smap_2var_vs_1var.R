#!/usr/bin/env Rscript
# =============================================================================
# smap_2var_vs_1var.R
#   二变量 S-map: 状态空间 [SIF(t), VPD_lag(t)] (embedded=TRUE), 读 ∂SIF/∂VPD,
#   控制了 SIF 自身状态。对 1-3tp 显著记录重算方向, 与现有单驱动版本对比。归一化。
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(rEDM); library(ggplot2); library(patchwork) })
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
BUF_R<-1000L; SIF_COL<-sprintf("SIF_buf%d",BUF_R); SIF_DT<-paste0(SIF_COL,"_dt")
STAT_DIR<-file.path(PROJ,"data_raw/hcsif/station_v3"); CACHE<-file.path(PROJ,"data_raw/hcsif/vpd_8day_cache.rds")
MIN_PTS<-40L; MIN_DAYS_WIN<-5L; SEED_BASE<-20260824L
CACHE_2V<-file.path(OUT,"sig1to3_coef_2var.rds")

ccm<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
ccm<-ccm[y_var==SIF_DT & x_var=="vpd_mean_dt"]; ccm[, sig:=p_surr<0.1&(rho-rho_min)>0]
full9<-ccm[,.N,by=meteo_stat][N==9,meteo_stat]
st<-ccm[meteo_stat%in%full9,.(nsig=sum(sig)),by=meteo_stat]
sel_st<-st[nsig>=1&nsig<=3,meteo_stat]
recs<-ccm[meteo_stat%in%sel_st & sig==TRUE,.(meteo_stat,tp,nsig=NA)]
recs<-merge(recs[,.(meteo_stat,tp)],st,by="meteo_stat")

if(file.exists(CACHE_2V)){
  cf2<-readRDS(CACHE_2V); cat("读取二变量缓存\n")
}else{
  fs<-list.files(STAT_DIR,pattern="^hcsif_station_[0-9]{4}\\.csv$",full.names=TRUE)
  sif<-rbindlist(lapply(fs,fread),fill=TRUE)[meteo_stat%in%sel_st,c("meteo_stat","year","doy","date",SIF_COL),with=FALSE]
  vpd<-readRDS(CACHE)[n_days>=MIN_DAYS_WIN & meteo_stat%in%sel_st]
  d<-merge(sif,vpd,by=c("meteo_stat","year","doy")); setorder(d,meteo_stat,year,doy)
  d[, idx8:=as.integer(round(as.numeric(as.Date(as.character(date),format="%Y%m%d")-as.Date("2000-01-01"))/8))]
  znorm<-function(x){s<-sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)return(rep(NA_real_,length(x)));(x-mean(x,na.rm=TRUE))/s}
  for(v in c(SIF_COL,"vpd_mean")) d[, paste0(v,"_dt"):=znorm(get(v)), by=meteo_stat]
  get_2v<-function(dd,tp_x){
    dd<-copy(dd); dd[, x_lag:=shift(vpd_mean_dt,n=tp_x,type="lag"), by=year]; dd<-dd[!is.na(idx8)]
    tmp<-dd[,.(idx8,y=get(SIF_DT),x=x_lag)]
    full<-tmp[data.table(idx8=seq.int(min(dd$idx8),max(dd$idx8))),on="idx8"]; setorder(full,idx8)
    if(full[!is.na(y)&!is.na(x),.N]<MIN_PTS) return(NULL)
    df<-data.frame(time=seq_len(nrow(full)),sif=full$y,vpd=full$x)
    set.seed(SEED_BASE+as.integer(dd$meteo_stat[1])+tp_x*1000L)
    # 二变量: 状态空间=[sif, vpd], E=2, embedded=TRUE, 读 ∂sif/∂vpd
    sm<-tryCatch(rEDM::SMap(dataFrame=df,E=2,theta=2,lib=paste("1",nrow(df)),pred=paste("1",nrow(df)),
                  columns=c("sif","vpd"),target="sif",embedded=TRUE),error=function(e)NULL)
    if(is.null(sm))return(NULL)
    co<-as.data.table(sm$coefficients); cc<-which(grepl("vpd",colnames(co),ignore.case=TRUE))[1]
    if(is.na(cc))return(NULL)
    tt<-co[["time"]]; cv<-co[[cc]]; ok<-!is.na(cv)&!is.nan(cv)&tt>=1&tt<=nrow(full); cv[ok]
  }
  cf2<-rbindlist(lapply(split(recs,seq_len(nrow(recs))),function(r){
    v<-get_2v(d[meteo_stat==r$meteo_stat],r$tp); if(is.null(v)||!length(v))return(NULL)
    data.table(meteo_stat=r$meteo_stat,tp=r$tp,nsig=r$nsig,coef=v)}),fill=TRUE)
  saveRDS(cf2,CACHE_2V)
}
cat("二变量: 系数点",nrow(cf2)," 记录",cf2[,uniqueN(paste(meteo_stat,tp))],"\n")

# 记录级方向: 单驱动(1v) vs 二变量(2v)
cf1<-readRDS(file.path(OUT,"sig1to3_coef_full.rds"))
dm<-function(x)fifelse(x>=0.75,"Promote",fifelse(x<=0.25,"Inhibit","Ambiguous"))
r1<-cf1[,.(mean1=mean(coef),med1=median(coef),fp1=mean(coef>0)),by=.(meteo_stat,tp,nsig)]
r2<-cf2[,.(mean2=mean(coef),med2=median(coef),fp2=mean(coef>0)),by=.(meteo_stat,tp)]
m<-merge(r1,r2,by=c("meteo_stat","tp"))
m[, `:=`(dir1=dm(fp1),dir2=dm(fp2),sign1=fifelse(mean1>0,"Promote","Inhibit"),sign2=fifelse(mean2>0,"Promote","Inhibit"))]
cat("\n配对记录:",nrow(m),"\n")
cat("平均系数相关 r =",round(cor(m$mean1,m$mean2),3),"\n")
cat("均值符号一致率 =",round(mean(m$sign1==m$sign2),3),"\n")
cat("主导(>=75%)方向一致率 =",round(mean(m$dir1==m$dir2),3),"\n")
cat("\n主导方向构成:\n单驱动:",paste(names(table(m$dir1)),table(m$dir1)),"\n二变量:",paste(names(table(m$dir2)),table(m$dir2)),"\n")
fwrite(m,file.path(OUT,"smap_1v_2v_records.csv"))

# 图
dircol<-c(Promote="#C1121F",Inhibit="#1A6FBF",Ambiguous="#9E9E9E")
pa<-ggplot(m,aes(mean1,mean2,color=dir2))+geom_abline(slope=1,linetype="dashed",color="grey50")+
  geom_hline(yintercept=0,color="grey80")+geom_vline(xintercept=0,color="grey80")+
  geom_point(size=1.2,alpha=.7)+scale_color_manual(values=dircol,name="2-var dominant dir")+
  labs(title="(a) Per-record mean S-map coefficient: 1-variable vs 2-variable",
       subtitle=sprintf("n=%d significant records; r=%.2f; sign-agreement=%.0f%%",nrow(m),cor(m$mean1,m$mean2),100*mean(m$sign1==m$sign2)),
       x="1-variable (VPD delay-embedding) mean coef",y="2-variable [SIF,VPD] mean coef")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
cc<-rbind(data.table(ver="1-variable",dir=m$dir1),data.table(ver="2-variable [SIF,VPD]",dir=m$dir2))[,.N,by=.(ver,dir)]
cc[, dir:=factor(dir,levels=c("Promote","Inhibit","Ambiguous"))]
pb<-ggplot(cc,aes(ver,N,fill=dir))+geom_col()+geom_text(aes(label=N),position=position_stack(vjust=.5),color="white",size=3.5)+
  scale_fill_manual(values=dircol,name="Dominant dir (>=75%)")+
  labs(title="(b) Direction category counts by S-map version",x=NULL,y="# significant records")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"smap_1v_vs_2v.png"),pa/pb,width=10,height=10,dpi=300)
cat("-> smap_1v_vs_2v.png\n")
