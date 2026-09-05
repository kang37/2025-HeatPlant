#!/usr/bin/env Rscript
# =============================================================================
# smap_timevarying_coef.R
#   因果方向= S-map 局部线性系数 ∂SIF/∂VPD, 每个时间点一个值; mean_coef 是其多年平均。
#   本脚本取 10 个显著 站点×tp, 还原逐时间点 S-map 系数, 看方向是否随年份翻转。
#   预处理=归一化(Ushio式), 与 ccm_hcsif_buf1000_norm_surr 一致。全英文标签。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(rEDM); library(ggplot2)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
BUF_R<-1000L; SIF_COL<-sprintf("SIF_buf%d",BUF_R)
STAT_DIR<-file.path(PROJ,"data_raw/hcsif/station_v3")
CACHE<-file.path(PROJ,"data_raw/hcsif/vpd_8day_cache.rds")
MIN_PTS<-40L; MIN_DAYS_WIN<-5L; SEED_BASE<-20260824L
set.seed(1)

# ---- 选 10 个显著 站×tp ----
ccm<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
ccm<-ccm[y_var==paste0(SIF_COL,"_dt") & x_var=="vpd_mean_dt"]
ccm[, sig:=p_surr<0.1 & (rho-rho_min)>0]
sigset<-ccm[sig==TRUE]
# 跨 mean_coef 谱选: 3强正 / 3强负 / 4近零(最可能翻转), 均要求rho较高以保证信号
setorder(sigset,-rho)
pick<-rbind(
  head(sigset[mean_coef>0.03],3),
  head(sigset[mean_coef< -0.03],3),
  head(sigset[abs(mean_coef)<0.02],4)
)[!duplicated(paste(meteo_stat,tp))]
sel<-pick[,.(meteo_stat,tp,E,mean_coef,rho,p_surr)]
cat("选中10对:\n"); print(sel)

# ---- 载入这些站的 SIF + VPD, 归一化 ----
ids<-unique(sel$meteo_stat)
fs<-list.files(STAT_DIR,pattern="^hcsif_station_[0-9]{4}\\.csv$",full.names=TRUE)
sif<-rbindlist(lapply(fs,fread),fill=TRUE)[meteo_stat %in% ids, c("meteo_stat","year","doy","date",SIF_COL),with=FALSE]
vpd<-readRDS(CACHE)[n_days>=MIN_DAYS_WIN & meteo_stat %in% ids]
d<-merge(sif,vpd,by=c("meteo_stat","year","doy")); setorder(d,meteo_stat,year,doy)
d[, idx8:=as.integer(round(as.numeric(as.Date(as.character(date),format="%Y%m%d")-as.Date("2000-01-01"))/8))]
znorm<-function(x){s<-sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)return(rep(NA_real_,length(x)));(x-mean(x,na.rm=TRUE))/s}
for(v in c(SIF_COL,"vpd_mean")) d[, paste0(v,"_dt"):=znorm(get(v)), by=meteo_stat]

# ---- 逐对: 还原时间点 S-map 系数 ----
get_coefs<-function(sid,tp_x,best_E){
  dd<-d[meteo_stat==sid]
  dd[, x_lag:=shift(vpd_mean_dt,n=tp_x,type="lag"), by=year]
  dd<-dd[!is.na(idx8)]
  tmp<-dd[,.(idx8,y=get(paste0(SIF_COL,"_dt")),x=x_lag)]
  full<-tmp[data.table(idx8=seq.int(min(dd$idx8),max(dd$idx8))),on="idx8"]; setorder(full,idx8)
  if(full[!is.na(y)&!is.na(x),.N]<MIN_PTS) return(NULL)
  df<-data.frame(time=seq_len(nrow(full)),sif=full$y,heat=full$x)
  set.seed(SEED_BASE+as.integer(sid)+tp_x*1000L)
  sm<-tryCatch(rEDM::SMap(dataFrame=df,E=best_E,theta=2,lib=paste("1",nrow(df)),
                pred=paste("1",nrow(df)),columns="heat",target="sif",embedded=FALSE),
               error=function(e)NULL)
  if(is.null(sm)) return(NULL)
  co<-as.data.table(sm$coefficients)
  cc<-which(grepl("heat",colnames(co),ignore.case=TRUE))[1]; if(is.na(cc))cc<-3L
  tt<-co[["time"]]                              # SMap的预测时间索引(对应df$time)
  coef<-co[[cc]]                                # ∂sif/∂heat(t-0): 同期偏导=因果方向
  ok<-!is.na(coef)&!is.nan(coef)&tt>=1&tt<=nrow(full)
  data.table(meteo_stat=sid,tp=tp_x,idx8=full$idx8[tt[ok]],coef=coef[ok])
}
res<-rbindlist(lapply(seq_len(nrow(sel)),function(i)
  get_coefs(sel$meteo_stat[i],sel$tp[i],sel$E[i])),fill=TRUE)
# idx8 -> 日期/年份
res[, date:=as.Date("2000-01-01")+idx8*8]
res[, year:=as.integer(format(date,"%Y"))]
res<-merge(res,sel[,.(meteo_stat,tp,mean_coef,rho)],by=c("meteo_stat","tp"))
res[, panel:=sprintf("stn %d, tp=%d\nmean=%.3f, rho=%.2f",meteo_stat,tp,mean_coef,rho)]

# 每对: 系数为正/负的时间比例, 是否翻转
flip<-res[,.(n=.N,frac_pos=mean(coef>0),frac_neg=mean(coef<0),
             mean_c=mean(coef),sd_c=sd(coef),
             flips=sum(diff(sign(coef))!=0,na.rm=TRUE)),by=.(meteo_stat,tp,mean_coef)]
cat("\n=== 每对: S-map系数随时间的符号构成与翻转次数 ===\n"); print(flip)
fwrite(res,file.path(OUT,"smap_timevarying_coef_series.csv"))
fwrite(flip,file.path(OUT,"smap_timevarying_coef_summary.csv"))

# ---- 图: 逐对 时间点S-map系数 时间序列 ----
p<-ggplot(res,aes(date,coef))+
  geom_hline(yintercept=0,color="grey40",linewidth=.4)+
  geom_hline(aes(yintercept=mean_coef),color="#2166AC",linetype="dashed",linewidth=.5)+
  geom_line(color="grey70",linewidth=.3)+
  geom_point(aes(color=coef>0),size=1)+
  scale_color_manual(values=c(`TRUE`="#C1121F",`FALSE`="#1A6FBF"),
                     labels=c(`TRUE`="Promote (dSIF/dVPD>0)","FALSE"="Inhibit (<0)"),name=NULL)+
  facet_wrap(~panel,scales="free",ncol=2)+
  labs(title="Time-varying S-map coefficient dSIF/dVPD for 10 significant station x tp",
       subtitle="Each point = one 8-day time step; dashed blue = multi-year mean (the 'direction' used so far). Sign flips across years are common.",
       x="Date",y="S-map coefficient (local dSIF/dVPD)")+
  theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),
    legend.position="top",strip.text=element_text(size=8))
ggsave(file.path(OUT,"smap_timevarying_coef.png"),p,width=12,height=13,dpi=300)
cat("\n-> smap_timevarying_coef.png\n")
