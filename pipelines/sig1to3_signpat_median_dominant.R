#!/usr/bin/env Rscript
# =============================================================================
# sig1to3_signpat_median_dominant.R
#   重做 sig1to3_lagpos_signpat 的两个版本: 因果方向不按 S-map 系数"均值"符号,
#   而按 (1) 中位数符号, (2) 主导符号(正或负>=75%, 否则 Ambiguous)。
#   需要每个显著 站×tp 的逐时间点 S-map 系数 → 重跑 S-map。归一化预处理。全英文。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(rEDM); library(ggplot2); library(patchwork)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
BUF_R<-1000L; SIF_COL<-sprintf("SIF_buf%d",BUF_R); SIF_DT<-paste0(SIF_COL,"_dt")
STAT_DIR<-file.path(PROJ,"data_raw/hcsif/station_v3"); CACHE<-file.path(PROJ,"data_raw/hcsif/vpd_8day_cache.rds")
MIN_PTS<-40L; MIN_DAYS_WIN<-5L; SEED_BASE<-20260824L
CACHE_COEF<-file.path(OUT,"sig1to3_coef_records.rds")

# ---- 显著记录(nsig 1-3 的站) ----
ccm<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
ccm<-ccm[y_var==SIF_DT & x_var=="vpd_mean_dt"]
ccm[, sig:=p_surr<0.1 & (rho-rho_min)>0]
full9<-ccm[,.N,by=meteo_stat][N==9,meteo_stat]
st<-ccm[meteo_stat%in%full9,.(nsig=sum(sig)),by=meteo_stat]
sel_st<-st[nsig>=1&nsig<=3,meteo_stat]
recs<-ccm[meteo_stat%in%sel_st & sig==TRUE,.(meteo_stat,tp,E,mean_coef)]
cat("站点(1-3 tp显著):",length(sel_st)," 显著记录:",nrow(recs),"\n")

if(file.exists(CACHE_COEF)){
  rc<-readRDS(CACHE_COEF); cat("读取缓存系数记录\n")
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
    dd<-copy(dd); dd[, x_lag:=shift(vpd_mean_dt,n=tp_x,type="lag"), by=year]
    dd<-dd[!is.na(idx8)]
    tmp<-dd[,.(idx8,y=get(SIF_DT),x=x_lag)]
    full<-tmp[data.table(idx8=seq.int(min(dd$idx8),max(dd$idx8))),on="idx8"]; setorder(full,idx8)
    if(full[!is.na(y)&!is.na(x),.N]<MIN_PTS) return(NULL)
    df<-data.frame(time=seq_len(nrow(full)),sif=full$y,heat=full$x)
    set.seed(SEED_BASE+as.integer(dd$meteo_stat[1])+tp_x*1000L)
    sm<-tryCatch(rEDM::SMap(dataFrame=df,E=best_E,theta=2,lib=paste("1",nrow(df)),
                  pred=paste("1",nrow(df)),columns="heat",target="sif",embedded=FALSE),error=function(e)NULL)
    if(is.null(sm))return(NULL)
    co<-as.data.table(sm$coefficients); cc<-which(grepl("heat",colnames(co),ignore.case=TRUE))[1]; if(is.na(cc))cc<-3L
    tt<-co[["time"]]; coef<-co[[cc]]; ok<-!is.na(coef)&!is.nan(coef)&tt>=1&tt<=nrow(full)
    coef[ok]
  }
  rc<-rbindlist(lapply(split(recs,seq_len(nrow(recs))),function(r){
    cf<-get_coefs(d[meteo_stat==r$meteo_stat],r$tp,r$E)
    if(is.null(cf)||!length(cf))return(NULL)
    data.table(meteo_stat=r$meteo_stat,tp=r$tp,mean_coef=r$mean_coef,
               n=length(cf),median_coef=median(cf),frac_pos=mean(cf>0))
  }),fill=TRUE)
  saveRDS(rc,CACHE_COEF)
}
cat("成功还原系数的记录:",nrow(rc),"/",nrow(recs),"\n")

# ---- 两种方向定义 ----
rc[, dir_mean   := fifelse(mean_coef>0,"Promote","Inhibit")]
rc[, dir_median := fifelse(median_coef>0,"Promote","Inhibit")]
rc[, dir_domin  := fifelse(frac_pos>=0.75,"Promote",
                    fifelse(frac_pos<=0.25,"Inhibit","Ambiguous"))]
cat("\n=== 三种方向定义的记录级构成 ===\n")
print(rc[,.(mean=sum(dir_mean=="Promote"),.N),by=.(dir_mean)][order(-N)])
cat("median: Promote=",sum(rc$dir_median=="Promote")," Inhibit=",sum(rc$dir_median=="Inhibit"),"\n")
cat("dominant75: Promote=",sum(rc$dir_domin=="Promote")," Inhibit=",sum(rc$dir_domin=="Inhibit"),
    " Ambiguous=",sum(rc$dir_domin=="Ambiguous"),"\n")
cat("均值与中位数方向不一致的记录:",sum(rc$dir_mean!=rc$dir_median),"\n")
fwrite(rc,file.path(OUT,"sig1to3_direction_defs.csv"))

# ---- 画图函数(与 sig1to3_lagpos_signpat 相同布局) ----
make_fig<-function(rc, dircol, title_tag, fname, cols){
  rc<-copy(rc); rc[, dir:=get(dircol)]
  # 每站 nsig 与 符号型
  info<-rc[,.(nsig=.N,
              signpat={u<-unique(dir)
                if("Ambiguous"%in%u && length(setdiff(u,"Ambiguous"))==0)"All ambiguous"
                else if(all(dir=="Promote"))"All promote"
                else if(all(dir=="Inhibit"))"All inhibit"
                else if("Ambiguous"%in%u)"Mixed w/ ambiguous"
                else "Sign-changing"}),by=meteo_stat]
  # (1) 显著tp落点
  tpc<-rc[,.N,by=.(tp,dir)]
  p1<-ggplot(tpc,aes(factor(tp),N,fill=dir))+geom_col()+
    scale_fill_manual(values=cols,name="Direction")+
    labs(title=sprintf("(a) Where significant lags fall — direction by %s",title_tag),
         subtitle=sprintf("%d significant (station,tp) records among 1-3 tp stations",nrow(rc)),
         x="Lag tp (each step = 8 days)",y="# significant records")+
    theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
  # (2) 符号型 x nsig
  splev<-intersect(c("All promote","All inhibit","Sign-changing","Mixed w/ ambiguous","All ambiguous"),unique(info$signpat))
  spcol<-c(`All promote`="#C1121F",`All inhibit`="#1A6FBF",`Sign-changing`="#F4A261",
           `Mixed w/ ambiguous`="#B39DDB",`All ambiguous`="#9E9E9E")
  sp<-info[,.N,by=.(nsig=factor(nsig),signpat=factor(signpat,levels=splev))]
  p2<-ggplot(sp,aes(nsig,N,fill=signpat))+geom_col()+
    geom_text(aes(label=N),position=position_stack(vjust=.5),size=3.3,color="white")+
    scale_fill_manual(values=spcol[splev],name="Sign pattern\n(sig tp only)")+
    labs(title=sprintf("(b) Per-station sign pattern — direction by %s",title_tag),
         subtitle="Split by number of significant tp",x="Number of significant tp",y="# stations")+
    theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
  ggsave(file.path(OUT,fname),p1/p2,width=10,height=9,dpi=300)
  cat("->",fname,"| signpat:",paste(sprintf("%s=%d",info[,.N,by=signpat]$signpat,info[,.N,by=signpat]$N),collapse=", "),"\n")
}
make_fig(rc,"dir_median","MEDIAN sign","sig1to3_lagpos_signpat_median.png",
         c(Promote="#C1121F",Inhibit="#1A6FBF"))
make_fig(rc,"dir_domin","DOMINANT sign (>=75%)","sig1to3_lagpos_signpat_dominant.png",
         c(Promote="#C1121F",Inhibit="#1A6FBF",Ambiguous="#9E9E9E"))
