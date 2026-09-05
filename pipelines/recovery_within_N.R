#!/usr/bin/env Rscript
# 固定视野 recovery: "N个8天周期内恢复了多少" —— 逐点残余 + 衰减模型两版本。
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
d<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
d<-d[y_var=="SIF_buf1000_dt"&x_var=="vpd_mean_dt"][order(meteo_stat,tp)]
SEL<-c(`A. Fast monotonic`=51814,`B. Rebound/overshoot`=51567,
       `C. Slow recovery`=53487,`D. Persistent`=54943)
sub<-d[meteo_stat%in%SEL]; sub[,panel:=factor(names(SEL)[match(meteo_stat,SEL)],levels=names(SEL))]

# 每站: 峰值, kappa, 逐点/模型 recovery 曲线
curves<-sub[,{
  b<-mean_coef;ab<-abs(b);i<-which.max(ab);tpstar<-tp[i];M<-ab[i]
  aft<-which(tp>=tpstar); fitd<-data.table(x=tp[aft]-tpstar,y=log(pmax(ab[aft],1e-4)))
  kap<- -coef(lm(y~x,fitd))[2]
  N<-0:(8-tpstar)
  rec_pt<-pmax(0,1-ab[aft]/M)               # 逐点残余
  rec_md<-1-exp(-kap*N)                      # 衰减模型
  rbind(data.table(N=N,rec=rec_pt,type="Point (1-|beta|/peak)"),
        data.table(N=N,rec=rec_md,type="Model (1-exp(-kappa N))"))
},by=.(meteo_stat,panel)]

# 汇总表: N=1,2,3 处的恢复比例(逐点)
tab<-sub[,{
  b<-mean_coef;ab<-abs(b);i<-which.max(ab);tpstar<-tp[i];M<-ab[i];aft<-which(tp>=tpstar)
  g<-function(n) if(tpstar+n<=8) round(1-ab[tp==tpstar+n]/M,2) else NA_real_
  .(tpstar=tpstar,Rec_1=g(1),Rec_2=g(2),Rec_3=g(3))
},by=.(meteo_stat,panel)]
cat("=== N个周期内恢复比例(逐点) ===\n"); print(tab)
fwrite(tab,file.path(OUT,"recovery_within_N_table.csv"))

p<-ggplot(curves,aes(N,rec,color=type))+
  geom_hline(yintercept=c(0,1),color="grey80",linewidth=.3)+
  geom_line(linewidth=.7)+geom_point(size=1.6)+
  geom_vline(xintercept=c(1,2,3),linetype="dotted",color="grey60")+
  scale_color_manual(values=c("Point (1-|beta|/peak)"="#C1121F","Model (1-exp(-kappa N))"="#2166AC"),name=NULL)+
  scale_y_continuous(limits=c(0,1.05))+scale_x_continuous(breaks=0:8)+
  facet_wrap(~panel,ncol=2)+
  labs(title="Recovery achieved within N periods (1 period = 8 days), counted from peak lag tp*",
       subtitle="Rec_N = fraction of the peak VPD->SIF effect dissipated by N periods after tp*. 1=fully recovered, 0=none. No censoring.",
       x="Horizon N (number of 8-day periods after peak)",y="Recovery fraction Rec_N")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),legend.position="top",
    strip.text=element_text(face="bold"),plot.subtitle=element_text(size=9,color="grey35"))
ggsave(file.path(OUT,"recovery_within_N.png"),p,width=12,height=8,dpi=300)
cat("-> recovery_within_N.png\n")
