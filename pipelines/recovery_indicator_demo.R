#!/usr/bin/env Rscript
# =============================================================================
# recovery_indicator_demo.R
#   以4个站点为例, 图解"考虑tp的recovery指标"如何从 beta(tp)=dSIF/dVPD 曲线构建。
#   beta(tp) 序列 = SIF对VPD脉冲的滞后响应(脉冲响应函数); recovery = 该响应随tp的衰减/回归。
#   全英文标签。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
ALPHA<-1/exp(1)   # 衰减阈值: 效应降到峰值的 1/e 视为"基本恢复"

d<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
d<-d[y_var=="SIF_buf1000_dt" & x_var=="vpd_mean_dt"][order(meteo_stat,tp)]
d[,sig:=p_surr<0.1&(rho-rho_min)>0]

SEL<-c(`A. Fast monotonic recovery`=51814,
       `B. Rebound / overshoot`=51567,
       `C. Slow recovery`=53487,
       `D. Persistent (censored)`=54943)
sub<-d[meteo_stat%in%SEL]
sub[,panel:=factor(names(SEL)[match(meteo_stat,SEL)],levels=names(SEL))]

# ---- 逐站计算指标 ----
calc<-function(x){
  b<-x$mean_coef; ab<-abs(b); tp<-x$tp
  i<-which.max(ab); tpstar<-tp[i]; M<-ab[i]
  aft<-which(tp>=tpstar)
  bel<-aft[ab[aft]<=ALPHA*M]
  Trec<-if(length(bel)) tp[bel[1]]-tpstar else NA_integer_
  tprec<-if(length(bel)) tp[bel[1]] else NA_integer_
  # 指数衰减率 kappa: log|beta| ~ (tp-tp*), tp>=tp*
  fitd<-data.table(x=tp[aft]-tpstar, y=log(pmax(ab[aft],1e-4)))
  kap<-if(nrow(fitd)>=3) -coef(lm(y~x,fitd))[2] else NA_real_
  s<-sign(b[aft]); reb<-any(diff(s)!=0)
  tflip<-if(reb) tp[aft][which(diff(s)!=0)[1]+1] else NA_integer_
  list(tpstar=tpstar,M=M,Trec=Trec,tprec=tprec,kappa=as.numeric(kap),
       rebound=reb,tflip=tflip,netsum=sum(b))
}
ann<-sub[,calc(.SD),by=.(meteo_stat,panel)]
ann[,`:=`(alphaM=ALPHA*M)]
cat("=== 4站指标 ===\n"); print(ann[,.(panel,meteo_stat,tpstar,M=round(M,3),
      alphaM=round(alphaM,3),Trec,kappa=round(kappa,2),rebound,netsum=round(netsum,3))])

# 指数拟合曲线(用于叠加): 幅度衰减 * 初始符号
fitcurve<-ann[,{
  ss<-sign(sub[meteo_stat==.BY$meteo_stat][tp==tpstar,mean_coef])
  tps<-seq(tpstar,8,.1)
  data.table(tp=tps, yfit=ss*M*exp(-kappa*(tps-tpstar)))
},by=.(meteo_stat,panel)]

# 注释文本
ann[,lab:=sprintf("peak |beta|=%.2f @ tp*=%d\nT_rec=%s (beta down to (1/e)*peak)\nkappa=%.2f /8d  half-life=%.1f steps\nrebound=%s   sum-beta=%.2f",
    M,tpstar,ifelse(is.na(Trec),"censored",as.character(Trec)),kappa,log(2)/kappa,
    ifelse(rebound,sprintf("yes @tp=%d",tflip),"no"),netsum)]

# ---- 图 ----
band<-ann[,.(panel,ymin=-alphaM,ymax=alphaM)]
p<-ggplot(sub,aes(tp,mean_coef))+
  geom_rect(data=band,aes(x=NULL,y=NULL,ymin=ymin,ymax=ymax,xmin=-Inf,xmax=Inf),
            fill="grey80",alpha=.45,inherit.aes=FALSE)+
  geom_hline(yintercept=0,color="grey45",linewidth=.4)+
  geom_line(data=fitcurve,aes(tp,yfit),color="#2166AC",linetype="dashed",linewidth=.5,inherit.aes=FALSE)+
  geom_line(color="grey55",linewidth=.5)+
  geom_point(aes(color=sig),size=2.6)+
  # 峰值tp*
  geom_point(data=merge(ann[,.(panel,meteo_stat,tpstar)],sub,by.x=c("panel","meteo_stat","tpstar"),
                        by.y=c("panel","meteo_stat","tp")),
             aes(tpstar,mean_coef),shape=8,size=4,color="black",inherit.aes=FALSE)+
  # 恢复点竖线
  geom_vline(data=ann[!is.na(tprec)],aes(xintercept=tprec),color="#C1121F",linetype="dotted",linewidth=.6)+
  # 恢复跨度箭头
  geom_segment(data=ann[!is.na(tprec)],aes(x=tpstar,xend=tprec,y=M*1.05,yend=M*1.05),
               arrow=arrow(length=unit(.15,"cm"),ends="both"),color="#C1121F",linewidth=.5,inherit.aes=FALSE)+
  geom_text(data=ann[!is.na(tprec)],aes(x=(tpstar+tprec)/2,y=M*1.18,label=sprintf("T_rec=%d",Trec)),
            color="#C1121F",size=3,inherit.aes=FALSE)+
  # 反弹标记
  geom_point(data=ann[rebound==TRUE&!is.na(tflip)],aes(tflip,0),shape=21,size=3.5,
             fill="#F4A261",color="black",inherit.aes=FALSE)+
  geom_text(data=ann,aes(x=8,y=-Inf,label=lab),hjust=1,vjust=-.15,size=2.7,color="grey20",inherit.aes=FALSE)+
  scale_color_manual(values=c(`TRUE`="#111111",`FALSE`="grey70"),
                     labels=c(`TRUE`="significant tp","FALSE"="not sig"),name=NULL)+
  scale_x_continuous(breaks=0:8)+
  facet_wrap(~panel,scales="free_y",ncol=2)+
  labs(title="Constructing a tp-aware recovery indicator from the lag-response beta(tp)",
       subtitle="beta(tp)=S-map dSIF/dVPD at lag tp (an impulse-response). Grey band=+/-(1/e)*peak; star=peak lag tp*; red=recovery point; orange=sign flip (rebound); blue dashed=exp-decay fit (kappa).",
       x="CCM lag tp (each step = 8 days)   (time since VPD pulse)",
       y="beta(tp) = dSIF/dVPD  (effect on SIF)")+
  theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),
    legend.position="top",strip.text=element_text(face="bold"),
    plot.subtitle=element_text(size=8.5,color="grey35"))
ggsave(file.path(OUT,"recovery_indicator_demo.png"),p,width=13,height=9,dpi=300)
cat("\n-> recovery_indicator_demo.png\n")
