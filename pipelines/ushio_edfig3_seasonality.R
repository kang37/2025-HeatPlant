#!/usr/bin/env Rscript
# =============================================================================
# ushio_edfig3_seasonality.R
#   复现 Ushio et al. 2018 (Nature) 处理季节性的思路(Extended Data Fig.3):
#   用"保季节替代零模型"检验 CCM —— 观测 rho 若不超过季节替代零模型, 视为无因果。
#   我们的替代=年块置换(保留DOY季节循环, 跨年份置换), 即 season-preserving null。
#   对两种预处理各做一遍: 归一化(norm) 与 线性去趋势(detrend)。
#   全英文标签。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
latest<-function(dir){f<-list.files(dir,pattern="vpd_.*rds$",full.names=TRUE);tail(sort(f[!grepl("_raw",f)]),1)}

load_one<-function(dir,proc){
  d<-as.data.table(readRDS(latest(dir)))
  d<-d[y_var=="SIF_buf1000_dt" & x_var=="vpd_mean_dt"]
  d[, `:=`(proc=proc, drho=rho-rho_min, excess=rho-rho_null95,
           sig=p_surr<0.1 & (rho-rho_min)>0)]
  d
}
dn<-load_one("data_proc/ccm_hcsif_buf1000_norm_surr","Normalize")
dd<-load_one("data_proc/ccm_hcsif_buf1000_surr","Detrend")
d<-rbind(dn,dd)
d[, proc:=factor(proc,levels=c("Normalize","Detrend"))]

cat("=== 各预处理: (站,tp)记录数 / 站点数 ===\n")
print(d[, .(records=.N, stations=uniqueN(meteo_stat)), by=proc])

# 每站峰值tp(观测rho最大处)——用于站点级检验
peak<-d[, .SD[which.max(rho)], by=.(proc,meteo_stat)]

# ---- 定量: 季节性贡献 ----
cat("\n=== 季节替代零模型能达到的CCM技巧(纯季节性) vs 观测 ===\n")
summ<-d[, .(
  rho_obs_med   = median(rho,na.rm=TRUE),
  rho_null_med  = median(rho_null95,na.rm=TRUE),   # 纯季节性可达的95%分位
  frac_below_null = mean(rho <= rho_null95, na.rm=TRUE),   # 观测被季节替代淹没的比例
  frac_sig      = mean(sig,na.rm=TRUE)
), by=proc]
print(summ)
cat("\n站点级(峰值tp)通过季节替代检验(p_surr<0.1 & drho>0)的站数:\n")
print(peak[, .(pass=sum(sig), total=.N), by=proc])
fwrite(summ, file.path(OUT,"edfig3_seasonality_summary.csv"))

# ---- 正确的季节性检验 ----
# 注意: rho_null95 是季节替代零模型的95%分位 = 显著性阈值(非零模型中心)。
# 季节性对CCM的影响 = 该阈值有多高(纯季节性能造出多强的交叉映射技巧)。
# 检验: 通过季节替代的显著链接是否超过假阳性率(p_surr<0.1 → 期望10%)?
cat("\n=== 季节性检验: 显著链接 vs 假阳性率(binomial, H0: 通过率=0.10) ===\n")
for(p in levels(d$proc)){
  x<-d[proc==p]; k<-sum(x$sig); n<-nrow(x)
  bt<-binom.test(k, n, 0.10, alternative="greater")
  cat(sprintf("%-10s 显著(站,tp)=%d/%d=%.1f%%  binom p=%.3g (H0=10%%)  → %s\n",
      p, k, n, 100*k/n, bt$p.value,
      ifelse(bt$p.value<0.05,"高于假阳性率","不超过假阳性率")))
}
cat("\n季节替代95%阈值 vs 观测rho 中位数(体现纯季节性可达的CCM技巧):\n")
print(d[, .(median_obs_rho=median(rho), median_seasonal_threshold=median(rho_null95),
            frac_obs_below_threshold=mean(rho<=rho_null95)), by=proc])

# =============================================================================
# 图 (ED Fig.3 精神: 观测 vs 季节替代零模型)
# =============================================================================
th<-theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"))

# (a) 观测rho vs 季节替代95% 散点 + 1:1线 (站点峰值tp)
pa<-ggplot(peak,aes(rho_null95,rho,color=sig))+
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey40")+
  geom_point(size=1.3,alpha=.7)+facet_wrap(~proc)+
  scale_color_manual(values=c(`TRUE`="#C1121F",`FALSE`="grey65"),
                     labels=c(`TRUE`="Passes surrogate","FALSE"="Not (<= seasonal null)"),name=NULL)+
  labs(title="(a) Observed CCM skill vs seasonal-surrogate null (per station, peak-rho tp)",
       subtitle="Points below the 1:1 line = seasonality alone matches/exceeds the observed skill (no genuine causation)",
       x="Seasonal-surrogate 95% rho (pure seasonality)",y="Observed rho")+th+
  theme(legend.position="top")

# (b) 观测rho 与 季节替代95%阈值 的分布对比 (所有站×tp)
dl<-melt(d[,.(proc,rho,rho_null95)],id.vars="proc",variable.name="which",value.name="val")
dl[, which:=factor(which,levels=c("rho","rho_null95"),labels=c("Observed rho","Seasonal-surrogate 95% threshold"))]
pb<-ggplot(dl,aes(val,fill=which))+geom_density(alpha=.5)+facet_wrap(~proc)+
  scale_fill_manual(values=c(`Observed rho`="#1A6FBF",`Seasonal-surrogate 95% threshold`="#F4A261"),name=NULL)+
  labs(title="(b) Observed skill vs the seasonal-null 95% threshold (all station x tp)",
       subtitle="The threshold is what pure seasonality can reach; most observed rho fall below it",
       x="Cross-map skill rho",y="Density")+th+theme(legend.position="top")

# (c) 超出季节零模型的余量 excess=rho-rho_null95 分布
pc<-ggplot(d,aes(excess,fill=proc))+geom_density(alpha=.5)+
  geom_vline(xintercept=0,linetype="dashed",color="grey40")+
  scale_fill_manual(values=c(Normalize="#4575B4",Detrend="#D73027"),name="Preprocessing")+
  labs(title="(c) Excess skill over the seasonal null (rho - surrogate95)",
       subtitle="Mass right of 0 = signal beyond seasonality; mass left = seasonality dominates",
       x="rho - seasonal-surrogate 95%",y="Density")+th+theme(legend.position="top")

# (d) 季节性剥离前后的"收敛站"数量: naive(仅drho>0) vs 季节替代校正后
bars<-rbind(
  d[, .(n=uniqueN(meteo_stat[drho>0])), by=proc][, crit:="Naive: converges (drho>0)"],
  d[, .(n=uniqueN(meteo_stat[sig])),   by=proc][, crit:="After seasonal-surrogate test"])
bars[, crit:=factor(crit,levels=c("Naive: converges (drho>0)","After seasonal-surrogate test"))]
pd<-ggplot(bars,aes(proc,n,fill=crit))+geom_col(position="dodge")+
  geom_text(aes(label=n),position=position_dodge(.9),vjust=-.3,size=3.5)+
  scale_fill_manual(values=c("Naive: converges (drho>0)"="#9E9E9E","After seasonal-surrogate test"="#C1121F"),name=NULL)+
  labs(title="(d) Stations retained: before vs after correcting for seasonality",
       subtitle="Seasonal-surrogate null removes the seasonality-driven false positives",
       x="Preprocessing",y="# stations (>=1 tp)")+th+theme(legend.position="top")

ggsave(file.path(OUT,"ushio_edfig3_seasonality.png"), (pa/pb)|(pc/pd),
       width=17, height=11, dpi=300)
cat("\n-> ushio_edfig3_seasonality.png\n")
