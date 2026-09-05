#!/usr/bin/env Rscript
# 合并 tp0-8(norm_surr) + tp9-12(补算), 在 tp0-12 上重做 宽松判据(p<0.1&drho>0) + 符号型分类。
# 输出 4类+mixed 分布, 并与 tp0-8 基线对比。
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify"); SIF<-"SIF_buf1000_dt"
latest<-function(dir){f<-list.files(dir,pattern="vpd_.*[0-9]\\.rds$",full.names=TRUE);tail(sort(f[!grepl("_raw",f)]),1)}
classify6<-function(c){s<-sign(c); if(any(is.na(s))||any(s==0))return("mixed")
  n<-sum(diff(s)!=0); if(n==0)return(if(s[1]>0)"always_promote" else "always_inhibit")
  if(n==1)return(if(s[1]<0)"inhibit_promote" else "promote_inhibit"); "mixed"}
lev<-c("always_promote","promote_inhibit","inhibit_promote","always_inhibit","mixed")

d08<-as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
d912<-as.data.table(readRDS(latest("data_proc/ccm_hcsif_buf1000_norm_tp9to12")))
common<-intersect(names(d08),names(d912))
d<-rbind(d08[,..common],d912[,..common])
d<-d[y_var==SIF & x_var=="vpd_mean_dt"][order(meteo_stat,tp)]
d[, sig:=p_surr<0.1 & (rho-rho_min)>0]

analyze<-function(dd,tpmax,tag){
  x<-dd[tp<=tpmax]; ntp<-tpmax+1
  full<-x[,.N,by=meteo_stat][N==ntp,meteo_stat]; x<-x[meteo_stat%in%full]
  keep<-x[,.(any=any(sig)),by=meteo_stat][any==TRUE,meteo_stat]
  cls<-x[meteo_stat%in%keep, .(clsA=classify6(mean_coef)), by=meteo_stat]
  tb<-cls[,.N,by=clsA]; setkey(tb,clsA)
  out<-setNames(tb[lev,on="clsA",N],lev); out[is.na(out)]<-0
  cat(sprintf("\n=== %s (tp 0-%d): 完整站=%d, 宽松保留=%d ===\n",tag,tpmax,length(full),length(keep)))
  print(data.frame(class=lev,n=as.integer(out)))
  data.table(scenario=tag,tpmax=tpmax,retained=length(keep),
             always_promote=out[1],promote_inhibit=out[2],inhibit_promote=out[3],
             always_inhibit=out[4],mixed=out[5])
}
r08<-analyze(d,8,"baseline tp0-8")
r12<-analyze(d,12,"extended tp0-12")
res<-rbind(r08,r12); fwrite(res,file.path(OUT,"tp12_class_distribution.csv")); print(res)

# 图: 两种tp下 分布对比
lab<-c(always_promote="Always promote",promote_inhibit="Promote->Inhibit",
       inhibit_promote="Inhibit->Promote",always_inhibit="Always inhibit",mixed="Mixed")
col<-setNames(c("#C1121F","#F4A261","#74C6E8","#1A6FBF","#9E9E9E"),lab)
m<-melt(res,id.vars=c("scenario","tpmax","retained"),measure.vars=lev,variable.name="class",value.name="n")
m[, class:=factor(lab[as.character(class)],levels=lab)]
m[, scen:=factor(scenario,levels=c("baseline tp0-8","extended tp0-12"),
    labels=c(sprintf("tp 0-8 (n=%d)",r08$retained),sprintf("tp 0-12 (n=%d)",r12$retained)))]
p<-ggplot(m,aes(scen,n,fill=class))+geom_col(position="dodge")+
  geom_text(aes(label=n),position=position_dodge(.9),vjust=-.3,size=3.5)+
  scale_fill_manual(values=col,name="Sign-type")+
  labs(title="VPD->SIF sign-type distribution: extending CCM lag tp 0-8 -> 0-12 (1km, normalize)",
       subtitle="Loose criteria p_surr<0.1 & drho>0; classification over the full sign sequence",
       x=NULL,y="# stations")+theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"tp12_class_distribution.png"),p,width=10,height=6,dpi=300)
cat("\n-> tp12_class_distribution.png\n")
