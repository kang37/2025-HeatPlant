#!/usr/bin/env Rscript
# =============================================================================
# loose_classify_drivers.R
#   宽松判据(p_surr<0.1 & drho>0)保留的 383 站 —— 多种分类 + 自变量驱动分析
#   数据: 归一化(Ushio式)+年块替代 CCM  data_proc/ccm_hcsif_buf1000_norm_surr
#   自变量: 20_drivers_all.R 的 BLK_B 全量(地形/气候/地表/城市/社经/土壤)
#   全英文标签(非交互 Rscript 中文缺字形)。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ranger); library(nnet)
  library(sf); library(rnaturalearth)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
set.seed(42)
winz <- function(x,p=.01){q<-quantile(x,c(p,1-p),na.rm=TRUE);pmin(pmax(x,q[1]),q[2])}

# ---- 1. CCM + 宽松判据 -------------------------------------------------------
d <- as.data.table(readRDS("data_proc/ccm_hcsif_buf1000_norm_surr/ccm_hcsif_buf1000_vpd_20260825_0404.rds"))
d <- d[y_var=="SIF_buf1000_dt" & x_var=="vpd_mean_dt"]
d[, drho := rho - rho_min]
d[, sig  := p_surr < 0.1 & drho > 0]
full9 <- d[, .N, by=meteo_stat][N==9, meteo_stat]
d <- d[meteo_stat %in% full9][order(meteo_stat, tp)]
keep <- d[, .(any=any(sig)), by=meteo_stat][any==TRUE, meteo_stat]
d <- d[meteo_stat %in% keep]
cat(sprintf("宽松保留站点: %d\n", length(keep)))

# ---- 2. 站点级分类特征 -------------------------------------------------------
# 分类只看"显著"的 tp(宽松判据下显著的滞后步)——因果证据所在
feat <- d[, {
  sg <- which(sig)
  i_peak <- sg[which.max(rho[sg])]          # 显著tp中rho最大的
  list(
    n_sig    = length(sg),
    peak_tp  = tp[i_peak],                   # 响应滞后(8天步)
    rho_peak = rho[i_peak],                  # 因果强度
    mc_peak  = mean_coef[i_peak],            # 峰值tp的带符号S-map系数
    drho_peak= drho[i_peak],
    mc_sig_mean = mean(mean_coef[sg]),       # 显著tp的平均带符号系数
    frac_pos = mean(mean_coef[sg] > 0),      # 显著tp里促进占比
    # 全9tp符号序列(用于历史式符号型分类)
    seq = list(mean_coef)
  )
}, by=meteo_stat]
feat[, `:=`(latitude = d$latitude[match(meteo_stat,d$meteo_stat)],
            longitude= d$longitude[match(meteo_stat,d$meteo_stat)])]

# 分类A: 符号型(历史5+1类, 基于全9tp符号序列)
classify6 <- function(c){s<-sign(c)
  if(any(is.na(s))||any(s==0))return("mixed")
  nch<-sum(diff(s)!=0)
  if(nch==0)return(if(s[1]>0)"always_promote" else "always_inhibit")
  if(nch==1)return(if(s[1]<0)"inhibit_promote" else "promote_inhibit"); "mixed"}
feat[, clsA := map_chr(seq, classify6)]

# 分类B: 峰值响应二分类(方向 = 显著峰值tp的S-map符号)
feat[, clsB := ifelse(mc_peak > 0, "Promote", "Inhibit")]

# 分类C: 显著tp主导方向(净符号)
feat[, clsC := fifelse(frac_pos > .5, "Promote-dominant",
                fifelse(frac_pos < .5, "Inhibit-dominant", "Balanced"))]

cat("\n=== 分类A 符号型(全9tp) ===\n"); print(feat[, .N, by=clsA][order(-N)])
cat("\n=== 分类B 峰值方向 ===\n");     print(feat[, .N, by=clsB][order(-N)])
cat("\n=== 分类C 主导方向 ===\n");     print(feat[, .N, by=clsC][order(-N)])

# ---- 3. 自变量(18/20个, 6块) ------------------------------------------------
GLC_YEAR<-2020L; NTL_YEAR<-2020L
dem <- fread("data_raw/covariates_1km/dem_cell_1km.csv")[,.(elev=mean(elev,na.rm=TRUE)),by=stat_id]
rs  <- rbindlist(lapply(2000:2022,function(y) fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv",y))))
rs  <- rs[,.(rsds=mean(rsds_mean,na.rm=TRUE)),by=.(stat_id,year)][,.(rsds_mean=mean(rsds),rsds_sd=sd(rsds)),by=stat_id]
g   <- fread(sprintf("data_raw/covariates_1km/glc_station_%d.csv",GLC_YEAR))
LCn <- function(...) rowSums(as.matrix(g[,sprintf("LC%02d",c(...)),with=FALSE]))
glc <- data.table(stat_id=g$stat_id, imperv=g$LC23, forest=LCn(4:13), grass=LCn(17),
                  crop=LCn(0:3), water=LCn(27))
nt  <- fread(sprintf("data_raw/covariates_1km/ntl_station_%d.csv",NTL_YEAR))[,.(stat_id,ntl=log1p(ntl_mean))]
old <- as.data.table(readRDS("data_proc/output_10y_built_up_05_01/station_covariates.rds"))
setnames(old,"meteo_stat_id","stat_id"); old[,stat_id:=as.integer(stat_id)]
old <- old[,.(stat_id,koppen_group,precip_mean,soil_clay,soil_sand,pa_built_10y,cgi_score,
              pop_10y,pgdp_10y,road_density,building_footprint,mean_height,building_vol_density)]
cov <- Reduce(function(a,b) merge(a,b,by="stat_id",all=TRUE), list(dem,rs,glc,nt,old))

# 全量18/20变量(6块) —— 用于"完整集"分析(样本受旧表覆盖限制)
PRED_FULL <- list(terrain="elev", climate=c("rsds_mean","rsds_sd","precip_mean"),
             landcover=c("imperv","forest","grass","crop","water"),
             urban=c("ntl","road_density","building_footprint","mean_height","building_vol_density","pa_built_10y"),
             socioecon=c("pop_10y","pgdp_10y","cgi_score"), soil=c("soil_clay","soil_sand"))
# 良好覆盖的"新变量"子集(地形+辐射+地表+灯光) —— 覆盖全部383站
PRED_A <- list(terrain="elev", climate=c("rsds_mean","rsds_sd"),
             landcover=c("imperv","forest","grass","crop","water"), urban="ntl")
vars      <- unlist(PRED_FULL, use.names=FALSE)
vars_a    <- unlist(PRED_A,    use.names=FALSE)
cat(sprintf("\n完整变量集: %d个; 良好覆盖A集: %d个\n", length(vars), length(vars_a)))

cov <- unique(cov, by="stat_id")   # 去重(旧表偶有重复stat_id)
dat <- merge(feat[, .(stat_id=meteo_stat, latitude, longitude, n_sig, peak_tp, rho_peak,
                      mc_peak, mc_sig_mean, frac_pos, clsA, clsB, clsC)],
             cov, by="stat_id")
for (v in vars) dat[, (v):=winz(as.numeric(get(v)))]
cat(sprintf("并入后: %d 站; A集完整=%d, 完整集完整=%d\n", nrow(dat),
            sum(complete.cases(dat[,c(vars_a),with=FALSE])),
            sum(complete.cases(dat[,c(vars),with=FALSE]))))
fwrite(dat[, setdiff(names(dat),"seq"), with=FALSE], file.path(OUT,"loose_stations_features_covariates.csv"))

# ---- 4. 无监督聚类(分类D): 用CCM特征发现自然类型 ----------------------------
fk <- dat[, .(peak_tp, rho_peak, mc_peak, n_sig, frac_pos, drho=NA)]
fk[, drho:=NULL]
Z  <- scale(as.matrix(fk))
ss <- sapply(2:6, function(k){ km<-kmeans(Z,centers=k,nstart=25); km$tot.withinss })
km4 <- kmeans(Z, centers=4, nstart=50)
dat[, clsD := factor(km4$cluster)]
cat("\n=== 分类D k-means(k=4) 各簇均值(CCM特征) ===\n")
prof <- dat[, .(n=.N, peak_tp=mean(peak_tp), rho_peak=mean(rho_peak),
                mc_peak=mean(mc_peak), n_sig=mean(n_sig), frac_pos=mean(frac_pos)), by=clsD][order(clsD)]
print(prof)
fwrite(prof, file.path(OUT,"clusterD_profiles.csv"))

# ---- 5. 驱动分析: 分类(A/D) ~ 自变量 ----------------------------------------
rf_importance <- function(yname, VS, tag){
  df <- dat[, c(yname, VS), with=FALSE]; setnames(df, yname, "y")
  df <- df[complete.cases(df)]; df[, y:=factor(y)]
  rf <- ranger(y~., df, num.trees=1500, importance="permutation", probability=TRUE, seed=42)
  imp <- data.table(variable=names(rf$variable.importance), importance=as.numeric(rf$variable.importance))[order(-importance)]
  fwrite(imp, file.path(OUT, sprintf("drivers_%s_rf_importance.csv",tag)))
  cat(sprintf("\n[%s] RF n=%d OOB-Brier=%.3f  top5: %s\n", tag, nrow(df), rf$prediction.error,
              paste(imp$variable[1:min(5,nrow(imp))],collapse=", ")))
  list(imp=imp, n=nrow(df), brier=rf$prediction.error)
}
# 主分析: 良好覆盖A集(n≈383)
rA  <- rf_importance("clsA", vars_a, "clsA_signtype_Aset")
rD  <- rf_importance("clsD", vars_a, "clsD_cluster_Aset")
# 稳健: 完整20变量(n≈146)
rAf <- rf_importance("clsA", vars,   "clsA_signtype_full")
rDf <- rf_importance("clsD", vars,   "clsD_cluster_full")

# ---- 6. 驱动分析: 连续因变量 ~ 自变量 --------------------------------------
cont_model <- function(yname, VS, BL, tag){
  df <- dat[, c(yname, VS), with=FALSE]; setnames(df, yname, "y")
  df <- df[complete.cases(df)]
  f  <- as.formula(paste("y ~", paste(VS, collapse="+")))
  m  <- lm(f, df); R2 <- summary(m)$r.squared
  z  <- copy(df); for(v in VS) set(z,j=v,value=as.numeric(scale(z[[v]]))); z[, y:=as.numeric(scale(y))]
  mz <- lm(f, z); co <- coef(summary(mz))
  std<- data.table(variable=rownames(co)[-1], beta=co[-1,1], se=co[-1,2], p=co[-1,4])[order(-abs(beta))]
  vp <- rbindlist(lapply(names(BL), function(b){
    rest<-setdiff(VS,BL[[b]]); if(!length(rest)) return(data.table(blk=b,uniq=R2))
    r2wo<-summary(lm(as.formula(paste("y~",paste(rest,collapse="+"))),df))$r.squared
    data.table(blk=b, uniq=R2-r2wo)
  }))[order(-uniq)]
  rf <- ranger(f, df, num.trees=1500, importance="permutation", seed=42)
  fwrite(std, file.path(OUT, sprintf("cont_%s_stdcoef.csv",tag)))
  fwrite(vp,  file.path(OUT, sprintf("cont_%s_varpart.csv",tag)))
  cat(sprintf("\n[连续:%s] n=%d lm-R2=%.3f rf-R2=%.3f\n top-beta: %s\n block-uniq: %s\n",
      tag,nrow(df),R2,rf$r.squared,
      paste(sprintf("%s=%.2f%s",std$variable[1:min(5,nrow(std))],std$beta[1:min(5,nrow(std))],
            ifelse(std$p[1:min(5,nrow(std))]<.05,"*","")),collapse=", "),
      paste(sprintf("%s=%.3f",vp$blk,vp$uniq),collapse=", ")))
  list(std=std, vp=vp, R2=R2, rfR2=rf$r.squared, tag=tag, n=nrow(df))
}
# 主分析: A集(n≈383)
cRho <- cont_model("rho_peak", vars_a, PRED_A, "rho_peak_Aset")
cLag <- cont_model("peak_tp",  vars_a, PRED_A, "peak_tp_Aset")
cCoef<- cont_model("mc_peak",  vars_a, PRED_A, "mc_peak_Aset")
# 稳健: 完整集(n≈146)
cRhoF<- cont_model("rho_peak", vars, PRED_FULL, "rho_peak_full")
cLagF<- cont_model("peak_tp",  vars, PRED_FULL, "peak_tp_full")
cCoefF<-cont_model("mc_peak",  vars, PRED_FULL, "mc_peak_full")

saveRDS(list(rA=rA,rD=rD,rAf=rAf,rDf=rDf,cRho=cRho,cLag=cLag,cCoef=cCoef,
             cRhoF=cRhoF,cLagF=cLagF,cCoefF=cCoefF,prof=prof), file.path(OUT,"models.rds"))
cat("\n=== 建模完成, 开始出图 ===\n")

# =============================================================================
# 图
# =============================================================================
theme_cn <- function() theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
lev6<-c("always_promote","promote_inhibit","inhibit_promote","always_inhibit","mixed")
lab6<-c("Always promote","Promote->Inhibit","Inhibit->Promote","Always inhibit","Mixed")
col6<-setNames(c("#C1121F","#F4A261","#74C6E8","#1A6FBF","#9E9E9E"),lab6)

# --- 图1: 分类A 中国地图(气候区底图) ---
bb<-c(72,136,17,54); k2g<-function(x)dplyr::case_when(x>=1&x<=4~"A",x>=5&x<=9~"B",x>=10&x<=17~"C",x>=18&x<=28~"D",TRUE~NA_character_)
world<-tryCatch(ne_countries(scale="medium",returnclass="sf"),error=function(e)NULL)
chn<-tryCatch(ne_countries(country=c("China","Taiwan"),scale="medium",returnclass="sf")|>st_union(),error=function(e)NULL)
prov<-tryCatch(ne_states(country="China",returnclass="sf"),error=function(e)NULL)
mapdat<-copy(dat); mapdat[, slab:=factor(lab6[match(clsA,lev6)],levels=lab6)]
pm<-ggplot()
if(!is.null(world)) pm<-pm+geom_sf(data=world,fill="grey90",color="grey78",linewidth=.15)
if(!is.null(prov))  pm<-pm+geom_sf(data=prov,fill=NA,color="grey60",linewidth=.2)
pm<-pm+geom_point(data=mapdat,aes(longitude,latitude,color=slab),size=1.9,alpha=.85)+
  scale_color_manual(values=col6,name="Sign-type (9-tp)")+
  coord_sf(xlim=bb[1:2],ylim=bb[3:4],expand=FALSE)+
  labs(title="Loose-criteria stations (n=383): VPD->SIF sign-type",
       subtitle="p_surr<0.1 & drho>0; Ushio-style normalization; S-map sign sequence over tp=0..8",
       x="Longitude",y="Latitude")+
  theme_minimal(base_size=12)+theme(panel.background=element_rect(fill="#EAF3FB",color=NA),
    plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,color="grey40"))
ggsave(file.path(OUT,"map_clsA_signtype.png"),pm,width=13,height=8.5,dpi=300)

# --- 图2: 分类D 聚类地图 ---
pmd<-ggplot()
if(!is.null(world)) pmd<-pmd+geom_sf(data=world,fill="grey90",color="grey78",linewidth=.15)
if(!is.null(prov))  pmd<-pmd+geom_sf(data=prov,fill=NA,color="grey60",linewidth=.2)
pmd<-pmd+geom_point(data=dat,aes(longitude,latitude,color=clsD),size=1.9,alpha=.85)+
  scale_color_brewer(palette="Set1",name="Data-driven cluster")+
  coord_sf(xlim=bb[1:2],ylim=bb[3:4],expand=FALSE)+
  labs(title="Data-driven CCM clusters (k-means, k=4) over China",
       subtitle="Features: peak_tp, rho_peak, mc_peak, n_sig, frac_pos (standardized)",
       x="Longitude",y="Latitude")+
  theme_minimal(base_size=12)+theme(panel.background=element_rect(fill="#EAF3FB",color=NA),
    plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,color="grey40"))
ggsave(file.path(OUT,"map_clsD_cluster.png"),pmd,width=13,height=8.5,dpi=300)

# --- 图3: 聚类特征剖面热图 ---
profz <- copy(prof); for(c in c("peak_tp","rho_peak","mc_peak","n_sig","frac_pos")) profz[[c]]<-as.numeric(scale(prof[[c]]))
plong<-melt(profz[, .(clsD,peak_tp,rho_peak,mc_peak,n_sig,frac_pos)],id.vars="clsD")
ph<-ggplot(plong,aes(variable,clsD,fill=value))+geom_tile(color="white")+
  geom_text(aes(label=sprintf("%.2f",value)),size=3.4)+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="z-score")+
  labs(title="Cluster profiles (z-scored CCM features)",
       subtitle=sprintf("k-means k=4; cell = cluster mean; cluster n: %s",paste(sprintf("C%s=%d",prof$clsD,prof$n),collapse=", ")),
       x="CCM feature",y="Cluster")+theme_cn()
ggsave(file.path(OUT,"clusterD_profile_heat.png"),ph,width=8.5,height=4.5,dpi=300)

# --- 图4: RF重要性(分类A + 分类D) ---
impdf<-rbind(rA$imp[,.(variable,importance,which="Sign-type (A)")],
             rD$imp[,.(variable,importance,which="Cluster (D)")])
impdf[, variable:=factor(variable, levels=rA$imp$variable)]
pimp<-ggplot(impdf,aes(importance,variable,fill=which))+geom_col(position="dodge")+
  scale_fill_manual(values=c("Sign-type (A)"="#4575B4","Cluster (D)"="#D73027"),name="Classification")+
  labs(title="What predicts the station class? (RF permutation importance)",
       subtitle=sprintf("Well-covered A-set (9 vars); n=%d; clsA OOB-Brier=%.3f, clsD OOB-Brier=%.3f",
                        rA$n,rA$brier,rD$brier),
       x="Permutation importance",y=NULL)+theme_cn()
ggsave(file.path(OUT,"drivers_class_rf_importance.png"),pimp,width=9,height=7,dpi=300)

# --- 图5: 连续因变量 标准化系数森林图 ---
cont_forest<-rbind(cRho$std[,.(variable,beta,se,p,y="Causal strength (rho_peak)")],
                   cLag$std[,.(variable,beta,se,p,y="Response lag (peak_tp)")],
                   cCoef$std[,.(variable,beta,se,p,y="Signed effect (mc_peak)")])
cont_forest[, variable:=factor(variable, levels=rev(vars_a))]
pcf<-ggplot(cont_forest,aes(beta,variable))+
  geom_vline(xintercept=0,linetype="dashed",color="grey55")+
  geom_errorbarh(aes(xmin=beta-1.96*se,xmax=beta+1.96*se),height=.3)+
  geom_point(aes(shape=p<.05),size=2)+scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide="none")+
  facet_wrap(~y,scales="free_x")+
  labs(title=sprintf("Continuous CCM outcomes ~ well-covered drivers (standardized OLS, n=%d)",cRho$n),
       subtitle="Filled = p<0.05; x = standardized beta (95% CI)",x="Standardized coefficient",y=NULL)+
  theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"drivers_continuous_forest.png"),pcf,width=13,height=6.5,dpi=300)

# --- 图6: 块方差分解 ---
vpall<-rbind(cRho$vp[,.(blk,uniq,y="rho_peak")],cLag$vp[,.(blk,uniq,y="peak_tp")],cCoef$vp[,.(blk,uniq,y="mc_peak")])
pvp<-ggplot(vpall,aes(reorder(blk,uniq),uniq,fill=y))+geom_col(position="dodge")+coord_flip()+
  scale_fill_brewer(palette="Dark2",name="Outcome")+
  labs(title="Unique variance explained by each driver block",
       subtitle="Hierarchical partitioning: R2(full) - R2(without block)",x=NULL,y="Unique R2")+theme_cn()
ggsave(file.path(OUT,"drivers_varpart_blocks.png"),pvp,width=9,height=5,dpi=300)

# --- 图7: 关键自变量在各符号型(clsA)间的梯度 boxplot ---
keyv<-c(elev="Elevation (m)",rsds_mean="Radiation mean",precip_mean="Precip (mm)",
        crop="Cropland frac",forest="Forest frac",imperv="Impervious frac")
bx<-melt(dat[!is.na(clsA), c("clsA",names(keyv)),with=FALSE],id.vars="clsA")
bx[, clsA:=factor(lab6[match(clsA,lev6)],levels=lab6)]
bx[, variable:=factor(keyv[as.character(variable)],levels=keyv)]
pbx<-ggplot(bx,aes(clsA,value,fill=clsA))+geom_boxplot(outlier.size=.5,alpha=.85)+
  facet_wrap(~variable,scales="free_y",ncol=3)+scale_fill_manual(values=col6,guide="none")+
  labs(title="Key drivers across VPD->SIF sign-types (n=383)",
       subtitle="Which environments host which causal pattern",x=NULL,y=NULL)+
  theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),
    axis.text.x=element_text(angle=30,hjust=1,size=8))
ggsave(file.path(OUT,"drivers_boxplot_by_signtype.png"),pbx,width=12,height=7,dpi=300)

# --- 图8: 数据驱动聚类的 Koppen 气候区构成 ---
zb<-dat[!is.na(koppen_group)&koppen_group%in%c("A","B","C","D"), .N, by=.(clsD,koppen_group)]
pz<-ggplot(zb,aes(clsD,N,fill=koppen_group))+geom_col(position="fill")+
  scale_fill_manual(values=c(A="#5A8F76",B="#EAD5A0",C="#A3B86C",D="#C2DFCD"),name="Koppen")+
  scale_y_continuous(labels=scales::percent)+
  labs(title="Climate-zone composition of data-driven clusters",
       subtitle="Fraction of each cluster in Koppen A/B/C/D",x="Cluster",y="Share")+theme_cn()
ggsave(file.path(OUT,"clusterD_by_koppen.png"),pz,width=8,height=5,dpi=300)

cat("\n完成. 输出目录:", OUT, "\n"); cat(list.files(OUT,pattern="png$"),sep="\n")
