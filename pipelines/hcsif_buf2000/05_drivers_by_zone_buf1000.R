#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 05_drivers_by_zone_buf1000.R  （仅 1000m 方案，分 Köppen 气候区）
#   分区回答：模式构成、以及2分类下转变速度的影响因素（城市聚类稳健SE）
#   A 区样本过少(5站)自动剔除；每个区×类事件<12 自动跳过。
#   全英文标签。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(sandwich); library(lmtest)
})
PROJ    <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
CCM_DIR <- file.path(PROJ, "data_proc/ccm_hcsif_buf1000")
COV_RDS <- file.path(PROJ, "data_proc/output_10y_built_up_05_01/station_covariates.rds")
OUT     <- file.path(PROJ, "data_proc/output_hcsif_buf1000")
SIF_Y   <- "SIF_buf1000_dt"; VPD_X <- "vpd_mean_dt"
theme_cn <- function() theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
winsorize <- function(x, p = 0.01) { q <- quantile(x, c(p,1-p), na.rm=TRUE); pmin(pmax(x,q[1]),q[2]) }
MIN_EVENTS <- 12L
zone_lab <- c(B="B (Arid)", C="C (Temperate)", D="D (Continental)")

CCM_RDS <- { fs <- list.files(CCM_DIR, pattern="^ccm_hcsif_buf1000_vpd_.*[0-9]\\.rds$", full.names=TRUE)
  tail(sort(fs[!grepl("_raw\\.rds$", fs)]), 1) }
ccm <- readRDS(CCM_RDS) %>% as.data.frame() %>% filter(y_var==SIF_Y, x_var==VPD_X)
N_TP <- length(unique(ccm$tp)); MAXTP <- max(ccm$tp)

find_transition <- function(coefs, to_neg=TRUE){cond<-if(to_neg)function(x)x<0 else function(x)x>0
  if(length(coefs)<2)return(NA_integer_);for(i in 2:length(coefs))if(!is.na(coefs[i])&&cond(coefs[i]))
    if(mean(sapply(coefs[i:length(coefs)],cond),na.rm=TRUE)>=0.5)return(as.integer(i-1L));NA_integer_}

pat <- ccm %>% arrange(meteo_stat,tp) %>% group_by(meteo_stat) %>%
  summarise(coef_seq=list(mean_coef),.groups="drop") %>%
  filter(map_lgl(coef_seq,~length(.x)==N_TP)) %>%
  mutate(tp0=map_lgl(coef_seq,~.x[[1]]>0),
    HTW=map_int(coef_seq,~find_transition(.x,TRUE)),
    ITW=map_int(coef_seq,~find_transition(.x,FALSE)),
    stype=case_when(tp0&is.na(HTW)~"always_promote",!tp0&is.na(ITW)~"always_inhibit",
      tp0&!is.na(HTW)~"promote_inhibit",!tp0&!is.na(ITW)~"inhibit_promote",TRUE~"always_inhibit"),
    class2=ifelse(stype%in%c("always_inhibit","inhibit_promote"),"Inhibit-first","Promote-first"),
    flipped=stype%in%c("inhibit_promote","promote_inhibit"),
    event_tp=case_when(stype=="inhibit_promote"~ITW,stype=="promote_inhibit"~HTW,TRUE~NA_integer_))

cov <- readRDS(COV_RDS) %>% rename(meteo_stat=meteo_stat_id) %>% mutate(meteo_stat=as.character(meteo_stat))
dat <- pat %>% mutate(meteo_stat=as.character(meteo_stat)) %>% left_join(cov,by="meteo_stat") %>%
  filter(!is.na(koppen_group), koppen_group %in% c("B","C","D")) %>%
  mutate(inv_w=winsorize(pa_built_10y),cgi_w=winsorize(cgi_score),precip_w=winsorize(precip_mean),
    pop_w=winsorize(pop_10y),road_w=winsorize(road_density),bvol_w=winsorize(building_vol_density),
    lon=longitude,lat=latitude,city=city_name)
pred_vars  <- c("inv_w","cgi_w","precip_w","pop_w","road_w","bvol_w","lon","lat")
var_labels <- c(inv_w="Investment",cgi_w="CGI",precip_w="Precipitation",pop_w="Population",
                road_w="Road density",bvol_w="Building vol.density",lon="Longitude",lat="Latitude")

# ===========================================================================
# 1. 模式构成 x 气候区
# ===========================================================================
comp <- dat %>% count(koppen_group, stype) %>% group_by(koppen_group) %>%
  mutate(frac=n/sum(n)) %>% ungroup()
write_csv(comp, file.path(OUT,"zone_pattern_composition.csv"))
st_lv <- c("always_promote","promote_inhibit","inhibit_promote","always_inhibit")
st_ll <- c("Always promote","Promote->Inhibit","Inhibit->Promote","Always inhibit")
st_col<- c("Always promote"="#C1121F","Promote->Inhibit"="#F4A261",
           "Inhibit->Promote"="#74C6E8","Always inhibit"="#1A6FBF")
p_comp <- comp %>% mutate(z=zone_lab[koppen_group],
    sl=factor(stype,levels=st_lv,labels=st_ll)) %>%
  ggplot(aes(z,frac,fill=sl))+geom_col(width=.65)+
  geom_text(aes(label=ifelse(frac>=.06,n,"")),position=position_stack(vjust=.5),size=3,color="white")+
  scale_fill_manual(values=st_col,name="Pattern")+
  scale_y_continuous(labels=scales::percent)+
  labs(title="Pattern composition by Köppen zone (buf1000)",
       subtitle="Bar label = station count; A zone (n=5) excluded",x=NULL,y="Share")+theme_cn()
ggsave(file.path(OUT,"zone_pattern_composition.png"),p_comp,width=9,height=5,dpi=300)
cat("-> zone_pattern_composition.png\n"); print(as.data.frame(comp))

# ===========================================================================
# 2. 分区 x 2分类：转变速度的离散时间风险模型（城市聚类）
# ===========================================================================
surv <- dat %>% mutate(time=ifelse(flipped,event_tp,MAXTP),event=as.integer(flipped))
pp <- surv %>% filter(!is.na(time),time>=1) %>% rowwise() %>%
  do({r<-.;tibble(meteo_stat=r$meteo_stat,koppen_group=r$koppen_group,class2=r$class2,city=r$city,
     tp=1:r$time,event=as.integer((1:r$time)==r$time & r$event==1)) %>% bind_cols(r[pred_vars])}) %>%
  ungroup()

haz_zone <- function(zone, cl) {
  d <- pp %>% filter(koppen_group==zone, class2==cl) %>%
    select(event,tp,city,all_of(pred_vars)) %>% drop_na()
  if (sum(d$event) < MIN_EVENTS) {
    cat(sprintf("[跳过] %s x %s: 事件=%d < %d\n", zone, cl, sum(d$event), MIN_EVENTS)); return(NULL) }
  form <- as.formula(paste("event ~ tp +", paste(pred_vars, collapse=" + ")))  # 分区样本小,用线性tp
  fit <- tryCatch(glm(form,data=d,family=binomial), error=function(e) NULL)
  if (is.null(fit)) { cat(sprintf("[跳过] %s x %s 拟合失败\n",zone,cl)); return(NULL) }
  ct <- tryCatch(lmtest::coeftest(fit,vcov=sandwich::vcovCL,cluster=~city), error=function(e) NULL)
  if (is.null(ct)) return(NULL)
  tibble(zone=zone,class=cl,variable=rownames(ct),coef=ct[,1],se=ct[,2],p_val=ct[,4],
         n_event=sum(d$event),n_row=nrow(d)) %>% filter(variable %in% pred_vars)
}
zones <- c("B","C","D")
res <- bind_rows(lapply(zones, function(z) bind_rows(haz_zone(z,"Inhibit-first"),
                                                     haz_zone(z,"Promote-first"))))
write_csv(res, file.path(OUT,"zone_hazard_coefs.csv"))
cat("\n各区×2分类 转变速度系数(城市聚类) 已存 zone_hazard_coefs.csv\n")
cat("\n投资(Investment)在各区的翻转风险系数:\n")
print(as.data.frame(res %>% filter(variable=="inv_w") %>%
        transmute(zone,class,coef=round(coef,4),p_val=round(p_val,4),n_event)))

# forest：分区 x 类
mkforest <- function(cl, tag) {
  sub <- res %>% filter(class==cl)
  if (!nrow(sub)) return(invisible())
  p <- sub %>% mutate(z=zone_lab[zone],
      vl=factor(var_labels[variable],levels=rev(var_labels[pred_vars])),sig=p_val<0.05) %>%
    ggplot(aes(coef,vl))+geom_vline(xintercept=0,linetype="dashed",color="grey50")+
    geom_errorbarh(aes(xmin=coef-1.96*se,xmax=coef+1.96*se),height=.25)+
    geom_point(aes(shape=sig),size=2.6)+scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide="none")+
    facet_wrap(~z)+
    labs(title=sprintf("Transition-speed drivers by zone — %s (buf1000)",cl),
         subtitle="Discrete-time hazard + city-clustered SE; coef>0 -> earlier transition; filled = p<0.05",
         x="Log-hazard coefficient (95% CI)",y=NULL)+theme_cn()
  ggsave(file.path(OUT,sprintf("zone_hazard_%s.png",tag)),p,width=13,height=4.5,dpi=300)
  cat(sprintf("-> zone_hazard_%s.png\n",tag))
}
mkforest("Inhibit-first","inhibit_first")
mkforest("Promote-first","promote_first")

# 投资效应跨区汇总图
inv <- res %>% filter(variable=="inv_w")
if (nrow(inv)) {
  p_inv <- inv %>% mutate(z=zone_lab[zone],sig=p_val<0.05) %>%
    ggplot(aes(coef,z,color=class))+geom_vline(xintercept=0,linetype="dashed",color="grey50")+
    geom_errorbarh(aes(xmin=coef-1.96*se,xmax=coef+1.96*se),height=.2,position=position_dodge(.5))+
    geom_point(aes(shape=sig),size=3,position=position_dodge(.5))+
    scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),guide="none")+
    scale_color_manual(values=c("Inhibit-first"="#3182BD","Promote-first"="#E6550D"),name="2-class")+
    labs(title="Investment effect on transition speed, by zone (buf1000)",
         subtitle="coef>0 -> higher investment brings EARLIER transition; filled = p<0.05 (city-clustered)",
         x="Investment log-hazard coefficient (95% CI)",y=NULL)+theme_cn()
  ggsave(file.path(OUT,"zone_investment_effect.png"),p_inv,width=9,height=4.5,dpi=300)
  cat("-> zone_investment_effect.png\n")
}
cat("\n完成\n")
