#!/usr/bin/env Rscript
# 为 1/2/3 km 各出一张 站点×tp 符号热力图(严格5类排序)，同 output_hcsif_buf1000 样式。
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(ggplot2)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
latest <- function(dir, pat){f<-list.files(dir,pattern=pat,full.names=TRUE);f<-f[!grepl("_raw",f)];tail(sort(f),1)}

BUFS <- list(
  list(r=1000, ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_v3"),"ccm_hcsif_buf1000_vpd_.*rds$")),
  list(r=2000, ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf2000"),"ccm_hcsif_buf2000_vpd_.*rds$")),
  list(r=3000, ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf3000"),"ccm_hcsif_buf3000_vpd_.*rds$")))

classify5 <- function(coefs){s<-sign(coefs)
  if(any(is.na(s))||any(s==0))return("other")
  nch<-sum(diff(s)!=0)
  if(nch==0)return(if(s[1]>0)"always_promote" else "always_inhibit")
  if(nch==1)return(if(s[1]<0)"inhibit_promote" else "promote_inhibit"); "other"}
trans_pos <- function(coefs){s<-sign(coefs)
  if(any(is.na(s))||any(s==0))return(NA_integer_)
  ch<-which(diff(s)!=0); if(length(ch)!=1)return(NA_integer_); as.integer(ch+1L)}

make_heat <- function(r, ccm_path) {
  SIF_Y <- sprintf("SIF_buf%d_dt", r); VPD_X <- "vpd_mean_dt"
  OUT <- file.path(PROJ, sprintf("data_proc/output_hcsif_buf%d", r))
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  ccm <- readRDS(ccm_path) %>% as.data.frame() %>% filter(y_var==SIF_Y, x_var==VPD_X)
  N_TP <- length(unique(ccm$tp))
  pat <- ccm %>% arrange(meteo_stat,tp) %>% group_by(meteo_stat) %>%
    summarise(coef_seq=list(mean_coef),.groups="drop") %>%
    filter(map_lgl(coef_seq,~length(.x)==N_TP)) %>%
    mutate(stype=map_chr(coef_seq,classify5), tpos=map_int(coef_seq,trans_pos))
  grp_rank <- c(always_promote=1L,promote_inhibit=2L,inhibit_promote=3L,always_inhibit=4L,other=5L)
  order_tbl <- pat %>% mutate(
    n_neg=map_int(coef_seq,~sum(.x<0,na.rm=TRUE)), g=grp_rank[stype],
    w=case_when(stype=="always_promote"~0,
                stype=="promote_inhibit"~ -as.numeric(tpos),
                stype=="inhibit_promote"~ as.numeric(tpos),
                stype=="always_inhibit"~0,
                stype=="other"~ as.numeric(n_neg))) %>%
    arrange(g,w) %>% mutate(station_order=row_number())
  heat_long <- ccm %>% select(meteo_stat,tp,mean_coef) %>%
    inner_join(order_tbl %>% select(meteo_stat,station_order,stype),by="meteo_stat") %>%
    mutate(sign=case_when(mean_coef>0~"Promote",mean_coef<0~"Inhibit",TRUE~NA_character_))
  gb <- order_tbl %>% group_by(stype) %>% summarise(ymax=max(station_order),.groups="drop") %>%
    arrange(ymax) %>% pull(ymax); gb <- gb[-length(gb)]+0.5
  glab <- order_tbl %>% group_by(stype) %>% summarise(y=mean(station_order),n=n(),.groups="drop") %>%
    mutate(lab=c(always_promote="Always promote",promote_inhibit="Promote->Inhibit",
                 inhibit_promote="Inhibit->Promote",always_inhibit="Always inhibit",
                 other="Other (mixed)")[stype])
  p <- ggplot(heat_long,aes(factor(tp),station_order,fill=sign))+geom_tile()+
    geom_hline(yintercept=gb,color="white",linewidth=.6)+
    scale_fill_manual(values=c("Promote"="#C1121F","Inhibit"="#1A6FBF"),na.value="grey85",name="VPD->SIF effect")+
    scale_y_continuous(expand=c(0,0),breaks=glab$y,labels=sprintf("%s\n(n=%d)",glab$lab,glab$n))+
    labs(title=sprintf("Sign of VPD->SIF causal effect across CCM lag — buffer %d m",r),
         subtitle=sprintf("SIF=HCSIF buf%d; each row=1 station (n=%d), each col=1 lag step (8 days);\nbottom->top: always promote -> promote-then-inhibit -> inhibit-then-promote -> always inhibit -> other",r,nrow(order_tbl)),
         x="CCM lag tp (each step = 8 days)",y=NULL)+
    theme_minimal(base_size=12)+
    theme(panel.grid=element_blank(),axis.text.y=element_text(size=10,lineheight=.9),
          plot.title=element_text(face="bold"))
  ggsave(file.path(OUT,"heatmap_station_tp_sign.png"),p,width=9,height=11,dpi=300)
  cat(sprintf("-> buf%d: %s | 5类分布: %s\n", r, file.path(OUT,"heatmap_station_tp_sign.png"),
              paste(sprintf("%s=%d",names(table(pat$stype)),as.integer(table(pat$stype))),collapse=", ")))
}
for (b in BUFS) make_heat(b$r, b$ccm)
