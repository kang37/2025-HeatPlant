#!/usr/bin/env Rscript
# 逐站 stype(严格5类) 在 1/2/3 km 缓冲区间的桑基/冲积图。
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggalluvial)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"
latest <- function(dir,pat){f<-list.files(dir,pattern=pat,full.names=TRUE);f<-f[!grepl("_raw",f)];tail(sort(f),1)}
OUT <- file.path(PROJ,"data_proc/output_hcsif_scale_compare"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

srcs <- list(
  `1km`=list(r=1000,ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf1000_v3"),"ccm_hcsif_buf1000_vpd_.*rds$")),
  `2km`=list(r=2000,ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf2000"),"ccm_hcsif_buf2000_vpd_.*rds$")),
  `3km`=list(r=3000,ccm=latest(file.path(PROJ,"data_proc/ccm_hcsif_buf3000"),"ccm_hcsif_buf3000_vpd_.*rds$")))

classify5 <- function(coefs){s<-sign(coefs)
  if(any(is.na(s))||any(s==0))return("other")
  nch<-sum(diff(s)!=0)
  if(nch==0)return(if(s[1]>0)"always_promote" else "always_inhibit")
  if(nch==1)return(if(s[1]<0)"inhibit_promote" else "promote_inhibit"); "other"}

stype_of <- function(r, ccm_path) {
  SIF_Y <- sprintf("SIF_buf%d_dt", r)
  ccm <- readRDS(ccm_path) %>% as.data.frame() %>% filter(y_var==SIF_Y, x_var=="vpd_mean_dt")
  N_TP <- length(unique(ccm$tp))
  ccm %>% arrange(meteo_stat,tp) %>% group_by(meteo_stat) %>%
    summarise(cs=list(mean_coef),.groups="drop") %>%
    filter(map_lgl(cs,~length(.x)==N_TP)) %>%
    mutate(stype=map_chr(cs,classify5)) %>% select(meteo_stat,stype)
}

tabs <- imap(srcs, ~ stype_of(.x$r, .x$ccm) %>% rename(!!.y := stype))
d <- reduce(tabs, inner_join, by="meteo_stat")   # 三尺度都可分类的站点
cat("三尺度均可分类的站点数:", nrow(d), "\n")

lev <- c("always_promote","promote_inhibit","inhibit_promote","always_inhibit","other")
lab <- c("Always promote","Promote->Inhibit","Inhibit->Promote","Always inhibit","Other (mixed)")
col <- setNames(c("#C1121F","#F4A261","#74C6E8","#1A6FBF","#9E9E9E"), lab)

long <- d %>%
  mutate(id=row_number()) %>%
  pivot_longer(c(`1km`,`2km`,`3km`), names_to="scale", values_to="stype") %>%
  mutate(scale=factor(scale,levels=c("1km","2km","3km")),
         stype=factor(lab[match(stype,lev)], levels=lab))

# 稳定性统计
same_all <- d %>% filter(`1km`==`2km` & `2km`==`3km`) %>% nrow()
cat(sprintf("三尺度类别完全一致: %d 站 (%.1f%%)\n", same_all, 100*same_all/nrow(d)))
write_csv <- data.table::fwrite
write_csv(as.data.table(d), file.path(OUT,"stype_by_scale_wide.csv"))

p <- ggplot(long, aes(x=scale, stratum=stype, alluvium=id, fill=stype)) +
  geom_flow(aes(fill=`stype`), stat="alluvium", alpha=0.55, lode.guidance="frontback", color=NA,
            decreasing=FALSE) +
  geom_stratum(width=0.42, color="white", decreasing=FALSE) +
  geom_text(stat="stratum", aes(label=after_stat(count)), size=3, decreasing=FALSE) +
  scale_fill_manual(values=col, name="Pattern (strict 5-class)") +
  labs(title="Station causal-type flow across buffer scales (buf1000_v3 -> 2km -> 3km)",
       subtitle=sprintf("Strict 5-class from VPD->SIF CCM sign sequence; n=%d stations classifiable at all three scales; %d (%.0f%%) unchanged",
                        nrow(d), same_all, 100*same_all/nrow(d)),
       x="Buffer radius", y="Number of stations") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"), panel.grid.major.x=element_blank())

ggsave(file.path(OUT,"sankey_stype_by_buffer.png"), p, width=11, height=8, dpi=300)
cat("-> ", file.path(OUT,"sankey_stype_by_buffer.png"), "\n")

# 1km->3km 转移矩阵(便于量化)
tm <- d %>% count(`1km`,`3km`) %>%
  mutate(`1km`=factor(lab[match(`1km`,lev)],levels=lab),
         `3km`=factor(lab[match(`3km`,lev)],levels=lab))
write_csv(as.data.table(tm), file.path(OUT,"transition_1km_to_3km.csv"))
cat("\n=== 1km -> 3km 转移(行=1km, 列=3km) ===\n")
print(as.data.frame(tidyr::pivot_wider(tm, names_from=`3km`, values_from=n, values_fill=0)))
