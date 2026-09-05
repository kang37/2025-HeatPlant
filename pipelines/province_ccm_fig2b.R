#!/usr/bin/env Rscript
# =============================================================================
# province_ccm_fig2b.R  —— 仿 s41467-024-48199-z Fig.2b
#   各省 站点 Δρ 的小提琴分布 + 用各省站点 p 值 Fisher 合并的省级/全国 meta 显著性。
#   (原文一个驱动=O3/AH/T 各一把小提琴; 本研究单驱动 VPD, 用"省"作分组。)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ<-"/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT<-file.path(PROJ,"data_proc/output_loose_classify")
pk<-fread(file.path(OUT,"province_ccm_stations.csv"))   # 由 province_ccm_fig2.R 产出

PFLOOR<-1/200                                            # N_SURR=199 → p 下限
pk[, p:=pmax(p_surr, PFLOOR)]
# Fisher 合并: chi2=-2*sum(ln p), df=2k
fisher_p<-function(p){k<-length(p); if(!k) return(NA_real_)
  pchisq(-2*sum(log(p)), df=2*k, lower.tail=FALSE)}

prov<-pk[, .(n=.N, nsig=sum(sig), med_drho=median(drho),
             p_fisher=fisher_p(p)), by=name_en]
prov[, p_meta_prov:=p.adjust(p_fisher, "BH")]           # 省间多重比较校正
prov[, sig_prov:=p_meta_prov < 0.05]
setorder(prov, -med_drho)
pk[, prov_f:=factor(name_en, levels=prov$name_en)]

# 全国 meta: 对全部站点 p 做 Fisher (原文口径: 合并出一个全国 p_meta)
p_meta_nat<-fisher_p(pk$p)
cat(sprintf("全国 meta 显著性 (Fisher over %d stations): p_meta = %.3e\n", nrow(pk), p_meta_nat))
cat(sprintf("省级 Fisher 合并 p<0.05 (BH校正) 的省: %d/%d\n", sum(prov$sig_prov), nrow(prov)))
print(prov[order(p_meta_prov)][1:12])
fwrite(prov, file.path(OUT,"province_meta_significance.csv"))

# 标签: 省名 + 显著星
prov[, lab:=sprintf("%s%s", name_en, ifelse(sig_prov," *",""))]
pk<-merge(pk, prov[,.(name_en,lab,sig_prov)], by="name_en")
pk[, lab:=factor(lab, levels=prov$lab)]

# ---- 图 2b: 省级 Δρ 小提琴 + 点 ----
p2b<-ggplot(pk, aes(lab, drho))+
  geom_hline(yintercept=0, color="grey55", linetype="dashed")+
  geom_violin(aes(fill=sig_prov), color="grey60", scale="width", width=.85, alpha=.5, linewidth=.3)+
  geom_jitter(aes(color=sig), width=.15, size=.7, alpha=.8)+
  scale_fill_manual(values=c(`TRUE`="#C1121F",`FALSE`="grey80"),
                    labels=c(`TRUE`="Province meta-sig (Fisher, BH p<0.05)","FALSE"="Not"),name="Province")+
  scale_color_manual(values=c(`TRUE`="#C1121F",`FALSE`="grey55"),
                     labels=c(`TRUE`="Station sig (p<0.1 & converge)","FALSE"="Not"),name="Station")+
  labs(title="(2b) Per-province distribution of convergence (delta-rho) and meta-significance",
       subtitle=sprintf("Violin = station-level drho within each province (ordered by median). Star = province Fisher-combined p<0.05 (BH).\nNationwide meta-significance (Fisher over %d stations): p_meta = %.2e   |   %d/%d provinces meta-significant",
                        nrow(pk), p_meta_nat, sum(prov$sig_prov), nrow(prov)),
       x=NULL, y="Convergence  delta-rho = rho - rho_min")+
  theme_bw(base_size=11)+
  theme(plot.title=element_text(face="bold"), legend.position="top",
        axis.text.x=element_text(angle=55,hjust=1,size=8),
        panel.grid.major.x=element_blank())
ggsave(file.path(OUT,"province_ccm_fig2b.png"), p2b, width=15, height=7.5, dpi=300)
cat("-> province_ccm_fig2b.png\n")
