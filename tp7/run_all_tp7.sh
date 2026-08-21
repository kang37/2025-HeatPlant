#!/bin/bash
# 依次跑 tp=7 截断版的全套分析, 每个脚本的控制台输出单独存一份日志
cd /Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant
mkdir -p tp7/logs
for s in 16_hcsif_figures_tp7 26_fig_rho_tp7 29_robust_se_tp7 32_final_conley_tp7 \
         34_zone_forest_tp7 35_importance_tp7 40_survival_fourtype_tp7; do
  echo "===== $s ====="
  Rscript tp7/$s.R > tp7/logs/$s.log 2>&1
  echo "exit=$? -> tp7/logs/$s.log"
  tail -3 tp7/logs/$s.log
done
echo "ALL DONE"
