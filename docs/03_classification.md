# 模块三：分类回归（当前状态：44-60 主线已删除，仅剩因果方向分类支线）

**2026-09-07 重大变更**：原 44-60(`44_class_data.R`~`60_s3_eval.R`)"18因素→抑制组/
促进组"driver 回归主线**已整体删除**（用户决定）。该主线用的 grp2 方向标签（单变量
S-map、固定 tp=0 符号，全部 904 站不筛因果显著性）已知有偏，且和"确认方向"支线的
134 站方向标签长期并存造成反复混淆（同一个"站点通过因果检验"的问题在对话里出现过
至少 4 次不同答案：107/117/134/345/408，参见 [00_cross_module_issues.md](00_cross_module_issues.md)）。
删除前的历史结果（AUC 0.79-0.88，8 个显著驱动变量：海拔/城镇化率/混交林/辐射均值/
水体/常绿阔叶林/落叶阔叶林/绿地投资）仍可在 git 历史里找回（`git log -- 44_class_data.R`
一类命令），但**不再是本项目的有效结论**，写论文/做后续分析不要引用。

## 现在唯一的因果确认标准

`pipelines/hcsif_buf1000/12_ccm_causal_confirmation.R`——项目里**唯一**的"CCM确认
VPD→SIF因果耦合"判据，取代之前所有其他版本（含旧134站判据、44-60主线宽松/严格判据）。
三条阈值/范围全部集中在脚本顶部，改标准不用碰逻辑：

1. tp∈[0,8] 范围内，存在某个 tp 同时满足 Δρ=rho-rho_min>0 且 p_surr<0.05
2. 在更大范围 tp∈[-8,8]（含负tp）里取全局 optimal_tp=argmax(rho)，要求 optimal_tp≥0
   （数据源必须是负tp扫描版本，否则这条会因为搜索范围本身不含负值而恒真，见
   [02_ccm.md](02_ccm.md) 2026-09-04 日志的教训）
3. 另外沿用旧脚本的地表覆盖断点站过滤（森林占比年际序列有结构性跳变的站，CCM单
   吸引子假设不成立）——这条不算严格意义上的第三条标准，是继承的数据质量把关，
   独立标注方便以后决定要不要保留。

**当前(2026-09-07)漏斗结果**（`data_proc/ccm_hcsif_buf1000_causal_confirmed/funnel.csv`）：

| 层级 | 剩余站数 |
|---|---|
| 0. 有CCM输出(904站) | 904 |
| 0b. 排除地表覆盖断点站 | 868 |
| 1. tp∈[0,8]内至少1个tp同时满足Δρ>0且p_surr<0.05 | 178 |
| 2. + 全局optimal_tp(tp∈[-8,8])≥0 | **117** |

最终站点清单：`data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv`
（含每站 optimal_tp/optimal_rho/optimal_drho/optimal_p_surr/nsig_0to8/nsig_all/
tp_sig_list）。**以后任何需要"CCM确认因果耦合站点清单"的分析，一律读这个文件，
不要自己重新拼一套判据**——这是本次改动的核心目的。

## 134站方向x tp分类（唯一保留的方向分类支线，脚本沿用历史命名）

脚本 `pipelines/hcsif_buf1000/06_smap_bivar_134.R`（二变量S-map方向，2026-09-07
已改为直接读上面 12 号脚本的统一站点清单，不再自己重建）、
`07_smap_robustness_theta_multivar.R`（30站稳健性抽查，基于旧134站清单，未重跑）、
`08_classify_134.R`（方向x tp分箱分类）、`09_fig_classify_134.R`（图，未重跑，
数据已过期）。产出在 `data_proc/smap_bivar_134/`。脚本/文件名仍叫"134"是历史包袱，
实际站数以 `classify_134.csv` 行数为准，不要按文件名猜数字。

**方法**：方向用二变量 S-map（状态空间`[SIF(t),VPD(t-optimal_tp)]`，E=2，theta=2）
逐时间点系数的中位数 + 主导符号占比≥75%阈值判定(Promote/Inhibit/Ambiguous)。

**2026-09-07 用新统一判据重跑的结果**：117 站，Promote 41 / Ambiguous 60 /
Inhibit 16（旧134站判据下是 Promote 29/Inhibit 22/Ambiguous 76，两套判据结果
不完全一样，属预期之中——判据本身变了）。

**已知局限（沿用旧结论，尚未重新验证是否在新117站上依然成立）**：07脚本发现这个
方向判断对状态空间维度不稳健（尤其Inhibit组，加一维温度后近半反转，见
[00_cross_module_issues.md](00_cross_module_issues.md)第1条）——07/09-11 号脚本
还没有用新的117站清单重跑，图和稳健性数字目前对应的是旧134站口径，**不要直接拿
旧数字去描述新的117站结果**。

**未决**：因为方向不稳健，这套分类结果目前**只适合探索性描述，不建议直接拿去做
驱动因素回归**——如果要往下做，需要先解决方向稳健性问题（见
[00_cross_module_issues.md](00_cross_module_issues.md)第1条），或者改用连续型
指标（见下）。

## Recovery/Resilience 交叉验证支线（`pipelines/recovery_*.R`）

不用离散方向标签，直接检验"因果方向/强度"和"独立测算的恢复力指标"是否一致，作为
方向分类结果的外部佐证。已改为统一读 12 号脚本的站点清单。

- `recovery_vs_134_direction.R`：Schwalm式事件-恢复构建(每季最大VPD窗口=事件,
  事件前2窗均值=基线, 回到基线=恢复, 季末未恢复=右删失) + 生存分析(Kaplan-Meier
  +站点聚类稳健SE的Cox回归), 按上面134站支线的方向标签(Promote/Ambiguous/
  Inhibit)分组。**2026-09-07新结果(117站,2001个station-year事件)**：log-rank
  p=0.0002；Cox(参照组=Promote)：Ambiguous和Inhibit都显著比Promote恢复更快
  (HR均≈1.20, p<0.01)——**Promote组是三组里恢复最慢的**，这个方向和用旧127站
  判据时的结果(同样Promote最慢，但Inhibit当时只是边缘显著p=0.097)方向一致、
  显著性更强，说明这个"Promote反而恢复慢"的模式对判据变动有一定稳健性。
  **已知局限**：direction_2v本身不稳健(见上)，此结果的可信度上限受此制约；事件
  定义朴素(每季单一最大VPD窗口)，未做参数敏感性检验。

- `recovery_ar1_pure.R`：不涉及方向标签，纯粹在SIF自身距平序列上比较两种连续型
  恢复力指标——纯AR(1)(Forzieri et al. 2022式, `sif_anom(t)~phi*sif_anom(t-1)`)
  vs 控制了VPD的ARX(De Keersmaecker et al. 2015式,
  `sif_anom(t)~phi*sif_anom(t-1)+beta*vpd_z(t)`)。**2026-09-07新结果(117站)**：
  两者几乎完全一致(Pearson r=0.999)——说明φ这个连续型指标对"要不要控制VPD"这种
  建模选择很稳健，和上面离散方向标签"判据一换、结论就翻"形成鲜明对比。**提示**：
  如果要往下做driver回归，连续型的φ可能比离散方向标签更适合当因变量。

- `recovery_schwalm_style_demo.R`：10站不筛因果的feasibility demo，用于论证"不
  先做CCM因果确认、直接对全部站点算recovery"这种(部分resilience文献常见的)做法
  在本数据上行不通——右删失率0-56.5%，事件对齐复合轨迹看不出恢复信号。这个结论
  不依赖具体判据版本，未重跑。

## 已删除的历史支线（不要再花时间找这些文件，git历史里能找到但不要用）

- 44-60 全部driver回归脚本及其grp2方向标签（见上，2026-09-07删除）。
- `recovery_vs_ccm_direction.R`、`recovery_arx_index.R`：用grp2做的recovery交叉
  验证，grp2已弃用，这两个脚本已删除（2026-09-05），结论已被上面用新判据的重跑
  结果取代。
