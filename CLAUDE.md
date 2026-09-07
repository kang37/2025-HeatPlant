# 2025-HeatPlant

中国 900+ 气象站（原始 924，协变量齐全后 898/904，视口径不同）2000–2022 年 5–9 月，
用 CCM（收敛交叉映射）+ S-map 研究 VPD 对 SIF（HCSIF，500m/8天合成）的因果作用及其影响因素。
R + data.table + ggplot2 + showtext。`data_raw/` 与 `data_proc/` 已 gitignore（靠 Dropbox
同步，git 未追踪的中间产物本机可能没有，需要时重跑脚本，不要假设文件存在）。

## 项目结构：三个先后模块

理想是顺序流水线，实际经常在下游模块里发现上游的方法问题，导致反复回补——见
[docs/00_cross_module_issues.md](docs/00_cross_module_issues.md)，开工前先看这份。

| 模块 | 脚本位置 | 状态文档 |
|---|---|---|
| 一、数据获取与协变量装配 | 根目录 `13_*.R` ~ `22_*.R` | [docs/01_data_acquisition.md](docs/01_data_acquisition.md) |
| 二、CCM 因果分析 + S-map | `pipelines/`（含 `hcsif_buf1000/`、`hcsif_buf2000/`、`hcsif_buf3000/`） | [docs/02_ccm.md](docs/02_ccm.md) |
| 三、因果确认标准 + 方向分类（原 44-60 driver 回归主线 2026-09-07 已删除） | `pipelines/hcsif_buf1000/12_ccm_causal_confirmation.R`（唯一CCM因果确认判据）+ `06/08_*_134.R`（方向分类）+ `pipelines/recovery_*.R`（恢复力交叉验证） | [docs/03_classification.md](docs/03_classification.md) |

`note/` 目录是论文/报告稿（qmd/html/docx），不是进度文档，不要往里面写状态记录。

## 工作流约定（重要）

**每完成一个任务后，把结果追加到对应模块的状态文档，不要只留在对话里。**
每份模块文档结构固定为：
- `## 状态速览`：当前对该模块的最新理解，直接覆盖更新，保持简短。
- `## 更新日志`：追加式，每条以日期开头，写清楚做了什么、发现了什么、留下什么坑，
  不要覆盖旧条目。

如果任务发现的问题涉及不止一个模块（比如 CCM 的方法缺陷影响分类模块的判据），
写入 [docs/00_cross_module_issues.md](docs/00_cross_module_issues.md) 而不是塞进单个模块文档。

## 环境坑（跨模块通用，新环境/新对话先看这里）

- Rscript 默认 locale 下中文标识符会解析失败，list 名和列名要加引号或用 ASCII。
- `setorder(dt, -abs(x))` 在 data.table 里不合法，用 `dt[order(-abs(x))]`。
- 子进程放在 `while read` 循环里必须加 `</dev/null`；管道接 `head`/`tail` 会因
  SIGPIPE 让上游误判失败。
- `pmax(0, 矩阵)` 会丢掉 `dim` 属性（属性只从第一参数复制），Conley 核矩阵要手动补
  `dim(K)`。
- xgboost 3.2.1.1：`xgb.cv` 最佳轮数在 `cv$early_stop$best_iteration`，不是
  `cv$best_iteration`；TreeSHAP 用 `predict(..., predcontrib=TRUE)`，交互用
  `predinteraction=TRUE`，都是精确解，不是近似。
- 已装包：survival, nnet, glmnet, ranger, data.table, sandwich, lmtest, car, MASS,
  xgboost。未装：lightgbm, shapviz, fastshap, DALEX, iml, cmprsk, mstate, brms,
  relaimpo（需要时先确认是否已补装）。
- GLC_FCS30D 瓦片命名 `Exxx`=西边界、`Nyy`=**北**边界（不是南边界）。
- curl 下载大文件不能加 `--retry`：断流重试会重复输出已传字节，导致 zlib 报错；
  改为自己重试 + 精确长度校验。
- rEDM 2.0.2 的 `SMap()$coefficients` 系数表列名是 `"Time"`(大写)，不是
  `"time"`；`co[["time"]]` 不会报错，静默返回 NULL，容易埋坑。系数表行序本来
  就与输入 `dataFrame` 的行序一致，直接用 `seq_len(nrow(co))` 做时间索引更稳妥。
- ggplot 画中文/CJK 文字，`font_add()`/`showtext_auto()` 前必须先
  `Sys.setlocale("LC_ALL", "en_US.UTF-8")`，否则不报错，只是图上每个中文字符
  静默渲染成一个占位圆点。
