---
name: run-benchmark
description: Run any of the 16 TiPhD benchmark algorithms end-to-end (single algorithm, single cancer, or all). Handles both A-class cell-level algorithms (scAB, SCIPAC, scPAS, ScPP, Scissor, PIPET, scSTAR2, SIDISH, scTREND, scSurv, SCAD, scDEAL, TiRank) and B-class deconvolution algorithms (MuSiC, scPER, Statescope) using the project's standardized 3-layer pipeline: `batch_run.py` (parallel per-pair run) → `evaluate.R` (permutation test + gold-standard comparison + metrics). Use when the user wants "run benchmark for X", "evaluate X on Y cancer", "跑 X 算法的评测", or invoke any specific algorithm by name. Triggers on "跑评测" / "运行算法评测" / "benchmark 流程" / "run the benchmark" / "evaluate algorithm".
---

# 通用 TiPhD 算法评测 Skill

> 适配本项目 **16 个算法**的统一 3 层流程：run（单配对）→ batch_run（批量调度）→ evaluate（置换+金标准+指标）。
> 不绑定单一算法；用户说哪个就跑哪个。

> **自包含版本**：所有脚本和第三方包在 `.claude/skills/<algo>/` 内。首次运行前请先调用 `setup-env` 与 `setup-data`。

---

## 〇、Prerequisites（首次运行必看）

1. **环境**——先调用 `setup-env` skill 装 R + Python（含 tiphd-torch + tiphd-stats 两个 conda env）。
2. **数据**——再调用 `setup-data` skill 从 ModelScope 下载数据集到 `<cwd>/data/`。

输入输出路径默认 `<cwd>/data` 与 `<cwd>/results`，可由 `TIPHD_DATA_DIR` / `TIPHD_OUT_DIR` 环境变量覆盖。

---

## 一、本 skill 做什么

收到用户指令后（指定或不指定算法 / 癌种），自动：

1. **解析用户意图**：从指令里识别（a）算法名、（b）癌种（可选）、（c）阶段（批量 / 评估 / 全跑）。
2. **校验环境**：R / Python / 数据目录 / 必备包（若缺失，触发 auto-recovery）。
3. **进入对应算法 skill 目录**（`.claude/skills/<algo>/`），跑第 2 层 + 第 3 层。
4. **汇报结果**：最终指标 CSV 路径、均值汇总、batch 日志摘要。

---

## 二、16 个算法速查表

### A 类（细胞级算法，13 个）— 评估用 `evaluate_cell_level`

| 算法 | 语言 | run 脚本 | plus_only | 自定义包 / 备注 |
|------|------|---------|:---:|------|
| scAB | R | `run_scAB_pair.R` | ✅ | scAB |
| SCIPAC | R | `run_SCIPAC_pair.R` |  | SCIPAC |
| scPAS | R | `run_scPAS_pair.R` |  | scPAS（自带 2000 次置换） |
| ScPP | R | `run_ScPP_pair.R` |  | （无，Cox + AUCell） |
| Scissor | R | `run_Scissor_pair.R` |  | Scissor（Rcpp 自定义包） |
| PIPET | R | `run_PIPET_pair.R` |  | PIPET（二分类） |
| scSTAR2 | R | `run_scSTAR2_pair.R` |  | （无，OPLS-DA + 本地 `src/`） |
| SIDISH | Python | `run_SIDISH_pair.py` | ✅ | SIDISH |
| scTREND | Python | `run_scTREND_pair.py` |  | sctrend |
| scSurv | Python | `run_scSurv_pair.py` |  | scsurv |
| SCAD | Python | `run_SCAD_pair.py` |  | （无，torch VAE） |
| scDEAL | Python | `run_scDEAL_pair.py` |  | （无，torch DaNN） |
| **TiRank** | Python | `run_TiRank_pair.py` |  | **TiRank（REO+NN+MMD，tiphd-stats conda env）** |

### B 类（反卷积/样本级算法，3 个）— 评估用 `evaluate_deconv`

| 算法 | run 脚本 | 自定义包 / 备注 |
|------|---------|------|
| MuSiC | `run_MuSiC_pair.R` | （无，加权 NNLS） |
| scPER | `run_scPER_pair.py` | （无，ADAE + xgboost；MAGIC 步可选 tiphd-stats env） |
| Statescope | `run_Statescope_pair.py` | （无，BLADE 贝叶斯反卷积） |

> A 类与 B 类的 run/evaluate 实现细节不同：B 类只做反卷积，置换在评估层；A 类 run 直接出 `Rank_Label`。

### 算法方法分类（四类）

| 分类 | 算法 |
|------|------|
| **deep learning** | SIDISH / scTREND / scSurv / SCAD / scDEAL / **TiRank** |
| **sparse linear model** | Scissor(L1) / scPAS(L0) / SCIPAC(elastic-net) |
| **classical statistical** | ScPP / scSTAR2 / PIPET / scAB |
| **deconvolution** | MuSiC / scPER / Statescope |

---

## 三、参数解析

从用户自然语言指令中识别：

| 字段 | 识别方式 | 默认值 |
|------|---------|--------|
| `algo` | 任意大小写匹配上表 16 个名字之一 | **必填**，否则用 `AskUserQuestion` 询问 |
| `cancer` | 5 个癌种名 `AML` / `CRC` / `HCC` / `LUAD` / `GC` 之一 | 全部 5 个 |
| `phase` | "批量" / "batch" → 步骤 2；"评估" / "evaluate" → 步骤 3；其他 → 全部 | 全部 |
| `out_suffix` | A=`.csv` / B=`_proportions.csv`，从 `batch_run.py` 一行配置读 | 自动 |
| `plus_only` | scAB / SIDISH → TRUE | 自动 |
| `lang` | R / Python，从 `batch_run.py` 一行配置读 | 自动 |

> 实际**不要**自己解析——直接看对应算法文件夹的 `batch_run.py` 和 `evaluate.R` 文件的三行配置即可，那是权威源。

---

## 四、标准作业流程

### 步骤 1 — 解析与校验

```bash
SKILL_ROOT="<path to .claude/skills>"   # 通常是 <project>/.claude/skills
algo="<from user>"
test -d "$SKILL_ROOT/$algo" && echo "OK algo skill exists" || { echo "FAIL unknown algo: $algo"; exit 1; }
test -d "${TIPHD_DATA_DIR:-$(pwd)/data}" && echo "OK data"
which Rscript && which python
```

> 若 algo 拼写不在 16 个里（如 "MUSIC" → 提示用 "music"），用 `AskUserQuestion` 让用户选。

### 步骤 2 — 批量调度

```bash
cd "$SKILL_ROOT/<algo>"
python batch_run.py
```

- **44 个配对**（5 癌种 × 各自的 (sc, bulk) 笛卡尔积）逐个独立进程。
- 已存在输出 → SKIP（断点续跑）。
- 单配对失败不中断，日志 `_batch_log.txt`。
- 跑完看 Summary：`Total: 44, Done: N, Skipped: M, Failed: K`。
- 若 `K > 0` 列出失败配对，让用户决定是否重跑。

### 步骤 3 — 评估

```bash
cd "$SKILL_ROOT/<algo>"
Rscript evaluate.R             # 5 癌种全跑
Rscript evaluate.R <cancer>    # 只跑 1 个癌种（如果用户指定）
```

最终产出 5 个 CSV（每癌种一份）：

```
$SKILL_ROOT/<algo>/<cancer>/<algo>_<cancer>_final_metrics.csv
```

列：`scRNA, Bulk, Sp, TP, FP, FN, Precision, Coverage, False_Rate`。

### 步骤 4 — 回报

- 列出 5 个最终指标 CSV 路径。
- 用 R 简单算每癌种 Precision / Coverage / False_Rate 的 mean（NA 略过）。
- 报告 batch Done / Skipped / Failed。

---

## 五、用户的可选参数

| 用户说法 | skill 行为 |
|---------|-----------|
| "跑 MuSiC" / "run MuSiC" / "MuSiC 评测" | 全部流程（5 癌种） |
| "跑 ScPP on CRC" / "ScPP 评测 CRC" | 步骤 2 跑全部 44 配对 + 步骤 3 只跑 CRC |
| "只跑评估 ScAB on LUAD" | 跳过 batch，只跑 `cd skills/scab && Rscript evaluate.R LUAD` |
| "只跑批量 MuSiC" | 只跑 `cd skills/music && python batch_run.py` |
| "TiRank 单配对 CRC GSE144735 GSE17536" | `cd skills/tirank && python run_TiRank_pair.py CRC GSE144735 GSE17536` |
| "跑全部算法" / "benchmark all" | 跑全部 16 个算法的批量+评估，**长任务**，先确认 |

---

## 六、算法特殊点

| 算法 | 特殊点 |
|------|--------|
| scAB | plus_only=TRUE（无 Rank-），Coverage 偏低是预期 |
| SIDISH | plus_only=TRUE（同上），Python 跑、**GPU 自动检测**（有 CUDA 版 torch 即用 GPU，无则回退 CPU），无 checkpoint，杀进程整对重跑 |
| scPAS | 算法内部**自带 2000 次置换**，evaluate 不再做细胞级置换（但仍跑 `permute_cell_level` —— 那是评估层的统一骨架，可接受） |
| Scissor | R 算法，需 Seurat ≥ 5.0；run 脚本内部将 Assay5 转为 Assay 以兼容 |
| scSTAR2 | source `scstar2/src/` 的 9 个 R 函数；包内含本地 patch（PLSconstruct/OPLSDA） |
| SCAD / scDEAL | torch 二分类，CPU 也行（自动 CPU） |
| MuSiC / scPER / Statescope | B 类，只做反卷积；样本级置换在 evaluate 完成 |
| MuSiC / scPER | B 类，只做反卷积；样本级置换在 evaluate 完成 |
| **Statescope** | B 类反卷积，**必须 tiphd-py310**（Python 3.10+）。耗时由细胞**类型数**主导（而非样本数），高类型数数据默认 Nrep=10 可能极慢，宜绕过批量超时直接跑；进度只有裸迭代号、AutoGeneS 阶段长时间静默均正常；无 checkpoint。`STATESCOPE_NREP=1` 仅用于显式标注的冒烟验证，不可当基准结果 |
| **TiRank** | **用 tiphd-stats conda env**（Python 3.9 + lifelines/optuna/leidenalg/python-igraph）。临床表必须恰好 2 列 [time, event]——run 脚本已按 config 的 tcol/scol 自动裁剪（自定义数据要自己保证两列）。GPU 需显式 `TIRANK_GPU=1` 且 env 内是 CUDA 版 torch（默认 yml 装的是 cpuonly）；CPU 能跑但慢 |
| AML | bulk 名是 `TCGA` / `wave12` / `wave34`（非 GSE 编号） |
| GC | `GSE183904` 137K 细胞，TiRank/SIDISH 可能数小时；Statescope 的主要瓶颈则是细胞类型数而非细胞数 |
| LUDA 文件名 | `gold_standard_all_LUDA.csv`（原始拼写，**不是** LUAD） |

---

## 七、注意事项（通用）

1. **不子采样**：所有 run 脚本读取完整细胞。
2. **相对路径**：所有路径走 `config.R` / `config.py`。
3. **幂等**：batch_run 跳过已有输出；evaluate 可重跑。
4. **失败兜底**：单配对 FAIL 不中断；但**全部 FAIL** 要怀疑环境问题。
5. **结果信任度**：plus-only 算法（scAB / SIDISH）Coverage 永远偏低，横向比较时只看 Precision。

---

## 八、失败兜底速查

| 现象 | 处理 |
|------|------|
| 算法 skill 文件夹不存在 | 拼写错误，给出最接近的 16 个名字让用户选 |
| `Rscript: command not found` | 装 R；或设 `export TIPHD_RSCRIPT=/abs/path/to/Rscript` |
| `library(SCIPAC) ... there is no package` | 调 `setup-env` skill |
| `ModuleNotFoundError: tirank` | 在 TiRank env 跑 `pip install -e third_party/python/TiRank` |
| batch 大量 FAIL（>30%） | 环境问题，先看 `_batch_log.txt` 第一条 FAIL 的 ERR |
| evaluate 报 "missing X output" | batch 没产出该配对，**先**跑 batch |
| TiRank 报 `0 Risk genes and 0 Protective genes` | **先查临床表，不是信号弱**：必须恰好 [time, event] 两列（run 脚本已自动裁剪 TiPhD 数据）；上游 `except: continue` 会吞掉真实 Cox 异常，可手动拟合一次 CoxPH 看报错。详见 tirank/SKILL.md |
| Statescope 长时间无输出 / 疑似卡死 | AutoGeneS 阶段静默、BLADE 只打裸迭代号、EM 各轮耗时不均都正常；看 GPU 活跃度/CPU 时间判断是否真挂。高细胞类型数数据默认配置可能极慢；只有冒烟验证才用 `STATESCOPE_NREP=1` 并标注非默认 |
| 进程 TIMEOUT | 默认 7200s（2 小时，可由环境变量 `TIPHD_PAIR_TIMEOUT` 调整）；超大配对（如 GC 137K 细胞、高细胞类型数的 Statescope）可直接运行 `run_*_pair.*` 绕过批量超时 |
| `modelscope` 缺失（setup-data） | `pip install modelscope` |
| 单个癌种 OOM / 读 h5ad 后无 traceback 被杀 | 多为 dense X 展开（GC 137K 细胞最甚）；TiRank/SIDISH 的 run 脚本已内置 backed→CSR 流式读取，确认用的是随 skill 附带的 run 脚本；再考虑加内存或限制线程 |
| 进程被杀后重跑从头开始 | 预期行为：SIDISH/TiRank/Statescope 均无 checkpoint；只有已产出 CSV 的配对会被跳过，未完成对整对重来（Statescope 连 AutoGeneS 也重跑），长任务放在可靠会话里启动 |

---

## 九、不在 skill 范围内

- 修改算法本身（公式、置换次数、超参）——属于改代码。
- 跨算法对比、画图——属于后处理。
- 装 R / Python 包——按 `setup-env` skill 走。
- 装环境变量（`TIPHD_RSCRIPT` 等）——用户自己设。
