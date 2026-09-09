---
name: setup-env
description: One-shot installer for the TiPhD benchmark R + Python environments (R 4.4/Seurat 5 + tiphd-torch/PyTorch + tiphd-stats/TiRank + optional tiphd-py310/Statescope). Detects R / Rscript / Python / conda, installs CRAN + Bioconductor + custom R packages, creates the tiphd-torch conda env (deep learning algos), tiphd-stats conda env (Python 3.9 + lifelines/optuna/leidenalg/python-igraph), and optional tiphd-py310 env (Statescope BLADE), pip-installs custom Python packages. Idempotent at every step. Use when the user says "setup env", "install TiPhD env", "装环境", "安装依赖", or before invoking any algorithm skill.
---

# setup-env Skill

> 一次性安装 TiPhD 评测所需 R + Python 环境（**多环境**）。**幂等**——可重复执行，已装的会跳过。

---

## 一、本 skill 做什么

收到用户指令后，自动完成 6 步安装：

1. **探测**——查 R / Rscript / Python / pip / conda / modelscope
2. **R 依赖**——装 CRAN + Bioconductor + 5 个自定义 R 包
3. **Conda tiphd-torch 环境**——从 `environment.yml` 建（深度学习算法用）
4. **Conda tiphd-stats 环境**——从 `environment-tirank.yml` 建（TiRank算法 + scPER MAGIC 用）
5. **Python 自定义包**——`pip install -e` SIDISH / scTREND / scSurv / TiRank（各装对应环境）
6. **验证**——`library()` + `import` 自检

---

## 二、关键路径

- 共享资源根：`.claude/skills/_toolkit/`
- 探测脚本：`_toolkit/detect_env.sh`
- R 安装脚本：`_toolkit/install_R_packages.R`
- Conda tiphd-torch 环境：`_toolkit/environment.yml`
- Conda tiphd-stats 环境：`_toolkit/environment-tirank.yml`（Python 3.9 + lifelines/optuna/leidenalg/python-igraph）
- Python 自定义包安装：`_toolkit/install_python_packages.sh`
- 第三方包源码：`_toolkit/third_party/R/` 与 `_toolkit/third_party/python/`

---

## 三、标准作业流程

### 步骤 1 — 探测环境

```bash
SKILL_DIR="<path to .claude/skills>"
bash "$SKILL_DIR/_toolkit/detect_env.sh"
```

- 缺 R / Rscript：提示装 R **4.4+**（数据为 Seurat v5 Assay5 格式，R 4.1 + Seurat 4 无法正确读取）；Windows 用户可设 `export TIPHD_RSCRIPT=/path/to/Rscript.exe`。
- R 版本过低（< 4.4）：**强烈建议新建 conda 环境** `conda create -n tiphd_r44 -c conda-forge r-base=4.4.2 r-seurat r-matrix r-nnls r-survival bioconductor-singlecellexperiment`。
- 缺 Seurat 或 Seurat < 5.0：数据 `.rds` 为 Assay5 格式，需 Seurat ≥ 5.0；Seurat 5 要求 Matrix ≥ 1.6.4（R 4.1 下 Matrix 最高 1.5.x，不兼容）。
- 缺 gcc（Windows）：提示装 Rtools44。
- 缺 conda：建议装 Miniconda；或跳过 conda 用 `python -m venv`。
- 缺 modelscope：提示 `pip install modelscope`（setup-data 用）。

### 版本兼容性矩阵

| R 版本 | Seurat 版本 | Matrix 版本 | 能读 v5 .rds? | 推荐 |
|--------|------------|------------|--------------|------|
| 4.4.x | 5.x | ≥ 1.6.4 | ✅ 完整支持 | **推荐** |
| 4.1.x | 4.x | 1.5.x | ❌ 缺 rownames/colnames | 不推荐 |
| 4.1.x | 5.x | ≥ 1.6.4 | ❌ R 4.1 无法装 Matrix ≥ 1.6.4 | 不可行 |
| 4.4.x | 4.x | 1.6.x | ⚠️ attr 兜底但可能缺 dimnames | 不推荐 |

### 工具 → 环境映射表

| 工具 | 语言 | 所需环境 | 关键依赖 | 数据格式 |
|------|------|---------|---------|---------|
| **MuSiC** | R | R 4.4 + Seurat 5 | Seurat, Matrix, SingleCellExperiment, nnls, survival | .rds (Assay5) |
| **scAB** | R | R 4.4 + Seurat 5 | Seurat, scAB(自定义) | .rds (Assay5) |
| **SCIPAC** | R | R 4.4 + Seurat 5 | Seurat, SCIPAC(自定义) | .rds (Assay5) |
| **Scissor** | R | R 4.4 + Seurat 5 | Seurat, Scissor(自定义, Rcpp) | .rds (Assay5) |
| **scPAS** | R | R 4.4 + Seurat 5 | Seurat, scPAS(自定义, Rcpp) | .rds (Assay5) |
| **ScPP** | R | R 4.4 + Seurat 5 | Seurat, Matrix, AUCell, survival | .rds (Assay5) |
| **scSTAR2** | R | R 4.4 + Seurat 5 | Seurat, Matrix, pls, MASS | .rds (Assay5) |
| **PIPET** | R | R 4.4 + Seurat 5 | Seurat, PIPET(自定义) | .rds (Assay5) |
| **SCAD** | Python | tiphd-torch (Python 3.9) | torch, anndata, scikit-learn | .h5ad |
| **SIDISH** | Python | tiphd-torch (Python 3.9) | torch, scanpy, anndata, SIDISH(自定义) | .h5ad |
| **scDEAL** | Python | tiphd-torch (Python 3.9) | torch, anndata, scikit-learn, DaNN(自定义) | .h5ad |
| **scSurv** | Python | tiphd-torch (Python 3.9) | torch, anndata, scsurv(自定义) | .h5ad |
| **scTREND** | Python | tiphd-torch (Python 3.9) | torch, scanpy, anndata, sctrend(自定义) | .h5ad |
| **scPER** | Python+R | tiphd-torch + R 4.4 | anndata(Python) + Seurat/Matrix(R subprocess) | .h5ad + .rds |
| **Statescope** | Python | **tiphd-py310 (Python 3.10+)** | torch, numba, anndata, autogenes, BLADE(自带) | .h5ad |
| **TiRank** | Python | tiphd-stats (Python 3.9) | lifelines, optuna, leidenalg, python-igraph, tirank(自定义) | .h5ad |

> **总结**：所有 R 类工具(8个)共享同一个 R 4.4 + Seurat 5 环境；深度学习类 Python 工具(5个)用 `tiphd-torch`；统计/图模型类 TiRank 用独立的 `tiphd-stats`（numpy/pandas 版本约束隔离）；scPER 需要 Python + R 双环境。
>
> **不需要为每个工具单独建环境或 skill**——setup-env 统一准备 3 个环境(R 4.4 / tiphd-torch / tiphd-stats)，各工具 skill 只需在 SKILL.md 中注明使用哪个环境即可。

### 步骤 2 — 装 R 依赖（与算法类型无关）

```bash
Rscript "$SKILL_DIR/_toolkit/install_R_packages.R"
```

- CRAN：Seurat / Matrix / survival / Rcpp / RcppEigen / AUCell / nnls / pls / dplyr / tidyr 等
- Bioconductor：`preprocessCore`
- 自定义 R 包（源码在 `_toolkit/third_party/R/`）：`scAB` `SCIPAC` `scPAS` `Scissor` `PIPET`
- **首次 5–30 分钟**（CRAN 下载 + scPAS/Scissor 的 Rcpp 编译）。

### 步骤 3 — 建 tiphd-torch Conda 环境（PyTorch 深度学习类工具）

适用工具：**SCAD / SIDISH / scDEAL / scSurv / scTREND**（VAE / GAN / 域自适应等深度学习算法）。

```bash
cd "$SKILL_DIR/_toolkit"
if conda env list | grep -qE "^\s*tiphd-torch\s"; then
  conda env update -f environment.yml --prune
else
  conda env create -f environment.yml
fi
conda activate tiphd-torch
```

- 装 PyTorch(CPU) + scanpy + anndata + numpy + pandas + scipy + scikit-learn。
- **Python 3.9**。
- 自定义包(SIDISH/scTREND/scSurv)在步骤5安装。

### 步骤 4 — 建 tiphd-stats Conda 环境（统计学习/图模型类工具）

适用工具：**TiRank**（Cox + 聚类 + 基因对图模型）+ **scPER 的 MAGIC 插补步骤**。

独立于 tiphd-torch 的原因：TiRank 依赖 `lifelines`/`optuna`/`leidenalg`，且要求 `numpy<2, pandas<2`，与深度学习环境的版本约束可能冲突。

```bash
cd "$SKILL_DIR/_toolkit"
if conda env list | grep -qE "^\s*tiphd-stats\s"; then
  conda env update -f environment-tirank.yml --prune
else
  conda env create -f environment-tirank.yml
fi
conda activate tiphd-stats
```

- 装 `lifelines`(Cox PH) `optuna`(超参搜索) `leidenalg`+`python-igraph`(聚类) —— TiRank 必须
- 装 `magic-impute` —— scPER 的 MAGIC 插补用
- **Python 3.9 + numpy<2 + pandas<2**（TiRank/magic-impute 兼容性约束）

### 步骤 4b — 建 tiphd-py310 Conda 环境（Statescope，可选）

适用工具：**Statescope**（BLADE 贝叶斯反卷积源码用了 Python 3.10+ 的类型语法 `X | None`）。

仅在运行 Statescope 时需要；与 tiphd-torch 独立。

```bash
conda create -y -n tiphd-py310 -c pytorch -c conda-forge \
  python=3.10 pytorch numpy pandas scikit-learn scanpy anndata \
  numba dill joblib seaborn statsmodels
conda activate tiphd-py310
pip install autogenes deap cachetools requests tqdm
```

- **Python 3.10+**；装 CUDA 版 PyTorch 可启用 BLADE 的 GPU 加速。

### 步骤 5 — 装 Python 自定义包

**tiphd-torch 环境**（跑 SIDISH / scTREND / scSurv / SCAD / scDEAL）：

```bash
conda activate tiphd-torch
bash "$SKILL_DIR/_toolkit/install_python_packages.sh"
```

依次 `pip install -e`：
- `_toolkit/third_party/python/SIDISH`
- `_toolkit/third_party/python/scTREND`
- `_toolkit/third_party/python/scsurv`（目录名 scSurv，import 名 scsurv）
- `_toolkit/third_party/python/scPER`（ADAE 部分；MAGIC 步骤需在 tiphd-stats 环境）

**tiphd-stats 环境**（跑 TiRank 算法本身 + scPER 的 MAGIC 步骤）：

```bash
conda activate tiphd-stats
pip install -e "$SKILL_DIR/_toolkit/third_party/python/TiRank"
```

### 步骤 6 — 验证

```bash
# R（任意环境都跑）
Rscript -e "suppressPackageStartupMessages({library(scAB); library(SCIPAC); library(scPAS); library(Scissor); library(PIPET)}); cat('R deps OK')"

# tiphd-torch 深度学习环境
conda activate tiphd-torch
python -c "from SIDISH.SIDISH import SIDISH; from sctrend import workflow; from scsurv.workflow import run_scSurv; import torch, scanpy, anndata; print('tiphd-torch OK')"

# tiphd-stats 统计/图模型环境
conda activate tiphd-stats
python -c "from tirank.Model import setup_seed; import lifelines, optuna, leidenalg, igraph; print('tiphd-stats OK')"
```

三者都打印 "OK" 即完成。

---

## 四、conda 环境与算法对应表

| 算法 | 用哪个 conda env |
|------|-----------------|
| MuSiC / scAB / SCIPAC / scPAS / ScPP / Scissor / PIPET / scSTAR2 | （无，仅 R 4.4 + Seurat 5） |
| SCAD / SIDISH / scDEAL / scSurv / scTREND | tiphd-torch |
| scPER | tiphd-torch + R 4.4（MAGIC 步可选 tiphd-stats） |
| **TiRank** | **tiphd-stats**（必须） |


---

## 五、用户的可选参数

| 用户说法 | skill 行为 |
|---------|-----------|
| "setup env" / "装环境" | 步骤 1 → 6（完整） |
| "只装 R 包" | 步骤 1 + 步骤 2 |
| "只装 tiphd-torch env" | 步骤 1 + 步骤 3 + 步骤 5（tiphd-torch 部分） |
| "只装 tiphd-stats env" | 步骤 1 + 步骤 4 + 步骤 5（tiphd-stats 部分） |
| "verify env" / "验证环境" | 步骤 1 + 步骤 6 |
| "Rscript 不在 PATH" | 提示设 `export TIPHD_RSCRIPT=/abs/path/to/Rscript` |

---

## 六、注意事项

1. **Windows 必须先装 Rtools44**，否则 scPAS/Scissor 的 C++ 代码无法编译。
2. **3 个环境**：R 4.4 + Seurat 5（系统 R 或 conda）/ tiphd-torch / tiphd-stats。
3. **conda activate 必须显式**——每个新 shell 都要重新激活；考虑加到 `~/.bashrc`。
4. **首次 install 慢**（5–60min），耐心等待。
5. **tiphd-stats 环境的 lifelines/optuna/leidenalg/python-igraph 必需**，否则 TiRank 跑不起来。
6. **不用 sudo**——所有包装在用户空间。

---

## 七、失败兜底

| 现象 | 处理 |
|------|------|
| `Rscript: command not found` | 装 R；或设 `export TIPHD_RSCRIPT=/abs/path/to/Rscript` |
| `cannot find function "sc.bulk.pca"` 等 | R 包没装好，重跑 `install_R_packages.R` |
| `Rcpp compilation failed` | 缺 Rtools44（Win）/ 缺 Xcode CLT（Mac） |
| `ERROR; return code from pthread_create() is 22` | 高核心数服务器上 preprocessCore 的 OpenMP 线程创建失败；`install_R_packages.R` 会自动检测并以 `--disable-threading` 重装，也可运行前设 `OMP_NUM_THREADS=1`（批量调度器 `run_batch` 默认已限制线程） |
| SIDISH `pip install -e` 解析依赖失败（要求 Python 3.10+） | `install_python_packages.sh` 已改用 `--no-deps` 安装包体并单独补装 pyro-ppl/torch_geometric/shap 等运行时依赖 |
| PIPET 安装时提示缺 survminer/DESeq2/tidyverse | 这些是包内**可视化/marker 提取**辅助功能的可选依赖，benchmark 核心不需要；已移到 `Suggests`，`setup-env` 安装的版本不要求它们 |
| Statescope 报 `unsupported operand type(s) for |` | BLADE 源码用了 Python 3.10+ 语法，需用 `tiphd-py310`（Python 3.10）环境运行 |
| `conda: command not found` | 装 Miniconda；或跳过 conda 用 `python -m venv` |
| `ModuleNotFoundError: SIDISH` / `tirank` | 跑对应 conda env 的 `install_python_packages.sh` |
| `No module named 'lifelines'` (tiphd-stats) | 在 tiphd-stats env 跑 `pip install lifelines optuna leidenalg python-igraph` |
| `pip install` 卡死 | 换 PyPI 镜像：`pip install -i https://pypi.tuna.tsinghua.edu.cn/simple` |

---

## 八、不在 skill 范围内

- 装 R 本身（用 https://cran.r-project.org/）。
- 装 Python 本身。
- 装 conda 本身。
- 装 Rtools44（Windows）。
- 装 Xcode CLT（Mac）。
