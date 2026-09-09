---
name: statescope
description: Run the Statescope benchmark pipeline. Statescope is a B-class (deconvolution) algorithm in Python using BLADE Bayesian latent-variable deconvolution to estimate cell-type proportions in bulk samples from a single-cell signature. Sample-level Cox permutation runs in the evaluate layer. Use when the user says "用Statescope跑评测", "Statescope评测", "run Statescope evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for Statescope.
---

# Statescope Benchmark Skill

Statescope performs Bayesian deconvolution (BLADE) of bulk RNA-seq using a single-cell-derived expression signature. It estimates per-sample cell-type fractions and, in the evaluate layer, tests association with survival.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | B (deconvolution) |
| Language | Python |
| Output | Sample x cell-type proportion matrix |
| Permutation | Sample-level (Cox PH) in evaluate.R |
| plus-only | No |
| GPU | Auto-detected (BLADE runs on CUDA if available) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.10 | BLADE source uses PEP-604 type unions (`X \| None`) |
| **torch** | >= 1.12 | Neural-network components (CUDA build for GPU) |
| **numba / dill / joblib** | any | JIT + parallel execution |
| **scanpy / anndata** | any | Reading `.h5ad` + log1p preprocessing |
| **numpy / pandas / scikit-learn / statsmodels** | any | Data handling |
| **autogenes / deap / cachetools** | any | Signature gene selection |
| **seaborn / matplotlib / requests / tqdm** | any | Utilities |

The Statescope/BLADE source is bundled in `third_party/python/Statescope` and loaded via `sys.path` by the run script (no pip install required).

> **Python version**: the BLADE source uses `pd.DataFrame | None` syntax which requires **Python 3.10 or later**. Use a dedicated Python 3.10 environment.
>
> **GPU**: BLADE is implemented in PyTorch and auto-detects CUDA (`torch.cuda.is_available()`). On CPU the Nrep Bayesian deconvolution is very slow; a CUDA GPU speeds it up substantially. Do not set `CUDA_VISIBLE_DEVICES=""` if you want GPU use.

### Recommended conda environment

```bash
# CPU-only
conda create -y -n tiphd-py310 -c pytorch -c conda-forge \
  python=3.10 pytorch cpuonly numpy pandas scikit-learn scanpy anndata \
  numba dill joblib seaborn statsmodels
# GPU (CUDA) — install a matching CUDA build of pytorch instead of cpuonly
conda activate tiphd-py310
pip install autogenes deap cachetools requests tqdm
```

## Prerequisites

1. **Environment** -- create the Python 3.10 environment above (or via `setup-env`).
2. **Data** -- run the `setup-data` skill.
3. **R** -- the evaluation layer (`evaluate.R`) needs R >= 4.4 with Seurat >= 5.0 and survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/Statescope/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/statescope

# Full benchmark (use the py310 interpreter)
/path/to/tiphd-py310/bin/python batch_run.py

# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
/path/to/tiphd-py310/bin/python run_Statescope_pair.py <cancer> <sc_name> <bulk_name>
```

## Data normalization note

BLADE expects **library-size-corrected linear, non-negative** bulk counts, while the TiPhD benchmark bulk data is log-normalized. The run script linearizes bulk input with `expm1()` (inverse of `log1p`) and clips negatives; the single-cell signature is expected to be log-normalized (anndata `.X`).

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `Statescope_<bulk>_<sc>_proportions.csv` | Sample x cell-type fraction matrix |
| `Statescope_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

BLADE runs Bayesian inference with repeated deconvolution (Nrep) and is computationally heavy. The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). Large single-cell datasets (especially GC, ~137K cells) may need longer or direct invocation:

```bash
python run_Statescope_pair.py GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'` | Python < 3.10; use the Python 3.10 environment |
| `ModuleNotFoundError: autogenes/deap/cachetools` | `pip install autogenes deap cachetools` in tiphd-py310 |
| `adata.X looks like raw counts; please normalise & log-transform` | The signature must be log-normalized `.h5ad` (TiPhD data already is) |
| `Bulk contains negative values` | Bulk must be non-negative; the run script clips and linearizes (expm1) automatically |
| `TIMEOUT (>7200s)` | BLADE is slow on large data; run the pair directly or reduce Nrep in the source |

## Scope

Does not modify algorithm parameters, install packages, or download data.
