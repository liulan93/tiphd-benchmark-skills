---
name: sctrend
description: Run the scTREND benchmark pipeline. scTREND is an A-class deep learning algorithm in Python combining a VAE latent representation with DeepCOLOR batch correction and a segmented risk model; cells in the top/bottom 20th percentile of the beta-z risk score are labeled Rank-/Rank+. Use when the user says "用scTREND跑评测", "scTREND评测", "run scTREND evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scTREND.
---

# scTREND Benchmark Skill

scTREND learns a VAE-based latent representation of single cells, aligns bulk and single-cell distributions with DeepCOLOR, and fits a segmented risk model. Each cell receives a `beta_z` risk score; cells beyond the top/bottom 20th percentile are labeled Rank-/Rank+.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | Python (PyTorch) |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (Cox PH) in evaluate.R |
| plus-only | No |
| GPU | Optional (runs on CPU; GPU accelerates training) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.9 | Runtime |
| **torch** | >= 1.12 | VAE + DeepCOLOR neural nets |
| **scanpy** | any | Single-cell preprocessing |
| **anndata** | >= 0.8 | Reading `.h5ad` |
| **numpy / pandas** | any | Data handling |
| **sctrend** | custom (in `third_party/python/scTREND`) | The algorithm itself |

The custom `sctrend` package is installed editable via `setup-env`.

### Recommended conda environment

```bash
conda env create -f skills/_toolkit/environment.yml   # tiphd-torch env
conda activate tiphd-torch
```

Or manually:

```bash
conda create -y -n tiphd-torch -c pytorch -c conda-forge \
  python=3.9 pytorch cpuonly numpy pandas scikit-learn scanpy anndata
conda activate tiphd-torch
pip install -e skills/_toolkit/third_party/python/scTREND
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (creates the `tiphd-torch` env).
2. **Data** -- run the `setup-data` skill.
3. **R** -- the evaluation layer (`evaluate.R`) needs R >= 4.4 with Seurat >= 5.0 and survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scTREND/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/sctrend

# Full benchmark
python batch_run.py
# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
python run_scTREND_pair.py <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scTREND_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scTREND_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

scTREND trains a VAE per pair, which is compute-intensive. Typical wall-clock time on CPU:

| Cancer | Cells | Typical runtime per pair |
|--------|-------|--------------------------|
| AML | ~36K | ~30 min |
| CRC | ~47K | ~30-60 min |
| HCC | ~31K | ~30-60 min |
| LUAD | ~24K | ~25-45 min |
| **GC** | **~137K** | **>2 h** |

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). Pairs that time out should be run directly (no batch timeout):

```bash
python run_scTREND_pair.py GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: sctrend` | Install editable: `pip install -e third_party/python/scTREND` in tiphd-torch |
| `ModuleNotFoundError: torch/scanpy` | Activate tiphd-torch: `conda activate tiphd-torch` |
| `TIMEOUT (>7200s)` | VAE training on large data; run the pair directly (see above) |
| CUDA errors | Runs on CPU by default (`CUDA_VISIBLE_DEVICES=""`); GPU is optional |
| `FileNotFoundError: .h5ad` | Run `setup-data` to download data |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
