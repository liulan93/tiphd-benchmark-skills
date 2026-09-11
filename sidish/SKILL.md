---
name: sidish
description: Run the SIDISH benchmark pipeline. SIDISH is an A-class plus-only deep learning algorithm in Python (VAE + DeepCox with SHAP interpretation) that identifies risk-associated cells; it outputs only Rank+ / Background labels. Use when the user says "用SIDISH跑评测", "SIDISH评测", "run SIDISH evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for SIDISH.
---

# SIDISH Benchmark Skill

SIDISH uses a Variational Autoencoder with a DeepCox (survival-oriented) model to identify cells associated with poor survival. It is a **plus-only** algorithm: only `Rank+` (risk-associated) cells are reported, all others are `Background`.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | Python (PyTorch) |
| Output | Cell-level labels: Rank+ / Background |
| Permutation | Cell-level (Cox PH / concordance) in evaluate.R |
| plus-only | Yes |
| GPU | Auto-detected (VAE + DeepCox run on CUDA if available) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.9 | Runtime |
| **torch / torchvision** | >= 1.12 | VAE + DeepCox |
| **pyro-ppl** | any | Probabilistic components |
| **torch_geometric** | any | Graph neural-network layers |
| **scanpy / anndata** | any | Single-cell preprocessing |
| **shap** | any | Feature interpretation |
| **lifelines** | any | Concordance / Cox utilities |
| **imbalanced-learn** | any | Resampling |
| **numpy / pandas / statsmodels / sklearn / seaborn / matplotlib / bioinfokit** | any | Data handling |
| **SIDISH** | custom (in `third_party/python/SIDISH`) | The algorithm itself |

> **Note**: SIDISH's `setup.py` pins exact dependency versions (some requiring Python 3.10+). The `setup-env` skill therefore installs SIDISH with `pip install -e <path> --no-deps` and then installs the runtime dependencies listed above separately, so it works on Python 3.9.
>
> **GPU**: the run script auto-detects CUDA (`device = "cuda" if torch.cuda.is_available() else "cpu"`). The VAE + DeepCox training is much faster on a GPU; on CPU a single large pair can take many hours. Install a CUDA build of PyTorch (not `cpuonly`) to use the GPU.

### Recommended conda environment

```bash
# CPU: use environment.yml (pytorch cpuonly)
conda env create -f skills/_toolkit/environment.yml   # tiphd-torch env
# GPU: install a CUDA-matching pytorch build instead, e.g.
#   conda install -c pytorch -c nvidia pytorch pytorch-cuda
conda activate tiphd-torch
pip install pyro-ppl torchvision torch_geometric shap bioinfokit imbalanced-learn lifelines
pip install -e skills/_toolkit/third_party/python/SIDISH --no-deps
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (handles the SIDISH `--no-deps` install and its runtime deps).
2. **Data** -- run the `setup-data` skill.
3. **R** -- the evaluation layer (`evaluate.R`) needs R >= 4.4 with Seurat >= 5.0 and survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/SIDISH/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/sidish

# Full benchmark
python batch_run.py
# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
python run_SIDISH_pair.py <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `SIDISH_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `SIDISH_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

SIDISH trains a VAE + DeepCox per pair (50 training epochs by default); runtimes scale with cell count and are much shorter on the auto-detected GPU. The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs are the slowest and should be run directly if they time out:

```bash
python run_SIDISH_pair.py GC GSE183904 GSETCGA
```

There is **no checkpoint**: killing a pair partway through training discards all epochs and the pair restarts from scratch (only pairs whose output CSV already exists are skipped). Launch long pairs from a resilient session. The run script streams dense `.h5ad` as CSR at load time, so very large dense matrices do not need to be expanded in RAM.

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: pyro / torch_geometric / shap / bioinfokit / lifelines` | Install the missing runtime dependency in tiphd-torch (see Dependencies) |
| `No matching distribution found for scanpy==1.10.4` | The pinned version needs Python 3.10+; install SIDISH with `--no-deps` and let `setup-env` provide deps |
| Training loss becomes `nan` | Learning rate too high for the data; the run script uses a conservative lr (3e-4) and reduced epochs to avoid this |
| concordance index error in evaluate.R | Survival arrays have edge cases for some pairs; expected for a subset of pairs |
| `TIMEOUT (>7200s)` | Large dataset on CPU; run on a GPU (auto-detected) or run the pair directly (see above) |
| `FileNotFoundError: .h5ad` | Run `setup-data` to download data |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
