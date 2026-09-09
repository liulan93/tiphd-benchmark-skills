---
name: scsurv
description: Run the scSurv benchmark pipeline. scSurv is an A-class survival analysis algorithm in Python using a VAE (Variational Autoencoder) that learns a latent representation of single-cell data and predicts survival-associated cell subpopulations. Cells with beta_z scores beyond +/-0.5 standard deviations are labeled as Rank+/Rank-. Use when the user says "用scSurv跑评测", "scSurv评测", "run scSurv evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scSurv.
---

# scSurv Benchmark Skill

scSurv uses a Variational Autoencoder to model single-cell transcriptomics with a survival-oriented latent space. It computes a per-cell risk score (beta_z) and identifies cells significantly associated with survival outcomes based on a threshold of +/-0.5 standard deviations from the mean.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | Python (PyTorch) |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (1000 Cox PH fits) |
| plus-only | No |
| GPU | Optional (CPU works) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.9 | Runtime |
| **torch** | >= 1.12 | VAE and neural network training |
| **anndata** | >= 0.8 | Reading `.h5ad` single-cell data |
| **numpy** | any | Numerical computation |
| **pandas** | any | Data handling |
| **lifelines** | any | Cox PH survival models |

scSurv custom package (`scsurv`) is installed from `third_party/python/scSurv`.

### Recommended conda environment

```bash
conda env create -f skills/_toolkit/environment.yml   # tiphd-torch env
conda activate tiphd-torch
pip install lifelines
```

Or manually:

```bash
conda create -y -n tiphd-torch -c pytorch -c conda-forge \
  python=3.9 pytorch cpuonly numpy pandas scikit-learn
conda activate tiphd-torch
pip install anndata lifelines
pip install -e skills/_toolkit/third_party/python/scSurv
```

## Prerequisites

1. **Environment** -- run `setup-env` skill.
2. **Data** -- run `setup-data` skill.
3. **R** -- evaluation layer requires R >= 4.4 with Seurat >= 5.0, survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scSurv/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scsurv

# Full benchmark
python batch_run.py
# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single cancer / pair
Rscript evaluate.R LUAD
python run_scSurv_pair.py AML GSE116256 TCGA
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scSurv_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scSurv_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Algorithm Details

1. **VAE encoder-decoder** learns a low-dimensional latent representation of single cells.
2. A **Cox PH model** is fit on bulk data, and the coefficient vector (beta) is projected onto the latent space.
3. Each cell receives a **beta_z score** (projection of its latent representation onto the survival direction).
4. Cells with `beta_z > mean + 0.5*std` are labeled Rank- (risk); cells with `beta_z < mean - 0.5*std` are labeled Rank+ (protective).

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 5-15 min |
| CRC | GSE132465 | ~47K | 10-20 min |
| HCC | GSE149614 | ~31K | 5-15 min |
| LUAD | GSE127465 | ~24K | 5-15 min |
| **GC** | **GSE183904** | **~137K** | **45-120 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). scSurv trains a VAE on all cells, so GC pairs (~137K cells) will likely exceed this. Run GC pairs directly or increase timeout:

```bash
python run_scSurv_pair.py GC GSE183904 GSETCGA
```

## Notes

- scSurv reads single-cell data from `.h5ad` files.
- The VAE training uses early stopping with a validation split.
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: scsurv` | Install: `pip install -e third_party/python/scSurv` |
| `ModuleNotFoundError: lifelines` | Install: `pip install lifelines` |
| `TIMEOUT (>7200s)` | Expected for GC pairs; see "Dataset Size and Timeout" section |
| `FileNotFoundError: .h5ad` | Run `setup-data` |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
