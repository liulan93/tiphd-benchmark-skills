---
name: scdeal
description: Run the scDEAL benchmark pipeline. scDEAL is an A-class deep learning algorithm in Python using Deep Adversarial Neural Networks with Maximum Mean Discrepancy (MMD) loss for domain adaptation. It trains a predictor on bulk data and transfers labels to single cells via domain-invariant features. Cells in top/bottom 20% of predictions are labeled Rank-/Rank+. Use when the user says "用scDEAL跑评测", "scDEAL评测", "run scDEAL evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scDEAL.
---

# scDEAL Benchmark Skill

scDEAL (single-cell Deep Adversarial Learning) uses a neural network with:
- **Feature extractor**: Encodes both bulk and single-cell expression into shared space
- **Label predictor**: Binary classifier (survival Dead/Alive) trained on bulk
- **Domain discriminator**: Distinguishes bulk vs single-cell, with MMD loss for distribution alignment

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
| **torch** | >= 1.12 | Neural network framework |
| **anndata** | >= 0.8 | Reading `.h5ad` single-cell data |
| **numpy** | any | Numerical computation |
| **pandas** | any | Data handling |
| **scikit-learn** | any | StandardScaler |

scDEAL custom modules (`models.py`, `DaNN/`) are bundled in `third_party/python/scDEAL`.

### Recommended conda environment

```bash
conda env create -f skills/_toolkit/environment.yml   # tiphd-torch env
conda activate tiphd-torch
```

Or manually:

```bash
conda create -y -n tiphd-torch -c pytorch -c conda-forge \
  python=3.9 pytorch cpuonly numpy pandas scikit-learn
conda activate tiphd-torch
pip install anndata
```

## Prerequisites

1. **Environment** -- run `setup-env` skill.
2. **Data** -- run `setup-data` skill.
3. **R** -- evaluation layer requires R >= 4.4 with Seurat >= 5.0, survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scDEAL/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scdeal

# Full benchmark
python batch_run.py
# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
python run_scDEAL_pair.py AML GSE116256 TCGA
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scDEAL_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scDEAL_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Model Architecture

The model uses a deep feature extractor with MMD (Maximum Mean Discrepancy) loss to align bulk and single-cell feature distributions. Training uses Adam optimizer with adversarial training between the predictor and domain discriminator.

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 5-15 min |
| CRC | GSE132465 | ~47K | 10-20 min |
| HCC | GSE149614 | ~31K | 5-15 min |
| LUAD | GSE127465 | ~24K | 5-15 min |
| **GC** | **GSE183904** | **~137K** | **45-120 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). scDEAL's MMD domain adaptation on ~137K GC cells will exceed this. Run GC pairs directly or increase timeout:

```bash
python run_scDEAL_pair.py GC GSE183904 GSETCGA
```

## Notes

- scDEAL reads single-cell data from `.h5ad` files.
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `IndexError: index out of bounds` in `searchsorted` | Gene names not sorted; handled by `pd.Index.get_indexer()` in current version |
| `ModuleNotFoundError: torch` | Activate tiphd-torch env |
| `TIMEOUT (>7200s)` | Expected for GC pairs; see "Dataset Size and Timeout" section |
| `FileNotFoundError: .h5ad` | Run `setup-data` |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
