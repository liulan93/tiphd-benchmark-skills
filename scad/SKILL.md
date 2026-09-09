---
name: scad
description: Run the SCAD benchmark pipeline. SCAD is an A-class deep learning algorithm in Python using a VAE (Variational Autoencoder) feature extractor with an adversarial domain classifier and binary predictor. It transfers bulk survival labels to single cells via adversarial transfer learning. Cells in top/bottom 20% of predicted logits are labeled Rank-/Rank+. Use when the user says "用SCAD跑评测", "SCAD评测", "run SCAD evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for SCAD.
---

# SCAD Benchmark Skill

SCAD (Single-Cell Adversarial Domain adaptation) uses a neural network with three components:
- **FX**: Feature extractor (VAE-style encoder)
- **MTLP**: Binary predictor (survival Dead/Alive)
- **Discriminator**: Adversarial domain classifier (bulk vs single-cell, with gradient reversal)

The model trains on bulk data and transfers phenotype labels to single cells through domain-invariant features.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | Python (PyTorch) |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (1000 Cox PH fits) |
| plus-only | No |
| GPU | Optional (CPU works, slower) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.9 | Runtime |
| **torch** | >= 1.12 | Neural network framework |
| **anndata** | >= 0.8 | Reading `.h5ad` single-cell data |
| **numpy** | any | Numerical computation |
| **pandas** | any | Data handling |
| **scikit-learn** | any | StandardScaler, metrics |

SCAD custom modules (`SCADmodules.py`) are bundled in `third_party/python/SCAD`.

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

1. **Environment** -- run `setup-env` skill (creates `tiphd-torch` conda env).
2. **Data** -- run `setup-data` skill (downloads TiPhD test data).
3. **R** -- the evaluation layer (`evaluate.R`) requires R >= 4.4 with Seurat >= 5.0, survival, dplyr (for reading `.rds` gold-standard files).

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/SCAD/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

### Full benchmark

```bash
cd skills/scad
python batch_run.py
# Evaluation requires R:
conda activate tiphd-r
Rscript evaluate.R
```

### Single cancer / pair

```bash
Rscript evaluate.R CRC
python run_SCAD_pair.py AML GSE116256 TCGA
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `SCAD_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `SCAD_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Model Architecture

| Component | Layers |
|-----------|--------|
| FX (encoder) | input_dim -> 512 -> 256 |
| MTLP (predictor) | 256 -> 64 -> 1 (sigmoid) |
| Discriminator | 256 -> 64 -> 1 (sigmoid, gradient reversal) |

Training uses Adam optimizer with learning rate 1e-4, batch size 64, and 50 epochs.

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 5-15 min |
| CRC | GSE132465 | ~47K | 10-20 min |
| HCC | GSE149614 | ~31K | 5-15 min |
| LUAD | GSE127465 | ~24K | 5-15 min |
| **GC** | **GSE183904** | **~137K** | **45-120 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). SCAD trains a VAE + adversarial network on all cells, so GC pairs (~137K cells) will almost certainly exceed this. Run GC pairs directly or increase timeout:

```bash
python run_SCAD_pair.py GC GSE183904 GSETCGA
```

## Notes

- SCAD reads single-cell data from `.h5ad` files (not `.rds`).
- Input data is standardized (StandardScaler) before training.
- The model uses CPU by default; set `device="cuda"` if GPU is available.
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `IndexError: index out of bounds` in `searchsorted` | Gene names not sorted; handled by `pd.Index.get_indexer()` in current version |
| `ModuleNotFoundError: torch` | Activate tiphd-torch: `conda activate tiphd-torch` |
| `CUDA out of memory` | Set `device="cpu"` in `run_SCAD_pair.py` |
| `TIMEOUT (>7200s)` | Expected for GC pairs (~137K cells, VAE training); see "Dataset Size and Timeout" section |
| `FileNotFoundError: .h5ad` | Run `setup-data` to download data |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
