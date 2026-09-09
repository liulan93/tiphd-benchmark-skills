---
name: scper
description: Run the scPER benchmark pipeline. scPER is a B-class (deconvolution) algorithm in Python that imputes expression with MAGIC, trains an adversarial denoising autoencoder (ADAE, TensorFlow/Keras), and estimates cell-type proportions with xgboost. It shells out to R for MAGIC imputation. Sample-level Cox permutation runs in the evaluate layer. Use when the user says "用scPER跑评测", "scPER评测", "run scPER evaluation", or wants the 3-stage pipeline for scPER.
---

# scPER Benchmark Skill

scPER performs bulk deconvolution in three stages: (1) MAGIC imputation + highly-variable-gene selection, (2) an adversarial denoising autoencoder (ADAE, implemented in TensorFlow/Keras), and (3) proportion estimation with xgboost. The run script orchestrates both Python and R (Rmagic) steps.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | B (deconvolution) |
| Language | Python (PyTorch/TensorFlow) + R (MAGIC) |
| Output | Sample x cell-type proportion matrix |
| Permutation | Sample-level (Cox PH) in evaluate.R |
| plus-only | No |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.9 | Runtime |
| **tensorflow / keras** | ~2.10 | ADAE model (TF 2.10 works on CPU) |
| **xgboost** | any | Proportion regression |
| **magic-impute** | any | Expression imputation (also via Rmagic) |
| **anndata / numpy / pandas** | any | Data handling |

### R Packages

| Package | Purpose |
|---------|---------|
| **Seurat / Matrix** | Single-cell object handling |
| **Rmagic / magic** | MAGIC imputation called via subprocess |

The scPER source is bundled in `third_party/python/scPER`.

> **Important**: scPER needs **both** the `tiphd-torch` Python environment and the R 4.4 + Seurat 5 environment because the MAGIC step is run via an R subprocess. TensorFlow 2.10 is recommended (newer versions changed the Keras API used by ADAE).

### Recommended conda environment

```bash
conda create -y -n tiphd-torch -c pytorch -c conda-forge \
  python=3.9 pytorch cpuonly numpy pandas scikit-learn anndata
conda activate tiphd-torch
pip install "tensorflow-cpu==2.10.0" xgboost magic-impute
# R side (separate R environment): Seurat, Matrix, plus Rmagic
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (sets up both Python and R).
2. **Data** -- run the `setup-data` skill.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scPER/<cancer>/` | `TIPHD_OUT_DIR` |
| Rscript binary | `Rscript` (from PATH) | `TIPHD_RSCRIPT` |

## Running the Pipeline

```bash
cd skills/scper

# Full benchmark (Python env active; Rscript on PATH)
python batch_run.py
# Evaluation (R)
Rscript evaluate.R

# Single pair
python run_scPER_pair.py <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scPER_<bulk>_<sc>_proportions.csv` | Sample x cell-type fraction matrix |
| `scPER_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Notes

- The ADAE uses a `Batch` column for batch correction; the run script uses patient/sample IDs as a batch proxy for single-dataset inputs.
- scPER writes intermediate files to a temporary working directory and resolves the original scripts' hard-coded paths automatically.
- GPU is not required (`CUDA_VISIBLE_DEVICES=""` forces CPU).

## Dataset Size and Timeout

scPER runs MAGIC imputation plus ADAE training, which is heavy for large datasets. The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs may need direct invocation:

```bash
python run_scPER_pair.py GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: tensorflow / xgboost / magic` | Install in tiphd-torch: `pip install tensorflow-cpu==2.10.0 xgboost magic-impute` |
| Rmagic / `Rscript` errors | Ensure R >= 4.4 with Seurat and Rmagic; set `TIPHD_RSCRIPT` |
| Keras API errors | Use TensorFlow 2.10.x; TF 2.16+ reorganized Keras namespaces |
| `TIMEOUT (>7200s)` | Heavy ADAE/MAGIC on large data; run the pair directly |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
