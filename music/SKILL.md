---
name: music
description: Run the MuSiC deconvolution benchmark pipeline. MuSiC is a B-class (deconvolution) algorithm in R using weighted non-negative least squares (NNLS). Sample-level Cox permutation (1000 iterations) and gold-standard comparison run in evaluate.R. Use when the user says "用MuSiC跑评测", "MuSiC评测", "run MuSiC evaluation", "run MuSiC", or wants to execute the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics summary) for MuSiC. Triggers on "跑评测" / "运行MuSiC" / "deconvolution benchmark".
---

# MuSiC Benchmark Skill

MuSiC is a **B-class deconvolution** algorithm that estimates cell type proportions in bulk RNA-seq samples using single-cell reference data with weighted NNLS. This skill runs the full benchmark: per-pair deconvolution -> Cox permutation testing -> gold-standard comparison -> metrics.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | B (deconvolution) |
| Language | R |
| Output | Sample x cell-type proportion matrix |
| Permutation | Sample-level (1000 Cox PH fits with shuffled survival) |
| plus-only | No (outputs both positive and negative significant labels) |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds` requires R >= 4.4) |
| **Seurat** | >= 5.0.0 | Reading single-cell `.rds` files (Assay5 format) |
| **Matrix** | >= 1.6.4 | Sparse matrix support (required by SeuratObject v5) |
| **SingleCellExperiment** | any | Building SCE object for `music_prop()` |
| **nnls** | any | Non-negative least squares solver |
| **survival** | any | Cox PH regression (evaluation layer) |
| **dplyr** | any | Data manipulation |

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-nnls r-survival r-dplyr r-data.table \
  bioconductor-singlecellexperiment
conda activate tiphd-r
```

> **Important**: The benchmark data uses Seurat v5 Assay5 objects. R 4.1 + Seurat 4 cannot read them correctly. You must use R >= 4.4 with Seurat >= 5.0.

## Prerequisites

Before running this skill, ensure:

1. **Environment ready** -- run the `setup-env` skill to install R packages and conda environments.
2. **Data downloaded** -- run the `setup-data` skill to download the TiPhD test dataset to `<cwd>/data/` from ModelScope.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output results | `<cwd>/results/MuSiC/<cancer>/` | `TIPHD_OUT_DIR` |
| Rscript binary | `Rscript` (from PATH) | `TIPHD_RSCRIPT` |

The shared configuration in `_toolkit/config.R` resolves all paths, preferring environment variables and falling back to `<cwd>/data` and `<cwd>/results`.

## Running the Pipeline

### Full benchmark (all 5 cancers, 44 pairs)

```bash
cd skills/music
python batch_run.py
Rscript evaluate.R
```

### Single cancer

```bash
cd skills/music
python batch_run.py            # runs all pairs (skips completed)
Rscript evaluate.R CRC         # evaluate one cancer
```

### Single pair

```bash
Rscript run_MuSiC_pair.R <cancer> <sc_name> <bulk_name>
```

### Re-evaluate only (batch already done)

```bash
Rscript evaluate.R [cancer]
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `MuSiC_<bulk>_<sc>_proportions.csv` | Per-pair proportion matrix (samples x cell types) |
| `MuSiC_<bulk>_<sc>.csv` | Per-cell-type permutation results (Obs_Coef, P_val, Rank_Label) |
| `MuSiC_<cancer>_final_metrics.csv` | Per-cancer summary: Precision, Coverage, False_Rate |

## Dataset Size and Timeout

The single-cell datasets vary in cell count by cancer type:

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 2-5 min |
| CRC | GSE132465 | ~47K | 5-10 min |
| HCC | GSE149614 | ~31K | 2-5 min |
| LUAD | GSE127465 | ~24K | 2-5 min |
| **GC** | **GSE183904** | **~137K** | **30-90 min** |

`batch_run.py` enforces a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var) per pair. GC pairs with large bulk datasets (e.g., GSE15459 with 192 samples, GSE26253 with 432, GSE66229 with 300, GSETCGA with 367) will likely exceed this limit.

**For GC pairs, run directly without the batch timeout:**

```bash
# Run a single GC pair directly (no timeout)
Rscript run_MuSiC_pair.R GC GSE183904 GSETCGA

# Or increase the timeout in batch_run.py:
# config.run_batch(..., timeout=7200)  # 2 hours
```

## Cancer-Specific Notes

- **AML bulk names**: `TCGA`, `wave12`, `wave34` (not GSE IDs).
- **Gold-standard file**: `gold_standard_all_LUDA.csv` (note: **LUDA**, not LUAD -- original dataset spelling).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `'RNA' is not an assay` | Seurat version too old; upgrade to Seurat >= 5.0 with R >= 4.4 |
| `Negative entry appears!` | Bulk expression contains negative values after log-normalization; `run_MuSiC_pair.R` clips these to 0 automatically |
| `Too few common genes: 0` | Gene names mismatch between scRNA and bulk; check that both use the same gene ID format |
| `TIMEOUT (>7200s)` | Expected for GC pairs (~137K cells); see "Dataset Size and Timeout" section above. Run directly: `Rscript run_MuSiC_pair.R GC GSE183904 <bulk>` or set `timeout=7200` in `batch_run.py` |
| `library(X) not found` | Run `setup-env` skill to install missing R packages |
| Data file not found | Run `setup-data` skill to download the dataset |

## Scope

This skill does **not**:
- Modify algorithm parameters, formulas, or permutation counts
- Install R/Python packages (use `setup-env`)
- Download data (use `setup-data`)
