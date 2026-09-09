---
name: scstar2
description: Run the scSTAR2 benchmark pipeline. scSTAR2 is an A-class binary classification algorithm in R using OPLS-DA. It trains a PLS-DA model on bulk data (Good vs Poor prognosis), projects single-cell data onto the model, and labels cells above/below the 20th percentile as Rank-/Rank+. Use when the user says "用scSTAR2跑评测", "scSTAR2评测", "run scSTAR2 evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scSTAR2.
---

# scSTAR2 Benchmark Skill

scSTAR2 uses OPLS-DA (Orthogonal Projections to Latent Structures Discriminant Analysis) to transfer bulk phenotype labels to single cells. A PLS-DA model trained on bulk samples (Good/Poor prognosis) projects single-cell transcriptomes onto a discriminant axis; cells in the top/bottom 20% are labeled Rank-/Rank+.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (1000 Cox PH fits with shuffled Rank_Label) |
| plus-only | No |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds`) |
| **Seurat** | >= 5.0.0 | Reading single-cell data |
| **Matrix** | >= 1.6.4 | Sparse matrix support |
| **pls** | any | PLS/OPLS regression |
| **MASS** | any | LDA/statistical utilities |
| **survival** | any | Cox PH regression (evaluation) |

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table r-pls r-mass
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run `setup-env` skill.
2. **Data** -- run `setup-data` skill to download TiPhD test data.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scSTAR2/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scstar2

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single cancer
Rscript evaluate.R LUAD

# Single pair
Rscript run_scSTAR2_pair.R CRC GSE144735 GSE14333
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scSTAR2_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scSTAR2_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 2-5 min |
| CRC | GSE132465 | ~47K | 5-15 min |
| HCC | GSE149614 | ~31K | 2-5 min |
| LUAD | GSE127465 | ~24K | 2-5 min |
| **GC** | **GSE183904** | **~137K** | **20-60 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). GC pairs may exceed this due to OPLS projection on ~137K cells. Run directly or increase timeout:

```bash
Rscript run_scSTAR2_pair.R GC GSE183904 GSETCGA
```

## Notes

- scSTAR2 uses the top 3000 variable genes from single-cell data for the PLS-DA model.
- Bulk binary labels are derived from survival status (1=Dead/Rank+, 0=Alive/Rank-).
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `replacement has 1 row, data has 0` | Empty cell type subset after projection; expected when no cells map to gold-standard types |
| `TIMEOUT (>7200s)` | Expected for GC pairs; see "Dataset Size and Timeout" section. Run directly or set `timeout=7200` |
| `'RNA' is not an assay` | Seurat >= 5.0 required with R >= 4.4 |
| `could not find function "plsr"` | Install pls: `install.packages("pls")` |

## Scope

Does not modify algorithm parameters, install packages, or download data.
