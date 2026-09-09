---
name: scipac
description: Run the SCIPAC benchmark pipeline. SCIPAC is an A-class Cox-based algorithm in R that identifies phenotype-associated cells by iteratively partitioning cells and fitting Cox models to find subpopulations associated with survival. Use when the user says "用SCIPAC跑评测", "SCIPAC评测", "run SCIPAC evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for SCIPAC.
---

# SCIPAC Benchmark Skill

SCIPAC (Single-Cell Identification of Phenotype-Associated Cells) uses an iterative partitioning approach: it recursively splits cell populations and fits Cox PH models to identify subpopulations significantly associated with bulk survival outcomes.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R |
| Output | Cell-level labels: Rank+ / Background |
| Permutation | Cell-level (1000 Cox PH fits) |
| plus-only | No |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5) |
| **Seurat** | >= 5.0.0 | Reading single-cell data |
| **Matrix** | >= 1.6.4 | Sparse matrix |
| **SCIPAC** | custom (in `third_party/R/SCIPAC`) | The algorithm itself |
| **glmnet** | any | Elastic-net regularization |
| **ordinalNet** | any | Ordinal regression |
| **survival** | any | Cox PH regression |
| **parallel** | any | Parallel computation |

SCIPAC is installed from source: `third_party/R/SCIPAC`.

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table r-glmnet
conda activate tiphd-r
# Install SCIPAC from source
Rscript -e 'install.packages("skills/_toolkit/third_party/R/SCIPAC", repos=NULL, type="source")'
```

Or use the `setup-env` skill which handles all custom R package installation.

## Prerequisites

1. **Environment** -- run `setup-env` skill (installs CRAN, Bioconductor, and 5 custom R packages).
2. **Data** -- run `setup-data` skill.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/SCIPAC/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scipac

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single cancer / pair
Rscript evaluate.R CRC
Rscript run_SCIPAC_pair.R AML GSE116256 TCGA
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `SCIPAC_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `SCIPAC_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 3-10 min |
| CRC | GSE132465 | ~47K | 5-20 min |
| HCC | GSE149614 | ~31K | 3-10 min |
| LUAD | GSE127465 | ~24K | 3-10 min |
| **GC** | **GSE183904** | **~137K** | **30-90 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). SCIPAC's iterative partitioning on ~137K GC cells with glmnet may exceed this. Run directly or increase timeout:

```bash
Rscript run_SCIPAC_pair.R GC GSE183904 GSETCGA
```

## Notes

- SCIPAC depends on `glmnet` and `ordinalNet` for its penalized Cox models.
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `there is no package called 'SCIPAC'` | Install from source: `install.packages("third_party/R/SCIPAC", repos=NULL)` |
| `there is no package called 'glmnet'` | Install: `install.packages("glmnet")` |
| `'RNA' is not an assay` | Seurat >= 5.0 required with R >= 4.4 |

## Scope

Does not modify algorithm parameters, install packages, or download data.
