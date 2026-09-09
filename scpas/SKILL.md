---
name: scpas
description: Run the scPAS benchmark pipeline. scPAS is an A-class Cox-based algorithm in R that identifies phenotype-associated subpopulations by combining bulk survival signatures with single-cell expression; it performs its own internal permutation testing (2000 iterations) with BH FDR correction. Cells are labeled Rank+/Rank-/Background. Use when the user says "用scPAS跑评测", "scPAS评测", "run scPAS evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scPAS.
---

# scPAS Benchmark Skill

scPAS (single-cell Phenotype-Associated Subpopulations) derives risk/protective feature signatures from bulk survival data, projects them onto single cells, and runs its own internal permutation testing (2000 iterations) with Benjamini-Hochberg FDR to call significant cells.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R (Rcpp) |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Internal (2000 iterations + BH FDR), plus cell-level in evaluate.R |
| plus-only | No |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds`) |
| **Seurat** | >= 5.0.0 | Single-cell object |
| **Matrix** | >= 1.6.4 | Sparse matrices |
| **preprocessCore** | any | Quantile normalization (Bioconductor) |
| **survival** | any | Cox PH regression |
| **Rcpp / RcppEigen** | any | Compiled C++ core |
| **scPAS** | custom (in `third_party/R/scPAS`) | The algorithm itself |

scPAS is installed from source via `setup-env` (includes Rcpp compilation; Windows requires Rtools44).

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table r-rcpp r-rcppeigen \
  bioconductor-preprocesscore
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill.
2. **Data** -- run the `setup-data` skill.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scPAS/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scpas

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single pair
Rscript run_scPAS_pair.R <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scPAS_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scPAS_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (GSE183904, ~137K cells) pairs are the slowest. If a pair times out, run it directly:

```bash
Rscript run_scPAS_pair.R GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `pthread_create() is 22` | Too many OpenMP threads; `setup-env` installs preprocessCore with `--disable-threading`. You can also set `OMP_NUM_THREADS=1` |
| `fatal error: scPAS.h: No such file or directory` | Missing placeholder header; the bundled source includes `inst/include/scPAS.h`. Reinstall from source |
| `Rcpp compilation failed` | Install Rtools44 (Windows) or Xcode CLT (macOS) |
| `there is no package called 'scPAS'` | Run `setup-env` to build the custom package |
| `'RNA' is not an assay` | Seurat >= 5.0 with R >= 4.4 required |

## Scope

Does not modify algorithm parameters, install packages, or download data.
