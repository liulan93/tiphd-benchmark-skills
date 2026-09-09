---
name: scissor
description: Run the Scissor benchmark pipeline. Scissor is an A-class algorithm in R that links bulk phenotype/survival data to single cells by leveraging a Seurat cell-similarity network and Cox regression, outputting Scissor+ / Scissor- / Background cells. Use when the user says "用Scissor跑评测", "Scissor评测", "run Scissor evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for Scissor.
---

# Scissor Benchmark Skill

Scissor identifies single cells most associated with bulk phenotypes (e.g. poor/good survival) by combining a Seurat-derived cell-neighbor network with a Cox/binary model over bulk samples. It labels cells as `Scissor+` (risk-associated), `Scissor-` (protective), or `Background`.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R (Rcpp) |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (Cox PH with shuffled labels) |
| plus-only | No |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds`) |
| **Seurat** | >= 5.0.0 | Single-cell object, neighbor network |
| **Matrix** | >= 1.6.4 | Sparse matrices |
| **preprocessCore** | any | Quantile normalization (Bioconductor) |
| **survival** | any | Cox PH regression |
| **Rcpp / RcppEigen** | any | Compiled C++ core |
| **Scissor** | custom (in `third_party/R/Scissor`) | The algorithm itself |

Scissor is installed from source via `setup-env` (includes Rcpp compilation; Windows requires Rtools44).

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table \
  r-rcpp r-rcppeigen bioconductor-preprocesscore
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill.
2. **Data** -- run the `setup-data` skill.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/Scissor/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scissor

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single pair
Rscript run_Scissor_pair.R <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `Scissor_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `Scissor_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

Scissor computes a dense cell-similarity network, so memory grows with cell count. The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs are the slowest; run them directly if they time out:

```bash
Rscript run_Scissor_pair.R GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `pthread_create() is 22` | Too many OpenMP threads; `setup-env` installs preprocessCore with `--disable-threading`. You can also set `OMP_NUM_THREADS=1` |
| `fatal error: Scissor.h: No such file or directory` | Missing placeholder header; the bundled source includes `inst/include/Scissor.h`. Reinstall from source |
| `Rcpp compilation failed` | Install Rtools44 (Windows) or Xcode CLT (macOS) |
| `there is no package called 'Scissor'` | Run `setup-env` to build the custom package |
| `'RNA' is not an assay` | Seurat >= 5.0 with R >= 4.4 required |

## Scope

Does not modify algorithm parameters, install packages, or download data.
