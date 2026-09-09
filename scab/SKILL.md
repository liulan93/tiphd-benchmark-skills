---
name: scab
description: Run the scAB benchmark pipeline. scAB is an A-class plus-only algorithm in R that identifies phenotype-associated cell subpopulations by adapting a bulk Cox regression to single-cell networks (single-cell adaptive boosting). It outputs only Rank+ (risk) / Background labels. Use when the user says "用scAB跑评测", "scAB评测", "run scAB evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for scAB.
---

# scAB Benchmark Skill

scAB (single-cell Adaptive Boosting) adapts a bulk survival association test to single cells using a cell-neighbor network (Seurat SNN graph). It iteratively selects K and identifies cells whose subpopulation profile is associated with survival. It is a **plus-only** algorithm: cells are labeled `Rank+` (risk-associated) or `Background` (no Rank- direction).

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R |
| Output | Cell-level labels: Rank+ / Background |
| Permutation | Cell-level (Cox PH with shuffled labels) |
| plus-only | Yes |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds`) |
| **Seurat** | >= 5.0.0 | Single-cell object, SNN graph |
| **Matrix** | >= 1.6.4 | Sparse matrices |
| **preprocessCore** | any | Quantile normalization (Bioconductor) |
| **survival** | any | Cox PH regression |
| **scAB** | custom (in `third_party/R/scAB`) | The algorithm itself |

scAB is installed from source via `setup-env`.

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table \
  bioconductor-preprocesscore
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (installs R packages including the custom `scAB` package).
2. **Data** -- run the `setup-data` skill to download the TiPhD test data.

> **Note**: scAB builds a dense cell-neighbor matrix and runs quantile normalization, which is memory- and thread-intensive. On high-core-count servers, set thread limits to avoid `pthread_create()` errors (see Troubleshooting).

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/scAB/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/scab

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single pair
Rscript run_scAB_pair.R <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `scAB_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `scAB_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Notes |
|--------|--------------|-------|-------|
| AML | GSE116256 | ~36K | Dense SNN matrix ~10 GiB |
| CRC | GSE132465 | ~47K | ~13 GiB |
| HCC | GSE149614 | ~31K | |
| LUAD | GSE127465 | ~24K | |
| **GC** | **GSE183904** | **~137K** | Very slow/memory-heavy |

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC pairs may be slow or hit memory limits; run them directly if needed:

```bash
Rscript run_scAB_pair.R GC GSE183904 GSE28541
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `pthread_create() is 22` | Too many OpenMP threads on a high-core machine; reinstall preprocessCore with `--disable-threading` (handled by `setup-env`) or set `OMP_NUM_THREADS=1` |
| `sparse->dense coercion: allocating vector of size ... GiB` | Large SNN matrix; needs ample RAM (>=32 GiB recommended for 36K+ cells) |
| `cannot allocate vector of ... Gb` | Insufficient memory; run smaller cancers first or on a higher-memory machine |
| `there is no package called 'scAB'` | Run `setup-env` to install the custom package from source |
| `'RNA' is not an assay` | Seurat >= 5.0 with R >= 4.4 required |

## Scope

Does not modify algorithm parameters, install packages, or download data.
