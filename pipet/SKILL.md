---
name: pipet
description: Run the PIPET benchmark pipeline. PIPET (Phenotypic Information Prediction) is an A-class binary classification algorithm in R that selects up/down-regulated markers from bulk Dead-vs-Alive log-fold-change and then classifies single cells by cosine distance with permutation-based significance. Use when the user says "用PIPET跑评测", "PIPET评测", "run PIPET evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for PIPET.
---

# PIPET Benchmark Skill

PIPET predicts phenotype-associated single cells from bulk-derived markers: it splits bulk samples by survival status (Dead/Alive), selects the top up-/down-regulated genes as markers, then scores each single cell against these markers using cosine distance with internal permutation testing to assign significance.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | R |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Internal (1000 permutations), plus cell-level in evaluate.R |
| plus-only | No |

## Dependencies

### R Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **R** | >= 4.4.0 | Runtime (Seurat v5 Assay5 `.rds`) |
| **Seurat** | >= 5.0.0 | Single-cell object |
| **Matrix** | >= 1.6.4 | Sparse matrices |
| **matrixStats** | any | Fast matrix statistics |
| **tibble / parallel / stats** | any | Core helpers |
| **PIPET** | custom (in `third_party/R/PIPET`) | The algorithm itself |

PIPET's optional visualization/marker-extraction helpers (ggplot2, cowplot, DESeq2, survminer, tidyverse) are **not** needed for the benchmark and are listed as `Suggests` so they do not block installation. PIPET is installed from source via `setup-env`.

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table r-matrixstats
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (installs the custom `PIPET` package; core imports only require Matrix/Seurat/matrixStats).
2. **Data** -- run the `setup-data` skill.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/PIPET/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/pipet

# Full benchmark
python batch_run.py
Rscript evaluate.R

# Single pair
Rscript run_PIPET_pair.R <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `PIPET_<bulk>_<sc>.csv` | Per-cell prediction, P-value, FDR, and Rank_Label |
| `PIPET_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs are the slowest. If a pair times out, run it directly:

```bash
Rscript run_PIPET_pair.R GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `lazy loading failed` / `there is no package called 'survminer'` | Installed from the full upstream package requiring visualization deps; the bundled version moves these to `Suggests`. Reinstall via `setup-env` |
| `there is no package called 'PIPET'` | Run `setup-env` to install from source |
| `'RNA' is not an assay` | Seurat >= 5.0 with R >= 4.4 required |

## Scope

Does not modify algorithm parameters, install packages, or download data.
