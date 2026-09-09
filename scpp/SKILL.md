---
name: scpp
description: Run the ScPP benchmark pipeline. ScPP is an A-class algorithm in R that identifies phenotype-associated cell subpopulations using univariate Cox PH on bulk data to select prognostic marker genes, then AUCell to score single cells. Cells above/below the 20th percentile are labeled Rank+/Rank-. Use when the user says "用ScPP跑评测", "ScPP评测", "run ScPP evaluation", "run ScPP", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for ScPP.
---

# ScPP Benchmark Skill

ScPP (Single-Cell Phenotype Predictor) identifies cell subpopulations associated with survival outcomes by: (1) running univariate Cox PH on each bulk gene to select prognostic markers, (2) using AUCell to score single cells by marker gene enrichment, (3) labeling top/bottom 20% cells as Rank+/Rank-.

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
| **AUCell** | any | Gene set enrichment scoring |
| **GSEABase** | any | Gene set infrastructure (AUCell dependency) |
| **survival** | any | Cox PH regression |

### Recommended conda environment

```bash
conda create -y -n tiphd-r -c conda-forge -c bioconda \
  r-base=4.4.2 r-seurat r-matrix r-survival r-dplyr r-data.table \
  bioconductor-aucell bioconductor-gseabase
conda activate tiphd-r
```

## Prerequisites

1. **Environment** -- run `setup-env` skill.
2. **Data** -- run `setup-data` skill to download TiPhD test data to `<cwd>/data/`.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/ScPP/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

### Full benchmark

```bash
cd skills/scpp
python batch_run.py
Rscript evaluate.R
```

### Single cancer / pair

```bash
Rscript evaluate.R CRC
Rscript run_ScPP_pair.R AML GSE116256 TCGA
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `ScPP_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type annotation |
| `ScPP_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Dataset Size and Timeout

| Cancer | scRNA dataset | Cells | Typical pair runtime |
|--------|--------------|-------|---------------------|
| AML | GSE116256 | ~36K | 2-5 min |
| CRC | GSE132465 | ~47K | 5-15 min |
| HCC | GSE149614 | ~31K | 2-5 min |
| LUAD | GSE127465 | ~24K | 2-5 min |
| **GC** | **GSE183904** | **~137K** | **20-60 min** |

`batch_run.py` uses a default **7200s (2 hour) timeout** (configurable via `TIPHD_PAIR_TIMEOUT` env var). GC pairs may exceed this due to AUCell scoring on ~137K cells. Run GC pairs directly or increase the timeout:

```bash
Rscript run_ScPP_pair.R GC GSE183904 GSETCGA
```

## Notes

- ScPP selects marker genes with Cox p-value < 0.05 from bulk data; if fewer than 5 genes pass, the pair is skipped ("marker genes insufficient").
- AML bulk names: `TCGA`, `wave12`, `wave34`.
- Gold-standard file: `gold_standard_all_LUDA.csv` (note: LUDA, not LUAD).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `could not find function "AUCell_buildRankings"` | Install AUCell: `BiocManager::install("AUCell")` |
| `marker genes insufficient (< 5)` | No significant prognostic genes in bulk; expected for some datasets |
| `TIMEOUT (>7200s)` | Expected for GC pairs; see "Dataset Size and Timeout" section. Run directly or set `timeout=7200` |
| `'RNA' is not an assay` | Seurat >= 5.0 required with R >= 4.4 |

## Scope

Does not modify algorithm parameters, install packages, or download data.
