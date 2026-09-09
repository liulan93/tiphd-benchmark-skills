---
name: tirank
description: Run the TiRank benchmark pipeline. TiRank is an A-class Cox algorithm in Python that uses relative-expression-ordering (REO) gene-pair features, a neural encoder (MLP/Transformer/DenseNet), and MMD domain adaptation to transfer bulk survival labels to single cells. Requires the tiphd-stats conda environment (Python 3.9 + lifelines/optuna/leidenalg/python-igraph). Use when the user says "用TiRank跑评测", "TiRank评测", "run TiRank evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for TiRank.
---

# TiRank Benchmark Skill

TiRank aligns bulk and single-cell data using relative-expression-ordering (REO) gene-pair features and a neural network trained with a multi-task loss (Cox survival / classification / regression) plus MMD domain adaptation. It predicts a per-cell risk score and labels cells Rank+/Rank-/Background.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | A (cell-level) |
| Language | Python |
| Output | Cell-level labels: Rank+ / Rank- / Background |
| Permutation | Cell-level (Cox PH) in evaluate.R |
| plus-only | No |
| GPU | Auto-detected (Optuna search + training run on CUDA if available) |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | 3.9 (recommended) | Runtime |
| **torch / timm** | any | Neural encoder (MLP/Transformer/DenseNet) |
| **scanpy / anndata** | any | Single-cell preprocessing and clustering |
| **lifelines** | any | Cox PH |
| **optuna** | any | Hyperparameter search |
| **leidenalg / python-igraph** | any | Leiden clustering |
| **imbalanced-learn** | any | SMOTE / resampling |
| **numpy<2 / pandas** | <2.0 | Version-pinned for compatibility |

The TiRank source is bundled in `third_party/python/TiRank`; `batch_run.py` adds it to `PYTHONPATH` automatically (no pip install needed).

> **Environment isolation**: TiRank requires `numpy<2, pandas<2` plus lifelines/optuna/leidenalg, which can conflict with the deep-learning environment. Use the dedicated **tiphd-stats** conda environment.

### Recommended conda environment

```bash
conda env create -f skills/_toolkit/environment-tirank.yml   # tiphd-stats env
conda activate tiphd-stats
pip install timm imbalanced-learn   # if not pulled in by the yml
```

Or manually:

```bash
conda create -y -n tiphd-stats -c conda-forge -c pytorch \
  python=3.9 pytorch cpuonly numpy<2 "pandas<2" scanpy anndata \
  lifelines optuna leidenalg python-igraph magic-impute
conda activate tiphd-stats
pip install timm imbalanced-learn
```

## Prerequisites

1. **Environment** -- run the `setup-env` skill (creates the `tiphd-stats` env).
2. **Data** -- run the `setup-data` skill.
3. **R** -- the evaluation layer (`evaluate.R`) needs R >= 4.4 with Seurat >= 5.0 and survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/TiRank/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/tirank

# Full benchmark with the tiphd-stats interpreter
conda activate tiphd-stats
python batch_run.py
# Evaluation (R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
python run_TiRank_pair.py <cancer> <sc_name> <bulk_name>
```

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `TiRank_<bulk>_<sc>.csv` | Per-cell Rank_Label with cell type |
| `TiRank_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Notes

- TiRank performs preprocessing (QC, normalization, Leiden clustering), REO gene-pair extraction, hyperparameter search (`n_trials=10` via Optuna), and prediction per pair.
- On CPU, hyperparameter search + training is slow; the official recommendation is a GPU. The run script auto-detects CUDA.
- The run script uses a relaxed gene-pair threshold (`p_value_threshold=0.2`, `top_var_genes=3000`) so that extraction succeeds on datasets with weak bulk-survival signal. If you still hit `A set of genes is empty` / `0 Risk genes`, loosen the p-value threshold further or increase `top_var_genes`.

## Dataset Size and Timeout

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs are the slowest. If a pair times out, run it directly:

```bash
python run_TiRank_pair.py GC GSE183904 GSETCGA
```

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: lifelines/optuna/leidenalg/igraph/timm/gseapy` | Install in tiphd-stats: `pip install lifelines optuna timm gseapy`; leidenalg/igraph via `conda install -c conda-forge python-igraph leidenalg` |
| `No module named 'tirank'` | The TiRank source is added to PYTHONPATH by batch_run.py; run via `batch_run.py` or set `PYTHONPATH` to `third_party/python/TiRank` |
| `A set of genes is empty` / `0 Risk genes` | Weak bulk-survival signal; loosen p-value threshold or increase `top_var_genes` |
| numpy/pandas version errors | Ensure `numpy<2, pandas<2` in tiphd-stats |
| `TIMEOUT (>7200s)` | CPU training is slow; run directly or use a GPU |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
