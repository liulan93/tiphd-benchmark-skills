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
| GPU | Opt-in: set `TIRANK_GPU=1` (and use an env with a CUDA build of PyTorch). Without it the run script hides CUDA and runs CPU. |

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

## Input contract — bulk clinical table MUST be exactly two columns

TiRank picks clinical columns **by position**, not by name: the gene-pair
extractor uses `clinical_data.columns[0:2]` while the Cox training split uses
`iloc[:, -2:]`. The pipeline therefore requires the bulk clinical CSV to be
exactly `[survival_time, event]` in that order. Several raw clinical tables in
the wild (and in the TiPhD data: HCC/CRC/GC) ship extra columns (age, gender,
stage, …), which silently produce wrong column pairs.

The run script performs **input alignment only**: it trims the clinical table
to the `tcol`/`scol` configured per cancer in `_toolkit/config.py`. No values
or model hyperparameters are changed. When running on custom data, either
pre-trim the clinical CSV to two columns yourself or register the time/event
column names in `config.py`.

A second trap: the upstream TiRank extractor wraps every per-gene Cox fit in
`try/except Exception: continue`, so any Cox failure (wrong columns, bad event
encoding, dtype issues) is swallowed and only surfaces much later as
`There are 0 Risk genes and 0 Protective genes`. If you see that error,
**first verify the clinical table** (manually reproduce one
`lifelines.CoxPHFitter.fit` to see the real exception) — do not assume weak
signal.

## Running the Pipeline

```bash
cd skills/tirank

# Full benchmark with the tiphd-stats interpreter (CPU)
conda activate tiphd-stats
python batch_run.py

# GPU: the env must contain a CUDA build of PyTorch (the default
# environment-tirank.yml installs the CPU-only build), then opt in:
TIRANK_GPU=1 python batch_run.py

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

- Per pair TiRank runs preprocessing (QC, normalization, Leiden clustering), REO gene-pair extraction, Optuna hyperparameter search (`n_trials=10`, 100 epochs each), and prediction. Runtime scales with cell count and varies widely across datasets.
- GPU is **opt-in**: the run script forces `CUDA_VISIBLE_DEVICES=""` unless `TIRANK_GPU=1` is set, and the conda env must actually contain a CUDA build of PyTorch (the shipped `environment-tirank.yml` installs `cpuonly`). GPU mainly accelerates the Optuna search and training.
- The run script uses a relaxed gene-pair threshold (`p_value_threshold=0.2`, `top_var_genes=3000`) so that extraction succeeds on datasets with weak bulk-survival signal. This is the skill's preset, not an upstream default — do not loosen it further without explicit user approval (it changes model inputs).
- Dense `.h5ad` inputs are streamed in backed mode and converted to CSR at load time. This changes memory layout only, not values, so very large dense matrices do not have to be expanded in RAM.
- Intermediate artifacts go to a per-pair temp directory that is deleted on exit. Do not remove `$TMPDIR/TiRank_*` while a pair is running — the prediction stage reads back a pickle written during tuning.

## Dataset Size and Timeout

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`). GC (~137K cells) pairs are the slowest. If a pair times out, run it directly (preferably with GPU enabled):

```bash
TIRANK_GPU=1 python run_TiRank_pair.py GC GSE183904 GSETCGA
```

There is no checkpoint: killing a pair discards all tuning/training progress and the pair restarts from scratch (finished pairs are skipped via their output CSV, so batch reruns only redo missing ones).

## Troubleshooting

| Error | Cause / Fix |
|-------|------------|
| `ModuleNotFoundError: lifelines/optuna/leidenalg/igraph/timm/gseapy` | Install in the TiRank env: `pip install lifelines optuna timm gseapy`; leidenalg/igraph via `conda install -c conda-forge python-igraph leidenalg` |
| `No module named 'tirank'` | The TiRank source is added to `sys.path` by the run script and needs no install; run via `batch_run.py`/`run_TiRank_pair.py` from the skill folder |
| `There are 0 Risk genes and 0 Protective genes` / `A set of genes is empty` | **First suspect the clinical table, not weak signal.** Cox fits failing inside the extractor's `try/except Exception: continue` are reported as "no significant genes". Checks, in order: (1) the bulk clinical CSV must be exactly `[time, event]` — the run script auto-trims the TiPhD files via `config.py` tcol/scol; custom data needs the same; (2) reproduce one `CoxPHFitter().fit(df, duration_col=…, event_col=…)` manually to see the swallowed exception; (3) only after the input checks out, discuss thresholds with the user (changing `p_value_threshold`/`top_var_genes` is a model-input change, not a free fix) |
| `TypeError`/`AttributeError` from pandas inside `GPextractor` (e.g. `DataFrame.append`) | pandas ≥ 2 removed `DataFrame.append`; the bundled TiRank source is patched to list-collect rows. On an unpatched TiRank copy, either stay on `pandas<2` or apply the equivalent list-collection fix |
| numpy/pandas version errors | The shipped env pins `numpy<2, pandas<2`; the patched bundled source is compatible with pandas 2.x — pick one env and keep it consistent |
| GPU not used despite `TIRANK_GPU=1` | The env's PyTorch must be a CUDA build (`python -c "import torch;print(torch.cuda.is_available())"`); the default `environment-tirank.yml` installs `cpuonly` — replace it with a CUDA-matching pytorch build |
| `TIMEOUT (>7200s)` | CPU training is slow; run the pair directly or use a GPU (`TIRANK_GPU=1` + CUDA pytorch in the env) |
| Process dies right after reading the `.h5ad`, no traceback | Likely an OS-level OOM kill from dense `X` expansion. The run script already streams dense HDF5 as CSR; check available RAM and that you are running the bundled run script |

## Scope

Does not modify model architecture, hyperparameters, install packages, or download data.
