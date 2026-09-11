---
name: statescope
description: Run the Statescope benchmark pipeline. Statescope is a B-class (deconvolution) algorithm in Python using BLADE Bayesian latent-variable deconvolution to estimate cell-type proportions in bulk samples from a single-cell signature. Sample-level Cox permutation runs in the evaluate layer. Use when the user says "用Statescope跑评测", "Statescope评测", "run Statescope evaluation", or wants the 3-stage pipeline (batch_run.py -> evaluate.R -> metrics) for Statescope.
---

# Statescope Benchmark Skill

Statescope performs Bayesian deconvolution (BLADE) of bulk RNA-seq using a single-cell-derived expression signature. It estimates per-sample cell-type fractions and, in the evaluate layer, tests association with survival.

## Algorithm Type

| Property | Value |
|----------|-------|
| Class | B (deconvolution) |
| Language | Python |
| Output | Sample x cell-type proportion matrix |
| Permutation | Sample-level (Cox PH) in evaluate.R |
| plus-only | No |
| GPU | Auto-detected (BLADE runs on CUDA if available). GPU is much faster than CPU but does **not** make very-high-cell-type configurations cheap — see Runtime and scaling |

## Dependencies

### Python Packages

| Package | Min Version | Purpose |
|---------|------------|---------|
| **Python** | >= 3.10 | BLADE source uses PEP-604 type unions (`X \| None`) |
| **torch** | >= 1.12 | Neural-network components (CUDA build for GPU) |
| **numba / dill / joblib** | any | JIT + parallel execution |
| **scanpy / anndata** | any | Reading `.h5ad` + log1p preprocessing |
| **numpy / pandas / scikit-learn / statsmodels** | any | Data handling |
| **autogenes / deap / cachetools** | any | Signature gene selection |
| **seaborn / matplotlib / requests / tqdm** | any | Utilities |

The Statescope/BLADE source is bundled in `third_party/python/Statescope` and loaded via `sys.path` by the run script (no pip install required).

> **Python version**: the BLADE source uses `pd.DataFrame | None` syntax which requires **Python 3.10 or later**. Use a dedicated Python 3.10 environment.
>
> **GPU**: BLADE is implemented in PyTorch and auto-detects CUDA (`torch.cuda.is_available()`). On CPU the Nrep Bayesian deconvolution is very slow; a CUDA GPU speeds it up substantially. Do not set `CUDA_VISIBLE_DEVICES=""` if you want GPU use.

### Recommended conda environment

```bash
# CPU-only
conda create -y -n tiphd-py310 -c pytorch -c conda-forge \
  python=3.10 pytorch cpuonly numpy pandas scikit-learn scanpy anndata \
  numba dill joblib seaborn statsmodels
# GPU (CUDA) — install a matching CUDA build of pytorch instead of cpuonly
conda activate tiphd-py310
pip install autogenes deap cachetools requests tqdm
```

## Prerequisites

1. **Environment** -- create the Python 3.10 environment above (or via `setup-env`).
2. **Data** -- run the `setup-data` skill.
3. **R** -- the evaluation layer (`evaluate.R`) needs R >= 4.4 with Seurat >= 5.0 and survival.

## Input/Output Paths

| Purpose | Default Path | Env Variable |
|---------|-------------|--------------|
| Input data | `<cwd>/data/` | `TIPHD_DATA_DIR` |
| Output | `<cwd>/results/Statescope/<cancer>/` | `TIPHD_OUT_DIR` |

## Running the Pipeline

```bash
cd skills/statescope

# Full benchmark (use the py310 interpreter)
/path/to/tiphd-py310/bin/python batch_run.py

# Evaluation (requires R)
conda activate tiphd-r
Rscript evaluate.R

# Single pair
/path/to/tiphd-py310/bin/python run_Statescope_pair.py <cancer> <sc_name> <bulk_name>
```

### Environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `TIPHD_DATA_DIR` / `TIPHD_OUT_DIR` | `./data` / `./results` | Standard data/result roots |
| `STATESCOPE_NREP` | `10` | Number of BLADE random initializations (the `Nrep` ensemble size). `Njob` is set automatically (≤4 workers, or 1 when `Nrep=1`). **The default 10 is the benchmark-faithful setting — do not change it for real results.** Setting `STATESCOPE_NREP=1` is appropriate only as an end-to-end smoke test that the pipeline runs; such an output is single-initialization and must be labelled as smoke, never quoted as a benchmark number |

`IterMax` (100 EM iterations) and every other BLADE hyperparameter stay at
upstream defaults; there is intentionally no env knob for them.

## Data normalization note

- **Bulk**: BLADE expects library-size-corrected, linear, non-negative counts, while the TiPhD benchmark bulk matrices are log-normalized. The run script clips tiny negatives and applies `expm1()` (inverse of `log1p`) before passing bulk to BLADE.
- **Single-cell signature**: the run script applies `sc.pp.log1p()` explicitly before building the signature. BLADE's `CreateSignature` uses a heuristic (`max < 50` plus an integer-ratio check) that can misclassify a normalized-but-float matrix as raw counts and abort with `looks like raw counts`; the explicit log1p makes the input unambiguously log-scale.
- A line `N/M gene expression differ between bulk and signature beyond ±1.96·SD` is an informational diagnostic from signature construction, not an error.

## Output Files

| File Pattern | Description |
|-------------|-------------|
| `Statescope_<bulk>_<sc>_proportions.csv` | Sample x cell-type fraction matrix |
| `Statescope_<cancer>_final_metrics.csv` | Per-cancer Precision/Coverage/False_Rate |

## Runtime and scaling — read before launching

BLADE runs Bayesian EM with an ensemble of `Nrep` random initializations (default 10), up to `IterMax=100` EM rounds each, stopping early on |objective change| < 1e-3. What governs cost:

- Cost scales primarily with the **number of annotated cell types** in the signature, not with bulk sample count: the parameter tensor involved in each L-BFGS evaluation has size Nsample × Ngene × Ncell, and an internal covariance term does O(Ncell²) work per EM round. Datasets with many cell types are therefore disproportionately expensive, and per-round cost grows as fitting proceeds (EM rounds are highly uneven — the first rounds can be quick while later ones take much longer).
- The `Nrep` reps run as threads in one process and share one CUDA context; they do not give linear GPU speed-up. The many small per-cell-type kernels can leave GPU utilization well below saturation even on a CUDA build.
- BLADE's heavy functions run on GPU automatically when a CUDA build of PyTorch is available, and are substantially slower on CPU.

Recommendations:

1. For real results, keep `Nrep=10`. On high-cell-type datasets, run pairs directly (the 7200 s batch timeout will trip) and plan for long runtimes.
2. If the goal is only to verify the pipeline works end-to-end (data → AutoGeneS → BLADE → output), `STATESCOPE_NREP=1` is a legitimate smoke option because reps are random-initialization replicas rather than a model hyperparameter — but label such output as a non-default smoke run and never quote it as a benchmark number. Do not truncate `IterMax` for convenience: that changes model fitting.
3. Prefer a CUDA build of PyTorch.

## Progress observability and no checkpoint

- The **AutoGeneS** marker-selection stage prints little or nothing while it runs; the log then shows the selected marker count and `Using GPU for computation` once BLADE starts. A long quiet period at this stage is normal.
- BLADE progress output is bare `print(i)` lines — one integer per EM round, with no label or timestamp, interleaved across reps. Each rep's sequence starts at `1`; the per-rep maximum is `IterMax-1 = 99`. Pair logs are opened in append mode by the batch driver, so numbers from an earlier killed attempt may precede the current run — judge progress from the tail of the log, not from the total number count.
- A long silence after the first iteration numbers means an `Optimize()` call is still running, not necessarily a hang — check process CPU time / GPU activity before restarting.
- **There is no checkpoint.** Killing the process discards everything (including the completed AutoGeneS stage); a restart reruns the whole pair from the beginning. Launch long pairs from a resilient session. The batch runner skips pairs whose output CSV already exists.

## Dataset Size and Timeout

The default per-pair timeout is **7200s** (configurable via `TIPHD_PAIR_TIMEOUT`); high-cell-type datasets at full defaults commonly exceed it, so run such pairs directly:

```bash
/path/to/tiphd-py310/bin/python run_Statescope_pair.py CRC GSE144735 GSE14333
# smoke only — NOT a benchmark result:
STATESCOPE_NREP=1 /path/to/tiphd-py310/bin/python run_Statescope_pair.py CRC GSE144735 GSE14333
```

## Troubleshooting

| Symptom | Cause / Fix |
|-------|------------|
| `TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'` | Python < 3.10; use the Python 3.10 environment |
| `ModuleNotFoundError: autogenes/deap/cachetools` | `pip install autogenes deap cachetools` in tiphd-py310 |
| `adata.X looks like raw counts; please normalise & log-transform` | The signature must be log-normalized. The run script applies `log1p` explicitly; if you bypass it, log-transform the `.h5ad` yourself |
| `Bulk contains negative values` | Bulk must be non-negative; the run script clips and linearizes (expm1) automatically |
| No log output for a long time after the AutoGeneS stage starts | Normal — marker selection is quiet until it finishes; not a hang |
| Only bare numbers (`1`, `2`, …) in the log, then a long silence | Normal BLADE progress output; EM rounds have uneven, growing cost. Check GPU activity and process CPU time to distinguish work from hang |
| `TIMEOUT (>7200s)` | Expected for high-cell-type datasets at full defaults; run the pair directly (see Runtime and scaling). Use `STATESCOPE_NREP=1` only for an explicitly-labelled smoke run — never silently for benchmark numbers |
| `Using GPU for computation` never appears / runs very slowly | Either no CUDA build of PyTorch in the env or `CUDA_VISIBLE_DEVICES` is hidden; check `torch.cuda.is_available()` |
| Pair restarted from AutoGeneS after a kill | Expected — no checkpoint exists; only completed pairs (output CSV present) are skipped |

## Scope

Does not modify algorithm parameters, install packages, or download data. The only sanctioned deviations from upstream defaults are `STATESCOPE_NREP` (default 10 preserved) for labelled smoke runs; everything else stays at defaults.
