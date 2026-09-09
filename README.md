# TiPhD Benchmark Skills

A benchmark suite for evaluating **16 single-cell phenotype-association and bulk deconvolution algorithms** on a standardized TiPhD test dataset (5 cancers: AML, CRC, HCC, LUAD, GC). Each algorithm is wrapped as a self-contained skill implementing a common three-stage pipeline.

## What is here

| Path | Purpose |
|------|---------|
| `setup-env/` | One-shot installer for R + Python environments |
| `setup-data/` | Download the TiPhD test dataset from ModelScope |
| `run-benchmark/` | Dispatcher skill to run any algorithm end-to-end |
| `_toolkit/` | Shared config, install scripts, third-party packages |
| `<algorithm>/` | One skill per algorithm (`batch_run.py`, `evaluate.R`, `run_*_pair.*`) |

## Algorithms

### A-class — cell-level (output Rank+/Rank-/Background per cell)

| Skill | Language | Method |
|-------|----------|--------|
| `music` → MuSiC | R | Weighted NNLS deconvolution * |
| `scpp` → ScPP | R | Cox markers + AUCell |
| `scstar2` → scSTAR2 | R | OPLS-DA projection |
| `scad` → SCAD | Python | VAE + adversarial domain adaptation |
| `scipac` → SCIPAC | R | Iterative partitioning + Cox |
| `scdeal` → scDEAL | Python | VAE + DaNN domain adaptation |
| `scsurv` → scSurv | Python | VAE + Cox latent space |
| `sctrend` → scTREND | Python | VAE + DeepCOLOR + segmented risk model |
| `scab` → scAB | R | Single-cell adaptive boosting (plus-only) |
| `scissor` → Scissor | R | Seurat network + Cox |
| `scpas` → scPAS | R | Cox signature (internal 2000 permutations) |
| `pipet` → PIPET | R | Marker cosine-distance classification |
| `sidish` → SIDISH | Python | VAE + DeepCox (plus-only, GPU) |
| `tirank` → TiRank | Python | REO gene pairs + NN + MMD (GPU) |

### B-class — sample-level deconvolution (output per-sample cell-type proportions)

| Skill | Language | Method |
|-------|----------|--------|
| `scper` → scPER | Python+R | MAGIC imputation + ADAE + xgboost |
| `statescope` → Statescope | Python | BLADE Bayesian deconvolution (GPU, Python 3.10+) |

\* MuSiC is grouped with deconvolution outputs in this benchmark.

## Installation

Clone the repository into your `.claude/skills/` directory. The repository root contains a `SKILL.md`, so the cloned folder itself is discovered as one skill (`tiphd-benchmark-skills`); that root skill documents and routes to the 16 algorithm skills and the bootstrap skills inside it.

```bash
cd ~/.claude/skills
git clone https://github.com/liulan93/tiphd-benchmark-skills.git
```

The resulting layout is `.claude/skills/tiphd-benchmark-skills/<algorithm>/SKILL.md`. If you prefer the algorithm folders themselves to be top-level skills, copy or move the repository contents (all subfolders + `_toolkit`) directly into `~/.claude/skills/` instead.

## Quick start

```bash
# 1. Install environments (R 4.4/Seurat 5, tiphd-torch, tiphd-stats, tiphd-py310)
#    follow setup-env/SKILL.md, or run _toolkit/install_*.{R,sh}

# 2. Download data to ./data
#    follow setup-data/SKILL.md (requires `pip install modelscope`)

# 3. Run an algorithm, e.g. MuSiC
cd music
python batch_run.py      # stage 2: run all 44 (scRNA, bulk) pairs
Rscript evaluate.R       # stage 3: permutation test + gold-standard metrics
```

Run from the project working directory; data defaults to `./data` and results to `./results`, overridable via `TIPHD_DATA_DIR` / `TIPHD_OUT_DIR`.

## Common pipeline

Every algorithm follows the same three stages:

1. **run** (`run_<ALGO>_pair.R/.py`) — one `(cancer, scRNA, bulk)` pair
2. **batch** (`batch_run.py`) — runs all pairs as independent subprocesses, skips existing outputs
3. **evaluate** (`evaluate.R`) — Cox permutation test + gold-standard comparison → `*_final_metrics.csv`

The batch per-pair timeout defaults to **7200 s**, configurable with `TIPHD_PAIR_TIMEOUT`. The GC dataset has ~137 K cells (3–6× the other cancers), so its pairs are much slower; large pairs can also be run directly without the batch timeout.

## Metrics

Each cancer produces Precision, Coverage (recall), and False Rate against a curated gold standard.

## Requirements

- **R ≥ 4.4 + Seurat ≥ 5** (data is stored as Seurat v5 Assay5 `.rds`; R 4.1/Seurat 4 cannot read it correctly)
- Python environments per the `setup-env` skill (3.9 for most tools, 3.10+ for Statescope)
- A CUDA GPU is optional but strongly recommended for the PyTorch tools (SIDISH, scTREND, scSurv, SCAD, scDEAL, TiRank, Statescope)
