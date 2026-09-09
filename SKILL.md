---
name: tiphd-benchmark-skills
description: TiPhD benchmark suite for 16 single-cell phenotype-association and bulk-deconvolution algorithms (MuSiC, ScPP, scSTAR2, SCAD, SCIPAC, scDEAL, scSurv, scTREND, scAB, Scissor, scPAS, PIPET, SIDISH, scPER, Statescope, TiRank) across 5 cancers (AML/CRC/HCC/LUAD/GC). Use when the user wants to set up environments, download the TiPhD test data, or run/evaluate any of these algorithms end-to-end (batch_run.py -> evaluate.R -> Precision/Coverage/False-Rate metrics). Triggers on "跑评测", "运行算法", "run benchmark", "TiPhD", an algorithm name, or installing/setting up the benchmark.
---

# TiPhD Benchmark Skills

This repository bundles **16 algorithm skills** plus the shared toolkit and two bootstrap skills (environment + data). It is a *meta-skill*: the individual algorithm folders are not auto-discovered by the agent when the repo is cloned as a nested directory, so this file tells you how to locate and invoke them by path.

## Layout (relative to this SKILL.md)

```
tiphd-benchmark-skills/
├── SKILL.md                  ← you are here
├── README.md
├── _toolkit/                 ← shared config (config.R/config.py), install scripts, third-party packages
├── setup-env/SKILL.md        ← install R + Python environments
├── setup-data/SKILL.md       ← download the TiPhD dataset from ModelScope
├── run-benchmark/SKILL.md    ← dispatcher for running any algorithm
└── <algorithm>/             ← one folder per algorithm (see table below)
    ├── SKILL.md              ← per-algorithm documentation
    ├── batch_run.py          ← run all (scRNA, bulk) pairs
    ├── evaluate.R            ← permutation test + gold-standard metrics
    └── run_*_pair.{R,py}     ← single-pair runner
```

Let `$ROOT` be the directory containing this file. Always read a sub-skill's `$ROOT/<algorithm>/SKILL.md` for its exact dependencies and commands.

## Algorithms

| Folder | Algorithm | Lang | Class |
|--------|-----------|------|-------|
| `music` | MuSiC | R | B (weighted NNLS deconvolution) |
| `scpp` | ScPP | R | A (Cox markers + AUCell) |
| `scstar2` | scSTAR2 | R | A (OPLS-DA) |
| `scad` | SCAD | Python | A (VAE + adversarial DA) |
| `scipac` | SCIPAC | R | A (iterative partition + Cox) |
| `scdeal` | scDEAL | Python | A (VAE + DaNN) |
| `scsurv` | scSurv | Python | A (VAE + Cox) |
| `sctrend` | scTREND | Python | A (VAE + DeepCOLOR) |
| `scab` | scAB | R | A (plus-only adaptive boosting) |
| `scissor` | Scissor | R | A (Seurat network + Cox) |
| `scpas` | scPAS | R | A (Cox signature) |
| `pipet` | PIPET | R | A (marker cosine classification) |
| `sidish` | SIDISH | Python | A (plus-only VAE + DeepCox, GPU) |
| `scper` | scPER | Python+R | B (MAGIC + ADAE + xgboost) |
| `statescope` | Statescope | Python | B (BLADE Bayesian; Python 3.10+, GPU) |
| `tirank` | TiRank | Python | A (REO gene pairs + NN + MMD, GPU) |

A-class tools output per-cell `Rank+`/`Rank-`/`Background` labels; B-class tools output per-sample cell-type proportions.

## Standard workflow

### 1. Resolve the root

Determine `$ROOT` from the location of this SKILL.md (the cloned repo directory). Run all commands from the user's project working directory; data/results default to `./data` and `./results` (env vars `TIPHD_DATA_DIR` / `TIPHD_OUT_DIR`).

### 2. Install environments — follow `setup-env/SKILL.md`

- R ≥ 4.4 + Seurat ≥ 5 (data is Seurat v5 Assay5 `.rds`)
- `tiphd-torch` conda env (Python 3.9 + PyTorch) — most Python tools
- `tiphd-stats` conda env (Python 3.9 + lifelines/optuna/leidenalg) — TiRank
- `tiphd-py310` conda env (Python 3.10) — Statescope only
- Install scripts: `$ROOT/_toolkit/install_R_packages.R`, `$ROOT/_toolkit/install_python_packages.sh`

### 3. Download data — follow `setup-data/SKILL.md`

Downloads the 5-cancer TiPhD dataset from ModelScope (`fansailing/TiPhD_test_data`) to `./data` (both `.rds` for R and `.h5ad` for Python tools).

### 4. Run an algorithm

Always read `$ROOT/<algorithm>/SKILL.md` first for the correct interpreter/environment. From the working directory:

```bash
cd $ROOT/<algorithm>
# stage 2 — all 44 pairs (idempotent; skips existing outputs)
python batch_run.py          # R tools use the R 4.4 env; Python tools use their conda env
# stage 3 — permutation test + gold-standard comparison
Rscript evaluate.R          # 5 cancers, one *_final_metrics.csv each
```

Single pair: `Rscript run_<ALGO>_pair.R <cancer> <sc_name> <bulk_name>` or the Python equivalent.

Refer to `run-benchmark/SKILL.md` for the canonical algorithm → environment map, the 44 pair list, and natural-language command routing.

## Important notes

- **Per-pair timeout** defaults to 7200 s (`TIPHD_PAIR_TIMEOUT`). The GC dataset has ~137 K cells (3–6× the others); its pairs are much slower and can be run directly without the batch timeout.
- **R/Seurat version is hard-required**: R ≥ 4.4 and Seurat ≥ 5 with Matrix ≥ 1.6.4; R 4.1/Seurat 4 cannot read the Assay5 data correctly.
- **GPU**: SIDISH, Statescope, and TiRank auto-detect CUDA and are much faster on a GPU. Install a CUDA build of PyTorch in their env. SCAD/scDEAL/scSurv/scTREND default to CPU.
- **Gold-standard filename quirk**: the file is `gold_standard_all_LUDA.csv` (original spelling, not `LUAD`) — do not rename.
- **AML bulk names** are `TCGA`, `wave12`, `wave34` (not GSE ids).

Do not modify algorithm parameters. If a pair fails, follow that algorithm's Troubleshooting table in its own SKILL.md.
