---
name: tiphd-benchmark-skills
description: TiPhD benchmark suite for 16 single-cell phenotype-association and bulk-deconvolution algorithms (MuSiC, ScPP, scSTAR2, SCAD, SCIPAC, scDEAL, scSurv, scTREND, scAB, Scissor, scPAS, PIPET, SIDISH, scPER, Statescope, TiRank) across 5 cancers (AML/CRC/HCC/LUAD/GC). Use for setting up environments, downloading the TiPhD test data, or running/evaluating any algorithm end-to-end (batch_run.py -> evaluate.R -> Precision/Coverage/False-Rate metrics), including running ALL tools at once. Natural-language entry point: "用xx skill跑评测", "运行所有评测工具", "run all benchmarks", "跑评测", "运行算法", an algorithm name, or installing/setting up the benchmark.
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
- `tiphd-torch` conda env (Python 3.9 + PyTorch) — most Python tools (incl. SIDISH; auto-uses CUDA when a CUDA build is installed)
- `tiphd-stats` conda env (Python 3.9 + lifelines/optuna/leidenalg) — TiRank. The shipped yml installs **CPU-only** PyTorch; for GPU install a CUDA-matching PyTorch build in this env and set `TIRANK_GPU=1`
- `tiphd-py310` conda env (Python 3.10) — Statescope only (use a CUDA PyTorch build for tractable runtimes)
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

### 5. Common natural-language requests

When the user phrases a request in natural language, map it to these commands. `$ROOT` is this repository directory.

| User says | Action |
|-----------|--------|
| "安装环境 / setup env / install dependencies" | Follow `setup-env/SKILL.md` (R env + conda envs) |
| "下载数据 / download the data" | Follow `setup-data/SKILL.md` (ModelScope → `./data`) |
| "用 \<algo\> 跑评测 / run \<algo\>" | `cd $ROOT/<algo>` then `python batch_run.py` with the correct interpreter, then `Rscript evaluate.R` |
| "只跑某个癌种 / run \<algo\> on \<cancer\>" | Run the pair runner for that cancer, then `Rscript evaluate.R <cancer>` |
| "运行所有评测工具 / run all algorithms / benchmark everything" | `bash $ROOT/run_all.sh` (batch) then `bash $ROOT/run_evaluate_all.sh` (evaluation) |
| "只跑 R 类 / Python 类工具" | `R_ONLY=1 bash run_all.sh` or `PY_ONLY=1 bash run_all.sh` |
| "汇总结果 / collect metrics" | `bash $ROOT/run_evaluate_all.sh`, then read each `results/<algo>/<cancer>/*_final_metrics.csv` |

Two convenience drivers are provided at the repository root:

- **`run_all.sh`** — runs every tool's `batch_run.py`, selecting the interpreter per tool (system `Rscript` for R tools; the `tiphd-torch`, `tiphd-stats`, `tiphd-py310` conda envs for Python tools). Accepts a list of tool names to restrict the run, or `R_ONLY=1` / `PY_ONLY=1`. It is idempotent (finished pairs are skipped) and prints a failure summary.
- **`run_evaluate_all.sh`** — runs `evaluate.R` for every tool and writes the per-cancer `*_final_metrics.csv`.

Activate the environments before calling them: an R ≥ 4.4 + Seurat 5 environment exposing `Rscript` on `PATH`, plus the three conda envs. Set `TIRANK_GPU=1` in the environment to run TiRank on GPU; Statescope/SIDISH auto-detect CUDA.

### 6. Algorithm → environment quick map

| Tool(s) | Interpreter / env |
|---------|-------------------|
| music, scpp, scstar2, scab, scipac, scissor, scpas, pipet | Rscript (R ≥ 4.4 + Seurat 5) |
| scad, scdeal, scsurv, sctrend, sidish, scper | `tiphd-torch` python (SIDISH auto-uses GPU) |
| tirank | `tiphd-stats` python (set `TIRANK_GPU=1` for GPU) |
| statescope | `tiphd-py310` python (Python 3.10+, auto-uses GPU) |

## Important notes

- **Per-pair timeout** defaults to 7200 s (`TIPHD_PAIR_TIMEOUT`). The GC dataset has ~137 K cells (3–6× the others); its pairs are much slower and can be run directly without the batch timeout.
- **R/Seurat version is hard-required**: R ≥ 4.4 and Seurat ≥ 5 with Matrix ≥ 1.6.4; R 4.1/Seurat 4 cannot read the Assay5 data correctly.
- **GPU**: SIDISH and Statescope auto-detect CUDA; TiRank is opt-in via `TIRANK_GPU=1` (its default env ships CPU-only PyTorch). Install a CUDA build of PyTorch in the relevant env. SCAD/scDEAL/scSurv/scTREND default to CPU.
- **TiRank clinical-table contract**: the bulk clinical CSV must be exactly `[time, event]` two columns — TiRank selects columns by position, and its per-gene Cox failures are swallowed by `except: continue`, surfacing much later as a misleading `0 Risk genes` error. The run script auto-trims the TiPhD tables via `_toolkit/config.py`; custom data needs the same treatment.
- **Statescope scaling**: BLADE cost grows with the number of annotated cell types (not bulk sample count); high-cell-type datasets at the default `Nrep=10` can be impractically slow even on a GPU. Its progress log is just bare EM iteration numbers, and it has no checkpoint. `STATESCOPE_NREP=1` exists only for explicitly-labelled end-to-end smoke runs; benchmark results always use the default 10.
- **No checkpoints**: SIDISH, TiRank, and Statescope pairs restart from scratch when killed; only pairs whose output CSV already exists are skipped. Launch long pairs from a resilient session.
- **Memory**: benchmark `.h5ad` files can store dense `X`; the Python run scripts stream these in backed mode as CSR (memory-layout change only, values untouched), so prefer the bundled run scripts over ad-hoc loading.
- **Gold-standard filename quirk**: the file is `gold_standard_all_LUDA.csv` (original spelling, not `LUAD`) — do not rename.
- **AML bulk names** are `TCGA`, `wave12`, `wave34` (not GSE ids).

Do not modify algorithm parameters or hyperparameters (input-alignment trimming and storage layout are not parameter changes). The sole sanctioned non-default knob is `STATESCOPE_NREP` for labelled smoke runs. If a pair fails, follow that algorithm's Troubleshooting table in its own SKILL.md.
