"""Statescope — single-pair runner (BLADE Bayesian deconvolution)
Usage: python run_Statescope_pair.py <cancer> <sc_name> <bulk_name>
Output: *_proportions.csv (samples x cell-type proportions); permutation testing is done in evaluate.R.

Calls the real Statescope Initialize_Statescope() + Deconvolution()
(BLADE Bayesian latent-variable deconvolution).
Dependencies: torch / numba / dill / joblib / anndata; the Statescope source is bundled under
  third_party/python/Statescope (from Statescope-master/src).
Note: BLADE expects Bulk to be linear (library-size corrected) counts; the project's bulk matrices
      are log-normalized, so they are converted back to a linear scale below (see input alignment).
"""
import sys, os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["OMP_NUM_THREADS"] = "4"
# GPU is not disabled here: BLADE auto-detects CUDA (uses GPU if available, otherwise CPU)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# Upstream Statescope source (BLADE + StateDiscovery, bundled under third_party)
SRC = os.path.join(config.CONFIG_DIR, "third_party", "python", "Statescope")
if os.path.isdir(SRC):
    sys.path.insert(0, SRC)
from Statescope.Statescope import Initialize_Statescope

import anndata

ALGO = "Statescope"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
prop_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}_proportions.csv")
if os.path.exists(prop_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== Statescope (BLADE) {cancer} | {sc_name} | {bulk_name} =====")

# ── 1. Load bulk (genes x samples, DataFrame) ───────────────────
bulk = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp",
                                f"{bulk_name}_exp.csv"), index_col=0)
# BLADE expects non-negative linear counts after library-size correction.
# The project bulk is log2-normalized (log2(CPM/10+1)-style), so convert it back to a linear scale.
bulk = bulk.clip(lower=0)              # clip tiny negative values (log-norm numeric error)
Bulk_df = np.expm1(bulk)               # inverse log1p: x -> exp(x)-1
Bulk_df = Bulk_df.clip(lower=0)
print(f"  Bulk: {Bulk_df.shape[0]} genes x {Bulk_df.shape[1]} samples (linearized)")

# ── 2. Load scRNA as the Signature (AnnData with cell types) ─────
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer,
                                       f"{sc_name}.h5ad"))
print(f"  scRNA: {adata.n_obs} cells x {adata.n_vars} genes")

# ── 2b. Signature log1p preprocessing (TiPhD h5ad is normalized but not log1p) ──
# CreateSignature.looks_logged uses two heuristics, max<50 and 99% integer values;
# normalized values that are fractional but not small are misclassified as raw and raise an error.
# Explicit log1p makes looks_logged pass (max becomes small; the data is accepted as log data).
import scanpy as sc
sc.pp.log1p(adata)
print(f"  scRNA after log1p: max={adata.X.max():.2f}")

# ── 3. Initialize_Statescope + BLADE Deconvolution ──────────────
try:
    ss = Initialize_Statescope(Bulk=Bulk_df, Signature=adata,
                               celltype_key=cc["ct_col"],
                               n_highly_variable=3000, Ncores=4)
    # [TiPhD benchmark] Default Nrep=10 (BLADE's multi-initialization ensemble; the benchmark keeps the upstream value).
    # BLADE cost grows quickly with the number of cell types in the signature, so references with many
    # cell types can be extremely slow under the default configuration.
    # STATESCOPE_NREP=1 is only for end-to-end smoke runs (rep is a count of initialization replicas, not a
    # model hyperparameter; fitting hyperparameters such as IterMax are untouched). Smoke results must be
    # explicitly marked as such by the caller and must not be treated as benchmark figures.
    nrep = int(os.environ.get("STATESCOPE_NREP", "10"))
    njob = min(4, nrep) if nrep > 1 else 1
    print(f"  BLADE Deconvolution: Nrep={nrep} Njob={njob} IterMax=100"
          + ("" if nrep == 10 else "  [NON-DEFAULT SMOKE CONFIG]"))
    ss.Deconvolution(Nrep=nrep, Njob=njob, IterMax=100)
    prop_df = ss.Fractions          # index=samples, columns=cell types
except Exception as e:
    print(f"  Statescope FAILED: {e}")
    sys.exit(1)

prop_df.index.name = None
prop_df.to_csv(prop_csv)
print(f"  Fractions: {prop_df.shape[0]} samples x {prop_df.shape[1]} cell types")
print(f"  Saved: {os.path.basename(prop_csv)}")
