"""SIDISH — single-pair runner (VAE + DeepCox)
Usage: python run_SIDISH_pair.py <cancer> <sc_name> <bulk_name>
Single-cell data is used in full (no downsampling). Outputs Rank+ / Background only."""
import sys, os, warnings, time
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

import torch
import scanpy as sc
import anndata
from SIDISH.SIDISH import SIDISH, preprocess

ALGO = "SIDISH"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]
# Auto-detect GPU: use CUDA when available, otherwise fall back to CPU
device = "cuda" if torch.cuda.is_available() else "cpu"

# Patient column per cancer (required by SIDISH)
patient_col = {"AML": "orig.ident", "CRC": "orig.ident", "HCC": "patient",
               "LUAD": "orig.ident", "GC": "orig.ident"}[cancer]

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== SIDISH {cancer} | {sc_name} | {bulk_name} =====")
t0 = time.time()

# ── 1. Load scRNA (full, no downsampling) ─────────────────────
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
# [TiPhD benchmark] Memory layout only (values unchanged): normalize on the sparse
# matrix first, and densify only after the <=2000-gene filter below, so the full dense
# matrix is never materialized on the largest references (SIDISH's internal copies
# would otherwise make OOM likely).
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
if patient_col not in adata.obs.columns:
    adata.obs[patient_col] = "batch1"
adata.obs["celltype_major"] = adata.obs[cc["ct_col"]].astype(str) \
    if cc["ct_col"] in adata.obs.columns else "Unknown"
print(f"  scRNA: {adata.n_obs} cells x {adata.n_vars} genes (full)")

# ── 2. Load bulk + clinical ───────────────────────────────────
bulk_exp = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)
clinical = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv"), index_col=0)
common = bulk_exp.columns.intersection(clinical.index)
bulk_exp, clinical = bulk_exp[common], clinical.loc[common]

bulk_t = bulk_exp.T.copy()
common_genes = np.intersect1d(adata.var_names.astype(str), bulk_t.columns.astype(str))
if len(common_genes) < 100:
    print("  ERROR: too few common genes!"); sys.exit(1)
if len(common_genes) > 2000:
    top_genes = bulk_t[common_genes].mean(axis=0).nlargest(2000).index.tolist()
else:
    top_genes = common_genes.tolist()
adata = adata[:, adata.var_names.astype(str).isin(top_genes)].copy()
bulk_t = bulk_t[adata.var_names.astype(str)]
# [TiPhD benchmark] Densify now (float32) that the matrix is filtered to <=2000 genes;
# numerically identical to densifying before normalization, only memory layout differs.
if hasattr(adata.X, "toarray"):
    adata.X = adata.X.toarray()
adata.X = adata.X.astype(np.float32)
print(f"  Final: {adata.n_obs} cells x {adata.n_vars} genes")

survival_df = clinical[[cc["tcol"], cc["scol"]]].copy()
survival_df.columns = ["Overall_survival_days", "Sample_Status"]
survival_df["Overall_survival_days"] = survival_df["Overall_survival_days"].astype(float).clip(lower=0.1)
survival_df["Sample_Status"] = survival_df["Sample_Status"].astype(int)
survival_df.index = survival_df.index.astype(str)

# ── 3. Preprocess (processed=True skips internal normalization) ─
adata_pp, bulk_pp = preprocess(
    adata, bulk_t, survival_df,
    patient_id=patient_col, celltype_name="celltype_major",
    processed=True, survival_="Overall_survival_days", status="Sample_Status")

# ── 4. Train SIDISH ───────────────────────────────────────────
# Lower learning rates (3e-4/2e-4) for training stability; avoid 1e-3, where the
# adversarial loss can diverge to NaN late in training.
sidish = SIDISH(adata_pp, bulk_pp, device=device, seed=42, use_spatial_graph=False)
sidish.init_Phase1(epochs=50, i_epochs=50, latent_size=32, layer_dims=[128, 64],
                   batch_size=min(64, adata_pp.n_obs), optimizer="Adam", lr=3e-4,
                   lr_3=2e-4, dropout=0.2, type="Normal")
sidish.init_Phase2(epochs=50, hidden=32, lr=3e-4, dropout=0.2,
                   test_size=0.3, batch_size_bulk=8)
path = os.path.join(out_dir, f"_sidish_{bulk_name}_{sc_name}")
os.makedirs(path, exist_ok=True)
adata_result = sidish.train(iterations=50, percentile=0.2, steepness=2.0,
                            path=path, num_workers=0, show=False)

# ── 5. Output (Rank+ / Background only) ───────────────────────
result_df = pd.DataFrame({"Cell_ID": adata_result.obs_names,
                          cc["ct_col"]: adata_result.obs[cc["ct_col"]].values
                          if cc["ct_col"] in adata_result.obs.columns
                          else adata_result.obs["celltype_major"].values})
result_df["Rank_Label"] = "Background"
if "SIDISH" in adata_result.obs.columns:
    result_df.loc[adata_result.obs["SIDISH"].astype(str).values == "h", "Rank_Label"] = "Rank+"
print(f"  Rank+: {(result_df['Rank_Label'] == 'Rank+').sum()}, "
      f"Bg: {(result_df['Rank_Label'] == 'Background').sum()}")
result_df.to_csv(out_csv, index=False)
print(f"  Saved: {os.path.basename(out_csv)} [{time.time()-t0:.1f}s]")
