"""SIDISH — 单配对运行器（VAE + DeepCox）
用法: python run_SIDISH_pair.py <cancer> <sc_name> <bulk_name>
单细胞数据不子采样，读取完整细胞。只输出 Rank+ / Background。"""
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
# 自动检测 GPU：有 CUDA 则用 GPU（VAE/DeepCox 可提速 10-50 倍），否则回退 CPU
device = "cuda" if torch.cuda.is_available() else "cpu"

# 各癌种的病人列名（SIDISH 需要）
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

# ── 1. 载入 scRNA（完整，不子采样） ────────────────────────────
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
if hasattr(adata.X, "toarray"):
    adata.X = adata.X.toarray()
adata.X = adata.X.astype(np.float32)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
if patient_col not in adata.obs.columns:
    adata.obs[patient_col] = "batch1"
adata.obs["celltype_major"] = adata.obs[cc["ct_col"]].astype(str) \
    if cc["ct_col"] in adata.obs.columns else "Unknown"
print(f"  scRNA: {adata.n_obs} cells x {adata.n_vars} genes (完整)")

# ── 2. 载入 bulk + clinical ───────────────────────────────────
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
print(f"  Final: {adata.n_obs} cells x {adata.n_vars} genes")

survival_df = clinical[[cc["tcol"], cc["scol"]]].copy()
survival_df.columns = ["Overall_survival_days", "Sample_Status"]
survival_df["Overall_survival_days"] = survival_df["Overall_survival_days"].astype(float).clip(lower=0.1)
survival_df["Sample_Status"] = survival_df["Sample_Status"].astype(int)
survival_df.index = survival_df.index.astype(str)

# ── 3. 预处理（processed=True，跳过内部归一化） ─────────────────
adata_pp, bulk_pp = preprocess(
    adata, bulk_t, survival_df,
    patient_id=patient_col, celltype_name="celltype_major",
    processed=True, survival_="Overall_survival_days", status="Sample_Status")

# ── 4. 训练 SIDISH（CPU 快速模式） ─────────────────────────────
# lr 降为 3e-4/2e-4 避免 NaN（之前 lr=1e-3 在 epoch 95 后 loss 爆炸）
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

# ── 5. 输出（仅 Rank+ / Background） ───────────────────────────
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
