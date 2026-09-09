"""scSurv — 单配对运行器（生存分析）
用法: python run_scSurv_pair.py <cancer> <sc_name> <bulk_name>
单细胞数据不子采样，读取完整细胞。"""
import sys, os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

import torch
import anndata
from scsurv.workflow import run_scSurv, post_process

ALGO = "scSurv"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]

bkey = {"AML": "orig.ident", "CRC": "orig.ident", "HCC": "sample",
        "LUAD": "orig.ident", "GC": "Sample"}[cancer]

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== scSurv {cancer} | {sc_name} | {bulk_name} =====")

# ── 1. 载入 scRNA（完整，不子采样），转换为整数 counts ─────────
sc_file = os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad")
sc_adata = anndata.read_h5ad(sc_file)
if sc_adata.raw is not None:
    sc_adata.X = sc_adata.raw.X.toarray() if hasattr(sc_adata.raw.X, "toarray") else np.array(sc_adata.raw.X)
elif hasattr(sc_adata.X, "toarray"):
    sc_adata.X = sc_adata.X.toarray()
else:
    sc_adata.X = np.array(sc_adata.X)
if sc_adata.X.dtype.kind == "f":
    x_sample = sc_adata.X[:min(100, sc_adata.X.shape[0])].flatten()
    if not np.allclose(x_sample, np.floor(x_sample)):
        sc_adata.X = np.maximum(0, np.floor(np.expm1(sc_adata.X)))
sc_adata.X = np.floor(sc_adata.X).astype(np.float32)
if bkey not in sc_adata.obs.columns:
    sc_adata.obs[bkey] = "batch1"
print(f"  scRNA: {sc_adata.n_obs} cells x {sc_adata.n_vars} genes (完整)")

# ── 2. 载入 bulk + clinical → bulk AnnData ────────────────────
bulk_exp = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)
clinical = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv"), index_col=0)
common = bulk_exp.columns.intersection(clinical.index)
bulk_exp, clinical = bulk_exp[common], clinical.loc[common]

bulk_adata = anndata.AnnData(
    X=np.maximum(0, bulk_exp.T.values).astype(np.float32),
    obs=clinical.copy(), var=pd.DataFrame(index=bulk_exp.index))
bulk_adata.obs[cc["tcol"]] = bulk_adata.obs[cc["tcol"]].astype(float)
bulk_adata.obs[cc["scol"]] = bulk_adata.obs[cc["scol"]].astype(float)
print(f"  Bulk: {bulk_adata.n_obs} samples x {bulk_adata.n_vars} genes")

# ── 3. 共同基因过滤（按 bulk 表达量均值选 top 5000） ────────────
common_genes = np.intersect1d(sc_adata.var_names, bulk_adata.var_names)
bulk_means = np.array(bulk_adata[:, common_genes].X.mean(axis=0)).flatten()
common_top = common_genes[np.argsort(bulk_means)[-5000:]]
sc_adata = sc_adata[:, common_top].copy()
bulk_adata = bulk_adata[:, common_top].copy()
sc_adata.X = np.maximum(0, sc_adata.X)
bulk_adata.X = np.maximum(0, bulk_adata.X)
print(f"  Filtered to {sc_adata.n_vars} genes")

# ── 4. 运行 scSurv ─────────────────────────────────────────────
param_path = os.path.join(out_dir, f"_checkpoint_{bulk_name}_{sc_name}.pt")
sc_adata, bulk_adata, _, _, scsurv_exp = run_scSurv(
    sc_adata, bulk_adata, param_save_path=param_path, epoch=100,
    batch_key=bkey, survival_time_label=cc["tcol"], survival_time_censor=cc["scol"])
sc_adata, bulk_adata, _ = post_process(scsurv_exp, sc_adata, bulk_adata, save_memory=True)

# ── 5. 输出（raw_beta_z 阈值） ─────────────────────────────────
result_df = pd.DataFrame({"Cell_ID": sc_adata.obs_names,
                          cc["ct_col"]: sc_adata.obs[cc["ct_col"]].values,
                          "raw_beta_z": sc_adata.obs["raw_beta_z"].values,
                          "beta_z": sc_adata.obs["beta_z"].values})
threshold = np.std(result_df["raw_beta_z"]) * 0.5
result_df["Rank_Label"] = np.where(result_df["raw_beta_z"] > threshold, "Rank+",
                           np.where(result_df["raw_beta_z"] < -threshold, "Rank-", "Background"))
print(f"  Rank+: {(result_df['Rank_Label'] == 'Rank+').sum()}, "
      f"Rank-: {(result_df['Rank_Label'] == 'Rank-').sum()}, "
      f"Bg: {(result_df['Rank_Label'] == 'Background').sum()}")
result_df.to_csv(out_csv, index=False)
print(f"  Saved: {os.path.basename(out_csv)}")

if os.path.exists(param_path):
    os.remove(param_path)
