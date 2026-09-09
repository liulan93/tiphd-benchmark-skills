"""scTREND — 单配对运行器（VAE + DeepCOLOR + 分段风险模型）
用法: python run_scTREND_pair.py <cancer> <sc_name> <bulk_name>
单细胞数据不子采样，读取完整细胞。"""
import sys, os, warnings, time
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

import torch
torch.set_num_threads(1)
import scanpy as sc
import anndata

# Windows 下强制 DataLoader 单进程（避免多进程崩溃）
import torch.utils.data
_orig = torch.utils.data.DataLoader.__init__
def _patched(self, *a, **k):
    k["num_workers"] = 0
    _orig(self, *a, **k)
torch.utils.data.DataLoader.__init__ = _patched

from sctrend import workflow

ALGO = "scTREND"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]

bkey = {"AML": "orig.ident", "CRC": "orig.ident", "HCC": "patient",
        "LUAD": "orig.ident", "GC": "orig.ident"}[cancer]

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== scTREND {cancer} | {sc_name} | {bulk_name} =====")
t0 = time.time()

# ── 1. 载入 scRNA（完整，不子采样） ────────────────────────────
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
if hasattr(adata.X, "toarray"):
    adata.X = np.round(adata.X.toarray()).astype(np.float32)
else:
    adata.X = np.round(np.asarray(adata.X)).astype(np.float32)
if bkey not in adata.obs.columns:
    adata.obs[bkey] = "batch1"
print(f"  scRNA: {adata.n_obs} cells x {adata.n_vars} genes (完整)")

# ── 2. 载入 bulk + clinical → bulk AnnData ────────────────────
bulk_exp = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)
clinical = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv"), index_col=0)
common = bulk_exp.columns.intersection(clinical.index)
bulk_exp, clinical = bulk_exp[common], clinical.loc[common]

bulk_X = np.maximum(0, np.round(2 ** bulk_exp.T.values - 1)).astype(np.float32)
bulk_adata = anndata.AnnData(
    X=bulk_X,
    obs=clinical.rename(columns={cc["tcol"]: "survival_times", cc["scol"]: "vital_status"}),
    var=pd.DataFrame(index=bulk_exp.index))
bulk_adata.obs["survival_times"] = bulk_adata.obs["survival_times"].astype(float).clip(lower=0.1)
bulk_adata.obs["vital_status"] = bulk_adata.obs["vital_status"].apply(lambda x: "Dead" if int(x) == 1 else "Alive")

# ── 3. 预处理 ─────────────────────────────────────────────────
sc_adata, bulk_adata = workflow.scTREND_preprocess(
    adata, bulk_adata, per=0.01, n_top_genes=5000,
    highly_variable="bulk", driver_genes=None)

# ── 4. 训练 ───────────────────────────────────────────────────
param_path = os.path.join(out_dir, f"_sctrend_{bulk_name}_{sc_name}.pt")
sc_adata, bulk_adata, _, _, _ = workflow.run_scTREND(
    sc_adata, bulk_adata, param_save_path=param_path, warm_path=None,
    epoch=100, batch_key=bkey, driver_genes=None, driver_bulk_adata=None,
    survival_time_label="survival_times", survival_time_censor="vital_status", edges=None)

# ── 5. 输出（上下 20% 分位） ───────────────────────────────────
scores = sc_adata.obs["beta_z_mean"].values if "beta_z_mean" in sc_adata.obs.columns \
    else np.zeros(sc_adata.n_obs)
probs = 0.2
q_hi, q_lo = np.quantile(scores, 1 - probs), np.quantile(scores, probs)
df = pd.DataFrame({"Cell_ID": sc_adata.obs_names,
                   cc["ct_col"]: sc_adata.obs[cc["ct_col"]].values})
df["Rank_Label"] = "Background"
df.loc[scores >= q_hi, "Rank_Label"] = "Rank+"
df.loc[scores <= q_lo, "Rank_Label"] = "Rank-"
print(f"  Rank+: {(df['Rank_Label'] == 'Rank+').sum()}, "
      f"Rank-: {(df['Rank_Label'] == 'Rank-').sum()}")
df.to_csv(out_csv, index=False)
print(f"  Saved: {os.path.basename(out_csv)} [{time.time()-t0:.1f}s]")
