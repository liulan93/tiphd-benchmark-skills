"""scDEAL — 单配对运行器（真·DAE + DaNN 域自适应二分类）
用法: python run_scDEAL_pair.py <cancer> <sc_name> <bulk_name>
输出: *_csv（Cell_ID + 细胞类型 + Rank_Label），置换检验在 evaluate.R 完成。

真实 scDEAL = 去噪自编码器(DAE/AE)预训练 + DaNN(MMD 域对齐) + 二分类预测器。
本脚本把 scDEAL 的域自适应框架用于生存二分类：bulk(Dead/Alive) 为 source，
scRNA 为 target，MMD 对齐两者 bottleneck 特征后预测单细胞二分类。
模型类源自 scDEAL-main/{models.py, DaNN/mmd.py}。
（原 scDEAL 是药物敏感性模型，此处按用户约定把标签换成 Dead/Alive。）
依赖：torch / scipy / sklearn / pandas / anndata。
"""
import sys, os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# scDEAL 真源码路径（已自包含到 third_party）
SCDEAL = os.path.join(config.CONFIG_DIR, "third_party", "python", "scDEAL")
sys.path.insert(0, SCDEAL)
from models import PretrainedPredictor, DaNN
import DaNN.mmd as mmd

import torch
import torch.nn as nn
import anndata
from sklearn.preprocessing import StandardScaler

ALGO = "scDEAL"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]
device = "cpu"

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== scDEAL (DaNN/MMD) {cancer} | {sc_name} | {bulk_name} =====")

# ── 1. 载入 bulk(source) 与 scRNA(target) ────────────────────────
bulk = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)
clin = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv"), index_col=0)
common = bulk.columns.intersection(clin.index)
bulk, clin = bulk[common], clin.loc[common]
y_src = clin[cc["scol"]].values.astype(np.float32)

adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
X_tar = np.asarray(adata.X.toarray()) if hasattr(adata.X, "toarray") else np.asarray(adata.X)

shared = np.intersect1d(adata.var_names, bulk.index)
X_src = bulk.loc[shared].T.values.astype(np.float32)
var_idx = pd.Index(adata.var_names).get_indexer(shared)
X_tar = X_tar[:, var_idx].astype(np.float32)

scaler = StandardScaler()
X_src = scaler.fit_transform(X_src); X_tar = scaler.transform(X_tar)
in_dim = len(shared)
print(f"  source: {X_src.shape}, target: {X_tar.shape}, genes {in_dim}")

# ── 2. 真实 scDEAL：PretrainedPredictor + DaNN(MMD) ──────────────
source_model = PretrainedPredictor(input_dim=in_dim, latent_dim=128, h_dims=[256],
                                   hidden_dims_predictor=[64], output_dim=1).to(device)
target_model = PretrainedPredictor(input_dim=in_dim, latent_dim=128, h_dims=[256],
                                   hidden_dims_predictor=[64], output_dim=1).to(device)
dann = DaNN(source_model, target_model).to(device)

Xs = torch.tensor(X_src).float().to(device)
Xt = torch.tensor(X_tar).float().to(device)
ys = torch.tensor(y_src.reshape(-1, 1)).float().to(device)

opt = torch.optim.Adam(dann.parameters(), lr=1e-3)
bce = nn.BCELoss()  # Predictor 输出已是 sigmoid 概率，用 BCELoss 而非 BCEWithLogitsLoss
for _ in range(200):
    dann.train()
    # MMD 要求源/目标 batch 等长：把 target(细胞) 子采样到 source(bulk) 大小
    idx = torch.randperm(Xt.shape[0])[:Xs.shape[0]]
    Xt_b = Xt[idx]
    y_pre, z_src, z_tar = dann(Xs, Xt_b)
    loss_c = bce(y_pre, ys)
    loss_mmd = mmd.mmd_loss(z_src, z_tar)
    loss = loss_c + 0.5 * loss_mmd
    opt.zero_grad(); loss.backward(); opt.step()

# ── 3. 预测单细胞 → Rank_Label（上下 20% 分位）──────────────────
dann.eval()
with torch.no_grad():
    scores = source_model.predictor(target_model.encode(Xt)).cpu().numpy().flatten()  # predictor 已含 sigmoid
q_hi, q_lo = np.quantile(scores, 0.8), np.quantile(scores, 0.2)
rank = np.where(scores >= q_hi, "Rank+", np.where(scores <= q_lo, "Rank-", "Background"))
df = pd.DataFrame({"Cell_ID": adata.obs_names, cc["ct_col"]: adata.obs[cc["ct_col"]].values, "Rank_Label": rank})
df.to_csv(out_csv, index=False)
print(f"  Rank+: {(rank=='Rank+').sum()}, Rank-: {(rank=='Rank-').sum()}, Bg: {(rank=='Background').sum()}")
print(f"  Saved: {os.path.basename(out_csv)}")
