"""SCAD — 单配对运行器（真·对抗迁移学习二分类）
用法: python run_SCAD_pair.py <cancer> <sc_name> <bulk_name>
输出: *_csv（Cell_ID + 细胞类型 + Rank_Label），置换检验在 evaluate.R 完成。

真实 SCAD = FX(特征提取器) + Discriminator(域判别器, 梯度反转对抗) + MTLP(二分类预测器)。
本脚本把 SCAD 的 AITL 框架用于生存二分类：bulk(Dead/Alive) 为 source，scRNA 为 target，
对抗对齐两域特征分布后预测单细胞二分类。模型类源自 SCAD-main/model/SCADmodules.py。
（原 SCAD 是药物反应模型，此处按用户约定把标签换成 Dead/Alive。）
依赖：torch / sklearn / pandas / anndata。
"""
import sys, os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# SCAD 真源码路径（模型类，已自包含到 third_party）
SCAD_DIR = os.path.join(config.CONFIG_DIR, "third_party", "python", "SCAD")
sys.path.insert(0, SCAD_DIR)
from SCADmodules import FX, MTLP, Discriminator, grad_reverse

import torch
import torch.nn as nn
import anndata
from sklearn.preprocessing import StandardScaler

ALGO = "SCAD"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]
device = "cpu"

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== SCAD (adversarial transfer) {cancer} | {sc_name} | {bulk_name} =====")

# ── 1. 载入 bulk(source, 带 Dead/Alive 标签) 与 scRNA(target) ────
bulk = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)
clin = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv"), index_col=0)
common = bulk.columns.intersection(clin.index)
bulk, clin = bulk[common], clin.loc[common]
y_src = clin[cc["scol"]].values.astype(np.float32)   # 0/1 (Alive/Dead)

adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
X_tar = np.asarray(adata.X.toarray()) if hasattr(adata.X, "toarray") else np.asarray(adata.X)

shared = np.intersect1d(adata.var_names, bulk.index)
X_src = bulk.loc[shared].T.values.astype(np.float32)
# np.searchsorted requires sorted var_names; use get_indexer for unsorted gene names
var_idx = pd.Index(adata.var_names).get_indexer(shared)
X_tar = X_tar[:, var_idx].astype(np.float32)

scaler = StandardScaler()
X_src = scaler.fit_transform(X_src)
X_tar = scaler.transform(X_tar)
print(f"  source: {X_src.shape} (bulk), target: {X_tar.shape} (scRNA), genes {len(shared)}")

# ── 2. 真实 SCAD 模型：FX + Discriminator + MTLP ─────────────────
in_dim = len(shared); h_dim = 512; z_dim = 256
fx = FX(dropout_rate=0.2, input_dim=in_dim, h_dim=h_dim, z_dim=z_dim).to(device)
pred = MTLP(dropout_rate=0.2, h_dim=h_dim, z_dim=64).to(device)          # 输入 = FX 的 h_dim 输出（FX 的 z_dim 参数未被使用）
disc = Discriminator(dropout_rate=0.2, h_dim=h_dim, z_dim=64).to(device) # Discriminator 内部已做 grad_reverse

Xs = torch.tensor(X_src).float().to(device)
Xt = torch.tensor(X_tar).float().to(device)
ys = torch.tensor(y_src.reshape(-1, 1)).float().to(device)

opt = torch.optim.Adam(list(fx.parameters()) + list(pred.parameters()) + list(disc.parameters()), lr=1e-3)
bce = nn.BCELoss()   # MTLP / Discriminator 输出均为 Sigmoid 概率

for _ in range(200):
    fx.train(); pred.train(); disc.train()
    zs, zt = fx(Xs), fx(Xt)
    pred_loss = bce(pred(zs), ys)
    # 对抗：域判别器区分 source(1) vs target(0)；grad_reverse 已内置于 Discriminator.forward
    disc_s, disc_t = disc(zs), disc(zt)
    disc_loss = bce(disc_s, torch.ones_like(disc_s)) + bce(disc_t, torch.zeros_like(disc_t))
    loss = pred_loss + 0.5 * disc_loss
    opt.zero_grad(); loss.backward(); opt.step()

# ── 3. 预测单细胞 → Rank_Label（上下 20% 分位）──────────────────
fx.eval(); pred.eval()
with torch.no_grad():
    scores = torch.sigmoid(pred(fx(Xt))).cpu().numpy().flatten()
q_hi, q_lo = np.quantile(scores, 0.8), np.quantile(scores, 0.2)
rank = np.where(scores >= q_hi, "Rank+", np.where(scores <= q_lo, "Rank-", "Background"))
df = pd.DataFrame({"Cell_ID": adata.obs_names, cc["ct_col"]: adata.obs[cc["ct_col"]].values, "Rank_Label": rank})
df.to_csv(out_csv, index=False)
print(f"  Rank+: {(rank=='Rank+').sum()}, Rank-: {(rank=='Rank-').sum()}, Bg: {(rank=='Background').sum()}")
print(f"  Saved: {os.path.basename(out_csv)}")
