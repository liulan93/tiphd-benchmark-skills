"""Statescope — 单配对运行器（真·BLADE 贝叶斯反卷积）
用法: python run_Statescope_pair.py <cancer> <sc_name> <bulk_name>
输出: *_proportions.csv（样本 × 细胞类型比例），置换检验在 evaluate.R 完成。

与旧版（普通 nnls）的区别：调用真实 Statescope 的 Initialize_Statescope()
+ Deconvolution()（BLADE 贝叶斯潜变量反卷积）。
依赖：torch / numba / dill / joblib / anndata；Statescope 源码已自包含在
  third_party/python/Statescope（源自 Statescope-master/src）。
注意：BLADE 期望 Bulk 为线性（library-size 校正）counts；本项目 bulk 为 log-norm，
     此处原样传入，尺度差异需在正式评估前确认。
"""
import sys, os, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["OMP_NUM_THREADS"] = "4"
# 不强制禁用 GPU：BLADE 会自动检测 CUDA（有 GPU 用 GPU，无则回退 CPU）

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# Statescope 真源码路径（BLADE + StateDiscovery，已自包含到 third_party）
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

# ── 1. 载入 bulk（基因 × 样本，DataFrame）───────────────────────
bulk = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp",
                                f"{bulk_name}_exp.csv"), index_col=0)
# BLADE 期望 library-size 校正后的线性 counts（非负）。
# 本项目 bulk 为 log2 归一化（log2(CPM/10+1) 类），此处还原为线性尺度。
bulk = bulk.clip(lower=0)              # 裁剪极小负值（log-norm 数值误差）
Bulk_df = np.expm1(bulk)               # log1p 逆变换：x -> exp(x)-1
Bulk_df = Bulk_df.clip(lower=0)
print(f"  Bulk: {Bulk_df.shape[0]} genes x {Bulk_df.shape[1]} samples (linearized)")

# ── 2. 载入 scRNA 作为 Signature（AnnData，含细胞类型）──────────
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer,
                                       f"{sc_name}.h5ad"))
print(f"  scRNA: {adata.n_obs} cells x {adata.n_vars} genes")

# ── 2b. Signature log1p 预处理（TiPhD h5ad 已是 normalized 但非 log1p） ──
# CreateSignature.looks_logged 用 max<50 + 99%int 两个启发式判断，
# 当归一化值有小数（浮点）但不算小时会被判为 raw 而报错。
# 显式 log1p 让 looks_logged 通过（max 自然变小、int 比例不变但作为 log 数据被认为 OK）。
import scanpy as sc
sc.pp.log1p(adata)
print(f"  scRNA after log1p: max={adata.X.max():.2f}")

# ── 3. Initialize_Statescope + BLADE Deconvolution ──────────────
try:
    ss = Initialize_Statescope(Bulk=Bulk_df, Signature=adata,
                               celltype_key=cc["ct_col"],
                               n_highly_variable=3000, Ncores=4)
    ss.Deconvolution(Nrep=10, Njob=4, IterMax=100)
    prop_df = ss.Fractions          # index=样本, columns=细胞类型
except Exception as e:
    print(f"  Statescope FAILED: {e}")
    sys.exit(1)

prop_df.index.name = None
prop_df.to_csv(prop_csv)
print(f"  Fractions: {prop_df.shape[0]} samples x {prop_df.shape[1]} cell types")
print(f"  Saved: {os.path.basename(prop_csv)}")
