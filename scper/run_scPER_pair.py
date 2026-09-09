"""scPER — 单配对运行器（真·对抗自编码器 ADAE + xgboost 比例估计）
用法: python run_scPER_pair.py <cancer> <sc_name> <bulk_name>
输出: *_proportions.csv（样本 × 细胞类型比例）。

真实 scPER 三步：MAGIC 插补+HVG → 对抗自编码器(ADAE, Keras/TF2.10) → xgboost 比例。
本脚本把三步串起来，调用 third_party/python/scPER 里的真实脚本/类。
依赖（重）：tensorflow==2.10.0, keras==2.10.0, magic-impute/Rmagic, xgboost。
注意：ADAE 需要标签里的 'Batch' 列做批次去混淆；本项目单数据集，用患者/样本 ID 作为
     Batch 代理（若 meta 里没有，则退化）。原始脚本含硬编码 /example、/results 路径，
     已由本脚本在临时目录下生成输入并用环境变量指定输出。
"""
import sys, os, warnings, subprocess, tempfile
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# scPER 真源码路径（ADAE + xgboost，已自包含到 third_party）
SCPER = os.path.join(config.CONFIG_DIR, "third_party", "python", "scPER")
RSCRIPT = config.RSCRIPT

ALGO = "scPER"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
prop_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}_proportions.csv")
if os.path.exists(prop_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== scPER (ADAE+xgboost) {cancer} | {sc_name} | {bulk_name} =====")

import anndata

# ── 1. 导出 scRNA / bulk / label 为 scPER 需要的 CSV ────────────
tmp = tempfile.mkdtemp(prefix="scPER_")
# 让 preprocess_data.R / main.py / proportion_estimate.R 把中间产物写到 tmp（替代原版硬编码 /example、/results）
os.environ["SCPER_OUT"] = tmp
# MAGIC 插补走 reticulate → tiphd-stats 环境里的 magic-impute（base 装不了）
os.environ["RETICULATE_PYTHON"] = config.MAGIC_PYTHON
adata = anndata.read_h5ad(os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad"))
bulk = pd.read_csv(os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv"), index_col=0)

# scRNA: 基因 × 细胞
if hasattr(adata.X, "toarray"):
    X = adata.X.toarray()
else:
    X = np.asarray(adata.X)
sc_df = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
sc_csv = os.path.join(tmp, "sc.csv"); sc_df.to_csv(sc_csv)
bulk_csv = os.path.join(tmp, "bulk.csv"); bulk.to_csv(bulk_csv)

# label: 行=细胞, 含 celltype 与 Batch（用患者/样本 ID 代理；缺则用细胞类型）
ct = adata.obs[cc["ct_col"]].astype(str)
batch_col = None
for c in ["Patient","patient","SampleName","Sample","sample","orig.ident"]:
    if c in adata.obs.columns and adata.obs[c].nunique() > 1:
        batch_col = c; break
batch = adata.obs[batch_col].astype(str) if batch_col else ct
label = pd.DataFrame({"Batch": batch.values, "celltype": ct.values}, index=adata.obs_names)
label_csv = os.path.join(tmp, "label.csv"); label.to_csv(label_csv)

print(f"  exported: sc {sc_df.shape}, bulk {bulk.shape}, label {label.shape} (Batch={batch_col or 'celltype'})")

# ── 2. 三步真实流水线 ───────────────────────────────────────────
def run(cmd, tag):
    print(f"  [{tag}] {' '.join(cmd[:6])} ...")
    r = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    if r.returncode != 0:
        print(f"  [{tag}] FAILED rc={r.returncode}\n{r.stderr[-1500:]}")
        return False
    return True

# 2a. 预处理（MAGIC 插补 + top5k HVG）
if not run([RSCRIPT, os.path.join(SCPER, "preprocess_data.R"), sc_csv, bulk_csv], "preprocess"):
    sys.exit(1)
# 2b. ADAE 训练（Keras，用 r-tensorflow 环境的 Python，因为 base 装不了 tf2.10）
if not run([config.SCPER_PYTHON, os.path.join(SCPER, "main.py"),
            os.path.join(tmp, "reference_top5k_imputation.csv"),
            label_csv, os.path.join(tmp, "bulk_5k_genes_matched.csv")], "ADAE"):
    sys.exit(1)
# 2c. xgboost 比例估计
if not run([RSCRIPT, os.path.join(SCPER, "proportion_estimate.R"),
            os.path.join(tmp, "ADAE_100_latents.tsv"),
            os.path.join(tmp, "bulk_100_latents.tsv"), label_csv], "proportion"):
    sys.exit(1)

# ── 3. 收集比例结果 ─────────────────────────────────────────────
res = os.path.join(tmp, "scPER_predicted_proportions.txt")
if not os.path.exists(res):
    print(f"  未找到结果文件 {res}")
    sys.exit(1)
prop = pd.read_csv(res, sep="\t", index_col=0)
prop.to_csv(prop_csv)
print(f"  Saved: {os.path.basename(prop_csv)} ({prop.shape})")
