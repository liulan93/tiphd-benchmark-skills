"""TiRank — 单配对运行器（真·基因对(REO) + 神经网络 + MMD 域自适应）
用法: python run_TiRank_pair.py <cancer> <sc_name> <bulk_name>
输出: {algo}_{bulk}_{sc}.csv（Cell_ID + 细胞类型列 + Rank_Label），置换检验在 evaluate.R 完成。

真实 TiRank = 相对表达序(REO)基因对特征 + 神经网络编码器(MLP/Transformer/DenseNet)
           + MMD 域自适应 + 多任务损失(Cox/分类/回归)。
本脚本用 Cox 生存模式：bulk(生存 time/status) 为 source，scRNA 为 target，
预测每个单细胞的 Rank_Score → Rank+/Rank-/Background。
真源码已自包含在 third_party/python/TiRank（源自 TiRank-main/tirank）。

依赖：torch / scanpy / anndata / lifelines / optuna / leidenalg / python-igraph / sklearn。
（TiRank 官方建议 GPU + Python3.9；本脚本在 CPU 上跑，需在装好上述依赖的环境里运行。）
"""
import sys, os, warnings, tempfile, pickle
import numpy as np
import pandas as pd
import torch

warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "4"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# TiRank 真源码路径（已自包含到 third_party）
TIRANK = os.path.join(config.CONFIG_DIR, "third_party", "python", "TiRank")
if os.path.isdir(TIRANK):
    sys.path.insert(0, TIRANK)

from tirank.Model import setup_seed, initial_model_para
from tirank.LoadData import load_sc_data, load_bulk_exp, load_bulk_clinical, check_bulk
from tirank.SCSTpreprocess import FilteringAnndata, Normalization, Logtransformation, Clustering, compute_similarity
from tirank.GPextractor import GenePairExtractor
from tirank.Dataloader import generate_val, PackData
from tirank.TrainPre import tune_hyperparameters, Predict

ALGO = "TiRank"
cancer, sc_name, bulk_name = sys.argv[1], sys.argv[2], sys.argv[3]
cc = config.CANCERS[cancer]
# 自动检测 GPU：有 CUDA 则用 GPU（Optuna 超参搜索 + 训练可大幅提速），否则回退 CPU
device = "cuda" if torch.cuda.is_available() else "cpu"

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== TiRank (Cox) {cancer} | {sc_name} | {bulk_name} =====")

setup_seed(619)
savePath = tempfile.mkdtemp(prefix="TiRank_")   # TiRank 中间产物工作目录
for sub in ("1_loaddata", "2_preprocessing", "3_Analysis"):
    os.makedirs(os.path.join(savePath, sub), exist_ok=True)

# ============================================================
# 1. 数据输入（格式要求见 DATA_INPUT.md，输入已按 TiRank 原生格式备好）
#     ① scRNA .h5ad：obs 含细胞类型列
#     ② bulk exp CSV：基因 × 样本（index_col=0）
#     ③ bulk clinical CSV：样本 × 2 列 [time, event]
# ============================================================
sc_h5ad = os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad")
bulk_exp_csv = os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv")
bulk_cli_csv = os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv")

scAnndata = load_sc_data(sc_h5ad, savePath)          # 读 h5ad 并存 anndata.pkl
bulkExp = load_bulk_exp(bulk_exp_csv)
bulkClinical = load_bulk_clinical(bulk_cli_csv)      # 已是 2 列 [time, event]
check_bulk(savePath, bulkExp, bulkClinical)          # 样本取交集并落盘

# ============================================================
# 2. scRNA 预处理（TiRank 自带 QC + 聚类 + 相似度）
#    FilteringAnndata 的 min_count 用 SC 口径(=3)，勿用 ST 默认 5000，
#    否则会误删大量 scRNA 细胞。
# ============================================================
scAnndata = FilteringAnndata(scAnndata, max_count=35000, min_count=3,
                             MT_propor=10, min_cell=1,
                             imgPath=os.path.join(savePath, "2_preprocessing"))
scAnndata = Normalization(scAnndata)
scAnndata = Logtransformation(scAnndata)
scAnndata = Clustering(scAnndata, infer_mode="SC", savePath=savePath)
compute_similarity(savePath=savePath, ann_data=scAnndata)
# 保存预处理后的 anndata 供 GPextractor.load_data() 读取
savePath_2 = os.path.join(savePath, "2_preprocessing")
os.makedirs(savePath_2, exist_ok=True)
with open(os.path.join(savePath_2, "scAnndata.pkl"), "wb") as f:
    pickle.dump(scAnndata, f)

# ============================================================
# 3. 临床划分 + 基因对(REO)提取
# ============================================================
mode = "Cox"
generate_val(savePath=savePath, validation_proportion=0.15, mode=mode)

gp = GenePairExtractor(savePath=savePath, analysis_mode=mode,
                       top_var_genes=3000, top_gene_pairs=1500,
                       p_value_threshold=0.2, max_cutoff=0.8, min_cutoff=-0.8)
gp.load_data()
gp.run_extraction()
gp.save_data()

# ============================================================
# 4. 训练 + 预测
# ============================================================
infer_mode = "SC"
PackData(savePath=savePath, mode=mode, infer_mode=infer_mode, batch_size=1024)
initial_model_para(savePath=savePath, nhead=2, nhid1=96, nhid2=8, n_output=32,
                   nlayers=3, n_pred=1, dropout=0.5, mode=mode,
                   encoder_type="MLP", infer_mode=infer_mode)
tune_hyperparameters(savePath=savePath, device=device, n_trials=10)
Predict(savePath=savePath, mode=mode, do_reject=True, tolerance=0.05, reject_mode="GMM")

# ============================================================
# 5. 收集结果 → 统一格式（Cell_ID + 细胞类型列 + Rank_Label）
#    TiRank 的 Predict 已输出 Rank_Label ∈ {Rank+, Rank-, Background}，与统一格式一致。
# ============================================================
pred_csv = os.path.join(savePath, "3_Analysis", "spot_predict_score.csv")
pred = pd.read_csv(pred_csv, index_col=0)
if cc["ct_col"] not in pred.columns:
    raise KeyError(f"TiRank 输出缺少细胞类型列 {cc['ct_col']}，请确认输入 h5ad 的 obs 含该列")

df = pd.DataFrame({
    "Cell_ID": pred.index.values,
    cc["ct_col"]: pred[cc["ct_col"]].values,
    "Rank_Label": pred["Rank_Label"].values,
})
df.to_csv(out_csv, index=False)
print(f"  Rank+: {(df['Rank_Label'] == 'Rank+').sum()}, "
      f"Rank-: {(df['Rank_Label'] == 'Rank-').sum()}, "
      f"Bg: {(df['Rank_Label'] == 'Background').sum()}")
print(f"  Saved: {os.path.basename(out_csv)}")
