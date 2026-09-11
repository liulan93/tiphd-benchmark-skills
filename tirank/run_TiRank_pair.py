"""TiRank — single-pair runner (gene-pair (REO) features + neural network + MMD domain adaptation)
Usage: python run_TiRank_pair.py <cancer> <sc_name> <bulk_name>
Output: {algo}_{bulk}_{sc}.csv (Cell_ID + cell-type column + Rank_Label); permutation testing is done in evaluate.R.

TiRank = relative-expression-ordering (REO) gene-pair features + neural-network encoder
         (MLP/Transformer/DenseNet) + MMD domain adaptation + multitask loss (Cox/classification/regression).
This runner uses Cox survival mode: bulk (survival time/status) is the source domain, scRNA the target;
each single cell receives a Rank_Score -> Rank+/Rank-/Background.
The upstream source is bundled under third_party/python/TiRank (from TiRank-main/tirank).

Dependencies: torch / scanpy / anndata / lifelines / optuna / leidenalg / python-igraph / sklearn.
(TiRank upstream recommends GPU + Python 3.9; the skill defaults to CPU — see TIRANK_GPU below —
and must run in an environment with the dependencies above installed.)
"""
import sys, os, warnings, tempfile, pickle, shutil, atexit
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import anndata

warnings.filterwarnings("ignore")
# The skill forces CPU by default; set TIRANK_GPU=1 to use CUDA (same fp32 model, numerically
# equivalent; CUDA substantially accelerates the optuna hyperparameter search and training).
if os.environ.get("TIRANK_GPU") != "1":
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ.setdefault("OMP_NUM_THREADS", "4")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

# Upstream TiRank source (bundled under third_party)
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
# Auto-detect GPU: use CUDA when available (greatly speeds up optuna search + training), else CPU
device = "cuda" if torch.cuda.is_available() else "cpu"

out_dir = os.path.join(config.OUT_BASE, ALGO, cancer)
os.makedirs(out_dir, exist_ok=True)
out_csv = os.path.join(out_dir, f"{ALGO}_{bulk_name}_{sc_name}.csv")
if os.path.exists(out_csv):
    print(f"SKIP {cancer}|{sc_name}|{bulk_name}")
    sys.exit(0)

print(f"===== TiRank (Cox) {cancer} | {sc_name} | {bulk_name} =====")

setup_seed(619)
savePath = tempfile.mkdtemp(prefix="TiRank_")   # working directory for TiRank intermediates
# [TiPhD benchmark] remove the temp directory on process exit (success or failure) to avoid disk fill
atexit.register(lambda: shutil.rmtree(savePath, ignore_errors=True))
for sub in ("1_loaddata", "2_preprocessing", "3_Analysis"):
    os.makedirs(os.path.join(savePath, sub), exist_ok=True)

# ============================================================
# 1. Data inputs (TiRank native format)
#     (1) scRNA .h5ad: obs contains the cell-type column
#     (2) bulk exp CSV: genes x samples (index_col=0)
#     (3) bulk clinical CSV: samples x exactly 2 columns [time, event]
#         TiRank selects columns positionally (GPextractor columns[0:2],
#         generate_val Cox branch iloc[:,-2:]); multi-column raw clinical
#         tables are trimmed below to exactly two columns per config
#         tcol/scol (input alignment only, no values changed)
# ============================================================
sc_h5ad = os.path.join(config.DATA_DIR, "scRNA", cancer, f"{sc_name}.h5ad")
bulk_exp_csv = os.path.join(config.DATA_DIR, "bulk", cancer, "exp", f"{bulk_name}_exp.csv")
bulk_cli_csv = os.path.join(config.DATA_DIR, "bulk", cancer, "clinical", f"{bulk_name}_clinical.csv")

# [TiPhD benchmark] Storage layout only: the benchmark h5ad files store X as dense arrays
# (a dense matrix for some cancer references can expand to tens of GB in memory), and QC/PCA
# creates several additional copies, so the dense path easily OOMs. Build CSR at load time:
#   * dense X: open the h5ad backed, read row chunks from HDF5 and assemble CSR from nonzeros,
#     peak memory ~ one chunk instead of the whole matrix;
#   * sparse X: read normally.
# Only the storage layout changes; no values are modified.
import scipy.sparse as _sp
import h5py as _h5py

def _load_sc_as_csr(path, chunk=4096):
    with _h5py.File(path, "r") as hf:
        x_is_dense = isinstance(hf["X"], _h5py.Dataset)
    if not x_is_dense:
        adata = sc.read_h5ad(path)                    # X already sparse
        if not _sp.issparse(adata.X):
            adata.X = _sp.csr_matrix(adata.X)
        return adata
    backed = sc.read_h5ad(path, backed="r")
    ds = backed.file["X"]
    n_obs, n_var = ds.shape
    rows_l, cols_l, vals_l = [], [], []
    for st in range(0, n_obs, chunk):
        blk = ds[st:st + chunk]                       # (chunk, n_var) dense block
        ri, ci = np.nonzero(blk)
        if len(ri):
            rows_l.append(ri + st); cols_l.append(ci); vals_l.append(blk[ri, ci])
        del blk
    Xcsr = _sp.csr_matrix(
        (np.concatenate(vals_l),
         (np.concatenate(rows_l), np.concatenate(cols_l))),
        shape=(n_obs, n_var)) if vals_l else _sp.csr_matrix((n_obs, n_var), dtype=ds.dtype)
    obsm = {k: backed.obsm[k][()].copy() for k in backed.obsm.keys()}
    adata = anndata.AnnData(X=Xcsr, obs=backed.obs.copy(), var=backed.var.copy(),
                            obsm=obsm)
    backed.file.close()
    return adata

scAnndata = _load_sc_as_csr(sc_h5ad)
bulkExp = load_bulk_exp(bulk_exp_csv)
bulkClinical = load_bulk_clinical(bulk_cli_csv)
# [TiPhD benchmark] Input alignment: TiRank selects clinical columns at fixed positions
# (GPextractor uses columns[0:2]; the Cox branch of generate_val uses iloc[:, -2:]),
# so the table must contain exactly the two columns [time, event]. Some bulk clinical
# files are raw multi-column tables (e.g. HCC GSE14520 ships Age/Gender/.../Cirrhosis,
# CRC GSE14333 ships 11 columns); trim here to config.CANCERS tcol/scol, as in the SIDISH
# runner. Only column selection — no values or model hyperparameters change.
bulkClinical = bulkClinical[[cc["tcol"], cc["scol"]]].copy()
check_bulk(savePath, bulkExp, bulkClinical)          # intersect samples and write to disk

# ============================================================
# 2. scRNA preprocessing (TiRank's built-in QC + clustering + similarity)
#    FilteringAnndata min_count uses the SC setting (=3); do NOT use the ST default
#    of 5000, which would wrongly drop large numbers of scRNA cells.
# ============================================================
scAnndata = FilteringAnndata(scAnndata, max_count=35000, min_count=3,
                             MT_propor=10, min_cell=1,
                             imgPath=os.path.join(savePath, "2_preprocessing"))
scAnndata = Normalization(scAnndata)
scAnndata = Logtransformation(scAnndata)
scAnndata = Clustering(scAnndata, infer_mode="SC", savePath=savePath)
compute_similarity(savePath=savePath, ann_data=scAnndata)
# Save the preprocessed anndata for GPextractor.load_data() to read
savePath_2 = os.path.join(savePath, "2_preprocessing")
os.makedirs(savePath_2, exist_ok=True)
# [TiPhD benchmark] Storage only: if X was densified, convert it back to CSR (values unchanged);
# GPextractor/PackData both accept sparse matrices, avoiding pickle bloat on large references.
import scipy.sparse as _sp
if not _sp.issparse(scAnndata.X):
    scAnndata.X = _sp.csr_matrix(scAnndata.X)
with open(os.path.join(savePath_2, "scAnndata.pkl"), "wb") as f:
    pickle.dump(scAnndata, f)

# ============================================================
# 3. Clinical split + gene-pair (REO) extraction
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
# 4. Training + prediction
# ============================================================
infer_mode = "SC"
PackData(savePath=savePath, mode=mode, infer_mode=infer_mode, batch_size=1024)
initial_model_para(savePath=savePath, nhead=2, nhid1=96, nhid2=8, n_output=32,
                   nlayers=3, n_pred=1, dropout=0.5, mode=mode,
                   encoder_type="MLP", infer_mode=infer_mode)
tune_hyperparameters(savePath=savePath, device=device, n_trials=10)
Predict(savePath=savePath, mode=mode, do_reject=True, tolerance=0.05, reject_mode="GMM")

# ============================================================
# 5. Collect results -> unified format (Cell_ID + cell-type column + Rank_Label)
#    TiRank's Predict already emits Rank_Label in {Rank+, Rank-, Background}, matching it.
# ============================================================
pred_csv = os.path.join(savePath, "3_Analysis", "spot_predict_score.csv")
pred = pd.read_csv(pred_csv, index_col=0)
if cc["ct_col"] not in pred.columns:
    raise KeyError(f"TiRank output is missing the cell-type column {cc['ct_col']}; ensure the input h5ad obs contains it")

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
