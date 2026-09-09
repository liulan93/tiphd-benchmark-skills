# ============================================================
# MuSiC — 单配对运行器（真·加权 NNLS 反卷积）
# 用法: Rscript run_MuSiC_pair.R <cancer> <sc_name> <bulk_name>
# 输出: *_proportions.csv（样本 × 细胞类型比例），置换检验在 evaluate.R 完成
#
# 与旧版（普通 nnls）的区别：调用真实 MuSiC 的 music_prop()，
# 用跨细胞/跨样本方差对基因加权（weighted-NNLS），而非等权 nnls。
# 真源码位于 ../third_party/R/MuSiC（源自 MuSiC-master）。
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "MuSiC"
cfg <- CANCERS[[cancer]]

# ── 载入真实 MuSiC 函数（music_prop / music_basis / music.basic）──
mus_dir <- file.path(SCRIPT_DIR, "..", "third_party", "R", "MuSiC")
for (f in c("analysis.R", "construct.R", "utils.R")) source(file.path(mus_dir, f))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(SingleCellExperiment); library(nnls)
})

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
prop_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, "_proportions.csv"))
if (file.exists(prop_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== MuSiC (weighted-NNLS) %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整）+ 细胞类型 / 样本标签 ─────────────────
sc <- load_scRNA_counts(cancer, sc_name)
meta <- sc$meta
counts <- sc$counts[rowSums(sc$counts) > 0, , drop = FALSE]

cell_types <- as.character(meta[[cfg$ct_col]])
names(cell_types) <- rownames(meta)

# 自动检测「样本/患者」列（MuSiC 需要跨样本方差）
cand <- c("Patient","patient","SampleName","SampleNameNew","Sample","sample",
          "Donor","donor","Subject","subject","Individual","individual","orig.ident")
sample_col <- NULL
for (cc in cand) if (cc %in% colnames(meta) && length(unique(meta[[cc]])) > 1) { sample_col <- cc; break }
if (is.null(sample_col)) {
  sample_ids <- colnames(counts); names(sample_ids) <- colnames(counts)
  cat("  [warn] 未找到多样本列，退化为每细胞独立（跨细胞方差）\n")
} else {
  sample_ids <- as.character(meta[[sample_col]]); names(sample_ids) <- rownames(meta)
  cat(sprintf("  样本列: %s（%d 个样本）\n", sample_col, length(unique(sample_ids))))
}

keep <- which(!is.na(cell_types) & cell_types != "")
counts <- counts[, keep, drop = FALSE]
cell_types <- cell_types[keep]
sample_ids <- sample_ids[keep]
rm(sc, meta); gc()
cat(sprintf("  scRNA: %d genes x %d cells, %d cell types\n",
            nrow(counts), ncol(counts), length(unique(cell_types))))

# ── 2. 构建 SingleCellExperiment（仅需 counts，MuSiC 内部自行 log 归一化）──
colnames(counts) <- names(cell_types)
sce <- SingleCellExperiment(assays = list(counts = counts))
sce$cell_type <- cell_types[colnames(sce)]
sce$sample_id <- sample_ids[colnames(sce)]

# ── 3. 载入 bulk（已 log-norm）──────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
bulk_exp <- bulk$exp
# 部分数据集经 log 归一化后含极小负值（log2(CPM<1) 导致），MuSiC 的 relative.ab() 要求非负，
# 将负值裁剪为 0（不改变主体分布，仅消除数值误差）。
neg_n <- sum(bulk_exp < 0, na.rm = TRUE)
if (neg_n > 0) {
  cat(sprintf("  [clip] bulk 含 %d 个负值，裁剪为 0\n", neg_n))
  bulk_exp[bulk_exp < 0] <- 0
}
cat(sprintf("  Bulk: %d genes x %d samples\n", nrow(bulk_exp), ncol(bulk_exp)))

# ── 4. music_prop：加权 NNLS 反卷积 ─────────────────────────────
shared <- intersect(rownames(counts), rownames(bulk_exp))
if (length(shared) < 100) stop(sprintf("Too few common genes: %d", length(shared)))

prop <- tryCatch(
  music_prop(bulk.mtx = bulk_exp[shared, , drop = FALSE],
             sc.sce = sce[shared, ],
             clusters = "cell_type",
             samples = "sample_id",
             verbose = FALSE),
  error = function(e) { cat(sprintf("  music_prop FAILED: %s\n", e$message)); return(NULL) })
if (is.null(prop)) quit(save = "no", status = 1)

prop_mat <- as.matrix(prop$Est.prop.weighted)  # 样本 × 细胞类型（Est.prop.weighted 行=样本，列=细胞类型）
prop_mat[is.na(prop_mat)] <- 0
cat(sprintf("  Deconvolution: %d samples x %d cell types (shared %d genes)\n",
            nrow(prop_mat), ncol(prop_mat), length(shared)))

write.csv(prop_mat, prop_csv)
cat(sprintf("  Saved: %s\n", basename(prop_csv)))
