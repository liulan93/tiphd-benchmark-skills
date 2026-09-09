# ============================================================
# SCIPAC — 单配对运行器（Cox）
# 用法: Rscript run_SCIPAC_pair.R <cancer> <sc_name> <bulk_name>
# 流程: 共同基因 + HVG 预处理 → sc.bulk.pca → seurat.ct → SCIPAC()
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "SCIPAC"
cfg <- CANCERS[[cancer]]

library(SCIPAC); library(Seurat)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== SCIPAC %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整） ──────────────────────────────────────
sc <- load_scRNA_counts(cancer, sc_name)
sc_counts <- sc$counts
cell_meta <- sc$meta
rm(sc); gc()
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_counts), ncol(sc_counts)))

# ── 2. 载入 bulk ──────────────────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
cat(sprintf("  Bulk: %d genes x %d samples\n", nrow(bulk$exp), ncol(bulk$exp)))

# ── 3. 预处理：共同基因 + HVG ──────────────────────────────────
custom_preprocess <- function(sc_counts, bulk_exp, hvg = 1000) {
  overlap_genes <- intersect(rownames(sc_counts), rownames(bulk_exp))
  sc_sub  <- sc_counts[overlap_genes, , drop = FALSE]
  bulk_sub <- bulk_exp[overlap_genes, , drop = FALSE]
  seurat_obj <- CreateSeuratObject(sc_sub)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = hvg, verbose = FALSE)
  hvg_genes <- VariableFeatures(seurat_obj)
  list(sc.dat.preprocessed   = as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[hvg_genes, , drop = FALSE]),
       bulk.dat.preprocessed = as.matrix(bulk_sub[hvg_genes, , drop = FALSE]))
}

preproc <- custom_preprocess(sc_counts, bulk$exp)
n_pc <- min(60, ncol(preproc$bulk.dat.preprocessed) - 1)
pca_res <- sc.bulk.pca(preproc$sc.dat.preprocessed, preproc$bulk.dat.preprocessed,
                       n.pc = n_pc, do.pca.sc = TRUE)
ct_res <- seurat.ct(pca_res$sc.dat.rot, res = 2.0)

# ── 4. SCIPAC (Cox) ────────────────────────────────────────────
sample_names <- rownames(pca_res$bulk.dat.rot)
times <- bulk$clinical[sample_names, cfg$tcol]; times[times <= 0] <- 0.1
y <- cbind(time = times, status = bulk$clinical[sample_names, cfg$scol])
rownames(y) <- sample_names
cat(sprintf("  Bulk: %d samples, %d events, K=%d\n", nrow(y), sum(y[, "status"]), ct_res$k))

scipac_res <- SCIPAC(pca_res$bulk.dat.rot, y, family = "cox", ct_res,
                     bt.size = 50, numCores = 7, nfold = 10)

# ── 5. 输出统一格式 ────────────────────────────────────────────
result_df <- data.frame(Cell_ID = rownames(scipac_res),
                        SCIPAC_sig = scipac_res$sig,
                        stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- cell_meta[result_df$Cell_ID, cfg$ct_col]
result_df$Rank_Label <- ifelse(result_df$SCIPAC_sig == "Sig.pos", "Rank+",
                               ifelse(result_df$SCIPAC_sig == "Sig.neg", "Rank-", "Background"))

cat(sprintf("  Rank+: %d, Rank-: %d, Bg: %d\n",
            sum(result_df$Rank_Label == "Rank+"),
            sum(result_df$Rank_Label == "Rank-"),
            sum(result_df$Rank_Label == "Background")))
write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
