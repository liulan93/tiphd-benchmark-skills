# ============================================================
# Scissor — 单配对运行器
# 用法: Rscript run_Scissor_pair.R <cancer> <sc_name> <bulk_name>
# 流程: Seurat 预处理(含网络) → Scissor() Cox → Rank+/Rank-/Background
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

# Prevent pthread_create errors on high-core machines
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "Scissor"
cfg <- CANCERS[[cancer]]

library(Scissor); library(Seurat)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== Scissor %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整）→ Seurat 预处理（Scissor 需要细胞网络） ──
sc <- load_scRNA_counts(cancer, sc_name)
sc_obj <- CreateSeuratObject(counts = sc$counts, meta.data = sc$meta, project = sc_name)
rm(sc); gc()
sc_obj <- NormalizeData(sc_obj, verbose = FALSE)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_obj <- ScaleData(sc_obj, verbose = FALSE)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(sc_obj), verbose = FALSE)
sc_obj <- FindNeighbors(sc_obj, dims = 1:30, verbose = FALSE)
# Seurat v5 → v4 Assay 兼容（Scissor 需要 v4 Assay）
sc_obj[["RNA"]] <- as(sc_obj[["RNA"]], "Assay")
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_obj), ncol(sc_obj)))

# ── 2. 载入 bulk + survival ───────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
surv <- build_survival(bulk$clinical, cfg$tcol, cfg$scol)
cat(sprintf("  Bulk: %d genes x %d samples, events: %d\n",
            nrow(bulk$exp), ncol(bulk$exp), sum(surv$status)))

# ── 3. 运行 Scissor (Cox) ──────────────────────────────────────
result <- tryCatch({
  Scissor(bulk_dataset = bulk$exp, sc_dataset = sc_obj,
          phenotype = surv, family = "cox", alpha = NULL, cutoff = 0.2,
          Save_file = file.path(out_dir, paste0("_inputs_", bulk_name, "_", sc_name, ".RData")))
}, error = function(e) { cat(sprintf("  Scissor FAILED: %s\n", e$message)); return(NULL) })

if (is.null(result)) { cat("  无结果\n"); quit(save = "no", status = 0) }
cat(sprintf("  Scissor+: %d, Scissor-: %d\n", length(result$Scissor_pos), length(result$Scissor_neg)))

# ── 4. 输出统一格式 ────────────────────────────────────────────
meta <- sc_obj@meta.data
result_df <- data.frame(Cell_ID = colnames(sc_obj), stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- meta[result_df$Cell_ID, cfg$ct_col]
result_df$Rank_Label <- "Background"
result_df$Rank_Label[result_df$Cell_ID %in% result$Scissor_pos] <- "Rank+"
result_df$Rank_Label[result_df$Cell_ID %in% result$Scissor_neg] <- "Rank-"

write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
