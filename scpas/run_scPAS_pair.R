# ============================================================
# scPAS — 单配对运行器（Cox，自带 2000 次置换检验）
# 用法: Rscript run_scPAS_pair.R <cancer> <sc_name> <bulk_name>
# 流程: run_Seurat 预处理 → scPAS() → scPAS+/scPAS-/0 → Rank_Label
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

# Prevent pthread_create errors on high-core machines (preprocessCore uses OpenMP)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "scPAS"
cfg <- CANCERS[[cancer]]

library(scPAS); library(Seurat)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== scPAS %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整）→ scPAS 预处理 ────────────────────────
sc <- load_scRNA_counts(cancer, sc_name)
sc_obj <- CreateSeuratObject(counts = sc$counts, meta.data = sc$meta, project = sc_name)
rm(sc); gc()

keep <- rownames(sc_obj)[rowSums(sc_obj@assays$RNA$counts > 1) >= 20]
sc_obj <- subset(sc_obj, features = keep)
sc_obj <- run_Seurat(sc_obj, verbose = TRUE)
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_obj), ncol(sc_obj)))

# ── 2. 载入 bulk + survival ───────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
surv <- build_survival(bulk$clinical, cfg$tcol, cfg$scol)

# 过滤低表达基因（>25% 样本为 0 的基因）
bulk_exp <- bulk$exp[apply(bulk$exp, 1, function(x) sum(x == 0) < 0.25 * ncol(bulk$exp)), , drop = FALSE]
cat(sprintf("  Bulk: %d genes x %d samples, events: %d\n",
            nrow(bulk_exp), ncol(bulk_exp), sum(surv$status)))

# ── 3. 运行 scPAS ──────────────────────────────────────────────
sc_result <- tryCatch({
  scPAS(bulk_dataset = bulk_exp, sc_dataset = sc_obj, phenotype = surv,
        assay = "RNA", imputation = TRUE, nfeature = 3000,
        alpha = NULL, network_class = "SC", family = "cox", FDR.threshold = 0.05)
}, error = function(e) { cat(sprintf("  scPAS FAILED: %s\n", e$message)); return(NULL) })

if (is.null(sc_result)) quit(save = "no", status = 0)

meta <- sc_result@meta.data
cat(sprintf("  scPAS+: %d, scPAS-: %d, 0: %d\n",
            sum(meta$scPAS == "scPAS+"), sum(meta$scPAS == "scPAS-"), sum(meta$scPAS == "0")))

# ── 4. 输出统一格式（scPAS 自带 2000 次置换 + FDR） ─────────────
result_df <- data.frame(
  Cell_ID = colnames(sc_result),
  scPAS_Label = meta$scPAS, scPAS_NRS = meta$scPAS_NRS,
  scPAS_Pvalue = meta$scPAS_Pvalue, scPAS_FDR = meta$scPAS_FDR,
  stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- meta[[cfg$ct_col]]
result_df$Rank_Label <- dplyr::case_when(
  result_df$scPAS_Label == "scPAS+" ~ "Rank+",
  result_df$scPAS_Label == "scPAS-" ~ "Rank-",
  TRUE ~ "Background")

write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
