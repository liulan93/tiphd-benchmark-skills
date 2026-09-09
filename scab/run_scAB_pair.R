# ============================================================
# scAB — 单配对运行器（只输出 Rank+ / Background，无 Rank-）
# 用法: Rscript run_scAB_pair.R <cancer> <sc_name> <bulk_name>
# 流程: run_seurat 预处理 → create_scAB → select_K → scAB → findSubset
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "scAB"
cfg <- CANCERS[[cancer]]

library(scAB); library(Seurat)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== scAB %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整）→ 基因过滤 → scAB 预处理 ──────────────
sc <- load_scRNA_counts(cancer, sc_name)
sc_obj <- CreateSeuratObject(counts = sc$counts, meta.data = sc$meta, project = sc_name)
rm(sc); gc()

keep <- rownames(sc_obj)[rowSums(sc_obj@assays$RNA$counts > 1) >= 20]
sc_obj <- subset(sc_obj, features = keep)
sc_obj <- run_seurat(sc_obj, verbose = TRUE)
sc_obj[["RNA"]] <- as(sc_obj[["RNA"]], "Assay")   # scAB 需要 v4 Assay
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_obj), ncol(sc_obj)))

# ── 2. 载入 bulk + survival ───────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
surv <- build_survival(bulk$clinical, cfg$tcol, cfg$scol)
cat(sprintf("  Bulk: %d genes x %d samples, events: %d\n",
            nrow(bulk$exp), ncol(bulk$exp), sum(surv$status)))

# ── 3. 运行 scAB ───────────────────────────────────────────────
scAB_data <- create_scAB(sc_obj, bulk$exp, surv, "survival")
K <- select_K(scAB_data, K_max = 10, repeat_times = 10, maxiter = 2000, verbose = TRUE)
cat(sprintf("  K = %d\n", K))
scAB_result <- scAB(Object = scAB_data, K = K, alpha = 0.005, alpha_2 = 0.005, maxiter = 2000)
sc_result <- findSubset(sc_obj, scAB_Object = scAB_result, tred = 2)

# ── 4. 输出统一格式（仅 Rank+ / Background） ────────────────────
meta <- sc_result@meta.data
cat(sprintf("  scAB+: %d, Other: %d\n",
            sum(meta$scAB_select == "scAB+ cells"), sum(meta$scAB_select == "Other cells")))

result_df <- data.frame(Cell_ID = colnames(sc_result), stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- meta[[cfg$ct_col]]
result_df$scAB_Label <- meta$scAB_select
result_df$Rank_Label <- ifelse(result_df$scAB_Label == "scAB+ cells", "Rank+", "Background")

write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
