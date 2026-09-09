# ============================================================
# PIPET — 单配对运行器（二分类 Dead/Alive）
# 用法: Rscript run_PIPET_pair.R <cancer> <sc_name> <bulk_name>
# 流程: bulk 按生存状态分 Dead/Alive → 上下调 marker → PIPET()
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "PIPET"
cfg <- CANCERS[[cancer]]

library(PIPET); library(Seurat)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== PIPET %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整，不子采样） ────────────────────────────
sc <- load_scRNA_counts(cancer, sc_name)
counts <- sc$counts
cell_meta <- sc$meta
rm(sc); gc()
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(counts), ncol(counts)))

# ── 2. 载入 bulk + clinical，构建 Dead/Alive 二分类 ────────────
bulk <- load_bulk(cancer, bulk_name)
status <- as.numeric(bulk$clinical[[cfg$scol]])
colData <- data.frame(status_class = factor(ifelse(status == 1, "Dead", "Alive"),
                                           levels = c("Alive", "Dead")),
                      row.names = rownames(bulk$clinical))
cat(sprintf("  Bulk: %d genes x %d samples, Dead=%d, Alive=%d\n",
            nrow(bulk$exp), ncol(bulk$exp), sum(status == 1), sum(status == 0)))

# ── 3. 手动构建 marker（Dead vs Alive 的 fold-change） ─────────
group_dead  <- which(colData$status_class == "Dead")
group_alive <- which(colData$status_class == "Alive")
dead_mean  <- rowMeans(bulk$exp[, group_dead, drop = FALSE])
alive_mean <- rowMeans(bulk$exp[, group_alive, drop = FALSE])
logFC <- dead_mean - alive_mean          # 正值 = Dead 高表达

n_top <- 100
up_idx   <- order(logFC, decreasing = TRUE)[1:n_top]
down_idx <- order(logFC, decreasing = FALSE)[1:n_top]
markers <- data.frame(
  genes = c(rownames(bulk$exp)[up_idx], rownames(bulk$exp)[down_idx]),
  class = c(rep("Dead", n_top), rep("Alive", n_top)),
  log2FoldChange = c(logFC[up_idx], logFC[down_idx]),
  stringsAsFactors = FALSE)
markers$class <- factor(markers$class, levels = c("Alive", "Dead"))

# ── 4. PIPET ───────────────────────────────────────────────────
pipet_res <- tryCatch({
  PIPET(SC_data = counts, markers = markers, gene_col = "genes", class_col = "class",
        nPerm = 1000, distance = "cosine", nCores = 1)
}, error = function(e) { cat(sprintf("  PIPET FAILED: %s\n", e$message)); return(NULL) })
if (is.null(pipet_res)) quit(save = "no", status = 0)

# ── 5. 输出统一格式 ────────────────────────────────────────────
result_df <- data.frame(Cell_ID = rownames(pipet_res),
                        PIPET_Prediction = pipet_res$prediction,
                        PIPET_Pvalue = pipet_res$Pvalue,
                        PIPET_FDR = pipet_res$FDR,
                        stringsAsFactors = FALSE)
result_df$Rank_Label <- ifelse(result_df$PIPET_FDR < 0.05,
                               ifelse(result_df$PIPET_Prediction == "Dead", "Rank+", "Rank-"),
                               "Background")
result_df[[cfg$ct_col]] <- cell_meta[result_df$Cell_ID, cfg$ct_col]

cat(sprintf("  Rank+: %d, Rank-: %d, Bg: %d\n",
            sum(result_df$Rank_Label == "Rank+"),
            sum(result_df$Rank_Label == "Rank-"),
            sum(result_df$Rank_Label == "Background")))
write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
