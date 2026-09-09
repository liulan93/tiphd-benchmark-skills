# ============================================================
# ScPP — 单配对运行器
# 用法: Rscript run_ScPP_pair.R <cancer> <sc_name> <bulk_name>
# 流程: 单变量 Cox 找预后 marker 基因 → AUCell 打分 → ScPP 选择
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "ScPP"
cfg <- CANCERS[[cancer]]

library(Seurat); library(Matrix); library(AUCell)

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== ScPP %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整细胞，不子采样） ────────────────────────
sc <- load_scRNA_counts(cancer, sc_name)
sc_obj <- CreateSeuratObject(counts = sc$counts, meta.data = sc$meta, project = sc_name)
sc_obj <- NormalizeData(sc_obj, verbose = FALSE)
rm(sc); gc()
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_obj), ncol(sc_obj)))

# ── 2. 载入 bulk + survival ───────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
surv <- build_survival(bulk$clinical, cfg$tcol, cfg$scol)
cat(sprintf("  Bulk: %d genes x %d samples, events: %d\n",
            nrow(bulk$exp), ncol(bulk$exp), sum(surv$status)))

# ── 3. 单变量 Cox 筛选 marker 基因 ─────────────────────────────
# bulk$exp: genes × samples (列名是样本ID), surv: samples × 2 (行名是样本ID)
shared <- intersect(colnames(bulk$exp), rownames(surv))
bulk_sub <- bulk$exp[, shared, drop = FALSE]
surv_sub <- surv[shared, , drop = FALSE]

bulk_t <- as.data.frame(t(bulk_sub))
colnames(bulk_t) <- make.names(colnames(bulk_t))
bulk_t$time <- surv_sub$time
bulk_t$status <- surv_sub$status

coef_vec <- c(); p_vec <- c()
for (g in make.names(rownames(bulk_sub))) {
  fit <- tryCatch(coxph(as.formula(paste("Surv(time, status) ~", g)), data = bulk_t),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)$coefficients
    coef_vec[g] <- s[1, "coef"]
    p_vec[g]   <- s[1, "Pr(>|z|)"]
  }
}
fdr <- p.adjust(p_vec, method = "fdr")
gene_pos <- names(which(fdr < 0.05 & coef_vec > 0 & abs(coef_vec) > 0.5))
gene_neg <- names(which(fdr < 0.05 & coef_vec < 0 & abs(coef_vec) > 0.5))
cat(sprintf("  gene_pos: %d, gene_neg: %d\n", length(gene_pos), length(gene_neg)))

if ((length(gene_pos) + length(gene_neg)) < 5) {
  cat("  marker 基因不足 (< 5)，跳过\n")
  quit(save = "no", status = 0)
}

# ── 4. AUCell 打分 ─────────────────────────────────────────────
geneList <- list(gene_pos = gene_pos, gene_neg = gene_neg)
cellrankings <- AUCell_buildRankings(sc_obj@assays$RNA$data, plotStats = FALSE)
cellAUC <- AUCell_calcAUC(geneList, cellrankings)

meta <- sc_obj@meta.data
meta$AUCup   <- as.numeric(getAUC(cellAUC)["gene_pos", ])
meta$AUCdown <- as.numeric(getAUC(cellAUC)["gene_neg", ])

# ── 5. ScPP 选择（上下 20% 分位） ──────────────────────────────
probs <- 0.2
ScPP_pos <- rownames(meta)[meta$AUCup   >= quantile(meta$AUCup,   probs = 1 - probs) &
                           meta$AUCdown <= quantile(meta$AUCdown, probs = probs)]
ScPP_neg <- rownames(meta)[meta$AUCup   <= quantile(meta$AUCup,   probs = probs) &
                           meta$AUCdown >= quantile(meta$AUCdown, probs = 1 - probs)]

# ── 6. 输出统一格式 ────────────────────────────────────────────
result_df <- data.frame(Cell_ID = colnames(sc_obj), stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- meta[result_df$Cell_ID, cfg$ct_col]
result_df$Rank_Label <- "Background"
result_df$Rank_Label[result_df$Cell_ID %in% ScPP_pos] <- "Rank+"
result_df$Rank_Label[result_df$Cell_ID %in% ScPP_neg] <- "Rank-"

cat(sprintf("  Rank+: %d, Rank-: %d, Bg: %d\n",
            sum(result_df$Rank_Label == "Rank+"),
            sum(result_df$Rank_Label == "Rank-"),
            sum(result_df$Rank_Label == "Background")))
write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
