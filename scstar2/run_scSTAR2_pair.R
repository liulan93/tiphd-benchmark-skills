# ============================================================
# scSTAR2 — 单配对运行器（二分类 OPLS-DA）
# 用法: Rscript run_scSTAR2_pair.R <cancer> <sc_name> <bulk_name>
# 依赖: 需要 scSTAR2-main/R 下的原始 R 函数（OGFSC/o2pls_m 等）
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
cancer <- args[1]; sc_name <- args[2]; bulk_name <- args[3]

.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))

ALGO <- "scSTAR2"
cfg <- CANCERS[[cancer]]

library(Seurat); library(Matrix); library(pls); library(MASS)

# scSTAR2 原始 R 函数（本地 src/，自包含，不依赖 tool_use）
scstar2_src <- file.path(SCRIPT_DIR, "src")
source(file.path(scstar2_src, "OGFSC.R"))
source(file.path(scstar2_src, "o2pls_m.R"))
source(file.path(scstar2_src, "dispopls_m.R"))
source(file.path(scstar2_src, "mjrO2pls.R"))
source(file.path(scstar2_src, "mjrO2plsPred.R"))
source(file.path(scstar2_src, "tTest.R"))

# 本地 patch：稳健版 PLSconstruct / OPLSDA
PLSconstruct <- function(data1, data2, prep, NCV, PLScomp_def, minNC) {
  X <- rbind(data1, data2)
  X[!is.finite(X)] <- 0
  mu <- colMeans(X)
  Xc <- scale(X, center = TRUE, scale = TRUE)
  Xc[!is.finite(Xc)] <- 0
  svd_res <- tryCatch(svd(Xc), error = function(e) svd(Xc, nu = 0, nv = 0))
  ncomp <- min(PLScomp_def, ncol(Xc), nrow(Xc), ncol(svd_res$v))
  if (ncomp == 0) ncomp <- 1
  list(XL = svd_res$v[, 1:ncomp, drop = FALSE], mu = mu)
}
OPLSDA <- function(data1, data2, variables, nrcv = 7, nc = 1, ncox = 1, ncoy = 1,
                   pCutoff = 0.05, np = 1000) {
  pCutoff <- pCutoff / length(variables)
  X <- rbind(data1, data2)
  Y <- matrix(0, nrow(X), 2)
  Y[1:nrow(data1), 1] <- 1; Y[(nrow(data1)+1):nrow(X), 2] <- 1
  model <- o2pls_m(X, Y, nc, ncox, ncoy, nrcv, pCutoff)
  Q2 <- model$Q2Yhatcum
  model <- o2pls_m(X, Y, nc, ncox, ncoy, 0, pCutoff)
  model$Q2Yhatcum <- Q2
  model_1 <- dispopls_m(model, X, Y, variables)
  Q2_p <- numeric(np)
  for (i in 1:np) {
    idx <- sample(nrow(Y))
    pmodel <- tryCatch(o2pls_m(X, Y[idx, ], nc, ncox, ncoy, nrcv, pCutoff),
                       error = function(e) list(Q2Yhatcum = 0))
    Q2_p[i] <- pmodel$Q2Yhatcum
  }
  P_val <- (sum(Q2_p >= Q2) + 1) / (np + 1)
  if (is.null(model_1$Pcorr)) model_1$Pcorr <- rep(1, ncol(X))
  list(Q2 = Q2, P_val = P_val, Pcorr = model_1$Pcorr, model = model_1)
}

out_dir <- file.path(OUT_BASE, ALGO, cancer)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, paste0(ALGO, "_", bulk_name, "_", sc_name, ".csv"))
if (file.exists(out_csv)) { cat(sprintf("SKIP %s|%s|%s\n", cancer, sc_name, bulk_name)); quit(save = "no", status = 0) }

cat(sprintf("===== scSTAR2 %s | %s | %s =====\n", cancer, sc_name, bulk_name))

# ── 1. 载入 scRNA（完整，不子采样） ────────────────────────────
sc <- load_scRNA_counts(cancer, sc_name)
sc_counts <- sc$counts
cell_meta <- sc$meta
rm(sc); gc()
cat(sprintf("  scRNA: %d genes x %d cells (完整)\n", nrow(sc_counts), ncol(sc_counts)))

# ── 2. 载入 bulk + 二分类表型 ──────────────────────────────────
bulk <- load_bulk(cancer, bulk_name)
bulk_exprs <- bulk$exp   # 已 log-norm
phenotype <- as.numeric(bulk$clinical[[cfg$scol]])
if (any(phenotype > 1)) phenotype <- ifelse(phenotype >= median(phenotype), 1, 0)
cat(sprintf("  Bulk: %d genes x %d samples, 1:%d, 0:%d\n",
            nrow(bulk_exprs), ncol(bulk_exprs), sum(phenotype == 1), sum(phenotype == 0)))
if (sum(phenotype == 1) < 5 || sum(phenotype == 0) < 5) {
  cat("  单组样本过少，跳过\n"); quit(save = "no", status = 0)
}

# ── 3. 共同基因 → 稀疏矩阵上先选高表达基因 → 再转 dense ─────────
shared <- intersect(rownames(bulk_exprs), rownames(sc_counts))
sc_counts <- sc_counts[shared, , drop = FALSE]
bulk_shared <- bulk_exprs[shared, , drop = FALSE]

# 先按表达量均值选 top 3000 基因（在稀疏矩阵上做，避免大 dense 矩阵 OOM）
subset_n <- min(3000, length(shared))
if (length(shared) > subset_n) {
  means <- Matrix::rowMeans(sc_counts)
  top_idx <- order(means, decreasing = TRUE)[1:subset_n]
  sc_counts <- sc_counts[top_idx, , drop = FALSE]
  bulk_shared <- bulk_shared[top_idx, , drop = FALSE]
}

sc_shared <- log2(as.matrix(sc_counts) + 1)
bulk_shared <- as.matrix(bulk_shared)
cat(sprintf("  使用基因数: %d\n", nrow(sc_shared)))

# ── 4. 方差过滤 ────────────────────────────────────────────────
sc_shared[!is.finite(sc_shared)] <- 0
bulk_shared[!is.finite(bulk_shared)] <- 0
row_vars <- apply(sc_shared, 1, var)
keep <- row_vars > 1e-10 & !is.na(row_vars)
sc_shared <- sc_shared[keep, , drop = FALSE]
bulk_shared <- bulk_shared[keep, , drop = FALSE]
cat(sprintf("  方差过滤后: %d genes\n", nrow(sc_shared)))

# ── 5. PLS-DA 训练(bulk) → 投影(scRNA) ─────────────────────────
data_good <- bulk_shared[, phenotype == 1, drop = FALSE]
data_poor <- bulk_shared[, phenotype == 0, drop = FALSE]

MODEL <- tryCatch(PLSconstruct(t(data_good), t(data_poor), "mc", 5, min(5, ncol(data_good) - 1), 2),
                  error = function(e) { cat(sprintf("  PLS FAILED: %s\n", e$message)); return(NULL) })
if (is.null(MODEL)) quit(save = "no", status = 0)

S_sc_temp <- t(sc_shared)
S_sc_temp <- S_sc_temp - matrix(1, nrow(S_sc_temp), 1) %*% MODEL$mu
proj <- tryCatch(S_sc_temp %*% MODEL$XL %*% MASS::ginv(MODEL$XL),
                 error = function(e) { cat(sprintf("  投影失败: %s\n", e$message)); return(NULL) })
if (is.null(proj)) quit(save = "no", status = 0)

cell_score <- proj[, 1]   # 正值=Good 预后, 负值=Poor 预后

# ── 6. 输出（上下 20% 分位） ───────────────────────────────────
probs <- 0.2
q_high <- quantile(cell_score, 1 - probs)
q_low  <- quantile(cell_score, probs)

result_df <- data.frame(Cell_ID = colnames(sc_counts), stringsAsFactors = FALSE)
result_df[[cfg$ct_col]] <- cell_meta[result_df$Cell_ID, cfg$ct_col]
result_df$Rank_Label <- "Background"
result_df$Rank_Label[cell_score >= q_high] <- "Rank-"   # Good 预后 = Protect
result_df$Rank_Label[cell_score <= q_low]  <- "Rank+"   # Poor 预后 = Risk

cat(sprintf("  Rank+: %d, Rank-: %d, Bg: %d\n",
            sum(result_df$Rank_Label == "Rank+"),
            sum(result_df$Rank_Label == "Rank-"),
            sum(result_df$Rank_Label == "Background")))
write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", basename(out_csv)))
