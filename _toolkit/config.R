# ============================================================
# TiPhD Benchmark — shared R configuration
# ------------------------------------------------------------
# 版本要求：R ≥ 4.4.0, Seurat ≥ 5.0.0, Matrix ≥ 1.6.4
# scRNA .rds 数据为 Seurat v5 Assay5 格式，R 4.1 + Seurat 4 无法正确读取。
# ------------------------------------------------------------
# 这是所有算法的 run_*.R 与 evaluate.R 共用的唯一配置源。
# 每个脚本开头都会 `source(../config.R)`，因此：
#   1) 数据路径、5 个癌种的参数、金标准名称映射只在这里维护一份；
#   2) 置换检验、金标准对比、指标计算的逻辑也只在这里实现一份。
#
# 约定（写入前请先了解）：
#   - 每个算法在 skills/<algo>/ 下有一个独立文件夹，内含
#     run_<algo>_pair.R（或 .py） / batch_run.py / evaluate.R。
#   - 运行脚本按 (cancer, sc_name, bulk_name) 三个命令行参数，单配对独立运行。
#   - 单细胞数据【不子采样】，读取完整细胞。
#   - 所有路径都从本文件位置推导（相对路径），不写死绝对路径。
# ============================================================

# ---- 定位 SCRIPT_DIR（仅供诊断） ----
.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
SCRIPT_DIR <- if (length(.f) == 0) normalizePath(getwd(), winslash = "/") else normalizePath(dirname(sub("^--file=", "", .f[1])), winslash = "/")

# ---- 数据与输出路径：环境变量优先，缺省 <cwd>/data 与 <cwd>/results ----
.tiphd_data <- Sys.getenv("TIPHD_DATA_DIR", unset = "")
.tiphd_out  <- Sys.getenv("TIPHD_OUT_DIR",  unset = "")
DATA_DIR <- normalizePath(
  if (nzchar(.tiphd_data)) .tiphd_data else file.path(getwd(), "data"),
  winslash = "/", mustWork = FALSE)
OUT_BASE <- normalizePath(
  if (nzchar(.tiphd_out)) .tiphd_out else file.path(getwd(), "results"),
  winslash = "/", mustWork = FALSE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# BASE_DIR 仅作向后兼容保留
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
})

# ============================================================
# 一、5 个癌种的统一配置（所有算法共用）
# ============================================================
CANCERS <- list(
  AML = list(
    ct_col = "CellType", tcol = "OS_time", scol = "OS_status",
    gold = "gold_standard_AML.csv",
    nmap = c("HSC-like" = "HSC-like", "Prog-like" = "Prog-like", "Mono-like" = "Mono-like"),
    scs = c("GSE116256"),
    bulks = c("TCGA", "wave12", "wave34")),
  CRC = list(
    ct_col = "Cell_subtype", tcol = "Time", scol = "Status",
    gold = "gold_standard_all_CRC.csv",
    nmap = c("T follicular helper cells (Tfh)" = "T follicular helper cells",
             "IgA+ Plasma cells" = "IgA+ Plasma",
             "SPP1+" = "SPP1+A",
             "CD8+ T cells" = "CD8+ T cells",
             "Regulatory T cells" = "Regulatory T cells",
             "Myofibroblasts" = "Myofibroblasts",
             "Stalk-like ECs" = "Stalk-like ECs",
             "T helper 17 cells" = "T helper 17 cells",
             "Tip-like ECs" = "Tip-like ECs"),
    scs = c("GSE144735", "GSE132465"),
    bulks = c("GSE14333", "GSE17536", "GSE33113", "GSE37892", "GSE39582")),
  HCC = list(
    ct_col = "subtype", tcol = "OS_time", scol = "OS_status",
    gold = "gold_standard_all_HCC.csv",
    nmap = c("CD4+ Tex" = "CD4_Tex", "CD4+ Th" = "CD4_Th", "Treg" = "Treg",
             "cDC1" = "cDC1_CLEC9A", "cDC2" = "cDC2_CD1C", "CD8+ Tex" = "CD8_Tex",
             "FOLR2+ TAMs (TAM1)" = "TAM_FOLR2", "IgA+ B cells" = "IgA+_Plasma_B_cells",
             "LAMP3+ DCs" = "LAMP3_DC", "NK cell" = "NK_cells",
             "PLVAP+ ECs" = "EC_PLVAP",
             "SPP1+ tumor-associated macrophages" = "TAM_SPP1"),
    scs = c("GSE149614", "GSE151530"),
    bulks = c("GSE116174", "GSE14520", "GSE76427")),
  LUAD = list(
    ct_col = "Cell_subtype", tcol = "OS_time", scol = "OS_status",
    gold = "gold_standard_all_LUDA.csv",
    nmap = c("CD8+ T cell" = "CD8+ T cell", "cDC1" = "cDC1", "cDC2" = "cDC2_CD1C",
             "Follicular B cells" = "Follicular B cells",
             "follicular helper T cell" = "follicular helper T cell",
             "IgG+ Plasma cells" = "IgG+ Plasma B cells", "Macro_ISG15" = "Macro_ISG15",
             "cDC_LAMP3 (Mature DCs)" = "cDC_LAMP3", "myCAF" = "myCAF",
             "Macro_PPARG" = "Macro_PPARG", "Macro_SPP1" = "Macro_SPP1", "Treg" = "Treg"),
    scs = c("GSE127465", "GSE131907", "GSE148071"),
    bulks = c("GSE31210", "GSE3141", "GSE37745", "GSE68465", "GSE72094")),
  GC = list(
    ct_col = "annotation", tcol = "Time", scol = "Status",
    gold = "gold_standard_all_GC.csv",
    nmap = c("eCAF" = "eCAF", "DC" = "DC",
             "exhausted CD8+ T cells" = "CD8-LAYN-exhausted",
             "iCAF" = "iCAF", "FOXP3+ treg cells" = "FOXP3-Treg",
             "Treg" = "FOXP3-Treg", "plasma cell" = "Plasma.cells", "CAFs" = "eCAF"),
    scs = c("GSE183904"),
    bulks = c("GSE15459", "GSE26253", "GSE26899", "GSE26901", "GSE28541",
              "GSE29272", "GSE34942", "GSE57303", "GSE66229", "GSETCGA"))
)

# 生成全部 (cancer, sc, bulk) 配对列表，供 batch_run 使用
PAIRS <- do.call(rbind, lapply(names(CANCERS), function(cn) {
  cfg <- CANCERS[[cn]]
  expand.grid(cancer = cn, sc = cfg$scs, bulk = cfg$bulks, stringsAsFactors = FALSE)
}))

# ============================================================
# 二、数据读取助手（相对路径，读取完整数据）
# ============================================================

# 读取 scRNA 的 counts 与 meta（R 侧统一读 .rds，完整细胞，不子采样）
load_scRNA_counts <- function(cancer, sc_name) {
  rds_file <- file.path(DATA_DIR, "scRNA", cancer, paste0(sc_name, ".rds"))
  if (!file.exists(rds_file)) stop("scRNA 文件不存在: ", rds_file)
  sc_obj <- readRDS(rds_file)
  # 兼容 Seurat v4/v5 .rds（数据可能存在多个位置）：
  #   - v4 Assay（v4 数据）：slot "@counts"
  #   - v5 Assay5（Seurat 5 正常读）：LayerData(obj, "RNA", "counts")
  #   - 兜底：Seurat 4 读 v5 .rds 时，Assay5 是空 S4（$ 不可用），数据在 attr(a, "layers")
  a <- sc_obj[["RNA"]]
  counts <- NULL
  if (isS4(a) && length(slotNames(a)) > 0 && "counts" %in% slotNames(a)) {
    # v4 Assay
    counts <- a@counts
  } else if (isS4(a) || is.null(a$counts)) {
    # 兜底：Seurat 4 读 v5 .rds（S4 无 slots，$ 不可用）→ attr/layers
    layers <- attr(a, "layers")
    if (!is.null(layers) && length(layers) > 0) counts <- layers[[1]]
  } else {
    # Seurat 5 正常路径：优先用 LayerData 拿到带 rownames 的矩阵
    if (requireNamespace("Seurat", quietly = TRUE) && packageVersion("Seurat") >= "5.0.0") {
      counts <- Seurat::LayerData(sc_obj, assay = "RNA", layer = "counts")
    } else {
      counts <- a$counts
    }
  }
  # 若 rownames 缺失（兜底路径下），从 features attr 或 meta 补全
  meta <- sc_obj@meta.data
  if (is.null(rownames(counts)) || length(rownames(counts)) == 0) {
    feat <- attr(a, "features")
    rn <- if (!is.null(feat)) rownames(feat) else NULL
    if (is.null(rn) || length(rn) != nrow(counts)) rn <- rownames(sc_obj)
    if (!is.null(rn) && length(rn) == nrow(counts)) rownames(counts) <- rn
  }
  # 若 colnames 缺失，从 meta rownames 补全
  if (is.null(colnames(counts)) || length(colnames(counts)) == 0) {
    cn <- colnames(sc_obj)
    if (is.null(cn) || length(cn) != ncol(counts)) cn <- rownames(meta)
    if (!is.null(cn) && length(cn) == ncol(counts)) colnames(counts) <- cn
  }
  if (is.null(counts)) stop("无法从 scRNA 对象提取 counts，Assay 类：", class(a))
  list(counts = counts, meta = meta)
}

# 读取 bulk 表达矩阵 + 临床信息，并按样本取交集
load_bulk <- function(cancer, bulk_name) {
  exp_file  <- file.path(DATA_DIR, "bulk", cancer, "exp", paste0(bulk_name, "_exp.csv"))
  clin_file <- file.path(DATA_DIR, "bulk", cancer, "clinical", paste0(bulk_name, "_clinical.csv"))
  bulk_exp <- as.matrix(read.csv(exp_file, row.names = 1))
  clinical <- read.csv(clin_file, row.names = 1)
  common <- intersect(colnames(bulk_exp), rownames(clinical))
  list(exp = bulk_exp[, common, drop = FALSE],
       clinical = clinical[common, , drop = FALSE])
}

# 从临床信息抽取 survival 两列 (time, status)
build_survival <- function(clinical, tcol, scol) {
  surv <- clinical[, c(tcol, scol), drop = FALSE]
  colnames(surv) <- c("time", "status")
  surv$time   <- as.numeric(surv$time)
  surv$status <- as.numeric(surv$status)
  surv$time[surv$time <= 0] <- 0.1
  surv
}

# 读取金标准，去重 + 规范 Phenotype 大小写
load_gold <- function(cancer) {
  gold <- read.csv(file.path(DATA_DIR, CANCERS[[cancer]]$gold), stringsAsFactors = FALSE)
  gold %>%
    distinct(Subtype, Phenotype) %>%
    mutate(Phenotype = tools::toTitleCase(tolower(Phenotype))) %>%
    filter(Phenotype %in% c("Protect", "Risk"), Subtype != "")
}

# ============================================================
# 三、置换检验（论文公式 3.1）
# ============================================================

# A 类（细胞级）：三分类 Rank+ / Rank- / Background
permute_cell_level <- function(df, ct_col, n_perm = 1000) {
  df_sub <- data.frame(
    CT  = factor(df[[ct_col]]),
    SCI = factor(df[["Rank_Label"]], levels = c("Rank+", "Rank-", "Background")))
  all_cts <- levels(df_sub$CT)

  obs_tbl <- table(df_sub$CT, df_sub$SCI)
  obs_stats <- as.data.frame.matrix(obs_tbl) %>%
    tibble::rownames_to_column(ct_col) %>%
    rename(Obs_Plus = `Rank+`, Obs_Minus = `Rank-`, Obs_BG = `Background`)

  perm_plus  <- matrix(0, nrow = length(all_cts), ncol = n_perm, dimnames = list(all_cts, NULL))
  perm_minus <- matrix(0, nrow = length(all_cts), ncol = n_perm, dimnames = list(all_cts, NULL))
  for (i in 1:n_perm) {
    shuffled <- sample(df_sub$SCI)
    tmp <- table(df_sub$CT, factor(shuffled, levels = levels(df_sub$SCI)))
    perm_plus[, i]  <- tmp[, "Rank+"]
    perm_minus[, i] <- tmp[, "Rank-"]
  }

  obs_stats %>%
    rowwise() %>%
    mutate(
      P_val_Plus  = (sum(perm_plus[.data[[ct_col]], ]  >= Obs_Plus)  + 1) / (n_perm + 1),
      P_val_Minus = (sum(perm_minus[.data[[ct_col]], ] >= Obs_Minus) + 1) / (n_perm + 1),
      Rank_Label  = case_when(
        P_val_Plus  < 0.05 ~ "Significantly_Positive",
        P_val_Minus < 0.05 ~ "Significantly_Negative",
        TRUE                ~ "Non_Significant")) %>%
    ungroup()
}

# A 类变体（scAB）：只有 Rank+ / Background，无 Rank-
permute_cell_level_plus <- function(df, ct_col, n_perm = 1000) {
  df_sub <- data.frame(
    CT  = factor(df[[ct_col]]),
    SCI = factor(df[["Rank_Label"]], levels = c("Rank+", "Background")))
  all_cts <- levels(df_sub$CT)

  obs_tbl <- table(df_sub$CT, df_sub$SCI)
  obs_stats <- as.data.frame.matrix(obs_tbl) %>%
    tibble::rownames_to_column(ct_col) %>%
    rename(Obs_Plus = `Rank+`, Obs_BG = `Background`)

  perm_plus <- matrix(0, nrow = length(all_cts), ncol = n_perm, dimnames = list(all_cts, NULL))
  for (i in 1:n_perm) {
    shuffled <- sample(df_sub$SCI)
    perm_plus[, i] <- table(df_sub$CT, factor(shuffled, levels = levels(df_sub$SCI)))[, "Rank+"]
  }

  obs_stats %>%
    rowwise() %>%
    mutate(
      P_val_Plus = (sum(perm_plus[.data[[ct_col]], ] >= Obs_Plus) + 1) / (n_perm + 1),
      Rank_Label = ifelse(P_val_Plus < 0.05, "Significantly_Positive", "Non_Significant")) %>%
    ungroup()
}

# B 类（样本级）：打乱样本 survival，对每种细胞类型的 proportion 做 Cox 置换
permute_sample_level <- function(prop_df, surv, n_perm = 1000) {
  results <- data.frame(stringsAsFactors = FALSE)
  for (ct in colnames(prop_df)) {
    prop <- prop_df[[ct]]
    if (sd(prop) < 1e-10) next

    cox_obs <- tryCatch(coxph(Surv(time, status) ~ prop, data = surv), error = function(e) NULL)
    if (is.null(cox_obs) || nrow(summary(cox_obs)$coefficients) == 0) next
    coef_obs <- summary(cox_obs)$coefficients["prop", "coef"]
    # 处理 NA：若 coef_obs 为 NA 则跳过该细胞类型
    if (is.na(coef_obs)) next

    null_coefs <- numeric(n_perm)
    for (i in 1:n_perm) {
      si <- sample(nrow(surv))
      surv_s <- data.frame(time = surv$time[si], status = surv$status[si])
      cox_null <- tryCatch(coxph(Surv(time, status) ~ prop, data = surv_s), error = function(e) NULL)
      null_coefs[i] <- if (!is.null(cox_null) && nrow(summary(cox_null)$coefficients) > 0)
        summary(cox_null)$coefficients["prop", "coef"] else 0
    }

    P_pos <- (sum(null_coefs >= coef_obs) + 1) / (n_perm + 1)
    P_neg <- (sum(null_coefs <= coef_obs) + 1) / (n_perm + 1)
    P_pos <- if (is.na(P_pos)) 1 else P_pos
    P_neg <- if (is.na(P_neg)) 1 else P_neg

    Rank_Label <- if (!is.na(coef_obs) && coef_obs > 0 && !is.na(P_pos) && P_pos < 0.05) "Significantly_Positive"
                  else if (!is.na(coef_obs) && coef_obs < 0 && !is.na(P_neg) && P_neg < 0.05) "Significantly_Negative"
                  else "Non_Significant"

    results <- rbind(results, data.frame(
      cell_type = ct, Obs_Coef = round(coef_obs, 6),
      P_val_Plus = round(P_pos, 6), P_val_Minus = round(P_neg, 6),
      Rank_Label = Rank_Label, stringsAsFactors = FALSE))
  }
  results
}

# ============================================================
# 四、金标准对比（TP / FP / FN）与指标计算
# ============================================================

# perm_result: 置换检验输出（首列为细胞类型名，含 Rank_Label 列）
# plus_only:   是否只允许 Risk 方向预测（scAB）
compare_gold <- function(perm_result, gold_unique, nmap, plus_only = FALSE) {
  ct_col <- colnames(perm_result)[1]
  gold_mapped <- gold_unique %>%
    mutate(mapped = dplyr::recode(Subtype, !!!nmap)) %>%
    filter(!is.na(mapped)) %>%
    distinct(mapped, Phenotype)          # 处理 GC 中 Treg/FOXP3 同名合并

  comparison <- perm_result %>%
    rename(cell_type = all_of(ct_col)) %>%
    inner_join(gold_mapped, by = c("cell_type" = "mapped")) %>%
    mutate(
      Prediction = if (plus_only) {
        ifelse(Rank_Label == "Significantly_Positive", "Risk", "Non_Sig")
      } else {
        case_when(Rank_Label == "Significantly_Positive" ~ "Risk",
                  Rank_Label == "Significantly_Negative" ~ "Protect",
                  TRUE ~ "Non_Sig")
      },
      Verdict = case_when(
        Prediction == "Non_Sig" ~ "FN",
        Prediction == Phenotype ~ "TP",
        Prediction != Phenotype ~ "FP",
        TRUE ~ "Unknown"))

  TP <- sum(comparison$Verdict == "TP")
  FP <- sum(comparison$Verdict == "FP")
  FN <- sum(comparison$Verdict == "FN")
  Sp <- nrow(comparison)

  data.frame(Sp = Sp, TP = TP, FP = FP, FN = FN,
             Precision  = ifelse(TP + FP > 0, TP / (TP + FP), NA),
             Coverage   = ifelse(Sp > 0, TP / Sp, NA),
             False_Rate = ifelse(TP + FP > 0, FP / (TP + FP), NA),
             stringsAsFactors = FALSE)
}

# ============================================================
# 五、评估驱动（A 类：细胞级）
#   用法：Rscript evaluate.R [cancer]
# ============================================================
evaluate_cell_level <- function(algo, cancer = NULL, plus_only = FALSE) {
  cancers <- if (is.null(cancer)) names(CANCERS) else cancer
  for (cn in cancers) {
    cfg <- CANCERS[[cn]]
    result_dir <- file.path(OUT_BASE, algo, cn)
    gold_unique <- load_gold(cn)

    cat(sprintf("\n========================================\n"))
    cat(sprintf("Evaluating: %s — %s\n", algo, cn))
    cat(sprintf("========================================\n"))

    all_metrics <- list()
    missing <- 0

    for (sc_name in cfg$scs) {
      for (bulk_name in cfg$bulks) {
        csv_file <- file.path(result_dir, paste0(algo, "_", bulk_name, "_", sc_name, ".csv"))
        if (!file.exists(csv_file)) { missing <- missing + 1; next }

        df <- read.csv(csv_file, stringsAsFactors = FALSE)
        perm_result <- if (plus_only) permute_cell_level_plus(df, cfg$ct_col)
                       else permute_cell_level(df, cfg$ct_col)

        perm_csv <- file.path(result_dir, paste0(algo, "_permutation_", bulk_name, "_", sc_name, ".csv"))
        write.csv(perm_result, perm_csv, row.names = FALSE)

        metrics <- compare_gold(perm_result, gold_unique, cfg$nmap, plus_only = plus_only)
        cat(sprintf("  %s | %s: Sp=%d, TP=%d, FP=%d, FN=%d → Prec=%.3f, Cov=%.3f, FR=%.3f\n",
                    sc_name, bulk_name, metrics$Sp, metrics$TP, metrics$FP, metrics$FN,
                    metrics$Precision, metrics$Coverage, metrics$False_Rate))

        all_metrics[[length(all_metrics) + 1]] <-
          data.frame(scRNA = sc_name, Bulk = bulk_name, metrics, stringsAsFactors = FALSE)
      }
    }

    if (length(all_metrics) > 0) {
      final_table <- bind_rows(all_metrics) %>%
        mutate(across(c(Precision, Coverage, False_Rate), ~ round(., 4)))
      write.csv(final_table, file.path(result_dir, paste0(algo, "_", cn, "_final_metrics.csv")),
                row.names = FALSE)
      cat(sprintf("✓ %s\n", file.path(result_dir, paste0(algo, "_", cn, "_final_metrics.csv"))))
    }
    if (missing > 0) cat(sprintf("  (缺失 %d 个配对输出)\n", missing))
  }
  cat("\nAll evaluations complete!\n")
}

# ============================================================
# 六、评估驱动（B 类：反卷积 / 样本级）
#   运行脚本只做反卷积并输出 *_proportions.csv；
#   这里读 proportions + clinical → 样本级置换 → 金标准对比。
# ============================================================
evaluate_deconv <- function(algo, cancer = NULL) {
  cancers <- if (is.null(cancer)) names(CANCERS) else cancer
  for (cn in cancers) {
    cfg <- CANCERS[[cn]]
    result_dir <- file.path(OUT_BASE, algo, cn)
    gold_unique <- load_gold(cn)

    cat(sprintf("\n========================================\n"))
    cat(sprintf("Evaluating: %s — %s\n", algo, cn))
    cat(sprintf("========================================\n"))

    all_metrics <- list()
    missing <- 0

    for (sc_name in cfg$scs) {
      for (bulk_name in cfg$bulks) {
        prop_csv <- file.path(result_dir, paste0(algo, "_", bulk_name, "_", sc_name, "_proportions.csv"))
        if (!file.exists(prop_csv)) { missing <- missing + 1; next }

        prop_df <- read.csv(prop_csv, row.names = 1, check.names = FALSE)
        bulk <- load_bulk(cn, bulk_name)
        surv <- build_survival(bulk$clinical, cfg$tcol, cfg$scol)

        perm_result <- permute_sample_level(prop_df, surv)

        out_csv <- file.path(result_dir, paste0(algo, "_", bulk_name, "_", sc_name, ".csv"))
        write.csv(perm_result, out_csv, row.names = FALSE)

        if (nrow(perm_result) == 0) next
        metrics <- compare_gold(perm_result, gold_unique, cfg$nmap)
        cat(sprintf("  %s | %s: Sp=%d, TP=%d, FP=%d, FN=%d → Prec=%.3f, Cov=%.3f, FR=%.3f\n",
                    sc_name, bulk_name, metrics$Sp, metrics$TP, metrics$FP, metrics$FN,
                    metrics$Precision, metrics$Coverage, metrics$False_Rate))

        all_metrics[[length(all_metrics) + 1]] <-
          data.frame(scRNA = sc_name, Bulk = bulk_name, metrics, stringsAsFactors = FALSE)
      }
    }

    if (length(all_metrics) > 0) {
      final_table <- bind_rows(all_metrics) %>%
        mutate(across(c(Precision, Coverage, False_Rate), ~ round(., 4)))
      write.csv(final_table, file.path(result_dir, paste0(algo, "_", cn, "_final_metrics.csv")),
                row.names = FALSE)
      cat(sprintf("✓ %s\n", file.path(result_dir, paste0(algo, "_", cn, "_final_metrics.csv"))))
    }
    if (missing > 0) cat(sprintf("  (缺失 %d 个配对 proportions)\n", missing))
  }
  cat("\nAll evaluations complete!\n")
}
