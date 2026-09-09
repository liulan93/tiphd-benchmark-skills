# ============================================================
# TiPhD Benchmark — R 依赖一键安装脚本
# 用法: Rscript install_R_packages.R
# 作用:
#   1) 安装 CRAN/Bioconductor 依赖（Seurat 等）
#   2) 从 third_party/R/ 源码安装自定义 R 包
# 注意:
#   - scPAS 和 Scissor 含 Rcpp C++ 代码，Windows 下需先装 Rtools
#   - 在高核心数服务器上，preprocessCore 默认的 OpenMP 多线程可能触发
#     `ERROR; return code from pthread_create() is 22`。此时需要用
#     --disable-threading 重装 preprocessCore（见下方逻辑）。
# ============================================================

# 定位本脚本所在目录（third_party 相对本目录）
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
HERE <- if (is.na(self)) getwd() else normalizePath(dirname(self), winslash = "/")

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 1. CRAN 依赖 ----
cran_pkgs <- c("Seurat", "Matrix", "MASS", "survival", "Rcpp", "RcppEigen",
               "diptest", "multimode", "progress", "pROC", "knitr", "rmarkdown",
               "dplyr", "tidyr", "tibble", "AUCell", "nnls", "pls",
               "GSEABase", "matrixStats", "glmnet", "ordinalNet", "preprocessCore")
# 上面 preprocessCore 是 Bioconductor 包，列出无害；BiocManager 安装逻辑见步骤 2
cran_pkgs <- setdiff(cran_pkgs, "preprocessCore")
missing <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) {
  cat("安装 CRAN 依赖:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
}

# ---- 2. Bioconductor 依赖 ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("preprocessCore", "SingleCellExperiment", "AUCell", "GSEABase")
for (bp in bioc_pkgs) {
  if (!requireNamespace(bp, quietly = TRUE)) {
    cat(sprintf("安装 Bioconductor 包 %s ...\n", bp))
    BiocManager::install(bp, ask = FALSE, update = FALSE)
  }
}

# ---- 2b. preprocessCore 线程兼容性修复 ----
# 在核心数很多（>100）的服务器上，preprocessCore 的 OpenMP 代码会尝试创建
# 与物理核心数等量的线程，可能超出进程的线程 ulimit，导致 quantile
# normalization 报 `pthread_create() is 22`。安装时禁用线程可避免该问题。
# 这里做一个运行时快速检测：在一个小矩阵上调用 normalize.quantiles，
# 若失败则以 --disable-threading 重装。
test_preprocessCore <- function() {
  ok <- tryCatch({
    library(preprocessCore)
    m <- matrix(rnorm(2000), nrow = 100)
    normalize.quantiles(m)
    TRUE
  }, error = function(e) {
    grepl("pthread_create", conditionMessage(e)) || TRUE
  })
  # 严格地：只要报错就重装禁用线程版本
  tryCatch({
    m <- matrix(rnorm(2000), nrow = 100)
    preprocessCore::normalize.quantiles(m)
    TRUE
  }, error = function(e) FALSE)
}
if (!test_preprocessCore()) {
  cat("检测到 preprocessCore 线程错误，以 --disable-threading 重装...\n")
  BiocManager::install("preprocessCore", ask = FALSE, update = FALSE, force = TRUE,
                        configure.args = c(preprocessCore = "--disable-threading"))
}

# ---- 3. 自定义 R 包（源码在 third_party/R/） ----
# scPAS / Scissor 的 RcppExports.cpp 引用 inst/include/<pkg>.h 占位头文件，
# 这些头文件已随源码提供（若自行从上游获取缺失，需创建空占位头文件）。
local_pkgs <- c("scAB", "SCIPAC", "scPAS", "Scissor", "PIPET")
for (p in local_pkgs) {
  pkg_dir <- file.path(HERE, "third_party", "R", p)
  if (!dir.exists(pkg_dir)) { cat("跳过（目录不存在）:", p, "\n"); next }
  if (requireNamespace(p, quietly = TRUE)) next
  cat(sprintf("安装本地包 %s ...\n", p))
  tryCatch(
    install.packages(pkg_dir, repos = NULL, type = "source"),
    error = function(e) cat(sprintf("  ! %s 安装失败: %s\n", p, conditionMessage(e)))
  )
}

cat("\n全部完成！可用 library(scAB) / library(scPAS) / ... 验证。\n")
cat("提示: 高核心服务器若仍遇到 pthread 错误，可在运行前设置\n")
cat("      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n")
