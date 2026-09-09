# MuSiC 评估（反卷积 / 样本级置换）
.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))
args <- commandArgs(trailingOnly = TRUE)
evaluate_deconv("MuSiC", cancer = if (length(args) > 0) args[1] else NULL)
