# scAB 评估（仅 Rank+，plus_only=TRUE）
.ca <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", .ca, value = TRUE)[1])), winslash = "/")
source(file.path(SCRIPT_DIR, "..", "_toolkit", "config.R"))
args <- commandArgs(trailingOnly = TRUE)
evaluate_cell_level("scAB", cancer = if (length(args) > 0) args[1] else NULL, plus_only = TRUE)
