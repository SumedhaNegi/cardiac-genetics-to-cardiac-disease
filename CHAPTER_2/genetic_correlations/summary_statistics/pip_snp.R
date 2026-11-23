#!/usr/bin/env Rscript
# =============================================================================
# Sort SBayesRC SNP-result tables by posterior inclusion probability (PIP)
# =============================================================================

# CODE MAP
# Orders SBayesRC SNP-result tables by the posterior inclusion probability field (historically column 10).
#
# The source workflow used Unix `sort -k10 -n` on column 10. This version uses
# R's order() so it is portable and does not mix shell commands into R code.
# Input files can be supplied as command-line arguments. If none are supplied,
# the historical filenames below are used when present in the working folder.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
input_files <- if (length(args)) args else c(
  "invnorm_RVEDV_LVEDV_ratio_cleaned.ma.snpRes",
  "invnorm_RVEF_cleaned.ma.snpRes",
  "invnorm_lvef_cleaned.ma.snpRes"
)
input_files <- input_files[file.exists(input_files)]
if (!length(input_files)) {
  stop("No input files found. Supply one or more SBayesRC .snpRes files.")
}

for (input_file in input_files) {
  # SNP-result files are tab-delimited and the historical PIP field is column 10.
  dat <- read.delim(input_file, header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(dat) < 10L) {
    warning("Skipping ", input_file, ": fewer than 10 columns.")
    next
  }
  pip <- suppressWarnings(as.numeric(dat[[10L]]))
  if (all(is.na(pip))) {
    warning("Skipping ", input_file, ": column 10 is not numeric PIP data.")
    next
  }
  # Ascending numeric sort matches the historical `sort -k10 -n` command.
  sorted <- dat[order(pip, na.last = TRUE), , drop = FALSE]
  output_file <- paste0(tools::file_path_sans_ext(input_file), ".sorted.tsv")
  write.table(sorted, output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
  message("Sorted PIP table written to: ", output_file)
}
