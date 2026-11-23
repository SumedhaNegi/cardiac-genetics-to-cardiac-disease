#!/usr/bin/env Rscript
# =============================================================================
# Historical GWAS summary-statistic formatting utilities
# =============================================================================

# CODE MAP
# Formats GWAS summary-statistic tables into SBayesRC/LDSC-style input columns; source column positions and sample sizes must be verified per phenotype.
#
# The supplied file combined PBS/AWK shell commands with several R fragments.
# This curated version expresses the visible transformations in R and exposes
# the source column selections explicitly. It does not infer missing GWAS
# metadata (sample size, allele frequency, genome build or effect direction).
#
# Usage:
#   Rscript Cleaning_sumstat.R awk_style input.tsv output.ma
#   Rscript Cleaning_sumstat.R sophie_style input.tsv output.ma [sample_N]
#
# `awk_style` preserves source columns 1,5,6,7,11,12,16,17 for rows whose first
# field begins with rs. `sophie_style` preserves columns 2,3,1,5,6,7,11,12,16,18
# and labels them as CHR, BP, SNP, A1, A2, FREQ, BETA, SE, P, N, as in the source.
# Confirm these positional mappings against the exact input format before use.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(paste(
    "Usage:",
    "Rscript Cleaning_sumstat.R <awk_style|sophie_style> <input.tsv> <output.ma> [sample_N]"
  ))
}
mode <- args[[1]]
input_file <- args[[2]]
output_file <- args[[3]]
if (!mode %in% c("awk_style", "sophie_style")) stop("Unknown mode: ", mode)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# Read all columns without modifying their names. The historical sources used
# tab-separated tables; character-first import avoids accidental ID coercion.
dat <- read.delim(input_file, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE, colClasses = "character")

if (mode == "awk_style") {
  selected_columns <- c(1L, 5L, 6L, 7L, 11L, 12L, 16L, 17L)
  if (ncol(dat) < max(selected_columns)) {
    stop("Input has ", ncol(dat), " columns; awk_style requires at least ",
         max(selected_columns), ".")
  }
  # Match the original AWK predicate: retain records with an rs-prefixed ID
  # in the first source column, then write the selected fields in source order.
  keep <- grepl("^rs", dat[[1L]], ignore.case = FALSE)
  out <- dat[keep, selected_columns, drop = FALSE]
} else {
  selected_columns <- c(2L, 3L, 1L, 5L, 6L, 7L, 11L, 12L, 16L, 18L)
  if (ncol(dat) < max(selected_columns)) {
    stop("Input has ", ncol(dat), " columns; sophie_style requires at least ",
         max(selected_columns), ".")
  }
  out <- dat[, selected_columns, drop = FALSE]
  names(out) <- c("CHR", "BP", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")
  out <- out[grepl("^rs", out$SNP), , drop = FALSE]
  # The historical follow-up replaced/added sample size with a fixed value.
  # Only do so when explicitly supplied, rather than silently assuming 41,135.
  if (length(args) >= 4L && nzchar(args[[4L]])) {
    sample_n <- suppressWarnings(as.numeric(args[[4L]]))
    if (!is.finite(sample_n) || sample_n <= 0) stop("sample_N must be positive.")
    out$N <- sample_n
  }
}

# Write a headered, tab-delimited output; row names are not part of the data.
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(out, output_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE, na = "")
message("Wrote ", nrow(out), " SNP records to ", output_file)
