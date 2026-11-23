# =============================================================================
# EXAMINER GUIDE
# Script: 03_summarise_sbayesrc.R
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Reads SBayesRC outputs and summarises annotation-specific posterior/heritability results for downstream plotting.
#

#!/usr/bin/env Rscript
# Combine SBayesRC per-annotation SNP-heritability enrichment result files.
# This script reads existing model outputs; it does not refit SBayesRC.
source("R/utils.R")
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "data/sbayesrc_results"
out_dir <- if (length(args) >= 2) args[2] else "results/tables"
pattern <- if (length(args) >= 3) args[3] else "AnnoPerSnpHsqEnrichment"
files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (!length(files)) stop("No SBayesRC enrichment files matched in: ", input_dir)
results <- bind_rows(lapply(files, function(f) {
  x <- fread(f, data.table = FALSE, check.names = FALSE)
  assert_columns(x, c("Annotation", "Enrich"), basename(f))
  x$trait_file <- basename(f)
  x
}))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(results, file.path(out_dir, "sbayesrc_enrichment_all_traits.tsv"))
message("Combined ", length(files), " result files.")
