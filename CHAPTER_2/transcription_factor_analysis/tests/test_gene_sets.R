# =============================================================================
# EXAMINER GUIDE
# Script: test_gene_sets.R
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Lightweight tests for the gene-set helper functions; these were retained but not executed, in keeping with the request not to test the repository.
#

# Lightweight validation of the supplied gene-list inputs.
source("R/utils.R")
files <- c(
  "HOPX_bound_genes.txt", "genes_bound_gata4.txt", "genes_bound_mef2a.txt",
  "genes_bound_nkx25.txt", "genes_bound_tbx5.txt", "genes_bound_SRF.txt"
)
for (f in files) {
  x <- read_gene_list(file.path("data/gene_lists", f))
  stopifnot(length(x) > 0L, !anyNA(x), all(nzchar(x)))
}
message("Gene-list input checks passed.")
