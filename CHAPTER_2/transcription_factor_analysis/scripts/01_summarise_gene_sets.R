# =============================================================================
# EXAMINER GUIDE
# Script: 01_summarise_gene_sets.R
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Summarises transcription-factor target gene sets and their overlap with supplied gene lists.
#

#!/usr/bin/env Rscript
# Summarise unique genes per transcription-factor gene set and report overlaps.
source("R/utils.R")
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "data/gene_lists"
out_dir <- if (length(args) >= 2) args[2] else "results/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tf_files <- c(
  HOPX = "HOPX_bound_genes.txt",
  GATA4 = "genes_bound_gata4.txt",
  MEF2A = "genes_bound_mef2a.txt",
  NKX2_5 = "genes_bound_nkx25.txt",
  TBX5 = "genes_bound_tbx5.txt",
  SRF = "genes_bound_SRF.txt"
)
gene_sets <- lapply(tf_files, function(f) read_gene_list(file.path(input_dir, f)))
names(gene_sets) <- names(tf_files)

summary <- tibble(
  transcription_factor = names(gene_sets),
  n_unique_genes = vapply(gene_sets, length, integer(1))
)
write_tsv(summary, file.path(out_dir, "tf_gene_set_sizes.tsv"))

# Pairwise Jaccard similarity: |A intersect B| / |A union B|.
jaccard <- expand.grid(tf_1 = names(gene_sets), tf_2 = names(gene_sets),
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(
    intersection_n = length(intersect(gene_sets[[tf_1]], gene_sets[[tf_2]])),
    union_n = length(union(gene_sets[[tf_1]], gene_sets[[tf_2]])),
    jaccard = ifelse(union_n == 0, NA_real_, intersection_n / union_n)
  ) |>
  ungroup()
write_tsv(jaccard, file.path(out_dir, "tf_pairwise_jaccard.tsv"))

# Export a long membership table for downstream enrichment and auditing.
membership <- bind_rows(lapply(names(gene_sets), function(tf) {
  tibble(transcription_factor = tf, gene_id = gene_sets[[tf]])
}))
write_tsv(membership, file.path(out_dir, "tf_gene_membership.tsv"))
message("Wrote gene-set sizes, pairwise Jaccard values and membership table to ", out_dir)
