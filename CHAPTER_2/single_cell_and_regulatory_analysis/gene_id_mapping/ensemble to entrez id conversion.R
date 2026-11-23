#!/usr/bin/env Rscript
# =============================================================================
# Ensembl-to-Entrez conversion for HOPX-bound gene lists
# =============================================================================

# CODE MAP
# Converts HOPX-bound human Ensembl gene identifiers to Entrez IDs for up- and down-regulated gene sets.
#
# Reads up- and down-regulated gene lists from the two named Excel worksheets,
# maps human Ensembl gene IDs to Entrez IDs, and writes separate tab-delimited
# outputs. This script assumes the first worksheet column contains Ensembl IDs.
# Package installation is intentionally not performed inside the analysis
# script; install readxl, AnnotationDbi and org.Hs.eg.db in the R environment.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly = TRUE)
input_xlsx <- if (length(args) >= 1L) args[[1]] else "Desktop/HOPX_bound_genes_TFs.xlsx"
out_dir <- if (length(args) >= 2L) args[[2]] else "Desktop/dry lab results"
if (!file.exists(input_xlsx)) stop("Excel input not found: ", input_xlsx)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the two source gene sets without assuming that the spreadsheet has headers.
upregulated <- read_xlsx(input_xlsx, sheet = "HOPX bound upregulated genes",
                         col_names = FALSE)
downregulated <- read_xlsx(input_xlsx, sheet = "HOPX bound downregulated genes",
                           col_names = FALSE)

convert_ensembl <- function(x, label) {
  if (ncol(x) < 1L) stop(label, " worksheet has no columns.")
  ids <- unique(trimws(as.character(x[[1]])))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  # Remove Ensembl version suffixes (e.g. ENSG... .4) for OrgDb matching while
  # retaining the original input identifier in the output table.
  lookup_ids <- sub("\\.[0-9]+$", "", ids)
  mapped <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = unique(lookup_ids),
                                  keytype = "ENSEMBL",
                                  columns = c("ENTREZID", "SYMBOL"))
  mapped$ENSEMBL_INPUT <- ids[match(mapped$ENSEMBL, lookup_ids)]
  mapped[, c("ENSEMBL_INPUT", "ENSEMBL", "ENTREZID", "SYMBOL")]
}

up_ids <- convert_ensembl(upregulated, "Upregulated")
down_ids <- convert_ensembl(downregulated, "Downregulated")
write.table(up_ids, file.path(out_dir, "hopx_upregulated_entrezid.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(down_ids, file.path(out_dir, "hopx_downregulated_entrezid.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote Ensembl-to-Entrez mappings to: ", normalizePath(out_dir))

# Optional Enrichr connectivity check retained as a separate, non-essential step.
# It is not required for identifier conversion and is not run automatically.
# if (requireNamespace("enrichR", quietly = TRUE)) {
#   dbs <- enrichR::listEnrichrDbs()
#   if (is.null(dbs)) message("Enrichr service is unavailable.") else print(utils::head(dbs))
# }
