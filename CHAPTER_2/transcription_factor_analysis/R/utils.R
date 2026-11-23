# =============================================================================
# EXAMINER GUIDE
# Script: utils.R
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Shared helper functions for reading gene-coordinate tables, validating required columns and supporting the transcription-factor annotation workflow.
#

# Shared utilities for the cardiac transcription-factor heritability workflow.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x

project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE))
  }
  normalizePath(".", mustWork = FALSE)
}

read_gene_list <- function(path, gene_col = NULL) {
  stopifnot(file.exists(path))
  x <- fread(path, header = TRUE, data.table = FALSE, check.names = FALSE)
  if (ncol(x) == 1L) {
    genes <- x[[1]]
  } else {
    gene_col <- gene_col %||% names(x)[1]
    if (!gene_col %in% names(x)) stop("Gene column '", gene_col, "' not found in ", path)
    genes <- x[[gene_col]]
  }
  genes <- str_trim(as.character(genes))
  sort(unique(genes[!is.na(genes) & nzchar(genes)]))
}

read_gene_coordinates <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  required <- c("GeneID", "chr", "start", "end")
  if (!all(required %in% names(x))) {
    stop("Coordinate file must contain: ", paste(required, collapse = ", "), " (", path, ")")
  }
  x <- x |>
    transmute(
      gene_id = str_to_lower(str_trim(as.character(GeneID))),
      chr = str_remove(as.character(chr), "^chr"),
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end))
    ) |>
    filter(!is.na(gene_id), gene_id != "", !is.na(chr), !is.na(start), !is.na(end)) |>
    mutate(
      # Gene annotations can be stored in transcriptional orientation (start > end).
      # Genomic intervals are normalised to ascending coordinates for overlap tests.
      interval_start = pmin(start, end),
      interval_end = pmax(start, end)
    ) |>
    distinct(gene_id, chr, interval_start, interval_end)
  x
}

assert_columns <- function(x, required, label = "input") {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path, na = "")
}
