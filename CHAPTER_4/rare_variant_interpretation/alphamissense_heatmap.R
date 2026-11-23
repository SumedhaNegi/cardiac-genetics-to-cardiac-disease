#!/usr/bin/env Rscript
# =============================================================================
# HOPX AlphaMissense amino-acid substitution heatmap
# =============================================================================

# CODE MAP
# Parses HOPX amino-acid substitution labels and plots AlphaMissense scores, with an optional pre-matched ClinVar classification overlay.
#
# The source fragment expected an object named `data` with V2 containing a
# protein substitution string (reference amino acid + residue position + mutant
# amino acid) and V3 containing the AlphaMissense score. A later plotting
# fragment expected a separately created `merged` object containing optional
# ClinVar Germline_classification annotations. The missing merge step cannot be
# reconstructed without the original ClinVar input and variant-matching rule.
#
# This cleaned script creates the AlphaMissense heatmap from a supplied TSV.
# Required columns: V2, V3. Optional column: Germline_classification.
# If classification is present, points are overlaid on the heatmap; this is a
# display annotation only and does not establish clinical pathogenicity.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1]] else "hopx_alphamissense_variants.tsv"
out_file <- if (length(args) >= 2L) args[[2]] else "hopx_alphamissense_heatmap.pdf"
if (!file.exists(input_file)) stop("Input variant table not found: ", input_file)

# Read the variant-level table and verify the columns used in the source code.
data <- read_tsv(input_file, show_col_types = FALSE, name_repair = "minimal")
required <- c("V2", "V3")
if (!all(required %in% names(data))) {
  stop("Input must contain columns named V2 (protein substitution) and V3 (score).")
}

# Split a substitution string into reference residue, position and alternate
# residue. The source convention assumes a one-letter amino-acid code at each
# end (for example R123W); malformed strings are excluded with a warning.
substitution <- as.character(data$V2)
valid <- !is.na(substitution) & grepl("^[A-Za-z][0-9]+[A-Za-z]$", substitution)
if (any(!valid)) warning(sum(!valid), " rows have unsupported V2 substitution strings and were excluded.")
data_split <- data[valid, , drop = FALSE] |>
  mutate(
    Reference_amino_acid = substr(as.character(V2), 1L, 1L),
    Residue_position = as.integer(sub("^[A-Za-z]([0-9]+)[A-Za-z]$", "\\1", as.character(V2))),
    mutated_amino_acid = substr(as.character(V2), nchar(as.character(V2)), nchar(as.character(V2))),
    Score = suppressWarnings(as.numeric(V3)),
    Reference_V2 = paste(Reference_amino_acid, Residue_position)
  ) |>
  filter(is.finite(Score))
if (nrow(data_split) == 0L) stop("No valid substitutions with numeric AlphaMissense scores were found.")

# Preserve the residue order from the input rather than alphabetically sorting.
data_split$Reference_V2 <- factor(data_split$Reference_V2,
                                  levels = unique(data_split$Reference_V2))

# Draw one tile per reference-position / alternate-residue observation. Multiple
# records for the same cell are retained as supplied; no score aggregation is
# introduced because the original scientific aggregation rule is unknown.
p <- ggplot(data_split, aes(x = Reference_V2, y = mutated_amino_acid, fill = Score)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       name = "AlphaMissense\nscore") +
  labs(title = "HOPX AlphaMissense amino-acid substitutions",
       x = "Reference amino acid and residue position",
       y = "Mutated amino acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Optional ClinVar overlay. It is only available if the input already contains
# a matched classification column; this script does not fabricate a merge.
if ("Germline_classification" %in% names(data_split)) {
  clinvar_rows <- data_split |>
    filter(!is.na(Germline_classification), nzchar(as.character(Germline_classification)))
  if (nrow(clinvar_rows) > 0L) {
    p <- p + geom_point(data = clinvar_rows,
                        aes(shape = Germline_classification),
                        colour = "black", size = 2.5,
                        inherit.aes = TRUE)
  }
}

ggsave(out_file, plot = p, width = 14, height = 8)
message("Heatmap written to: ", normalizePath(out_file, mustWork = FALSE))
