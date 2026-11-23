#!/usr/bin/env Rscript
# =============================================================================
# CardioExcyte raw-export processing
# =============================================================================

# CODE MAP
# Parses CardioExcyte instrument exports, joins wells to the experiment-specific condition map, and creates tidy measurement tables and exploratory sweep plots.
#
# Purpose
#   Convert the tab-delimited CardioExcyte export into a tidy table with one
#   row per well, sweep and measurement parameter, then attach experimental
#   condition labels and produce exploratory plots.
#
# Input
#   run1.txt (or a path supplied as the first command-line argument).
#   The export layout assumed here stores sweep, parameter and elapsed-time
#   metadata in rows 5, 6 and 7, respectively; well measurements begin at row 10.
#
# Output
#   cardioexcyte_tidy.tsv and two exploratory PDF figures in the output folder.
#
# Important
#   The well-to-condition map below is experiment-specific. Confirm it against
#   the plate map before using the output in a thesis or downstream analysis.
#   This script deliberately does not install packages or submit HPC jobs.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# Read input/output paths from the command line, with convenient local defaults.
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "run1.txt"
out_dir <- if (length(args) >= 2) args[[2]] else "results/cardioexcyte"
if (!file.exists(input_file)) stop("CardioExcyte export not found: ", input_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the export as text first. Instrument exports can contain mixed metadata
# and numeric measurement fields, so conversion is performed after row parsing.
raw <- read.delim(input_file, header = FALSE, sep = "\t", fill = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character")
if (nrow(raw) < 10L || ncol(raw) < 2L) {
  stop("Unexpected CardioExcyte export layout: expected >=10 rows and >=2 columns.")
}

# Extract the three metadata rows. Their first cell is a row label, not data.
metadata_rows <- list(
  sweep = unlist(raw[5, -1, drop = FALSE], use.names = FALSE),
  parameter = unlist(raw[6, -1, drop = FALSE], use.names = FALSE),
  time = unlist(raw[7, -1, drop = FALSE], use.names = FALSE)
)
metadata <- tibble(
  col_index = seq_along(metadata_rows$sweep),
  Sweep = str_remove(as.character(metadata_rows$sweep), "^Sweep\\s*"),
  Parameter = as.character(metadata_rows$parameter),
  Time_of_Sweep = suppressWarnings(as.numeric(metadata_rows$time))
)

# Measurement rows start after the first nine export/header rows. Retain only
# plausible well identifiers (A1-H12 style, including zero-padded exports).
measurement_block <- raw[-seq_len(9L), , drop = FALSE]
well_ids <- str_trim(as.character(measurement_block[[1]]))
valid_well <- !is.na(well_ids) & grepl("^[A-H]0?[1-9]$|^[A-H]1[0-2]$", well_ids)
if (!any(valid_well)) {
  stop("No well identifiers found in column 1 after the first nine rows.")
}
well_ids <- str_replace(well_ids[valid_well], "^([A-H])0([1-9])$", "\\1\\2")
measurement_values <- measurement_block[valid_well, -1, drop = FALSE]
measurement_values[] <- lapply(measurement_values, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

# Ensure the export's measurement columns align with its metadata columns.
if (ncol(measurement_values) != nrow(metadata)) {
  stop("Measurement columns (", ncol(measurement_values),
       ") do not match metadata columns (", nrow(metadata),"). Check export rows 5-7.")
}
colnames(measurement_values) <- as.character(metadata$col_index)
measurement_values$Well_ID <- well_ids

# Convert the wide instrument matrix to tidy long form and attach sweep metadata.
values_long <- measurement_values |>
  pivot_longer(cols = -Well_ID, names_to = "col_index", values_to = "Value") |>
  mutate(col_index = as.integer(col_index)) |>
  left_join(metadata, by = "col_index") |>
  select(Well_ID, Sweep, Parameter, Time_of_Sweep, Value) |>
  filter(!is.na(Parameter), nzchar(Parameter))

# EXPERIMENT-SPECIFIC PLATE MAP: verify these assignments against the source
# plate layout. Each comma-separated pair is treated as replicate wells.
conditions_df <- tribble(
  ~Condition, ~Wells,
  "ZTT", "A1, A2", "VDW", "A7, A8",
  "BFP", "B1, B2", "YHH", "B7, B8",
  "GRP", "C1, C2", "BKH", "C7, C8",
  "UNR", "D1, D2", "YCP", "D7, D8",
  "RYU", "E1, E2", "JNJ", "E7, E8",
  "MTG", "F1, F2"
)
well_condition_map <- conditions_df |>
  separate_rows(Wells, sep = ",\\s*") |>
  rename(Well_ID = Wells)

# An inner join intentionally drops wells not represented in the plate map.
values_with_condition <- values_long |>
  inner_join(well_condition_map, by = "Well_ID")
if (nrow(values_with_condition) == 0L) {
  stop("No measurements matched the configured well-to-condition map.")
}

# Save the tidy data so later analysis does not depend on an interactive R session.
write_tsv(values_with_condition, file.path(out_dir, "cardioexcyte_tidy.tsv"))

# Plot the selected sweep as a condition comparison, faceted by parameter.
# The historical source called this Sweep 12 and therefore uses identifier 012.
sweep12_data <- values_with_condition |>
  filter(as.character(Sweep) %in% c("012", "12"))
if (nrow(sweep12_data) > 0L) {
  p_sweep12 <- ggplot(sweep12_data, aes(x = Condition, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    facet_wrap(~Parameter, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Sweep 12: values by condition and parameter",
         x = "Condition", y = "Measurement")
  ggsave(file.path(out_dir, "sweep12_by_condition.pdf"), p_sweep12,
         width = 12, height = 8)
} else {
  warning("Sweep 012/12 was not present; the selected-sweep plot was skipped.")
}

# Plot each well's trajectory across sweeps. Lines connect repeated measurements
# from the same well and condition; facets keep measurement scales separate.
values_with_condition <- values_with_condition |>
  mutate(Sweep_numeric = suppressWarnings(as.numeric(Sweep)))
p_all <- ggplot(values_with_condition,
                aes(x = Sweep_numeric, y = Value, colour = Condition,
                    group = interaction(Condition, Well_ID))) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_minimal() +
  labs(title = "CardioExcyte parameter trends across sweeps",
       x = "Sweep", y = "Measurement") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "all_sweeps_by_condition.pdf"), p_all,
       width = 12, height = 8)

message("Tidy data and plots written to: ", normalizePath(out_dir))
