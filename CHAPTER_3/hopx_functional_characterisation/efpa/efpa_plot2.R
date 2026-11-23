# =============================================================================
# EXAMINER GUIDE
# Script: efpa_plot2.R
# Purpose: Plots extracellular field-potential action (EFPA) mean-beat waveforms from CardioExcyte-derived CSV files.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# Processes EFPA measurements and generates exploratory/summary plots for the HOPX functional-characterisation experiments.
#

# ============================================================
# EFPA MEAN-BEAT WAVEFORM PLOTTING
#
# Expected folder contents:
#   Baseline_Sweep006_mean_beats.csv
#   PostNA_Sweep009_mean_beats.csv
#   Baseline_Sweep006_parameters.csv
#   PostNA_Sweep009_parameters.csv
#
# The script automatically identifies mean-beat CSVs using the
# presence of a Time column. Parameter-summary CSVs are ignored.
#
# Outputs:
#   1. Individual well traces
#   2. Baseline vs post-NA overlays for each well
#   3. Variant mean traces with SEM
#   4. Faceted variant-screen plots
#   5. Separate TTN/PB10.5 plots
#   6. Combined tidy waveform CSV
#   7. PDF and SVG files
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "svglite"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(svglite)


# ============================================================
# 2. USER SETTINGS
# ============================================================

# Folder containing all EFPA CSV outputs
input_directory <- "Downloads/baseline_results /"

# Folder where plots and combined data will be saved
output_directory <- file.path(
  input_directory,
  "Mean_beat_plots"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# Plot time units:
# EFPA mean-beat files may use seconds.
# TRUE converts time from seconds to milliseconds.
convert_time_to_ms <- TRUE

# Signal scaling:
# The EFPA data may be expressed in raw current units.
# TRUE multiplies signal by 10^6 and displays micro-units.
scale_signal_to_micro <- TRUE

# Display replicate traces behind the condition mean?
show_replicates <- TRUE

# Add SEM ribbon around each variant mean?
show_sem_ribbon <- TRUE

# Apply baseline centring to each mean beat?
# This subtracts the mean signal from the earliest pre-beat region.
baseline_centre_traces <- TRUE

# Proportion of the beginning of each waveform used for centring.
# For example, 0.10 uses the first 10% of the waveform.
baseline_proportion <- 0.10

# Optional time limits for plots in milliseconds.
# Use NULL to show the complete EFPA mean-beat window.
plot_time_limits_ms <- NULL

# Example:
# plot_time_limits_ms <- c(-200, 1000)


# ============================================================
# 3. WELL-TO-CONDITION LOOKUP
# ============================================================

well_lookup <- tibble::tribble(
  ~Well,  ~Variant,     ~Experiment,       ~Replicate,
  "E1",   "WT",         "Variant screen",  1,
  "E2",   "WT",         "Variant screen",  2,
  "D1",   "P9R",        "Variant screen",  1,
  "D2",   "P9R",        "Variant screen",  2,
  "F1",   "D12Y",       "Variant screen",  1,
  "F2",   "D12Y",       "Variant screen",  2,
  "B7",   "V24I",       "Variant screen",  1,
  "B8",   "V24I",       "Variant screen",  2,
  "A1",   "K26R",       "Variant screen",  1,
  "A2",   "K26R",       "Variant screen",  2,
  "C1",   "S30F",       "Variant screen",  1,
  "C2",   "S30F",       "Variant screen",  2,
  "D7",   "A38V",       "Variant screen",  1,
  "D8",   "A38V",       "Variant screen",  2,
  "C7",   "R54H",       "Variant screen",  1,
  "C8",   "R54H",       "Variant screen",  2,
  "E7",   "W58*",       "Variant screen",  1,
  "E8",   "W58*",       "Variant screen",  2,
  "A7",   "S70=",       "Variant screen",  1,
  "A8",   "S70=",       "Variant screen",  2,
  "B1",   "KD",         "Variant screen",  1,
  "B2",   "KD",         "Variant screen",  2,
  "F7",   "PB10.5",     "TTN experiment",  1,
  "F8",   "PB10.5",     "TTN experiment",  2,
  "F9",   "TTN",        "TTN experiment",  1,
  "F10",  "TTN",        "TTN experiment",  2,
  "F11",  "TTN+D12Y",   "TTN experiment",  1,
  "F12",  "TTN+D12Y",   "TTN experiment",  2
)

variant_levels <- c(
  "WT",
  "P9R",
  "D12Y",
  "V24I",
  "K26R",
  "S30F",
  "A38V",
  "R54H",
  "W58*",
  "S70=",
  "KD",
  "PB10.5",
  "TTN",
  "TTN+D12Y"
)

timepoint_levels <- c(
  "Baseline",
  "Post-NA"
)


# ============================================================
# 4. FIND CSV FILES
# ============================================================

all_csv_files <- list.files(
  path = input_directory,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(all_csv_files) == 0) {
  stop(
    "No CSV files were found in:\n",
    input_directory
  )
}

cat(
  "CSV files found:\n",
  paste(basename(all_csv_files), collapse = "\n"),
  "\n\n",
  sep = ""
)


# ============================================================
# 5. IDENTIFY MEAN-BEAT FILES
# ============================================================

has_time_column <- function(file_path) {
  
  header_names <- tryCatch(
    {
      names(
        readr::read_delim(
          file = file_path,
          delim = NULL,
          n_max = 1,
          show_col_types = FALSE,
          progress = FALSE,
          name_repair = "minimal"
        )
      )
    },
    error = function(e) {
      character(0)
    }
  )
  
  any(
    tolower(trimws(header_names)) == "time"
  )
}

mean_beat_file_status <- vapply(
  all_csv_files,
  has_time_column,
  logical(1)
)

mean_beat_files <- all_csv_files[
  mean_beat_file_status
]

if (length(mean_beat_files) == 0) {
  stop(
    "No mean-beat CSV files with a Time column were found."
  )
}

cat(
  "Mean-beat files identified:\n",
  paste(basename(mean_beat_files), collapse = "\n"),
  "\n\n",
  sep = ""
)


# ============================================================
# 6. IDENTIFY TIME POINT FROM FILENAME
# ============================================================

identify_timepoint <- function(file_path) {
  
  file_name <- tolower(
    basename(file_path)
  )
  
  baseline_match <- str_detect(
    file_name,
    "baseline|sweep[_ -]?0*6([^0-9]|$)|sweep[_ -]?006"
  )
  
  post_na_match <- str_detect(
    file_name,
    "post.?na|norad|30.?min|sweep[_ -]?0*9([^0-9]|$)|sweep[_ -]?009"
  )
  
  if (baseline_match && !post_na_match) {
    return("Baseline")
  }
  
  if (post_na_match && !baseline_match) {
    return("Post-NA")
  }
  
  return(NA_character_)
}

file_metadata <- tibble(
  File = mean_beat_files,
  File_name = basename(mean_beat_files),
  Timepoint = vapply(
    mean_beat_files,
    identify_timepoint,
    character(1)
  )
)

print(file_metadata)

if (any(is.na(file_metadata$Timepoint))) {
  
  unidentified_files <- file_metadata$File_name[
    is.na(file_metadata$Timepoint)
  ]
  
  stop(
    "The time point could not be identified for:\n",
    paste(unidentified_files, collapse = "\n"),
    "\n\nRename the files so they contain Baseline/Sweep006 ",
    "or PostNA/Sweep009."
  )
}

if (any(duplicated(file_metadata$Timepoint))) {
  warning(
    "More than one mean-beat CSV was identified for the same time point. ",
    "All matching files will be included. Check that this is intentional."
  )
}


# ============================================================
# 7. READ ONE MEAN-BEAT FILE
# ============================================================

read_mean_beat_file <- function(
    file_path,
    timepoint_label
) {
  
  raw_data <- readr::read_delim(
    file = file_path,
    delim = NULL,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  
  names(raw_data) <- trimws(
    names(raw_data)
  )
  
  time_column <- names(raw_data)[
    tolower(names(raw_data)) == "time"
  ][1]
  
  if (is.na(time_column)) {
    stop(
      "No Time column was found in:\n",
      basename(file_path)
    )
  }
  
  raw_data <- raw_data %>%
    rename(Time_raw = all_of(time_column))
  
  available_well_columns <- intersect(
    well_lookup$Well,
    names(raw_data)
  )
  
  if (length(available_well_columns) == 0) {
    stop(
      "No recognised well columns were found in:\n",
      basename(file_path)
    )
  }
  
  missing_wells <- setdiff(
    well_lookup$Well,
    available_well_columns
  )
  
  if (length(missing_wells) > 0) {
    message(
      "These wells were absent from ",
      basename(file_path),
      ": ",
      paste(missing_wells, collapse = ", ")
    )
  }
  
  long_data <- raw_data %>%
    select(
      Time_raw,
      all_of(available_well_columns)
    ) %>%
    mutate(
      Time_raw = suppressWarnings(
        as.numeric(Time_raw)
      )
    ) %>%
    pivot_longer(
      cols = all_of(available_well_columns),
      names_to = "Well",
      values_to = "Signal_raw"
    ) %>%
    mutate(
      Signal_raw = suppressWarnings(
        as.numeric(Signal_raw)
      )
    ) %>%
    filter(
      !is.na(Time_raw),
      !is.na(Signal_raw)
    ) %>%
    left_join(
      well_lookup,
      by = "Well"
    ) %>%
    mutate(
      Timepoint = timepoint_label,
      Source_file = basename(file_path)
    )
  
  long_data
}


# ============================================================
# 8. COMBINE BASELINE AND POST-NA FILES
# ============================================================

mean_beat_long <- purrr::map2_dfr(
  .x = file_metadata$File,
  .y = file_metadata$Timepoint,
  .f = read_mean_beat_file
)

if (nrow(mean_beat_long) == 0) {
  stop("No waveform data were successfully imported.")
}


# ============================================================
# 9. FORMAT TIME AND SIGNAL
# ============================================================

mean_beat_long <- mean_beat_long %>%
  mutate(
    Time = if (convert_time_to_ms) {
      Time_raw * 1000
    } else {
      Time_raw
    },
    Signal = if (scale_signal_to_micro) {
      Signal_raw * 1e6
    } else {
      Signal_raw
    },
    Timepoint = factor(
      Timepoint,
      levels = timepoint_levels
    ),
    Variant = factor(
      Variant,
      levels = variant_levels
    )
  )

time_axis_label <- if (convert_time_to_ms) {
  "Time relative to beat alignment (ms)"
} else {
  "Time relative to beat alignment (s)"
}

signal_axis_label <- if (scale_signal_to_micro) {
  "Mean EFP signal (×10⁻⁶)"
} else {
  "Mean EFP signal"
}


# ============================================================
# 10. BASELINE-CENTRE EACH MEAN BEAT
# ============================================================

if (baseline_centre_traces) {
  
  baseline_values <- mean_beat_long %>%
    group_by(
      Source_file,
      Timepoint,
      Well
    ) %>%
    arrange(
      Time,
      .by_group = TRUE
    ) %>%
    mutate(
      Baseline_n = max(
        1,
        floor(n() * baseline_proportion)
      )
    ) %>%
    summarise(
      Signal_baseline = mean(
        Signal[
          seq_len(
            first(Baseline_n)
          )
        ],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  mean_beat_long <- mean_beat_long %>%
    left_join(
      baseline_values,
      by = c(
        "Source_file",
        "Timepoint",
        "Well"
      )
    ) %>%
    mutate(
      Signal = Signal - Signal_baseline
    )
}


# ============================================================
# 11. QUALITY CHECKS
# ============================================================

qc_summary <- mean_beat_long %>%
  group_by(
    Timepoint,
    Variant,
    Well
  ) %>%
  summarise(
    Number_of_points = n(),
    Time_min = min(Time, na.rm = TRUE),
    Time_max = max(Time, na.rm = TRUE),
    Signal_min = min(Signal, na.rm = TRUE),
    Signal_max = max(Signal, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  qc_summary,
  file.path(
    output_directory,
    "Mean_beat_waveform_QC.csv"
  ),
  row.names = FALSE
)

print(qc_summary)


# ============================================================
# 12. CONDITION-LEVEL MEAN AND SEM
# ============================================================

variant_summary <- mean_beat_long %>%
  group_by(
    Experiment,
    Variant,
    Timepoint,
    Time
  ) %>%
  summarise(
    Mean_signal = mean(
      Signal,
      na.rm = TRUE
    ),
    SD_signal = sd(
      Signal,
      na.rm = TRUE
    ),
    N_wells = sum(
      !is.na(Signal)
    ),
    SEM_signal = ifelse(
      N_wells > 1,
      SD_signal / sqrt(N_wells),
      NA_real_
    ),
    .groups = "drop"
  )


# ============================================================
# 13. COMMON THEME
# ============================================================

theme_efpa <- theme_classic(
  base_size = 11,
  base_family = "Arial"
) +
  theme(
    axis.title = element_text(
      face = "plain"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.title = element_blank(),
    legend.position = "top",
    panel.spacing = grid::unit(
      1,
      "lines"
    )
  )


# ============================================================
# 14. INDIVIDUAL WELL TRACES
# ============================================================

p_individual <- ggplot(
  mean_beat_long,
  aes(
    x = Time,
    y = Signal,
    group = interaction(
      Well,
      Timepoint
    ),
    linetype = Timepoint
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.55,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Variant,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    x = time_axis_label,
    y = signal_axis_label,
    title = "EFPA mean-beat waveforms",
    subtitle = "Individual well traces at baseline and after noradrenaline"
  ) +
  theme_efpa

if (!is.null(plot_time_limits_ms)) {
  p_individual <- p_individual +
    coord_cartesian(
      xlim = plot_time_limits_ms
    )
}

print(p_individual)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_individual_wells.pdf"
  ),
  plot = p_individual,
  width = 12,
  height = 10,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_individual_wells.svg"
  ),
  plot = p_individual,
  width = 12,
  height = 10,
  device = svglite
)


# ============================================================
# 15. VARIANT MEANS WITH REPLICATE TRACES AND SEM
# ============================================================

p_variant_mean <- ggplot()

if (show_replicates) {
  
  p_variant_mean <- p_variant_mean +
    geom_line(
      data = mean_beat_long,
      aes(
        x = Time,
        y = Signal,
        group = interaction(
          Timepoint,
          Well
        ),
        linetype = Timepoint
      ),
      linewidth = 0.35,
      alpha = 0.25
    )
}

if (show_sem_ribbon) {
  
  p_variant_mean <- p_variant_mean +
    geom_ribbon(
      data = variant_summary,
      aes(
        x = Time,
        ymin = Mean_signal - SEM_signal,
        ymax = Mean_signal + SEM_signal,
        group = Timepoint
      ),
      alpha = 0.15,
      inherit.aes = FALSE
    )
}

p_variant_mean <- p_variant_mean +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    data = variant_summary,
    aes(
      x = Time,
      y = Mean_signal,
      group = Timepoint,
      linetype = Timepoint
    ),
    linewidth = 1
  ) +
  facet_wrap(
    ~ Variant,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    x = time_axis_label,
    y = signal_axis_label,
    title = "Mean EFP waveform response to noradrenaline",
    subtitle = "Thin traces show individual wells; thick traces show variant means"
  ) +
  theme_efpa

if (!is.null(plot_time_limits_ms)) {
  p_variant_mean <- p_variant_mean +
    coord_cartesian(
      xlim = plot_time_limits_ms
    )
}

print(p_variant_mean)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_variant_means_with_SEM.pdf"
  ),
  plot = p_variant_mean,
  width = 12,
  height = 10,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_variant_means_with_SEM.svg"
  ),
  plot = p_variant_mean,
  width = 12,
  height = 10,
  device = svglite
)


# ============================================================
# 16. VARIANT-SCREEN FIGURE
# ============================================================

variant_screen_data <- mean_beat_long %>%
  filter(
    Experiment == "Variant screen"
  )

variant_screen_summary <- variant_summary %>%
  filter(
    Experiment == "Variant screen"
  )

p_variant_screen <- ggplot()

if (show_replicates) {
  
  p_variant_screen <- p_variant_screen +
    geom_line(
      data = variant_screen_data,
      aes(
        x = Time,
        y = Signal,
        group = interaction(
          Well,
          Timepoint
        ),
        linetype = Timepoint
      ),
      linewidth = 0.3,
      alpha = 0.25
    )
}

if (show_sem_ribbon) {
  
  p_variant_screen <- p_variant_screen +
    geom_ribbon(
      data = variant_screen_summary,
      aes(
        x = Time,
        ymin = Mean_signal - SEM_signal,
        ymax = Mean_signal + SEM_signal,
        group = Timepoint
      ),
      alpha = 0.15
    )
}

p_variant_screen <- p_variant_screen +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    data = variant_screen_summary,
    aes(
      x = Time,
      y = Mean_signal,
      group = Timepoint,
      linetype = Timepoint
    ),
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ Variant,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    x = time_axis_label,
    y = signal_axis_label,
    title = "HOPX variant mean-beat waveforms",
    subtitle = "Baseline compared with the acute noradrenaline response"
  ) +
  theme_efpa

if (!is.null(plot_time_limits_ms)) {
  p_variant_screen <- p_variant_screen +
    coord_cartesian(
      xlim = plot_time_limits_ms
    )
}

print(p_variant_screen)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_HOPX_variant_screen.pdf"
  ),
  plot = p_variant_screen,
  width = 11,
  height = 9,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_HOPX_variant_screen.svg"
  ),
  plot = p_variant_screen,
  width = 11,
  height = 9,
  device = svglite
)


# ============================================================
# 17. TTN EXPERIMENT FIGURE
# ============================================================

ttn_data <- mean_beat_long %>%
  filter(
    Experiment == "TTN experiment"
  )

ttn_summary <- variant_summary %>%
  filter(
    Experiment == "TTN experiment"
  )

p_ttn <- ggplot()

if (show_replicates) {
  
  p_ttn <- p_ttn +
    geom_line(
      data = ttn_data,
      aes(
        x = Time,
        y = Signal,
        group = interaction(
          Well,
          Timepoint
        ),
        linetype = Timepoint
      ),
      linewidth = 0.35,
      alpha = 0.3
    )
}

if (show_sem_ribbon) {
  
  p_ttn <- p_ttn +
    geom_ribbon(
      data = ttn_summary,
      aes(
        x = Time,
        ymin = Mean_signal - SEM_signal,
        ymax = Mean_signal + SEM_signal,
        group = Timepoint
      ),
      alpha = 0.15
    )
}

p_ttn <- p_ttn +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_line(
    data = ttn_summary,
    aes(
      x = Time,
      y = Mean_signal,
      group = Timepoint,
      linetype = Timepoint
    ),
    linewidth = 1
  ) +
  facet_wrap(
    ~ Variant,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    x = time_axis_label,
    y = signal_axis_label,
    title = "Mean-beat waveforms in the TTN experiment",
    subtitle = "PB10.5, TTN and TTN+D12Y before and after noradrenaline"
  ) +
  theme_efpa

if (!is.null(plot_time_limits_ms)) {
  p_ttn <- p_ttn +
    coord_cartesian(
      xlim = plot_time_limits_ms
    )
}

print(p_ttn)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_TTN_experiment.pdf"
  ),
  plot = p_ttn,
  width = 11,
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_directory,
    "Mean_beats_TTN_experiment.svg"
  ),
  plot = p_ttn,
  width = 11,
  height = 4.5,
  device = svglite
)


# ============================================================
# 18. ONE PDF PER VARIANT
# ============================================================

variants_present <- mean_beat_long %>%
  filter(
    !is.na(Variant)
  ) %>%
  distinct(
    Variant
  ) %>%
  pull(
    Variant
  ) %>%
  as.character()

individual_variant_directory <- file.path(
  output_directory,
  "Individual_variants"
)

dir.create(
  individual_variant_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

for (current_variant in variants_present) {
  
  current_data <- mean_beat_long %>%
    filter(
      as.character(Variant) == current_variant
    )
  
  current_summary <- variant_summary %>%
    filter(
      as.character(Variant) == current_variant
    )
  
  p_current <- ggplot()
  
  if (show_replicates) {
    
    p_current <- p_current +
      geom_line(
        data = current_data,
        aes(
          x = Time,
          y = Signal,
          group = interaction(
            Well,
            Timepoint
          ),
          linetype = Timepoint
        ),
        linewidth = 0.45,
        alpha = 0.35
      )
  }
  
  if (show_sem_ribbon) {
    
    p_current <- p_current +
      geom_ribbon(
        data = current_summary,
        aes(
          x = Time,
          ymin = Mean_signal - SEM_signal,
          ymax = Mean_signal + SEM_signal,
          group = Timepoint
        ),
        alpha = 0.15
      )
  }
  
  p_current <- p_current +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.25,
      linetype = "dashed"
    ) +
    geom_line(
      data = current_summary,
      aes(
        x = Time,
        y = Mean_signal,
        group = Timepoint,
        linetype = Timepoint
      ),
      linewidth = 1.1
    ) +
    labs(
      title = paste0(
        current_variant,
        " mean-beat waveform"
      ),
      subtitle = "Baseline compared with post-noradrenaline",
      x = time_axis_label,
      y = signal_axis_label
    ) +
    theme_efpa
  
  if (!is.null(plot_time_limits_ms)) {
    p_current <- p_current +
      coord_cartesian(
        xlim = plot_time_limits_ms
      )
  }
  
  safe_variant_name <- str_replace_all(
    current_variant,
    "[^A-Za-z0-9]+",
    "_"
  )
  
  ggsave(
    filename = file.path(
      individual_variant_directory,
      paste0(
        "Mean_beat_",
        safe_variant_name,
        ".pdf"
      )
    ),
    plot = p_current,
    width = 7,
    height = 5,
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(
      individual_variant_directory,
      paste0(
        "Mean_beat_",
        safe_variant_name,
        ".svg"
      )
    ),
    plot = p_current,
    width = 7,
    height = 5,
    device = svglite
  )
}


# ============================================================
# 19. EXPORT COMBINED DATA
# ============================================================

write.csv(
  mean_beat_long %>%
    select(
      Source_file,
      Experiment,
      Timepoint,
      Variant,
      Well,
      Replicate,
      Time,
      Signal
    ),
  file.path(
    output_directory,
    "EFPA_mean_beats_combined_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  variant_summary,
  file.path(
    output_directory,
    "EFPA_mean_beats_variant_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. FINAL MESSAGE
# ============================================================

cat(
  "\nMean-beat plotting complete.\n\n",
  "Output folder:\n",
  output_directory,
  "\n\nCreated:\n",
  "- individual-well waveform figure\n",
  "- variant mean ± SEM figure\n",
  "- HOPX variant-screen figure\n",
  "- TTN experiment figure\n",
  "- one PDF and SVG per variant\n",
  "- combined long-format CSV\n",
  "- condition-level summary CSV\n",
  sep = ""
)
