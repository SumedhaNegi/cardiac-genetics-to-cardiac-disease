# =============================================================================
# EXAMINER GUIDE
# Script: cardioexcyte_thesis_figures.R
# Purpose: Builds thesis CardioExcyte extracellular field-potential figures, including baseline normalisation, post-noradrenaline windows, ANOVA and Dunnett comparisons.
# Inputs: An in-memory data.frame named mydata with the columns documented below.
# Outputs: CSV audit/statistics tables and PDF/PNG thesis figures in a generated output directory.
# Methodological note: Review the detailed inline section headings for data
# filtering, baseline handling, aggregation, statistical model and plotting.
# Reproducibility: Requires original instrument exports and the stated R packages.
# Review note: No testing or statistical re-analysis was performed during curation.
# =============================================================================

# CODE MAP
# Processes CardioExcyte sweeps for thesis figures, including baseline/acute windows, per-well normalisation, replicate aggregation, statistical comparisons and plots; review experimental-unit assumptions.
#

# ============================================================
# FINAL CARDIOEXCYTE EFP THESIS FIGURES
#
# Input:
#   mydata is already loaded
#
# Parameters:
#   Beat rate (BPM)
#   Extracellular field-potential amplitude (V)
#   Beat-rate regularity CoV (s)
#   Field-potential duration, abs1 (s)
#   Downstroke velocity (V/s)
#
# Data processing:
#   - Sweeps 001 and 002 excluded
#   - Missing measurements from prior outlier removal allowed
#   - Baseline calculated separately for each well and parameter
#   - Post-NA values expressed as fold-change from that
#     well-specific baseline
#   - The two wells are averaged within each variant and sweep
#
# Time windows:
#   Baseline: sweeps 003–006
#   Acute: sweeps 007–010
#   24 h: up to five sweeps closest to 24 h post-NA
#         within 20–28 h post-NA
#   48 h: up to five sweeps closest to 48 h post-NA
#         within 44–52 h post-NA
#
# Statistics:
#   Separate one-way ANOVA for each parameter and time window
#
#       Analysis_Value ~ Variant
#
#   Dunnett-adjusted post hoc comparisons:
#       each variant versus WT within the same time window
# ============================================================


# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "readr",
  "tibble",
  "purrr",
  "emmeans"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(tibble)
library(purrr)
library(emmeans)


# ============================================================
# 2. CHECK INPUT OBJECT
# ============================================================

if (!exists("mydata")) {
  stop(
    "The object 'mydata' is not loaded."
  )
}

if (!is.data.frame(mydata)) {
  stop(
    "'mydata' must be a data frame or tibble."
  )
}

required_columns <- c(
  "Variant",
  "Condition",
  "Well_ID",
  "Sweep",
  "Parameter",
  "Time_Hours",
  "Replicate_Value"
)

missing_columns <- setdiff(
  required_columns,
  names(mydata)
)

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing:\n",
    paste(
      missing_columns,
      collapse = "\n"
    )
  )
}


# ============================================================
# 3. CREATE OUTPUT DIRECTORY
# ============================================================

output_directory <- file.path(
  getwd(),
  "CardioExcyte_final_EFP_thesis_figures"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. DEFINE VARIANTS AND PARAMETERS
# ============================================================

variant_order <- c(
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
  "KD"
)

condition_lookup <- c(
  "RYU" = "WT",
  "UNR" = "P9R",
  "MTG" = "D12Y",
  "YHH" = "V24I",
  "ZTT" = "K26R",
  "GRP" = "S30F",
  "YCP" = "A38V",
  "BKH" = "R54H",
  "JNJ" = "W58*",
  "VDW" = "S70=",
  "BFP" = "KD"
)

parameter_order <- c(
  "Beat Rate",
  "Amplitude",
  "BR Regularity",
  "FP Duration",
  "Downstroke Velocity"
)

parameter_labels <- c(
  "Beat Rate" = "Beat rate (BPM)",
  "Amplitude" = "Extracellular field-potential amplitude (V)",
  "BR Regularity" = "Beat-rate regularity CoV (s)",
  "FP Duration" = "Field-potential duration, abs1 (s)",
  "Downstroke Velocity" = "Downstroke velocity (V/s)"
)


# ============================================================
# 5. REBUILD DATA FROM THE RAW COLUMNS IN mydata
#
# Existing calculated columns are ignored:
#   Baseline_Mean
#   Relative_Response
#   Time_Window
#   Plot_Value
# ============================================================

cardio_raw <- mydata |>
  transmute(
    Variant_original = as.character(Variant),

    Condition = as.character(Condition),

    Well_ID = as.character(Well_ID),

    Sweep_number = readr::parse_number(
      as.character(Sweep)
    ),

    Parameter_original = as.character(Parameter),

    Time_Hours = as.numeric(Time_Hours),

    Replicate_Value = as.numeric(
      Replicate_Value
    )
  ) |>
  mutate(
    Variant_from_condition = unname(
      condition_lookup[Condition]
    ),

    Variant = case_when(
      Variant_original %in% variant_order ~
        Variant_original,

      Variant_from_condition %in% variant_order ~
        Variant_from_condition,

      TRUE ~ NA_character_
    ),

    Parameter = case_when(
      str_detect(
        Parameter_original,
        regex(
          "^Beat Rate$|Beat Rate|BPM",
          ignore_case = TRUE
        )
      ) ~ "Beat Rate",

      str_detect(
        Parameter_original,
        regex(
          "^Amplitude$|Amplitude",
          ignore_case = TRUE
        )
      ) ~ "Amplitude",

      str_detect(
        Parameter_original,
        regex(
          "Regularity|CoV",
          ignore_case = TRUE
        )
      ) ~ "BR Regularity",

      str_detect(
        Parameter_original,
        regex(
          "FP Duration|Field Potential",
          ignore_case = TRUE
        )
      ) ~ "FP Duration",

      str_detect(
        Parameter_original,
        regex(
          "Downstroke",
          ignore_case = TRUE
        )
      ) ~ "Downstroke Velocity",

      TRUE ~ NA_character_
    )
  ) |>
  filter(
    Variant %in% variant_order,
    Parameter %in% parameter_order,
    is.finite(Sweep_number),
    is.finite(Time_Hours),
    is.finite(Replicate_Value),
    !is.na(Well_ID),
    !Sweep_number %in% c(
      1,
      2
    )
  ) |>
  mutate(
    Variant = factor(
      Variant,
      levels = variant_order
    ),

    Parameter = factor(
      Parameter,
      levels = parameter_order
    ),

    Sweep = sprintf(
      "%03d",
      Sweep_number
    )
  )


# ============================================================
# 6. CHECK VARIANT MAPPING
# ============================================================

variant_mapping_check <- cardio_raw |>
  distinct(
    Condition,
    Variant_original,
    Variant_from_condition,
    Variant
  ) |>
  arrange(
    Variant
  )

cat(
  "\n============================================================\n"
)

cat(
  "VARIANT MAPPING\n"
)

cat(
  "============================================================\n\n"
)

print(
  variant_mapping_check,
  n = Inf
)

variant_mapping_mismatches <- cardio_raw |>
  filter(
    !is.na(Variant_original),
    !is.na(Variant_from_condition),
    Variant_original != Variant_from_condition
  ) |>
  distinct(
    Condition,
    Variant_original,
    Variant_from_condition
  )

if (nrow(variant_mapping_mismatches) > 0) {
  print(
    variant_mapping_mismatches,
    n = Inf
  )

  stop(
    "The Variant and Condition columns contain inconsistent mappings."
  )
}


# ============================================================
# 7. CHECK EACH WELL BELONGS TO ONE VARIANT
# ============================================================

well_variant_check <- cardio_raw |>
  distinct(
    Well_ID,
    Variant
  ) |>
  count(
    Well_ID,
    name = "Number_of_variants"
  ) |>
  filter(
    Number_of_variants != 1
  )

if (nrow(well_variant_check) > 0) {
  print(
    well_variant_check,
    n = Inf
  )

  stop(
    "At least one well is assigned to more than one variant."
  )
}

cat(
  "\nPASS: Every well is assigned to one variant.\n"
)


# ============================================================
# 8. CREATE ONE VALUE PER WELL X SWEEP X PARAMETER
#
# Any duplicate source rows are averaged.
# ============================================================

well_sweep_data <- cardio_raw |>
  group_by(
    Variant,
    Condition,
    Well_ID,
    Sweep_number,
    Sweep,
    Parameter
  ) |>
  summarise(
    Time_Hours = median(
      Time_Hours,
      na.rm = TRUE
    ),

    Replicate_Value = mean(
      Replicate_Value,
      na.rm = TRUE
    ),

    Number_of_source_rows = n(),

    .groups = "drop"
  )


# ============================================================
# 9. REPORT DUPLICATE SOURCE ROWS
# ============================================================

duplicate_source_rows <- well_sweep_data |>
  filter(
    Number_of_source_rows > 1
  )

if (nrow(duplicate_source_rows) > 0) {
  warning(
    paste0(
      "Some well x sweep x parameter combinations contained ",
      "multiple source rows. These rows were averaged."
    )
  )

  write.csv(
    duplicate_source_rows,
    file.path(
      output_directory,
      "Duplicate_source_rows.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# 10. CHECK NUMBER OF WELLS PER VARIANT
# ============================================================

variant_well_check <- well_sweep_data |>
  distinct(
    Variant,
    Well_ID
  ) |>
  count(
    Variant,
    name = "Number_of_wells"
  )

cat(
  "\n============================================================\n"
)

cat(
  "NUMBER OF WELLS PER VARIANT\n"
)

cat(
  "============================================================\n\n"
)

print(
  variant_well_check,
  n = Inf
)

if (any(variant_well_check$Number_of_wells != 2)) {
  warning(
    paste0(
      "At least one variant does not contain exactly two wells. ",
      "Inspect the printed table."
    )
  )
}


# ============================================================
# 11. CREATE SWEEP-TIME LOOKUP
# ============================================================

sweep_time_lookup <- well_sweep_data |>
  group_by(
    Sweep_number,
    Sweep
  ) |>
  summarise(
    Median_Time_Hours = median(
      Time_Hours,
      na.rm = TRUE
    ),

    Minimum_Time_Hours = min(
      Time_Hours,
      na.rm = TRUE
    ),

    Maximum_Time_Hours = max(
      Time_Hours,
      na.rm = TRUE
    ),

    Time_Range = Maximum_Time_Hours -
      Minimum_Time_Hours,

    Number_of_rows = n(),

    .groups = "drop"
  ) |>
  arrange(
    Sweep_number
  )


# ============================================================
# 12. ESTIMATE NORADRENALINE ADDITION TIME
#
# Sweep 007 was identified as the first rapid beat-rate
# response. NA was therefore added between sweeps 006 and 007.
# ============================================================

sweep_006_time <- sweep_time_lookup |>
  filter(
    Sweep_number == 6
  ) |>
  pull(
    Median_Time_Hours
  )

sweep_007_time <- sweep_time_lookup |>
  filter(
    Sweep_number == 7
  ) |>
  pull(
    Median_Time_Hours
  )

if (
  length(sweep_006_time) != 1 ||
  length(sweep_007_time) != 1
) {
  stop(
    "Sweeps 006 and 007 must both be present."
  )
}

estimated_na_time <- mean(
  c(
    sweep_006_time,
    sweep_007_time
  )
)

sweep_time_lookup <- sweep_time_lookup |>
  mutate(
    Time_Post_NA = Median_Time_Hours -
      estimated_na_time,

    Distance_From_24h = abs(
      Time_Post_NA - 24
    ),

    Distance_From_48h = abs(
      Time_Post_NA - 48
    )
  )

cat(
  "\n============================================================\n"
)

cat(
  "NORADRENALINE TIMING\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Sweep 006: ",
  round(
    sweep_006_time,
    4
  ),
  " h\n",
  sep = ""
)

cat(
  "Sweep 007: ",
  round(
    sweep_007_time,
    4
  ),
  " h\n",
  sep = ""
)

cat(
  "Estimated NA addition: ",
  round(
    estimated_na_time,
    4
  ),
  " h after recording began\n",
  sep = ""
)


# ============================================================
# 13. DEFINE BASELINE AND ACUTE SWEEPS
# ============================================================

baseline_sweeps <- c(
  3,
  4,
  5,
  6
)

acute_sweeps <- c(
  7,
  8,
  9,
  10
)

missing_baseline_sweeps <- setdiff(
  baseline_sweeps,
  sweep_time_lookup$Sweep_number
)

missing_acute_sweeps <- setdiff(
  acute_sweeps,
  sweep_time_lookup$Sweep_number
)

if (length(missing_baseline_sweeps) > 0) {
  stop(
    "Missing baseline sweeps: ",
    paste(
      sprintf(
        "%03d",
        missing_baseline_sweeps
      ),
      collapse = ", "
    )
  )
}

if (length(missing_acute_sweeps) > 0) {
  stop(
    "Missing acute sweeps: ",
    paste(
      sprintf(
        "%03d",
        missing_acute_sweeps
      ),
      collapse = ", "
    )
  )
}


# ============================================================
# 14. SELECT UP TO FIVE SWEEPS NEAREST 24 HOURS POST-NA
#
# Candidate interval:
#   20–28 hours post-NA
# ============================================================

candidate_24h_sweeps <- sweep_time_lookup |>
  filter(
    Time_Post_NA >= 20,
    Time_Post_NA <= 28
  ) |>
  arrange(
    Distance_From_24h,
    Sweep_number
  )

sweeps_24h <- candidate_24h_sweeps |>
  slice_head(
    n = 5
  ) |>
  arrange(
    Sweep_number
  ) |>
  pull(
    Sweep_number
  )

if (length(sweeps_24h) == 0) {
  warning(
    paste0(
      "No sweeps were found within 20–28 hours post-NA. ",
      "The 24-hour window will be omitted."
    )
  )
}


# ============================================================
# 15. SELECT UP TO FIVE SWEEPS NEAREST 48 HOURS POST-NA
#
# Candidate interval:
#   44–52 hours post-NA
# ============================================================

candidate_48h_sweeps <- sweep_time_lookup |>
  filter(
    Time_Post_NA >= 44,
    Time_Post_NA <= 52
  ) |>
  arrange(
    Distance_From_48h,
    Sweep_number
  )

sweeps_48h <- candidate_48h_sweeps |>
  slice_head(
    n = 5
  ) |>
  arrange(
    Sweep_number
  ) |>
  pull(
    Sweep_number
  )

if (length(sweeps_48h) == 0) {
  warning(
    paste0(
      "No sweeps were found within 44–52 hours post-NA. ",
      "The 48-hour window will be omitted."
    )
  )
}


# ============================================================
# 16. PRINT THE SELECTED SWEEPS
# ============================================================

selected_windows <- bind_rows(
  tibble(
    Time_Window = "Baseline",
    Sweep_number = baseline_sweeps
  ),

  tibble(
    Time_Window = "30 min post-NA",
    Sweep_number = acute_sweeps
  ),

  tibble(
    Time_Window = "24 h post-NA",
    Sweep_number = sweeps_24h
  ),

  tibble(
    Time_Window = "48 h post-NA",
    Sweep_number = sweeps_48h
  )
) |>
  left_join(
    sweep_time_lookup |>
      select(
        Sweep_number,
        Median_Time_Hours,
        Time_Post_NA,
        Distance_From_24h,
        Distance_From_48h
      ),
    by = "Sweep_number"
  ) |>
  mutate(
    Sweep = sprintf(
      "%03d",
      Sweep_number
    )
  )

cat(
  "\n============================================================\n"
)

cat(
  "SELECTED ANALYSIS SWEEPS\n"
)

cat(
  "============================================================\n\n"
)

print(
  selected_windows,
  n = Inf,
  width = Inf
)

write.csv(
  selected_windows,
  file.path(
    output_directory,
    "Selected_analysis_sweeps.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. DEFINE AVAILABLE WINDOW ORDER
# ============================================================

available_window_levels <- c(
  "Baseline",
  "30 min post-NA"
)

if (length(sweeps_24h) > 0) {
  available_window_levels <- c(
    available_window_levels,
    "24 h post-NA"
  )
}

if (length(sweeps_48h) > 0) {
  available_window_levels <- c(
    available_window_levels,
    "48 h post-NA"
  )
}


# ============================================================
# 18. ASSIGN TIME WINDOWS
# ============================================================

analysis_data <- well_sweep_data |>
  mutate(
    Time_Window = case_when(
      Sweep_number %in% baseline_sweeps ~
        "Baseline",

      Sweep_number %in% acute_sweeps ~
        "30 min post-NA",

      Sweep_number %in% sweeps_24h ~
        "24 h post-NA",

      Sweep_number %in% sweeps_48h ~
        "48 h post-NA",

      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(Time_Window)
  ) |>
  mutate(
    Time_Window = factor(
      Time_Window,
      levels = available_window_levels
    ),

    Sweep_factor = factor(
      Sweep,
      levels = sprintf(
        "%03d",
        c(
          baseline_sweeps,
          acute_sweeps,
          sweeps_24h,
          sweeps_48h
        )
      )
    )
  )


# ============================================================
# 19. CALCULATE THE RETAINED BASELINE MEAN FOR EACH WELL
#
# Missing baseline sweeps are allowed because outliers were
# removed before this analysis.
# ============================================================

well_baseline <- analysis_data |>
  filter(
    Time_Window == "Baseline"
  ) |>
  group_by(
    Variant,
    Condition,
    Well_ID,
    Parameter
  ) |>
  summarise(
    Baseline_Mean = mean(
      Replicate_Value,
      na.rm = TRUE
    ),

    Number_of_baseline_sweeps = n_distinct(
      Sweep_number
    ),

    Present_baseline_sweeps = paste(
      sort(
        unique(
          Sweep
        )
      ),
      collapse = ", "
    ),

    .groups = "drop"
  )


# ============================================================
# 20. REPORT INCOMPLETE BUT VALID BASELINES
# ============================================================

incomplete_baselines <- well_baseline |>
  filter(
    Number_of_baseline_sweeps <
      length(baseline_sweeps),
    Number_of_baseline_sweeps >= 1,
    is.finite(Baseline_Mean),
    Baseline_Mean != 0
  )

cat(
  "\n============================================================\n"
)

cat(
  "BASELINE COVERAGE AFTER OUTLIER REMOVAL\n"
)

cat(
  "============================================================\n\n"
)

if (nrow(incomplete_baselines) > 0) {
  cat(
    "Some baseline means use fewer than four retained sweeps.\n\n"
  )

  print(
    incomplete_baselines,
    n = Inf,
    width = Inf
  )

  write.csv(
    incomplete_baselines,
    file.path(
      output_directory,
      "Incomplete_but_valid_baseline_means.csv"
    ),
    row.names = FALSE
  )
} else {
  cat(
    "All baseline means use four retained sweeps.\n"
  )
}


# ============================================================
# 21. CHECK FOR GENUINELY INVALID BASELINES
# ============================================================

invalid_baselines <- well_baseline |>
  filter(
    Number_of_baseline_sweeps < 1 |
      !is.finite(Baseline_Mean) |
      Baseline_Mean == 0
  )

if (nrow(invalid_baselines) > 0) {
  print(
    invalid_baselines,
    n = Inf,
    width = Inf
  )

  stop(
    paste0(
      "At least one well and parameter has no usable baseline ",
      "or has a zero/non-finite baseline mean."
    )
  )
}

cat(
  "\nPASS: Every well and parameter has a usable baseline.\n"
)


# ============================================================
# 22. NORMALISE EACH WELL TO ITS OWN BASELINE
# ============================================================

analysis_data <- analysis_data |>
  left_join(
    well_baseline |>
      select(
        Variant,
        Condition,
        Well_ID,
        Parameter,
        Baseline_Mean,
        Number_of_baseline_sweeps
      ),
    by = c(
      "Variant",
      "Condition",
      "Well_ID",
      "Parameter"
    )
  ) |>
  mutate(
    Relative_Fold = case_when(
      is.finite(Baseline_Mean) &
        Baseline_Mean != 0 ~
        Replicate_Value / Baseline_Mean,

      TRUE ~ NA_real_
    ),

    Percentage_Change = case_when(
      is.finite(Relative_Fold) ~
        100 * (
          Relative_Fold - 1
        ),

      TRUE ~ NA_real_
    )
  )


# ============================================================
# 23. CHECK BASELINE NORMALISATION
# ============================================================

baseline_normalisation_check <- analysis_data |>
  filter(
    Time_Window == "Baseline"
  ) |>
  group_by(
    Variant,
    Well_ID,
    Parameter
  ) |>
  summarise(
    Mean_Baseline_Fold = mean(
      Relative_Fold,
      na.rm = TRUE
    ),

    Number_of_Retained_Sweeps = n_distinct(
      Sweep_number
    ),

    .groups = "drop"
  )

maximum_baseline_deviation <- max(
  abs(
    baseline_normalisation_check$Mean_Baseline_Fold -
      1
  ),
  na.rm = TRUE
)

if (
  !is.finite(maximum_baseline_deviation) ||
  maximum_baseline_deviation > 1e-10
) {
  print(
    baseline_normalisation_check,
    n = Inf,
    width = Inf
  )

  stop(
    "The baseline normalisation check failed."
  )
}

cat(
  "PASS: Retained baseline values average to one within each well.\n"
)


# ============================================================
# 24. AVERAGE THE TWO WELLS WITHIN EACH VARIANT AND SWEEP
#
# Each resulting row is one:
#   variant x parameter x sweep
# ============================================================

variant_sweep_data <- analysis_data |>
  group_by(
    Variant,
    Parameter,
    Time_Window,
    Sweep_number,
    Sweep_factor
  ) |>
  summarise(
    Absolute_Value = mean(
      Replicate_Value,
      na.rm = TRUE
    ),

    Relative_Fold = mean(
      Relative_Fold,
      na.rm = TRUE
    ),

    Percentage_Change = mean(
      Percentage_Change,
      na.rm = TRUE
    ),

    Number_of_wells = n_distinct(
      Well_ID
    ),

    Wells_included = paste(
      sort(
        unique(
          Well_ID
        )
      ),
      collapse = ", "
    ),

    .groups = "drop"
  )


# ============================================================
# 25. REPORT SWEEP VALUES WITH FEWER THAN TWO WELLS
#
# This can occur if one well observation was removed as an
# outlier. The remaining observation is retained and reported.
# ============================================================

incomplete_two_well_averages <- variant_sweep_data |>
  filter(
    Number_of_wells < 2
  )

if (nrow(incomplete_two_well_averages) > 0) {
  warning(
    paste0(
      "Some variant-by-sweep values contain only one retained ",
      "well after outlier removal."
    )
  )

  write.csv(
    incomplete_two_well_averages,
    file.path(
      output_directory,
      "Sweep_values_with_fewer_than_two_wells.csv"
    ),
    row.names = FALSE
  )
} else {
  cat(
    "PASS: Every sweep-level value is the average of two wells.\n"
  )
}


# ============================================================
# 26. DEFINE THE VALUE USED FOR PLOTTING AND STATISTICS
#
# Baseline:
#   absolute EFP measurement
#
# Post-NA:
#   fold-change from the well-specific baseline
# ============================================================

variant_sweep_data <- variant_sweep_data |>
  mutate(
    Analysis_Value = case_when(
      Time_Window == "Baseline" ~
        Absolute_Value,

      Time_Window != "Baseline" ~
        Relative_Fold,

      TRUE ~ NA_real_
    )
  ) |>
  filter(
    is.finite(Analysis_Value)
  )


# ============================================================
# 27. CALCULATE MEAN AND SEM ACROSS SWEEPS
# ============================================================

figure_summary <- variant_sweep_data |>
  group_by(
    Variant,
    Parameter,
    Time_Window
  ) |>
  summarise(
    Mean = mean(
      Analysis_Value,
      na.rm = TRUE
    ),

    SD = sd(
      Analysis_Value,
      na.rm = TRUE
    ),

    N_Sweeps = sum(
      is.finite(Analysis_Value)
    ),

    SEM = case_when(
      N_Sweeps > 1 ~
        SD / sqrt(N_Sweeps),

      TRUE ~ 0
    ),

    Mean_Number_of_Wells = mean(
      Number_of_wells,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ============================================================
# 28. EMPTY STATISTICS TABLES
# ============================================================

empty_anova_table <- function() {

  tibble(
    Parameter = character(),
    Time_Window = character(),
    Term = character(),
    Df = numeric(),
    Sum_Sq = numeric(),
    Mean_Sq = numeric(),
    F_Value = numeric(),
    P_Value = numeric()
  )
}

empty_dunnett_table <- function() {

  tibble(
    Parameter = character(),
    Time_Window = character(),
    Contrast = character(),
    Variant = character(),
    Estimate = numeric(),
    SE = numeric(),
    Degrees_Freedom = numeric(),
    Test_Statistic = numeric(),
    P_Value = numeric()
  )
}


# ============================================================
# 29. ONE-WAY ANOVA FUNCTION
#
# Separate ANOVA for each parameter and time window:
#
#   Analysis_Value ~ Variant
#
# Post hoc:
#   every variant compared with WT using Dunnett adjustment
# ============================================================

run_one_way_anova <- function(
    parameter_name,
    window_name
) {

  data_subset <- variant_sweep_data |>
    filter(
      Parameter == parameter_name,
      Time_Window == window_name,
      is.finite(Analysis_Value)
    ) |>
    mutate(
      Variant = relevel(
        factor(
          Variant,
          levels = variant_order
        ),
        ref = "WT"
      )
    ) |>
    droplevels()

  if (
    nrow(data_subset) == 0 ||
    n_distinct(data_subset$Variant) < 2 ||
    !"WT" %in% levels(data_subset$Variant)
  ) {
    return(
      list(
        anova = empty_anova_table(),
        dunnett = empty_dunnett_table()
      )
    )
  }

  fitted_model <- aov(
    Analysis_Value ~ Variant,
    data = data_subset
  )

  raw_anova <- as.data.frame(
    summary(fitted_model)[[1]]
  ) |>
    rownames_to_column(
      "Term"
    )

  anova_table <- raw_anova |>
    transmute(
      Parameter = parameter_name,

      Time_Window = window_name,

      Term = Term,

      Df = Df,

      Sum_Sq = `Sum Sq`,

      Mean_Sq = `Mean Sq`,

      F_Value = `F value`,

      P_Value = `Pr(>F)`
    )

  estimated_means <- emmeans::emmeans(
    fitted_model,
    ~Variant
  )

  comparisons <- emmeans::contrast(
    estimated_means,
    method = "trt.vs.ctrl",
    ref = "WT",
    adjust = "dunnett"
  )

  comparison_table <- as.data.frame(
    summary(
      comparisons,
      infer = TRUE
    )
  )

  statistic_column <- intersect(
    c(
      "t.ratio",
      "z.ratio"
    ),
    names(comparison_table)
  )

  if (length(statistic_column) == 0) {
    comparison_table$Test_Statistic <- NA_real_
  } else {
    comparison_table$Test_Statistic <- comparison_table[[statistic_column[1]]]
  }

  dunnett_table <- comparison_table |>
    transmute(
      Parameter = parameter_name,

      Time_Window = window_name,

      Contrast = contrast,

      Variant = str_remove(
        contrast,
        " - WT$"
      ),

      Estimate = estimate,

      SE = SE,

      Degrees_Freedom = df,

      Test_Statistic = Test_Statistic,

      P_Value = p.value
    )

  list(
    anova = anova_table,
    dunnett = dunnett_table
  )
}


# ============================================================
# 30. RUN ALL ONE-WAY ANOVAS
# ============================================================

analysis_grid <- expand_grid(
  Parameter = parameter_order,
  Time_Window = available_window_levels
)

anova_results <- purrr::map2(
  analysis_grid$Parameter,
  analysis_grid$Time_Window,
  run_one_way_anova
)

anova_table <- bind_rows(
  lapply(
    anova_results,
    function(result) {
      result$anova
    }
  )
)

dunnett_table <- bind_rows(
  lapply(
    anova_results,
    function(result) {
      result$dunnett
    }
  )
) |>
  mutate(
    Variant = factor(
      Variant,
      levels = variant_order
    ),

    Parameter = factor(
      Parameter,
      levels = parameter_order
    ),

    Time_Window = factor(
      Time_Window,
      levels = available_window_levels
    ),

    Significance = case_when(
      P_Value < 0.0001 ~ "****",
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),

    Plot_Significance = case_when(
      P_Value < 0.0001 ~ "****",
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write.csv(
  anova_table,
  file.path(
    output_directory,
    "One_way_ANOVA_by_parameter_and_window.csv"
  ),
  row.names = FALSE
)

write.csv(
  dunnett_table,
  file.path(
    output_directory,
    "Dunnett_each_variant_vs_WT.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 31. CALCULATE SIGNIFICANCE-LABEL POSITIONS
# ============================================================

panel_ranges <- variant_sweep_data |>
  group_by(
    Parameter,
    Time_Window
  ) |>
  summarise(
    Panel_Minimum = min(
      Analysis_Value,
      na.rm = TRUE
    ),

    Panel_Maximum = max(
      Analysis_Value,
      na.rm = TRUE
    ),

    Panel_Range = Panel_Maximum -
      Panel_Minimum,

    .groups = "drop"
  ) |>
  mutate(
    Panel_Range = case_when(
      is.finite(Panel_Range) &
        Panel_Range > 0 ~
        Panel_Range,

      TRUE ~
        pmax(
          abs(Panel_Maximum) * 0.1,
          0.1
        )
    )
  )

variant_maximum <- variant_sweep_data |>
  group_by(
    Variant,
    Parameter,
    Time_Window
  ) |>
  summarise(
    Variant_Maximum = max(
      Analysis_Value,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

significance_labels <- dunnett_table |>
  filter(
    Plot_Significance != ""
  ) |>
  left_join(
    variant_maximum,
    by = c(
      "Variant",
      "Parameter",
      "Time_Window"
    )
  ) |>
  left_join(
    panel_ranges,
    by = c(
      "Parameter",
      "Time_Window"
    )
  ) |>
  mutate(
    Label_Y = Variant_Maximum +
      0.08 * Panel_Range
  )


# ============================================================
# 32. THESIS FIGURE THEME
# ============================================================

thesis_theme <- theme_classic(
  base_size = 12
) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      colour = "black"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    ),

    strip.text = element_text(
      face = "bold"
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),

    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),

    legend.position = "bottom",

    legend.title = element_text(
      face = "bold"
    ),

    panel.spacing = grid::unit(
      0.9,
      "lines"
    )
  )


# ============================================================
# 33. FIGURE 1: BASELINE ABSOLUTE EFP PHENOTYPE
# ============================================================

baseline_plot_data <- variant_sweep_data |>
  filter(
    Time_Window == "Baseline"
  )

baseline_plot_summary <- figure_summary |>
  filter(
    Time_Window == "Baseline"
  )

baseline_significance <- significance_labels |>
  filter(
    Time_Window == "Baseline"
  )

baseline_plot <- ggplot(
  baseline_plot_summary,
  aes(
    x = Variant,
    y = Mean
  )
) +
  geom_col(
    width = 0.68,
    fill = "grey82",
    colour = "black",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    width = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    data = baseline_plot_data,
    aes(
      x = Variant,
      y = Analysis_Value,
      colour = Sweep_factor
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.09,
      height = 0,
      seed = 123
    ),
    size = 2.3,
    alpha = 0.85
  ) +
  geom_text(
    data = baseline_significance,
    aes(
      x = Variant,
      y = Label_Y,
      label = Plot_Significance
    ),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  facet_wrap(
    ~Parameter,
    scales = "free_y",
    ncol = 2,
    labeller = as_labeller(
      parameter_labels
    )
  ) +
  scale_colour_viridis_d(
    option = "D",
    end = 0.9,
    drop = TRUE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.18
      )
    )
  ) +
  labs(
    title = "Baseline extracellular field-potential phenotype",
    subtitle = paste(
      "Each point represents the average of retained wells within one sweep;",
      "one-way ANOVA with Dunnett comparisons against WT"
    ),
    x = NULL,
    y = "Absolute CardioExcyte measurement",
    colour = "Sweep"
  ) +
  thesis_theme

print(
  baseline_plot
)

ggsave(
  filename = file.path(
    output_directory,
    "Figure_1_baseline_EFP_vs_WT.pdf"
  ),
  plot = baseline_plot,
  width = 12,
  height = 10
)

ggsave(
  filename = file.path(
    output_directory,
    "Figure_1_baseline_EFP_vs_WT.png"
  ),
  plot = baseline_plot,
  width = 12,
  height = 10,
  dpi = 600
)


# ============================================================
# 34. FIGURE 2: POST-NORADRENALINE EFP RESPONSE
# ============================================================

post_na_plot_data <- variant_sweep_data |>
  filter(
    Time_Window != "Baseline"
  )

post_na_plot_summary <- figure_summary |>
  filter(
    Time_Window != "Baseline"
  )

post_na_significance <- significance_labels |>
  filter(
    Time_Window != "Baseline"
  )

post_na_plot <- ggplot(
  post_na_plot_summary,
  aes(
    x = Variant,
    y = Mean
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey35"
  ) +
  geom_col(
    width = 0.68,
    fill = "grey82",
    colour = "black",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    width = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    data = post_na_plot_data,
    aes(
      x = Variant,
      y = Analysis_Value,
      colour = Sweep_factor
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.09,
      height = 0,
      seed = 123
    ),
    size = 2.1,
    alpha = 0.85
  ) +
  geom_text(
    data = post_na_significance,
    aes(
      x = Variant,
      y = Label_Y,
      label = Plot_Significance
    ),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_grid(
    Parameter ~ Time_Window,
    scales = "free_y",
    labeller = labeller(
      Parameter = as_labeller(
        parameter_labels
      )
    )
  ) +
  scale_colour_viridis_d(
    option = "D",
    end = 0.9,
    drop = TRUE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.18
      )
    )
  ) +
  labs(
    title = "Extracellular field-potential response following noradrenaline",
    subtitle = paste(
      "Retained wells were averaged within each sweep;",
      "each variant was compared with WT in the same time window"
    ),
    x = NULL,
    y = "Fold change relative to baseline",
    colour = "Sweep"
  ) +
  thesis_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),

    strip.text = element_text(
      face = "bold",
      size = 9
    )
  )

print(
  post_na_plot
)

ggsave(
  filename = file.path(
    output_directory,
    "Figure_2_post_NA_EFP_response_vs_WT.pdf"
  ),
  plot = post_na_plot,
  width = 16,
  height = 14
)

ggsave(
  filename = file.path(
    output_directory,
    "Figure_2_post_NA_EFP_response_vs_WT.png"
  ),
  plot = post_na_plot,
  width = 16,
  height = 14,
  dpi = 600
)


# ============================================================
# 35. CREATE ONE POST-NA FIGURE FOR EACH PARAMETER
# ============================================================

make_parameter_plot <- function(parameter_name) {

  parameter_data <- post_na_plot_data |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_summary <- post_na_plot_summary |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_significance <- post_na_significance |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_plot <- ggplot(
    parameter_summary,
    aes(
      x = Variant,
      y = Mean
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.7,
      colour = "grey35"
    ) +
    geom_col(
      width = 0.68,
      fill = "grey82",
      colour = "black",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(
        ymin = Mean - SEM,
        ymax = Mean + SEM
      ),
      width = 0.16,
      linewidth = 0.65
    ) +
    geom_point(
      data = parameter_data,
      aes(
        x = Variant,
        y = Analysis_Value,
        colour = Sweep_factor
      ),
      inherit.aes = FALSE,
      position = position_jitter(
        width = 0.09,
        height = 0,
        seed = 123
      ),
      size = 2.3,
      alpha = 0.85
    ) +
    geom_text(
      data = parameter_significance,
      aes(
        x = Variant,
        y = Label_Y,
        label = Plot_Significance
      ),
      inherit.aes = FALSE,
      size = 4.2,
      fontface = "bold"
    ) +
    facet_wrap(
      ~Time_Window,
      scales = "free_y",
      nrow = 1
    ) +
    scale_colour_viridis_d(
      option = "D",
      end = 0.9,
      drop = TRUE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.04,
          0.18
        )
      )
    ) +
    labs(
      title = paste(
        parameter_labels[[parameter_name]],
        "following noradrenaline"
      ),
      subtitle = paste(
        "Each point represents the average of retained wells within one sweep;",
        "Dunnett-adjusted comparisons against WT"
      ),
      x = NULL,
      y = "Fold change relative to baseline",
      colour = "Sweep"
    ) +
    thesis_theme

  safe_parameter_name <- str_replace_all(
    parameter_name,
    "[^A-Za-z0-9]+",
    "_"
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "Figure_",
        safe_parameter_name,
        "_post_NA_vs_WT.pdf"
      )
    ),
    plot = parameter_plot,
    width = 15,
    height = 6.5
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "Figure_",
        safe_parameter_name,
        "_post_NA_vs_WT.png"
      )
    ),
    plot = parameter_plot,
    width = 15,
    height = 6.5,
    dpi = 600
  )

  return(
    parameter_plot
  )
}

parameter_plots <- lapply(
  parameter_order,
  make_parameter_plot
)

names(parameter_plots) <- parameter_order


# ============================================================
# 36. CREATE ONE BASELINE FIGURE FOR EACH PARAMETER
# ============================================================

make_baseline_parameter_plot <- function(parameter_name) {

  parameter_data <- baseline_plot_data |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_summary <- baseline_plot_summary |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_significance <- baseline_significance |>
    filter(
      Parameter == parameter_name
    ) |>
    droplevels()

  parameter_plot <- ggplot(
    parameter_summary,
    aes(
      x = Variant,
      y = Mean
    )
  ) +
    geom_col(
      width = 0.68,
      fill = "grey82",
      colour = "black",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(
        ymin = Mean - SEM,
        ymax = Mean + SEM
      ),
      width = 0.16,
      linewidth = 0.65
    ) +
    geom_point(
      data = parameter_data,
      aes(
        x = Variant,
        y = Analysis_Value,
        colour = Sweep_factor
      ),
      inherit.aes = FALSE,
      position = position_jitter(
        width = 0.09,
        height = 0,
        seed = 123
      ),
      size = 2.3,
      alpha = 0.85
    ) +
    geom_text(
      data = parameter_significance,
      aes(
        x = Variant,
        y = Label_Y,
        label = Plot_Significance
      ),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold"
    ) +
    scale_colour_viridis_d(
      option = "D",
      end = 0.9,
      drop = TRUE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.04,
          0.18
        )
      )
    ) +
    labs(
      title = paste(
        "Baseline",
        parameter_labels[[parameter_name]]
      ),
      subtitle = paste(
        "Each point represents the average of retained wells within one sweep;",
        "Dunnett-adjusted comparisons against WT"
      ),
      x = NULL,
      y = parameter_labels[[parameter_name]],
      colour = "Sweep"
    ) +
    thesis_theme

  safe_parameter_name <- str_replace_all(
    parameter_name,
    "[^A-Za-z0-9]+",
    "_"
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "Baseline_",
        safe_parameter_name,
        "_vs_WT.pdf"
      )
    ),
    plot = parameter_plot,
    width = 10,
    height = 7
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        "Baseline_",
        safe_parameter_name,
        "_vs_WT.png"
      )
    ),
    plot = parameter_plot,
    width = 10,
    height = 7,
    dpi = 600
  )

  return(
    parameter_plot
  )
}

baseline_parameter_plots <- lapply(
  parameter_order,
  make_baseline_parameter_plot
)

names(baseline_parameter_plots) <- parameter_order


# ============================================================
# 37. CREATE SWEEP-TIMING CHECK FIGURE
# ============================================================

time_plot <- ggplot(
  sweep_time_lookup,
  aes(
    x = Sweep_number,
    y = Time_Post_NA
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 24,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = 48,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_point(
    data = sweep_time_lookup |>
      filter(
        Sweep_number %in% c(
          baseline_sweeps,
          acute_sweeps,
          sweeps_24h,
          sweeps_48h
        )
      ),
    size = 3.8,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  scale_x_continuous(
    breaks = sweep_time_lookup$Sweep_number
  ) +
  labs(
    title = "CardioExcyte sweep timing relative to noradrenaline",
    subtitle = paste0(
      "Noradrenaline addition estimated at ",
      round(
        estimated_na_time,
        3
      ),
      " recorded hours"
    ),
    x = "Sweep",
    y = "Hours post-noradrenaline"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

print(
  time_plot
)

ggsave(
  filename = file.path(
    output_directory,
    "Sweep_timing_relative_to_NA.pdf"
  ),
  plot = time_plot,
  width = 11,
  height = 7
)

ggsave(
  filename = file.path(
    output_directory,
    "Sweep_timing_relative_to_NA.png"
  ),
  plot = time_plot,
  width = 11,
  height = 7,
  dpi = 600
)


# ============================================================
# 38. SAVE PROCESSED DATA TABLES
# ============================================================

write.csv(
  sweep_time_lookup,
  file.path(
    output_directory,
    "Sweep_time_lookup_post_NA.csv"
  ),
  row.names = FALSE
)

write.csv(
  well_baseline,
  file.path(
    output_directory,
    "Retained_baseline_means_by_well.csv"
  ),
  row.names = FALSE
)

write.csv(
  analysis_data,
  file.path(
    output_directory,
    "Well_level_normalised_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  variant_sweep_data,
  file.path(
    output_directory,
    "Variant_sweep_data_after_well_averaging.csv"
  ),
  row.names = FALSE
)

write.csv(
  figure_summary,
  file.path(
    output_directory,
    "Figure_mean_and_SEM.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 39. SAVE COMPLETE RESULTS OBJECT
# ============================================================

cardioexcyte_final_results <- list(
  estimated_na_time = estimated_na_time,

  baseline_sweeps = baseline_sweeps,

  acute_sweeps = acute_sweeps,

  sweeps_24h = sweeps_24h,

  sweeps_48h = sweeps_48h,

  selected_windows = selected_windows,

  sweep_time_lookup = sweep_time_lookup,

  well_baseline = well_baseline,

  incomplete_baselines = incomplete_baselines,

  well_level_data = analysis_data,

  variant_sweep_data = variant_sweep_data,

  figure_summary = figure_summary,

  anova_table = anova_table,

  dunnett_table = dunnett_table,

  baseline_plot = baseline_plot,

  post_na_plot = post_na_plot,

  parameter_plots = parameter_plots,

  baseline_parameter_plots = baseline_parameter_plots
)

saveRDS(
  cardioexcyte_final_results,
  file.path(
    output_directory,
    "CardioExcyte_final_EFP_thesis_results.rds"
  )
)


# ============================================================
# 40. FINAL SUMMARY
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FINAL CARDIOEXCYTE EFP ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Baseline sweeps: ",
  paste(
    sprintf(
      "%03d",
      baseline_sweeps
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Acute sweeps: ",
  paste(
    sprintf(
      "%03d",
      acute_sweeps
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "24-hour sweeps: ",
  ifelse(
    length(sweeps_24h) > 0,
    paste(
      sprintf(
        "%03d",
        sweeps_24h
      ),
      collapse = ", "
    ),
    "No sweeps found"
  ),
  "\n",
  sep = ""
)

cat(
  "48-hour sweeps: ",
  ifelse(
    length(sweeps_48h) > 0,
    paste(
      sprintf(
        "%03d",
        sweeps_48h
      ),
      collapse = ", "
    ),
    "No sweeps found"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Correct parameter interpretation:\n"
)

cat(
  "- Beat rate: BPM\n"
)

cat(
  "- Amplitude: extracellular field-potential amplitude, not force\n"
)

cat(
  "- Beat-rate regularity: beat-to-beat timing variability\n"
)

cat(
  "- Field-potential duration: extracellular waveform duration\n"
)

cat(
  "- Downstroke velocity: descending EFP waveform velocity\n\n"
)

cat(
  "Missing values caused by outlier removal were allowed.\n"
)

cat(
  "The retained baseline mean was calculated separately for each well.\n"
)

cat(
  "The two wells were averaged within each variant and sweep.\n"
)

cat(
  paste0(
    "Separate one-way ANOVAs were run for each parameter and ",
    "time window.\n"
  )
)

cat(
  paste0(
    "Dunnett-adjusted comparisons tested every variant against ",
    "WT within the same time window.\n\n"
  )
)

cat(
  "Results saved to:\n",
  normalizePath(output_directory),
  "\n",
  sep = ""
)
