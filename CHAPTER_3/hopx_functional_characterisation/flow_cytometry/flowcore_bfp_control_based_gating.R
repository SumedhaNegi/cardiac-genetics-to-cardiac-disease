# =============================================================================
# EXAMINER GUIDE
# Script: flowcore_bfp_control_based_gating.R
# Purpose: Flow-cytometry/BFP analysis using control-based gating; the source documents the gating sequence and downstream plots.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# Imports FCS files, applies control-based fluorescence gating and derives flow-cytometry populations and plots.
#

# ============================================================
# CONTROL-BASED FACS GATING AND FINAL BFP OVERLAY
#
# Gating sequence:
# P1 = main cells: FSC-A vs SSC-A
# P3 = singlets: FSC-A vs FSC-H
# P4 = viable cells: FSC-A vs R670/30-A
# P5 = BFP positive: V450/50-A
#
# Final output:
# - Highest BFP-positive sample from A1–A8
# - Filled control distribution
# - One PDF only
# ============================================================


# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install(
    "flowCore",
    ask = FALSE,
    update = FALSE
  )
}

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "purrr",
  "stringr",
  "tibble"
)

missing_packages <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(flowCore)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)


# ============================================================
# 2. SETTINGS
# ============================================================

fcs_directory <- paste0(
  "/Users/uqsnegi/Downloads/",
  "clinvar_virus_titration"
)

output_directory <- file.path(
  fcs_directory,
  "Control_based_gating"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# Exact channel names from the FCS files
fsc_a_channel <- "FSC-A"
ssc_a_channel <- "SSC-A"
fsc_h_channel <- "FSC-H"
viability_channel <- "R670/30-A"
bfp_channel <- "V450/50-A"

# Filename text used to identify the negative control
control_pattern <- "Control"

# Exclude unwanted samples
exclude_pattern <- "water"

# P5 threshold:
# 0.5% of the gated control is permitted above the gate
# Equivalent to the 99.5th percentile
control_false_positive_percent <- 0.5

# Density-curve smoothing
density_smoothing <- 1.5

# Number of events displayed while drawing gates
maximum_display_points <- 25000

# Maximum number of events used per density curve
maximum_plot_events <- 50000

# Use generic font to avoid macOS PDF font errors
figure_font <- "sans"

# Final PDF filename
final_pdf <- file.path(
  output_directory,
  "Highest_BFP_positive_A1_to_A8_overlay.pdf"
)


# ============================================================
# 3. FIND FCS FILES
# ============================================================

fcs_files <- list.files(
  path = fcs_directory,
  pattern = "\\.fcs$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

fcs_files <- fcs_files[
  !grepl(
    exclude_pattern,
    basename(fcs_files),
    ignore.case = TRUE
  )
]

if (length(fcs_files) == 0) {
  stop(
    "No FCS files were found in:\n",
    fcs_directory
  )
}

message(
  "Number of FCS files found: ",
  length(fcs_files)
)


# ============================================================
# 4. CREATE SAMPLE NAMES
# ============================================================

clean_sample_name <- function(file_path) {

  basename(file_path) |>
    tools::file_path_sans_ext() |>
    stringr::str_remove(
      stringr::regex(
        "^Cardios[_-]?",
        ignore_case = TRUE
      )
    ) |>
    stringr::str_squish()
}

sample_names <- vapply(
  fcs_files,
  clean_sample_name,
  character(1)
)

sample_names <- make.unique(
  sample_names,
  sep = "_"
)


# ============================================================
# 5. READ FCS FILES
# ============================================================

message("Reading FCS files...")

flow_frames <- lapply(
  fcs_files,
  function(file_path) {

    tryCatch(
      flowCore::read.FCS(
        filename = file_path,
        transformation = FALSE,
        truncate_max_range = FALSE
      ),
      error = function(e) {

        warning(
          "Could not read:\n",
          file_path,
          "\n",
          conditionMessage(e)
        )

        NULL
      }
    )
  }
)

names(flow_frames) <- sample_names

successful_files <- !vapply(
  flow_frames,
  is.null,
  logical(1)
)

flow_frames <- flow_frames[
  successful_files
]

if (length(flow_frames) == 0) {
  stop("None of the FCS files could be read.")
}

message(
  "Successfully read ",
  length(flow_frames),
  " FCS files."
)


# ============================================================
# 6. IDENTIFY CONTROL SAMPLE
# ============================================================

control_matches <- grep(
  pattern = control_pattern,
  x = names(flow_frames),
  ignore.case = TRUE
)

if (length(control_matches) == 0) {

  stop(
    "No control sample was found.\n\n",
    "Available samples:\n",
    paste(
      names(flow_frames),
      collapse = "\n"
    )
  )
}

control_sample <- names(flow_frames)[
  control_matches[1]
]

message(
  "Control sample: ",
  control_sample
)


# ============================================================
# 7. CHECK REQUIRED CHANNELS
# ============================================================

required_channels <- c(
  fsc_a_channel,
  ssc_a_channel,
  fsc_h_channel,
  viability_channel,
  bfp_channel
)

common_channels <- Reduce(
  intersect,
  lapply(
    flow_frames,
    function(flow_frame) {
      colnames(
        flowCore::exprs(flow_frame)
      )
    }
  )
)

missing_channels <- setdiff(
  required_channels,
  common_channels
)

if (length(missing_channels) > 0) {

  print(common_channels)

  stop(
    "These required channels were not found:\n",
    paste(
      missing_channels,
      collapse = "\n"
    )
  )
}


# ============================================================
# 8. APPLY A COMMON LOGICLE TRANSFORMATION
# ============================================================

control_untransformed <- flow_frames[
  [control_sample]
]

common_logicle <- flowCore::estimateLogicle(
  control_untransformed,
  channels = c(
    viability_channel,
    bfp_channel
  )
)

transformed_frames <- lapply(
  flow_frames,
  function(flow_frame) {

    flowCore::transform(
      flow_frame,
      common_logicle
    )
  }
)

names(transformed_frames) <- names(
  flow_frames
)

control_frame <- transformed_frames[
  [control_sample]
]

stopifnot(
  methods::is(
    control_frame,
    "flowFrame"
  )
)


# ============================================================
# 9. INTERACTIVE PLOT WINDOW
# ============================================================

open_gate_window <- function(
    window_title,
    width = 7,
    height = 7
) {

  if (
    identical(
      Sys.info()[["sysname"]],
      "Darwin"
    )
  ) {

    tryCatch(
      grDevices::quartz(
        title = window_title,
        width = width,
        height = height
      ),
      error = function(e) {

        grDevices::dev.new(
          width = width,
          height = height
        )
      }
    )

  } else {

    grDevices::dev.new(
      width = width,
      height = height
    )
  }
}


# ============================================================
# 10. PREPARE DATA FOR INTERACTIVE GATING
# ============================================================

get_display_data <- function(
    flow_frame,
    x_channel,
    y_channel,
    maximum_points = maximum_display_points
) {

  plot_matrix <- flowCore::exprs(
    flow_frame
  )[
    ,
    c(
      x_channel,
      y_channel
    ),
    drop = FALSE
  ]

  keep_rows <- (
    is.finite(plot_matrix[, 1]) &
      is.finite(plot_matrix[, 2])
  )

  plot_matrix <- plot_matrix[
    keep_rows,
    ,
    drop = FALSE
  ]

  if (nrow(plot_matrix) > maximum_points) {

    set.seed(123)

    selected_rows <- sample(
      seq_len(
        nrow(plot_matrix)
      ),
      maximum_points
    )

    plot_matrix <- plot_matrix[
      selected_rows,
      ,
      drop = FALSE
    ]
  }

  plot_matrix
}


# ============================================================
# 11. FUNCTION TO DRAW A POLYGON GATE
# ============================================================

draw_polygon_gate <- function(
    flow_frame,
    x_channel,
    y_channel,
    gate_name,
    plot_title,
    instructions
) {

  plot_matrix <- get_display_data(
    flow_frame = flow_frame,
    x_channel = x_channel,
    y_channel = y_channel
  )

  open_gate_window(
    plot_title
  )

  current_device <- grDevices::dev.cur()

  plot(
    plot_matrix[, 1],
    plot_matrix[, 2],
    pch = 16,
    cex = 0.25,
    col = grDevices::adjustcolor(
      "#2166AC",
      alpha.f = 0.25
    ),
    xlab = x_channel,
    ylab = y_channel,
    main = plot_title
  )

  mtext(
    instructions,
    side = 3,
    line = 0.2,
    cex = 0.75,
    col = "red"
  )

  message("\n", gate_name)
  message(instructions)
  message(
    "Click around the population, ",
    "then press Esc."
  )

  selected_points <- graphics::locator(
    n = 100,
    type = "l",
    col = "red",
    lwd = 2
  )

  if (
    is.null(selected_points) ||
      length(selected_points$x) < 3
  ) {

    if (
      grDevices::dev.cur() ==
        current_device
    ) {
      grDevices::dev.off()
    }

    stop(
      gate_name,
      " requires at least three clicked points."
    )
  }

  gate_coordinates <- cbind(
    selected_points$x,
    selected_points$y
  )

  colnames(gate_coordinates) <- c(
    x_channel,
    y_channel
  )

  polygon(
    selected_points$x,
    selected_points$y,
    border = "red",
    lwd = 2
  )

  Sys.sleep(0.5)

  if (
    grDevices::dev.cur() ==
      current_device
  ) {
    grDevices::dev.off()
  }

  flowCore::polygonGate(
    filterId = gate_name,
    .gate = gate_coordinates
  )
}


# ============================================================
# 12. FUNCTION TO SELECT THE P4 THRESHOLD
# ============================================================

select_horizontal_threshold <- function(
    flow_frame,
    x_channel,
    y_channel,
    plot_title,
    instructions
) {

  plot_matrix <- get_display_data(
    flow_frame = flow_frame,
    x_channel = x_channel,
    y_channel = y_channel
  )

  open_gate_window(
    plot_title
  )

  current_device <- grDevices::dev.cur()

  plot(
    plot_matrix[, 1],
    plot_matrix[, 2],
    pch = 16,
    cex = 0.25,
    col = grDevices::adjustcolor(
      "#2166AC",
      alpha.f = 0.25
    ),
    xlab = x_channel,
    ylab = y_channel,
    main = plot_title
  )

  mtext(
    instructions,
    side = 3,
    line = 0.2,
    cex = 0.75,
    col = "red"
  )

  message("\nP4 viability gate")
  message(instructions)
  message(
    "Click once between the low-R670 ",
    "and high-R670 populations."
  )

  selected_point <- graphics::locator(
    n = 1,
    type = "p",
    pch = 4,
    col = "red",
    lwd = 2,
    cex = 1.5
  )

  if (is.null(selected_point)) {

    if (
      grDevices::dev.cur() ==
        current_device
    ) {
      grDevices::dev.off()
    }

    stop(
      "No P4 threshold was selected."
    )
  }

  threshold <- as.numeric(
    selected_point$y
  )

  abline(
    h = threshold,
    col = "red",
    lwd = 2,
    lty = 2
  )

  Sys.sleep(0.5)

  if (
    grDevices::dev.cur() ==
      current_device
  ) {
    grDevices::dev.off()
  }

  threshold
}


# ============================================================
# 13. RECTANGLE-GATE HELPER
# ============================================================

make_rectangle_gate <- function(
    channel,
    lower,
    upper,
    gate_name
) {

  gate_arguments <- c(
    list(
      filterId = gate_name
    ),
    stats::setNames(
      list(
        c(
          lower,
          upper
        )
      ),
      channel
    )
  )

  do.call(
    flowCore::rectangleGate,
    gate_arguments
  )
}


# ============================================================
# 14. DRAW P1 ON THE CONTROL
# ============================================================

P1_gate <- draw_polygon_gate(
  flow_frame = control_frame,
  x_channel = fsc_a_channel,
  y_channel = ssc_a_channel,
  gate_name = "P1",
  plot_title = paste0(
    control_sample,
    ": P1 main cell population"
  ),
  instructions = paste0(
    "Select the broad main cell cloud; ",
    "exclude debris near the bottom-left"
  )
)

control_P1 <- flowCore::Subset(
  control_frame,
  P1_gate
)

message(
  "P1 retained ",
  nrow(
    flowCore::exprs(control_P1)
  ),
  " control events."
)


# ============================================================
# 15. DRAW P3 ON THE CONTROL
# ============================================================

P3_gate <- draw_polygon_gate(
  flow_frame = control_P1,
  x_channel = fsc_a_channel,
  y_channel = fsc_h_channel,
  gate_name = "P3",
  plot_title = paste0(
    control_sample,
    ": P3 singlets"
  ),
  instructions = paste0(
    "Draw a narrow polygon around ",
    "the diagonal singlet band"
  )
)

control_P3 <- flowCore::Subset(
  control_P1,
  P3_gate
)

message(
  "P3 retained ",
  nrow(
    flowCore::exprs(control_P3)
  ),
  " control events."
)


# ============================================================
# 16. SET P4 ON THE CONTROL
# ============================================================

P4_threshold <- select_horizontal_threshold(
  flow_frame = control_P3,
  x_channel = fsc_a_channel,
  y_channel = viability_channel,
  plot_title = paste0(
    control_sample,
    ": P4 viability gate"
  ),
  instructions = paste0(
    "Click between the low-R670 viable ",
    "population and high-R670 dead population"
  )
)

P4_gate <- make_rectangle_gate(
  channel = viability_channel,
  lower = -Inf,
  upper = P4_threshold,
  gate_name = "P4"
)

control_P4 <- flowCore::Subset(
  control_P3,
  P4_gate
)

message(
  "P4 threshold: ",
  round(
    P4_threshold,
    4
  )
)

message(
  "P4 retained ",
  nrow(
    flowCore::exprs(control_P4)
  ),
  " viable control events."
)


# ============================================================
# 17. SET P5 AUTOMATICALLY FROM THE CONTROL
# ============================================================

control_bfp_values <- flowCore::exprs(
  control_P4
)[
  ,
  bfp_channel
]

control_bfp_values <- as.numeric(
  control_bfp_values
)

control_bfp_values <- control_bfp_values[
  is.finite(control_bfp_values)
]

if (length(control_bfp_values) < 20) {

  stop(
    "Too few P4-gated control events ",
    "remained to calculate P5."
  )
}

control_negative_fraction <- (
  1 -
    control_false_positive_percent / 100
)

P5_threshold <- as.numeric(
  stats::quantile(
    control_bfp_values,
    probs = control_negative_fraction,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
)

P5_gate <- make_rectangle_gate(
  channel = bfp_channel,
  lower = P5_threshold,
  upper = Inf,
  gate_name = "P5"
)

control_P5 <- flowCore::Subset(
  control_P4,
  P5_gate
)

control_positive_percent <- (
  100 *
    nrow(
      flowCore::exprs(control_P5)
    ) /
    nrow(
      flowCore::exprs(control_P4)
    )
)

cat("\n====================================\n")
cat(
  "P5 threshold:",
  round(P5_threshold, 6),
  "\n"
)
cat(
  "Control BFP-positive:",
  round(control_positive_percent, 3),
  "%\n"
)
cat("====================================\n\n")


# ============================================================
# 18. SAVE THE CONTROL-BASED GATES
# ============================================================

gate_information <- list(
  P1_gate = P1_gate,
  P3_gate = P3_gate,
  P4_gate = P4_gate,
  P5_gate = P5_gate,
  P4_threshold = P4_threshold,
  P5_threshold = P5_threshold,
  control_sample = control_sample,
  control_false_positive_percent =
    control_false_positive_percent,
  logicle_transformation =
    common_logicle,
  channels = list(
    FSC_A = fsc_a_channel,
    SSC_A = ssc_a_channel,
    FSC_H = fsc_h_channel,
    viability = viability_channel,
    BFP = bfp_channel
  )
)

saveRDS(
  gate_information,
  file.path(
    output_directory,
    "Control_based_gates.rds"
  )
)


# ============================================================
# 19. APPLY THE SAME GATES TO EVERY SAMPLE
# ============================================================

apply_control_gates <- function(
    flow_frame
) {

  P1_population <- flowCore::Subset(
    flow_frame,
    P1_gate
  )

  P3_population <- flowCore::Subset(
    P1_population,
    P3_gate
  )

  P4_population <- flowCore::Subset(
    P3_population,
    P4_gate
  )

  P5_population <- flowCore::Subset(
    P4_population,
    P5_gate
  )

  list(
    All = flow_frame,
    P1 = P1_population,
    P3 = P3_population,
    P4 = P4_population,
    P5 = P5_population
  )
}

message(
  "Applying control-based gates ",
  "to all samples..."
)

gated_samples <- lapply(
  transformed_frames,
  apply_control_gates
)

names(gated_samples) <- names(
  transformed_frames
)


# ============================================================
# 20. CREATE GATING SUMMARY
# ============================================================

count_events <- function(flow_frame) {

  nrow(
    flowCore::exprs(flow_frame)
  )
}

safe_percent <- function(
    numerator,
    denominator
) {

  if (denominator == 0) {
    return(NA_real_)
  }

  100 * numerator / denominator
}

gating_summary <- purrr::imap_dfr(
  gated_samples,
  function(populations, sample_name) {

    n_all <- count_events(
      populations$All
    )

    n_P1 <- count_events(
      populations$P1
    )

    n_P3 <- count_events(
      populations$P3
    )

    n_P4 <- count_events(
      populations$P4
    )

    n_P5 <- count_events(
      populations$P5
    )

    tibble::tibble(
      Sample = sample_name,

      All_events = n_all,

      P1_events = n_P1,
      P1_percent_total = safe_percent(
        n_P1,
        n_all
      ),

      P3_events = n_P3,
      P3_percent_parent = safe_percent(
        n_P3,
        n_P1
      ),

      P4_events = n_P4,
      P4_percent_parent = safe_percent(
        n_P4,
        n_P3
      ),

      P5_events = n_P5,
      BFP_positive_percent = safe_percent(
        n_P5,
        n_P4
      )
    )
  }
)

print(
  gating_summary,
  n = Inf
)

write.csv(
  gating_summary,
  file.path(
    output_directory,
    "Control_based_gating_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SELECT HIGHEST BFP-POSITIVE SAMPLE FROM A1–A8
# ============================================================

selected_samples <- gating_summary |>
  dplyr::mutate(
    Tube = stringr::str_extract(
      Sample,
      "^A[1-8](?=_)"
    )
  ) |>
  dplyr::filter(
    !is.na(Tube),
    !is.na(BFP_positive_percent)
  ) |>
  dplyr::group_by(Tube) |>
  dplyr::slice_max(
    order_by = BFP_positive_percent,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Tube_number = as.numeric(
      stringr::str_remove(
        Tube,
        "A"
      )
    )
  ) |>
  dplyr::arrange(
    Tube_number
  )

cat(
  "\nHighest BFP-positive sample ",
  "selected for each tube:\n\n"
)

print(
  selected_samples |>
    dplyr::select(
      Tube,
      Sample,
      BFP_positive_percent
    ),
  n = Inf
)

missing_tubes <- setdiff(
  paste0(
    "A",
    1:8
  ),
  selected_samples$Tube
)

if (length(missing_tubes) > 0) {

  warning(
    "No selected sample was found for: ",
    paste(
      missing_tubes,
      collapse = ", "
    )
  )
}


# ============================================================
# 22. EXTRACT P4-GATED CONTROL VALUES
# ============================================================

set.seed(123)

control_values <- flowCore::exprs(
  gated_samples[
    [control_sample]
  ]$P4
)[
  ,
  bfp_channel
]

control_values <- as.numeric(
  control_values
)

control_values <- control_values[
  is.finite(control_values)
]

if (
  length(control_values) >
    maximum_plot_events
) {

  control_values <- sample(
    control_values,
    maximum_plot_events
  )
}

if (length(control_values) < 10) {

  stop(
    "Too few P4-gated control events ",
    "were available for plotting."
  )
}

control_plot_data <- tibble::tibble(
  BFP = control_values
)


# ============================================================
# 23. EXTRACT THE EIGHT SELECTED SAMPLE DISTRIBUTIONS
# ============================================================

sample_plot_data <- purrr::map2_dfr(
  selected_samples$Sample,
  selected_samples$Tube,
  function(sample_name, tube_name) {

    if (
      !sample_name %in%
        names(gated_samples)
    ) {

      stop(
        "Sample missing from gated_samples: ",
        sample_name
      )
    }

    values <- flowCore::exprs(
      gated_samples[
        [sample_name]
      ]$P4
    )[
      ,
      bfp_channel
    ]

    values <- as.numeric(values)

    values <- values[
      is.finite(values)
    ]

    if (
      length(values) >
        maximum_plot_events
    ) {

      values <- sample(
        values,
        maximum_plot_events
      )
    }

    tibble::tibble(
      BFP = values,
      Tube = tube_name,
      Sample = sample_name
    )
  }
)

sample_plot_data$Tube <- factor(
  sample_plot_data$Tube,
  levels = paste0(
    "A",
    1:8
  )
)


# ============================================================
# 24. CALCULATE COMMON AXIS LIMITS
# ============================================================

all_plot_values <- c(
  control_plot_data$BFP,
  sample_plot_data$BFP
)

all_plot_values <- all_plot_values[
  is.finite(all_plot_values)
]

x_limits <- as.numeric(
  stats::quantile(
    all_plot_values,
    probs = c(
      0.001,
      0.999
    ),
    na.rm = TRUE,
    names = FALSE
  )
)

x_min <- x_limits[1]
x_max <- x_limits[2]

stopifnot(
  length(x_limits) == 2,
  is.finite(x_min),
  is.finite(x_max),
  x_min < x_max
)


# ============================================================
# 25. CREATE FINAL OVERLAY
# ============================================================

final_overlay <- ggplot2::ggplot() +

  # Filled negative-control distribution
  ggplot2::geom_density(
    data = control_plot_data,
    ggplot2::aes(
      x = BFP,
      y = after_stat(scaled)
    ),
    adjust = density_smoothing,
    fill = "grey75",
    colour = "grey35",
    alpha = 0.70,
    linewidth = 1,
    trim = TRUE,
    n = 1024
  ) +

  # Highest A1–A8 samples
  ggplot2::geom_density(
    data = sample_plot_data,
    ggplot2::aes(
      x = BFP,
      y = after_stat(scaled),
      colour = Tube,
      group = Tube
    ),
    adjust = density_smoothing,
    fill = NA,
    linewidth = 1.35,
    trim = TRUE,
    n = 1024
  ) +

  # Fixed P5 threshold
  ggplot2::geom_vline(
    xintercept = P5_threshold,
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black"
  ) +

  ggplot2::scale_colour_brewer(
    palette = "Dark2",
    drop = FALSE
  ) +

  ggplot2::coord_cartesian(
    xlim = c(
      x_min,
      x_max
    ),
    ylim = c(
      0,
      1.05
    ),
    expand = FALSE
  ) +

  ggplot2::labs(
    x = paste0(
      bfp_channel,
      " fluorescence"
    ),
    y = "Normalized density",
    colour = NULL
  ) +

  ggplot2::theme_classic(
    base_size = 12,
    base_family = figure_font
  ) +

  ggplot2::theme(
    axis.text = ggplot2::element_text(
      colour = "black"
    ),
    axis.title = ggplot2::element_text(
      colour = "black"
    ),
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 10
    ),
    legend.key.width = grid::unit(
      1.2,
      "cm"
    ),
    plot.margin = ggplot2::margin(
      10,
      12,
      10,
      10
    )
  )

print(
  final_overlay
)


# ============================================================
# 26. SAVE PDF ONLY
# ============================================================

ggplot2::ggsave(
  filename = final_pdf,
  plot = final_overlay,
  device = grDevices::cairo_pdf,
  width = 180,
  height = 115,
  units = "mm",
  bg = "white"
)

cat(
  "\nFinal PDF saved to:\n",
  final_pdf,
  "\n"
)

cat(
  "\nGating summary saved to:\n",
  file.path(
    output_directory,
    "Control_based_gating_summary.csv"
  ),
  "\n"
)

cat(
  "\nGates saved to:\n",
  file.path(
    output_directory,
    "Control_based_gates.rds"
  ),
  "\n"
)
