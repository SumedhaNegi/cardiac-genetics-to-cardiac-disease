# =============================================================================
# EXAMINER GUIDE
# Script: Heritability_plots_M1_report.R
# Purpose: Visualise per-SNP heritability enrichment estimates and standard
#          errors from SBayesRC annotation result summaries.
# Provenance: Reconstructed from the plotting section of the historical script.
# Input: A tab-delimited, already-extracted results table with columns:
#        trait, annot, Per_SNP_enrichment, error.
# Output: Overall, LV, RV and outflow-tract trait plots as PDF files.
# Important: The original source referenced undefined objects case1, c1 and h1
#            and did not provide a reliable parser for .parSetRes logs. This
#            script therefore requires a prepared table; it does not guess the
#            undocumented raw log format or calculate new model estimates.
# =============================================================================

suppressPackageStartupMessages({library(ggplot2); library(data.table)})

# Supply the prepared table path as the first command-line argument or set
# HERITABILITY_RESULTS. The expected columns are listed above.
args <- commandArgs(trailingOnly = TRUE)
results_file <- if (length(args)) args[[1]] else Sys.getenv("HERITABILITY_RESULTS", "")
output_dir <- Sys.getenv("HERITABILITY_PLOT_DIR", "visualization")

if (!nzchar(results_file)) {
  message("No input table specified. Set HERITABILITY_RESULTS or run: ")
  message("Rscript Heritability_plots_M1_report.R path/to/prepared_per_snp_results.tsv")
  message("Expected columns: trait, annot, Per_SNP_enrichment, error")
} else {
  if (!file.exists(results_file)) stop("Results file not found: ", results_file)
  dat <- data.table::fread(results_file, data.table = FALSE)
  required <- c("trait", "annot", "Per_SNP_enrichment", "error")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Input table missing columns: ", paste(missing, collapse = ", "))
  dat$trait <- as.character(dat$trait)
  dat$annot <- as.character(dat$annot)
  dat$Per_SNP_enrichment <- as.numeric(dat$Per_SNP_enrichment)
  dat$error <- as.numeric(dat$error)
  dat <- dat[is.finite(dat$Per_SNP_enrichment) & is.finite(dat$error), , drop = FALSE]
  if (!nrow(dat)) stop("No valid enrichment estimates remain in input table.")
  dat$trait <- gsub("\\|\\./invnorm_|\\.ma.parSetRes:|\\./|_cleaned", "", dat$trait)
  dat$trait <- gsub("^.*invnorm_", "", dat$trait)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # The original plots compare the median enrichment for cardiac enhancer and
  # HOPX enhancer annotations. Use annotation labels to identify those groups.
  cardiac_rows <- grepl("cardiac", dat$annot, ignore.case = TRUE)
  hopx_rows <- grepl("hopx", dat$annot, ignore.case = TRUE)
  cardiac_median <- if (any(cardiac_rows)) median(dat$Per_SNP_enrichment[cardiac_rows], na.rm = TRUE) else NA_real_
  hopx_median <- if (any(hopx_rows)) median(dat$Per_SNP_enrichment[hopx_rows], na.rm = TRUE) else NA_real_

  make_plot <- function(d, title, filename, width = 8, height = 6) {
    if (!nrow(d)) {
      warning("No rows available for plot: ", title)
      return(invisible(NULL))
    }
    p <- ggplot(d, aes(x = Per_SNP_enrichment, y = reorder(trait, Per_SNP_enrichment), colour = annot)) +
      geom_point() +
      geom_errorbar(aes(xmin = Per_SNP_enrichment - error,
                        xmax = Per_SNP_enrichment + error), width = 0.25) +
      theme_classic() +
      labs(title = title, x = "Per-SNP heritability enrichment", y = "Trait", colour = "Annotation")
    if (is.finite(cardiac_median)) p <- p + geom_vline(xintercept = cardiac_median, colour = "grey50", alpha = 0.6)
    if (is.finite(hopx_median)) p <- p + geom_vline(xintercept = hopx_median, colour = "purple", alpha = 0.5)
    ggsave(file.path(output_dir, filename), plot = p, width = width, height = height)
    invisible(p)
  }

  # Broad phenotypic groupings are preserved from the source. These are
  # substring matches; inspect the selected traits in the saved input table.
  lv_terms <- c("bsa", "lv_indexed", "lvef", "svi", "ratio")
  rv_terms <- c("RVEDV_indexed", "ratio", "RVEF", "RVESV_indexed", "RV_SV_indexed")
  oft_terms <- c("aa_diam_indexed", "pa_ao")
  matches_any <- function(x, terms) Reduce(`|`, lapply(terms, function(z) grepl(z, x, ignore.case = TRUE)))

  make_plot(dat, "Heritability enrichment across all traits", "enhancer_all_traits.pdf", 8, 8)
  make_plot(dat[matches_any(dat$trait, lv_terms), ], "Left-ventricular traits", "lv_traits_enhancer.pdf", 8, 6)
  make_plot(dat[matches_any(dat$trait, rv_terms), ], "Right-ventricular traits", "rv_traits_enhancer.pdf", 8, 5)
  make_plot(dat[matches_any(dat$trait, oft_terms), ], "Outflow-tract / great-vessel traits", "oft_traits_enhancer.pdf", 8, 4)
  data.table::fwrite(dat, file.path(output_dir, "per_snp_results_used_for_plots.tsv"), sep = "\t")
  message("Plots and the exact plotting input table saved under: ", normalizePath(output_dir))
}

# Historical HPC commands (qsub/module load/grep) have been moved to a separate
# shell workflow. Raw .parSetRes extraction remains dependent on the precise
# SBayesRC output format and must be reconciled with the original report.
