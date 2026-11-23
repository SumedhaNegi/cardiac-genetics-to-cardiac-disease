# =============================================================================
# EXAMINER GUIDE
# Script: corr3.R
# Purpose: Plot a selected genetic-correlation matrix and count significant
#          correlations for each displayed row trait.
# Provenance: Reconstructed from the surviving body of the supplied script.
# Inputs: A long-format LDSC rg table with columns p1, p2, rg and q; two ordered
#         trait vectors. q is assumed to contain FDR-adjusted p-values.
# Outputs: PDF and PNG correlation plots, plus an optional count table.
# Important: The original function header and setup were missing. This version
#            makes those assumptions explicit rather than claiming recovery of
#            the original code verbatim.
# =============================================================================

suppressPackageStartupMessages({
  library(reshape2)
  library(corrplot)
})

#' Plot a rectangular genetic-correlation matrix and row-wise significance counts.
#'
#' @param rg Long table containing p1, p2, rg, q.
#' @param trait1 Traits to display as rows, in desired order.
#' @param trait2 Traits to display as columns, in desired order.
#' @param order corrplot ordering; use "original" to preserve supplied order.
#' @param output_file Optional output path. Extension controls PDF/PNG device.
#' @return Invisibly returns correlation matrix, FDR matrix and significance counts.
corrplot_nsquare_with_barplot <- function(rg, trait1, trait2,
                                          order = "original",
                                          output_file = NULL) {
  required <- c("p1", "p2", "rg", "q")
  missing <- setdiff(required, names(rg))
  if (length(missing)) stop("rg is missing required columns: ", paste(missing, collapse = ", "))
  if (!length(trait1) || !length(trait2)) stop("trait1 and trait2 must not be empty.")

  rg <- as.data.frame(rg)
  rg$p1 <- as.character(rg$p1); rg$p2 <- as.character(rg$p2)
  rg$rg <- suppressWarnings(as.numeric(rg$rg))
  rg$q <- suppressWarnings(as.numeric(rg$q))
  rg <- rg[is.finite(rg$rg) & is.finite(rg$q), , drop = FALSE]
  if (!nrow(rg)) stop("No finite rg/q rows remain after numeric conversion.")

  # dcast can aggregate duplicate trait pairs; mean reproduces the surviving source.
  rg_corr <- aggregate(rg$rg, list(p1 = rg$p1, p2 = rg$p2), mean, na.rm = TRUE)
  names(rg_corr)[3] <- "rg"
  rg_q <- aggregate(rg$q, list(p1 = rg$p1, p2 = rg$p2), mean, na.rm = TRUE)
  names(rg_q)[3] <- "q"

  make_matrix <- function(d, value) {
    wide <- reshape2::dcast(d, p1 ~ p2, value.var = value)
    m <- as.matrix(wide[, -1, drop = FALSE])
    rownames(m) <- wide$p1
    m
  }
  corr <- make_matrix(rg_corr, "rg")
  qmat <- make_matrix(rg_q, "q")

  # Align both matrices to the requested trait order; missing pairs remain NA.
  align <- function(m, rows, cols) {
    out <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                  dimnames = list(rows, cols))
    rr <- intersect(rows, rownames(m)); cc <- intersect(cols, colnames(m))
    if (length(rr) && length(cc)) out[rr, cc] <- m[rr, cc, drop = FALSE]
    out
  }
  corr <- align(corr, trait1, trait2)
  qmat <- align(qmat, trait1, trait2)
  corr <- pmax(-1, pmin(1, corr))
  # Missing rg is neutral for display; missing q is non-significant.
  corr[is.na(corr)] <- 0
  qmat[is.na(qmat)] <- 1
  if (identical(trait1, trait2) && nrow(corr) == ncol(corr)) {
    diag(corr) <- 1
    diag(qmat) <- -1
  }

  sig_counts <- rowSums(qmat < 0.05, na.rm = TRUE)
  draw_plots <- function() {
    graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
    par(mar = c(1, 1, 3, 1))
    corrplot::corrplot(corr, method = "square", order = order,
                       p.mat = qmat, sig.level = 0.05, insig = "label_sig",
                       pch = "*", pch.cex = 1.1, na.label = "square",
                       na.label.col = "grey30", tl.cex = 0.7,
                       title = "Genetic correlations (FDR-adjusted significance)")
    par(mar = c(9, 4, 2, 1))
    graphics::barplot(sig_counts,
      main = "Significant correlations per row trait (FDR < 0.05)",
      xlab = "Trait", ylab = "Number of significant correlations",
      las = 2, col = "skyblue")
  }
  if (is.null(output_file)) {
    draw_plots()
  } else {
    ext <- tolower(tools::file_ext(output_file))
    if (!ext %in% c("pdf", "png")) stop("output_file must end in .pdf or .png")
    if (ext == "pdf") grDevices::pdf(output_file, width = 16, height = 12, family = "Helvetica")
    if (ext == "png") grDevices::png(output_file, width = 16, height = 12, units = "in", res = 300)
    draw_plots()
    grDevices::dev.off()
  }
  invisible(list(correlation = corr, fdr = qmat, significant_counts = sig_counts))
}

# Example (edit paths and object names to match the original LDSC results):
# rg <- read.delim("rg_results.tsv")
# traitlist1 <- read.delim("traitlist_rows.tsv")
# traitlist2 <- read.delim("traitlist_columns.tsv")
# corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT,
#                               output_file = "corrplot2_rg_disease_with_barplot.pdf")
