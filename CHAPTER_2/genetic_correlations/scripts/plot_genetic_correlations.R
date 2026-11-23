# =============================================================================
# EXAMINER GUIDE
# Script: plot_genetic_correlations.R
# Purpose: Supports summary-statistic preparation, LDSC genetic-correlation analysis or visualisation; see chapter workflow documentation.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Reads LDSC pairwise results and trait metadata, validates pair consistency, builds a symmetric rg/q matrix, plots cross-category correlations and performs Ward.D2 clustering.
#

#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Chapter 2 | Genetic-correlation visualisation and hierarchical clustering
# -----------------------------------------------------------------------------
# Purpose:
#   1. Plot quantitative cardiac traits against cardiovascular disease traits.
#   2. Plot the all-trait genetic-correlation matrix in hierarchical-cluster order.
#   3. Generate a Ward.D2 dendrogram using the thesis-reported distance d = 1-rg.
#
# Inputs:
#   results/rg_results_original.tsv  LDSC pairwise results (includes rg, p and q)
#   metadata/traitlist.tsv           trait-to-category annotation
#
# Important provenance note:
#   This is a cleaned, consolidated plotting script based on the plotting code
#   in the reviewed repository's Annotated_scripts/ldsc_corrplot.R. The archived
#   plotting script contains several successive drafts and exploratory code.
#   The Ward.D2 code was present in that plotting script, not in the small
#   standalone R scripts in the Genetic_correlation archive.
#
#   The script uses the q-values already stored in the LDSC results table. It
#   does NOT recalculate FDR. Repeated records for the same unordered trait pair
#   are collapsed only after checking that their rg, p and q values agree.
#   Missing pairwise correlations are not silently imputed for clustering.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(corrplot)
})

args <- commandArgs(trailingOnly = TRUE)
rg_file <- if (length(args) >= 1) args[[1]] else
  "CHAPTER_2/genetic_correlations/results/rg_results_original.tsv"
trait_file <- if (length(args) >= 2) args[[2]] else
  "CHAPTER_2/genetic_correlations/metadata/traitlist.tsv"
out_dir <- if (length(args) >= 3) args[[3]] else
  "CHAPTER_2/genetic_correlations/figures/generated"

if (!file.exists(rg_file)) stop("Cannot find LDSC results file: ", rg_file)
if (!file.exists(trait_file)) stop("Cannot find trait annotation file: ", trait_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rg <- read.delim(rg_file, check.names = FALSE, stringsAsFactors = FALSE)
traits <- read.delim(trait_file, check.names = FALSE, stringsAsFactors = FALSE,
                     fileEncoding = "UTF-8")
required_rg <- c("p1", "p2", "rg", "p", "q")
required_traits <- c("CATEGORY", "TRAIT")
if (!all(required_rg %in% names(rg))) {
  stop("LDSC results must contain: ", paste(required_rg, collapse = ", "))
}
if (!all(required_traits %in% names(traits))) {
  stop("Trait annotation must contain: ", paste(required_traits, collapse = ", "))
}

rg <- rg[complete.cases(rg[, required_rg]), required_rg]
rg$rg <- as.numeric(rg$rg)
rg$p <- as.numeric(rg$p)
rg$q <- as.numeric(rg$q)
if (any(!is.finite(rg$rg)) || any(rg$rg < -1 | rg$rg > 1)) {
  stop("Invalid or out-of-range genetic correlations found in results.")
}

# Canonicalise pair order so reverse-direction rows are recognised as repeats.
pair_lo <- ifelse(rg$p1 <= rg$p2, rg$p1, rg$p2)
pair_hi <- ifelse(rg$p1 <= rg$p2, rg$p2, rg$p1)
rg$pair_id <- paste(pair_lo, pair_hi, sep = "\r")

# Verify that repeated records are equivalent before retaining one record.
by_pair <- split(rg, rg$pair_id)
consistent <- vapply(by_pair, function(x) {
  all(vapply(c("rg", "p", "q"), function(v) {
    length(unique(signif(x[[v]], 10))) == 1
  }, logical(1)))
}, logical(1))
if (!all(consistent)) {
  bad <- names(consistent)[!consistent]
  stop(length(bad), " trait pairs have conflicting rg/p/q records. Inspect the source results before plotting.")
}
rg_unique <- do.call(rbind, lapply(by_pair, function(x) x[1, , drop = FALSE]))
rg_unique$p1 <- pair_lo[match(rg_unique$pair_id, rg$pair_id)]
rg_unique$p2 <- pair_hi[match(rg_unique$pair_id, rg$pair_id)]

# Build a symmetric matrix over the annotated traits. No unobserved correlations
# are set to zero: clustering stops if the requested matrix is incomplete.
trait_names <- unique(traits$TRAIT)
make_matrix <- function(value, trait_names) {
  m <- matrix(NA_real_, nrow = length(trait_names), ncol = length(trait_names),
              dimnames = list(trait_names, trait_names))
  for (i in seq_len(nrow(rg_unique))) {
    a <- rg_unique$p1[[i]]
    b <- rg_unique$p2[[i]]
    if (a %in% trait_names && b %in% trait_names) {
      m[a, b] <- rg_unique[[value]][[i]]
      m[b, a] <- rg_unique[[value]][[i]]
    }
  }
  diag(m) <- if (value == "rg") 1 else 0
  m
}

rg_matrix <- make_matrix("rg", trait_names)
q_matrix <- make_matrix("q", trait_names)

# Figure 2.3-style plot: quantitative cardiac traits versus disease traits.
# Categories are taken from the archived trait annotation, not inferred here.
quant_categories <- c("BP", "Heart rate", "contractility", "MFD",
                      "CONtractility", "LVWT", "LV", "RV",
                      "RV_LV_ratio", "OFT", "Comorbidity",
                      "Fractal", "LV_MRI", "OF_MRI", "RV_MRI",
                      "V_MRI", "Wall thickness", "LV volume",
                      "Fractal slice", "LV wall thickness", "RV volume",
                      "Ventricular", "Contraction", "OF", "RA")
disease_categories <- c("Disease", "disease", "CVD", "Other", "PM")
quant_traits <- traits$TRAIT[traits$CATEGORY %in% quant_categories]
disease_traits <- traits$TRAIT[traits$CATEGORY %in% disease_categories]
quant_traits <- quant_traits[quant_traits %in% rownames(rg_matrix)]
disease_traits <- disease_traits[disease_traits %in% colnames(rg_matrix)]

if (length(quant_traits) > 0 && length(disease_traits) > 0) {
  cross_rg <- rg_matrix[quant_traits, disease_traits, drop = FALSE]
  cross_q <- q_matrix[quant_traits, disease_traits, drop = FALSE]
  pdf(file.path(out_dir, "quantitative_traits_vs_disease.pdf"), width = 16, height = 9)
  corrplot(cross_rg, method = "square", is.corr = FALSE,
           col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200),
           cl.lim = c(-1, 1), p.mat = cross_q, sig.level = 0.05,
           insig = "pch", pch = "*", pch.cex = 0.8,
           tl.cex = 0.65, tl.col = "black", na.label = "")
  dev.off()
} else {
  warning("Could not identify both quantitative and disease trait categories; disease plot skipped.")
}

# All-trait clustering, following the thesis description d = 1 - rg.
# Ward.D2 is applied to the dissimilarity object derived from this matrix.
if (anyNA(rg_matrix)) {
  missing_n <- sum(is.na(rg_matrix[upper.tri(rg_matrix)]))
  stop("Cannot perform complete all-trait clustering: ", missing_n,
       " off-diagonal trait pairs are missing from the supplied results.")
}
if (!isTRUE(all.equal(rg_matrix, t(rg_matrix), tolerance = 1e-10))) {
  stop("Genetic-correlation matrix is not symmetric.")
}
distance_matrix <- 1 - rg_matrix
diag(distance_matrix) <- 0
hc <- hclust(as.dist(distance_matrix), method = "ward.D2")

pdf(file.path(out_dir, "all_trait_wardD2_dendrogram.pdf"), width = 14, height = 10)
plot(hc, main = "Hierarchical clustering of traits (Ward.D2)",
     xlab = "", sub = "Distance = 1 - genetic correlation", cex = 0.55)
dev.off()

ordered_traits <- hc$labels[hc$order]
ordered_rg <- rg_matrix[ordered_traits, ordered_traits, drop = FALSE]
ordered_q <- q_matrix[ordered_traits, ordered_traits, drop = FALSE]
pdf(file.path(out_dir, "all_trait_genetic_correlation_matrix.pdf"), width = 14, height = 12)
corrplot(ordered_rg, method = "square", is.corr = FALSE,
         col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200),
         cl.lim = c(-1, 1), p.mat = ordered_q, sig.level = 0.05,
         insig = "pch", pch = "*", pch.cex = 0.45,
         tl.cex = 0.45, tl.col = "black", tl.srt = 45, na.label = "")
dev.off()

write.table(data.frame(trait = ordered_traits, cluster_order = seq_along(ordered_traits)),
            file.path(out_dir, "wardD2_trait_order.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Plots and Ward.D2 trait order written to: ", normalizePath(out_dir))
