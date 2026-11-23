# =============================================================================
# EXAMINER GUIDE
# Script: DamID_enhancer_heritability.R
# Purpose: Create binary SNP annotations for cardiac enhancer and HOPX-enriched
#          enhancer intervals for downstream SBayesRC heritability enrichment.
# Source method retained: BED intervals are tested against UK Biobank/HM3 BIM
#          SNP positions; an annotation and its complement are written.
# Inputs: BED files (chromosome, start, end, optional extra columns) and a
#         standard PLINK BIM file (CHR, SNP, CM, BP, A1, A2).
# Outputs: Tab-delimited annotation matrices keyed by SNP ID.
# Important: Genome build and coordinate convention must match between BED,
#            BIM and LD reference. BED is 0-based half-open; BIM BP is 1-based.
#            This script converts BED starts to 1-based before overlap testing.
#            Optional WT/differential-accessibility files are only included if
#            their paths are supplied below. No biological result is inferred.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---- User configuration -----------------------------------------------------
# Set paths here or provide environment variables with the same names.
BIM_FILE <- Sys.getenv("BIM_FILE", "ukbEURu_v3_hm3_maf01.bim")
CARDIAC_BED <- Sys.getenv("CARDIAC_ENHANCER_BED", "cardiac_enhancer_hg38.bed")
HOPX_BED <- Sys.getenv("HOPX_ENHANCER_BED", "hopx_wt_ratio1_average_peaks2_enhancers.bed")
# These source objects were referenced but never defined in the supplied file.
# Supply BED paths if the corresponding WT/combined annotation is required.
WT_BED <- Sys.getenv("WT_ENHANCER_BED", "")
DIFF_BED <- Sys.getenv("DIFF_ACCESSIBILITY_BED", "")
OUTPUT_DIR <- Sys.getenv("ANNOTATION_OUTPUT_DIR", ".")

# ---- Input readers and validation ------------------------------------------
read_bed_intervals <- function(path, label) {
  if (!nzchar(path) || !file.exists(path)) stop(label, " BED file not found: ", path)
  x <- data.table::fread(path, header = FALSE, sep = "\t", data.table = FALSE)
  if (ncol(x) < 3L) stop(label, " BED file must have at least 3 columns.")
  x <- x[, 1:3, drop = FALSE]
  names(x) <- c("chr", "start0", "end0")
  x$chr <- sub("^chr", "", as.character(x$chr), ignore.case = TRUE)
  x$chr <- toupper(x$chr)
  x$start0 <- suppressWarnings(as.numeric(x$start0))
  x$end0 <- suppressWarnings(as.numeric(x$end0))
  x <- x[!is.na(x$chr) & is.finite(x$start0) & is.finite(x$end0) &
           x$start0 >= 0 & x$end0 > x$start0, , drop = FALSE]
  if (!nrow(x)) stop(label, " BED file contains no valid intervals.")
  # BED is 0-based, half-open. Convert to 1-based closed coordinates.
  x$start <- x$start0 + 1
  x$end <- x$end0
  x[, c("chr", "start", "end"), drop = FALSE]
}

read_plink_bim <- function(path) {
  if (!file.exists(path)) stop("PLINK BIM file not found: ", path)
  x <- data.table::fread(path, header = FALSE, sep = "\t", data.table = FALSE)
  if (ncol(x) < 6L) stop("BIM file must have six standard PLINK columns.")
  x <- x[, 1:6, drop = FALSE]
  names(x) <- c("chr", "SNP", "CM", "BP", "A1", "A2")
  x$chr <- toupper(sub("^chr", "", as.character(x$chr), ignore.case = TRUE))
  x$BP <- suppressWarnings(as.numeric(x$BP))
  if (anyNA(x$BP)) stop("BIM BP column contains non-numeric positions.")
  if (anyDuplicated(x$SNP)) warning("BIM contains duplicated SNP IDs; rows are retained.")
  x
}

# Efficient interval overlap using a chromosome-wise sorted interval scan.
# This avoids the original sapply loop over >1 million SNPs x all intervals.
flag_overlaps <- function(snps, intervals) {
  # GenomicRanges uses interval indexing, avoiding the original nested scan
  # over every SNP and every enhancer interval.
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop("Install Bioconductor packages GenomicRanges and IRanges for interval overlap.")
  }
  snp_gr <- GenomicRanges::GRanges(
    seqnames = snps$chr,
    ranges = IRanges::IRanges(start = snps$BP, end = snps$BP)
  )
  interval_gr <- GenomicRanges::GRanges(
    seqnames = intervals$chr,
    ranges = IRanges::IRanges(start = intervals$start, end = intervals$end)
  )
  as.logical(GenomicRanges::overlapsAny(snp_gr, interval_gr))
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
snps <- read_plink_bim(BIM_FILE)
intervals <- list(
  cardiac_enhancer = read_bed_intervals(CARDIAC_BED, "Cardiac enhancer"),
  hopx_enriched = read_bed_intervals(HOPX_BED, "HOPX-enriched enhancer")
)
if (nzchar(WT_BED)) intervals$wt_enriched <- read_bed_intervals(WT_BED, "WT enhancer")
if (nzchar(DIFF_BED)) intervals$differential_accessibility <- read_bed_intervals(DIFF_BED, "Differential accessibility")

# One SNP row per BIM record; columns are binary annotation indicators.
annotation <- data.frame(SNP = snps$SNP, stringsAsFactors = FALSE)
for (nm in names(intervals)) {
  annotation[[nm]] <- as.integer(flag_overlaps(snps, intervals[[nm]]))
  annotation[[paste0(nm, "_opposite")]] <- 1L - annotation[[nm]]
  message(nm, ": ", sum(annotation[[nm]]), " overlapping SNPs of ", nrow(snps))
}

# Save one combined file and one file per annotation, retaining the SNP key.
combined_file <- file.path(OUTPUT_DIR, "hm3_enhancer_annotations.txt")
data.table::fwrite(annotation, combined_file, sep = "\t", quote = FALSE)
for (nm in names(intervals)) {
  cols <- c("SNP", nm, paste0(nm, "_opposite"))
  out <- file.path(OUTPUT_DIR, paste0("hm3_", nm, "_snps.txt"))
  data.table::fwrite(annotation[, cols, drop = FALSE], out, sep = "\t", quote = FALSE)
}

# The historical file also wrote an awk/sort/uniq command here. It is not
# required for annotation generation and is intentionally not run in R.
message("Annotation files written to: ", normalizePath(OUTPUT_DIR))

# SBayesRC execution is deliberately separate from this R annotation script.
# Use the matching SLURM template and confirm that BIM, BED and LD reference
# all use the same genome build before submitting a model.
