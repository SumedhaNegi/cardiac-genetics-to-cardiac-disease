#!/usr/bin/env Rscript
# =============================================================================
# Cardiac enhancer interval annotation with SNP positions
# =============================================================================

# CODE MAP
# Expands cardiac enhancer intervals by 50 kb and annotates SNPs that overlap the expanded intervals.
#
# The historical analysis expanded each cardiac-enhancer BED interval by
# 50 kb on either side and labelled UK Biobank SNPs that overlapped any window.
# Inputs are supplied as command-line arguments to avoid machine-specific paths:
#   Rscript cardiac_enhancer_annpt.R enhancers_hg19.bed snp.info output.tsv
# BED coordinates are assumed to be 0-based, half-open. GRanges uses 1-based,
# closed coordinates, so the BED start is incremented by one on import.
# Confirm the genome assembly of BOTH files before interpreting overlaps.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript cardiac_enhancer_annpt.R <enhancers.bed> <snp.info> <output.tsv>")
}
enhancer_file <- args[[1]]
snp_file <- args[[2]]
output_file <- args[[3]]
if (!file.exists(enhancer_file)) stop("Enhancer BED file not found: ", enhancer_file)
if (!file.exists(snp_file)) stop("SNP reference file not found: ", snp_file)

# Read BED without a header. The first three columns are chromosome, start, end.
enhancers <- fread(enhancer_file, header = FALSE, select = 1:3)
if (ncol(enhancers) != 3L || nrow(enhancers) == 0L) {
  stop("Enhancer BED must contain at least three columns and one interval.")
}
setnames(enhancers, c("chr", "start_bed", "end_bed"))
enhancers[, `:=`(
  chr = sub("^chr", "", as.character(chr)),
  start = as.integer(start_bed) - 50000L + 1L,
  end = as.integer(end_bed) + 50000L
)]
enhancers <- enhancers[!is.na(chr) & !is.na(start) & !is.na(end) & end >= start]
if (nrow(enhancers) == 0L) stop("No valid enhancer intervals remain after parsing.")

# Historical SNP.info layout: columns 1, 2, 5, 6 and 7 are chromosome,
# rsID, physical position, allele 1 and allele 2, respectively.
snps <- fread(snp_file, header = FALSE)
if (ncol(snps) < 7L) stop("SNP reference file has fewer than seven columns.")
snp2 <- snps[, .(Chromosome = as.character(V1), rsid = as.character(V2),
                  PhysPos = as.integer(V5), A1 = as.character(V6),
                  A2 = as.character(V7))]
snp2[, Chromosome := sub("^chr", "", Chromosome)]
snp2 <- snp2[!is.na(Chromosome) & !is.na(PhysPos)]

# Convert both coordinate sets to GRanges and mark each SNP once if it overlaps
# at least one expanded enhancer interval.
snp_ranges <- GRanges(seqnames = snp2$Chromosome,
                       ranges = IRanges(start = snp2$PhysPos, width = 1L))
enhancer_ranges <- GRanges(seqnames = enhancers$chr,
                           ranges = IRanges(start = enhancers$start,
                                            end = enhancers$end))
snp2[, enhancer := as.integer(overlapsAny(snp_ranges, enhancer_ranges,
                                          ignore.strand = TRUE))]
fwrite(snp2, output_file, sep = "\t", na = "NA")
message(sum(snp2$enhancer), " of ", nrow(snp2),
        " SNP rows overlap a cardiac-enhancer interval expanded by 50 kb.")
message("Output written to: ", normalizePath(output_file, mustWork = FALSE))
