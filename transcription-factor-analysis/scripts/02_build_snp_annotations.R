#!/usr/bin/env Rscript
# Intersect SNP positions with TF-associated gene intervals using GenomicRanges.
# Required input: SNP reference table with chromosome, position and rsid columns.
# The source archive's SNP membership matrices do not contain genomic positions;
# provide the matching SNP reference table separately.
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(readr)
})
source("R/utils.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/02_build_snp_annotations.R <snp_info.tsv> <genome_build_label> [output_dir]")
}
snp_file <- args[1]
genome_build <- args[2]
out_dir <- if (length(args) >= 3) args[3] else "results/annotations"
coordinate_mode <- if (length(args) >= 4) args[4] else "window"
if (!coordinate_mode %in% c("window", "no_window")) {
  stop("coordinate_mode must be either 'window' or 'no_window'.")
}
if (!file.exists(snp_file)) stop("SNP reference file not found: ", snp_file)
if (!nzchar(genome_build)) stop("A genome-build label is required.")

snps <- fread(snp_file, data.table = FALSE, check.names = FALSE)
# Map columns explicitly via config/analysis_config.yml; defaults below match the
# UKB snp.info layout recorded in the historical analysis notes.
cfg <- yaml::read_yaml("config/analysis_config.yml")
sc <- cfg$snp_columns
assert_columns(snps, unlist(sc), "SNP reference table")
snps <- snps |>
  transmute(
    rsid = as.character(.data[[sc$rsid]]),
    chr = str_remove(as.character(.data[[sc$chromosome]]), "^chr"),
    position = as.numeric(.data[[sc$position]])
  ) |>
  filter(!is.na(rsid), !is.na(chr), !is.na(position)) |>
  distinct(rsid, .keep_all = TRUE)

snp_gr <- GRanges(seqnames = snps$chr, ranges = IRanges(snps$position, snps$position))
coord_files <- c(
  HOPX = if (coordinate_mode == "window") "HOPX_bound_genes_location.txt" else "HOPX_bound_genesnowindow_location.txt",
  GATA4 = if (coordinate_mode == "window") "genes_bound_gata4_location.txt" else "genes_bound_gata4nowindow_location.txt",
  MEF2A = if (coordinate_mode == "window") "genes_bound_mef2a_location.txt" else "genes_bound_mef2anowindow_location.txt",
  NKX2_5 = if (coordinate_mode == "window") "genes_bound_nkx25_location.txt" else "genes_bound_nkx25nowindow_location.txt",
  TBX5 = if (coordinate_mode == "window") "genes_bound_tbx5_location.txt" else "genes_bound_tbx5nowindow_location.txt",
  SRF = if (coordinate_mode == "window") "genes_bound_SRF_location.txt" else "genes_bound_SRFnowindow_location.txt"
)
annotation_list <- list()
for (tf in names(coord_files)) {
  path <- file.path("data/gene_coordinates", coord_files[[tf]])
  coords <- read_gene_coordinates(path)
  gr <- GRanges(seqnames = coords$chr,
                ranges = IRanges(coords$interval_start, coords$interval_end))
  hits <- findOverlaps(snp_gr, gr, ignore.strand = TRUE)
  membership <- rep(0L, length(snps))
  membership[unique(queryHits(hits))] <- 1L
  annotation_list[[tf]] <- membership
  message(tf, ": ", sum(membership), " SNPs overlap at least one interval")
}
annotation <- tibble(rsid = snps$rsid, !!!setNames(annotation_list, names(annotation_list)))
attr(annotation, "genome_build") <- genome_build
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(annotation, file.path(out_dir, paste0("tf_snp_annotation_", coordinate_mode, "_", genome_build, ".tsv")))
write_tsv(tibble(
  annotation = names(annotation_list),
  n_snps = vapply(annotation_list, sum, integer(1)),
  n_total_snps = nrow(snps),
  genome_build = genome_build,
  coordinate_mode = coordinate_mode
), file.path(out_dir, paste0("tf_snp_annotation_qc_", coordinate_mode, "_", genome_build, ".tsv")))
