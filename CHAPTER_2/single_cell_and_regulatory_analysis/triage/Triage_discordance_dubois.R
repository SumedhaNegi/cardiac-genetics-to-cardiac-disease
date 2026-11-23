# =============================================================================
# EXAMINER GUIDE
# Script: Triage_discordance_dubois.R
# Purpose: Apply the TRIAGE repressive-score weighting to a mouse cardiac
#          single-cell expression matrix, average by annotated cell type, rank
#          genes, and optionally build gene-set SNP annotations for SBayesRC.
# Source method: Weighted expression = expression * TRIAGE repressive score;
#                average weighted expression by cell type; rank genes within
#                each cell type; select top-N genes; map genes to geneMatrix
#                coordinates expanded by 50 kb; test SNP overlap.
# Input: Seurat RDS; TRIAGE repressive score and gene-ID mapping tables;
#        optional geneMatrix.tsv and SNP-info table for annotation stage.
# Outputs: TRIAGE-weighted matrix, cell-type means, ranked gene lists, top-N
#          lists, and optional binary SNP annotation tables.
# Important: Mouse-to-human conversion code in the historical file was not
#            used downstream. This reconstruction does not silently convert
#            species. The exact cell annotation column and mapping-file schemas
#            must match the source data. No scripts have been run or tested.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# ---- Configuration (override via environment variables) --------------------
INPUT_SEURAT <- Sys.getenv("TRIAGE_SEURAT_RDS", "GSE205950_TP_merge_cardiac.rds")
TRIAGE_DIR <- Sys.getenv("TRIAGE_REFERENCE_DIR", "/QRISdata/Q5149/TRIAGE1")
OUTPUT_DIR <- Sys.getenv("TRIAGE_OUTPUT_DIR", "triage_discordance_output")
CELLTYPE_COLUMN <- Sys.getenv("TRIAGE_CELLTYPE_COLUMN", "cell.type.anno.man")
SPECIES <- Sys.getenv("TRIAGE_SPECIES", "Mouse")
GENE_ID_MODE <- tolower(Sys.getenv("TRIAGE_GENE_ID_MODE", "ensembl")) # ensembl or symbol
LOG_INPUT <- tolower(Sys.getenv("TRIAGE_LOG_INPUT", "false")) %in% c("true", "1", "yes")
TOP_N <- as.integer(strsplit(Sys.getenv("TRIAGE_TOP_N", "25,50,100,250,500"), ",")[[1]])
GENE_MATRIX_FILE <- Sys.getenv("TRIAGE_GENE_MATRIX", "")
SNP_INFO_FILE <- Sys.getenv("TRIAGE_SNP_INFO", "")
FLANK_BP <- as.integer(Sys.getenv("TRIAGE_GENE_FLANK_BP", "50000"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ---- TRIAGE reference tables ------------------------------------------------
# The source used these exact filenames. The original code used global <<-
# assignment and repeated row scans; this version validates columns and uses
# vectorised matching while retaining the same expression-times-score method.
load_triage_reference <- function(reference_dir, species = "Mouse", gene_id_mode = "ensembl") {
  score_file <- file.path(reference_dir, "repressive_hg19_cont_prop.txt")
  mapping_file <- file.path(reference_dir, "ensembl_id_mapping.txt")
  if (!file.exists(score_file)) stop("TRIAGE score file not found: ", score_file)
  if (!file.exists(mapping_file)) stop("TRIAGE gene mapping file not found: ", mapping_file)
  score <- data.table::fread(score_file, header = FALSE, data.table = FALSE)
  mapping <- data.table::fread(mapping_file, header = TRUE, data.table = FALSE)
  if (ncol(score) < 2L) stop("TRIAGE score table needs gene and score columns.")
  if (ncol(score) < 2L) stop("Could not identify gene/score columns in TRIAGE score file.")
  score_genes <- as.character(score[[1]])
  score_values <- suppressWarnings(as.numeric(score[[2]]))
  if (anyDuplicated(score_genes)) warning("Duplicate genes in TRIAGE score table; first occurrence is used.")
  keep <- !duplicated(score_genes) & !is.na(score_genes) & is.finite(score_values)
  score <- data.frame(gene = score_genes[keep], score = score_values[keep])

  # The original constructed these mapping headers from species and gene ID mode.
  species_col <- if (identical(species, "Human")) "Gene" else paste(species, "gene")
  id_col <- if (identical(gene_id_mode, "ensembl")) "name" else "stable ID"
  mapping_col <- paste(species_col, id_col)
  if (!mapping_col %in% names(mapping)) {
    stop("Mapping table does not contain expected column '", mapping_col,
         "'. Available columns: ", paste(names(mapping), collapse = ", "))
  }
  # The historical function selected column 2 as the mapped HUGO/ortholog name.
  mapped_gene_col <- names(mapping)[2]
  list(score = score, mapping = mapping, input_id_col = mapping_col,
       mapped_gene_col = mapped_gene_col)
}

apply_triage_weights <- function(expression, reference, log_input = FALSE) {
  expression <- as(expression, "dgCMatrix")
  input_ids <- rownames(expression)
  if (is.null(input_ids)) stop("Expression matrix must have gene identifiers as row names.")
  mapping <- reference$mapping
  idx <- match(input_ids, as.character(mapping[[reference$input_id_col]]))
  mapped <- rep(NA_character_, length(idx))
  ok <- !is.na(idx)
  mapped[ok] <- as.character(mapping[[reference$mapped_gene_col]][idx[ok]])
  score_idx <- match(mapped, reference$score$gene)
  weights <- numeric(length(input_ids))
  has_score <- !is.na(score_idx)
  weights[has_score] <- reference$score$score[score_idx[has_score]]
  # Original code assigned zero weight to unmappable genes.
  if (log_input) expression <- log1p(expression)
  weighted <- Matrix::Diagonal(x = weights) %*% expression
  rownames(weighted) <- input_ids
  attr(weighted, "triage_mapping_summary") <- data.frame(
    input_genes = length(input_ids), mapped_to_reference = sum(ok),
    assigned_nonzero_weight = sum(weights != 0), stringsAsFactors = FALSE)
  weighted
}

# ---- Seurat expression extraction and cell-type averaging ------------------
if (!file.exists(INPUT_SEURAT)) stop("Seurat RDS not found: ", INPUT_SEURAT)
suppressPackageStartupMessages(library(Seurat))
seurat_obj <- readRDS(INPUT_SEURAT)
if (!inherits(seurat_obj, "Seurat")) stop("Input RDS is not a Seurat object.")
if (!"RNA" %in% names(seurat_obj@assays)) stop("Seurat object has no RNA assay.")
meta <- seurat_obj@meta.data
if (!CELLTYPE_COLUMN %in% names(meta)) {
  stop("Cell-type metadata column '", CELLTYPE_COLUMN, "' not found. Available: ",
       paste(names(meta), collapse = ", "))
}

# Match the historical use of RNA normalized data (not raw counts). Seurat v4/v5
# accessors differ; GetAssayData handles the common versions with layer/slot.
expr <- tryCatch(
  Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
  error = function(e) Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data")
)
if (!identical(colnames(expr), rownames(meta))) {
  common <- intersect(colnames(expr), rownames(meta))
  if (!length(common)) stop("No matching cell identifiers between RNA assay and metadata.")
  expr <- expr[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
}
cell_type <- as.character(meta[[CELLTYPE_COLUMN]])
valid_cells <- !is.na(cell_type) & nzchar(cell_type)
expr <- expr[, valid_cells, drop = FALSE]
cell_type <- cell_type[valid_cells]
if (!ncol(expr)) stop("No cells remain after removing missing cell-type labels.")

reference <- load_triage_reference(TRIAGE_DIR, SPECIES, GENE_ID_MODE)
triage_expr <- apply_triage_weights(expr, reference, LOG_INPUT)
write.csv(attr(triage_expr, "triage_mapping_summary"),
          file.path(OUTPUT_DIR, "triage_mapping_summary.csv"), row.names = FALSE)

# Sparse indicator matrix: columns correspond to cell types; each cell belongs
# to exactly one group. Multiplication calculates the mean expression per type.
cell_levels <- unique(cell_type)
grouping <- Matrix::sparse.model.matrix(~ 0 + factor(cell_type, levels = cell_levels))
colnames(grouping) <- cell_levels
means <- triage_expr %*% grouping
cell_counts <- as.numeric(Matrix::colSums(grouping))
means <- means %*% Matrix::Diagonal(x = 1 / cell_counts)
colnames(means) <- cell_levels
saveRDS(triage_expr, file.path(OUTPUT_DIR, "triage_weighted_expression.rds"))
saveRDS(means, file.path(OUTPUT_DIR, "mean_triage_expression_by_cell_type.rds"))

# Rank genes independently within each cell type by decreasing mean weighted
# expression, preserving the gene IDs used in the original expression matrix.
ranked <- lapply(seq_len(ncol(means)), function(j) {
  vals <- as.numeric(means[, j])
  ord <- order(vals, decreasing = TRUE, na.last = NA)
  data.frame(cell_type = colnames(means)[j], gene = rownames(means)[ord],
             triage_discordance_score = vals[ord], stringsAsFactors = FALSE)
})
names(ranked) <- colnames(means)
rank_dir <- file.path(OUTPUT_DIR, "ordered_genes_files")
dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
for (nm in names(ranked)) {
  data.table::fwrite(ranked[[nm]], file.path(rank_dir,
    paste0("ordered_genes_", safe_name(nm), ".txt")), sep = "\t")
}

# Generate top-N subsets. The historical script used 25, 50, 250 and 500;
# top 100 is included because it is discussed elsewhere in the source workflow.
TOP_N <- unique(TOP_N[is.finite(TOP_N) & TOP_N > 0])
for (n in TOP_N) {
  dest <- file.path(OUTPUT_DIR, paste0("top_", n))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(ranked)) {
    top <- utils::head(ranked[[nm]], n)
    data.table::fwrite(top, file.path(dest,
      paste0("top_", n, "_", safe_name(nm), ".txt")), sep = "\t")
  }
}

# ---- Optional gene-set-to-SNP overlap annotations --------------------------
# This reproduces the source's geneMatrix coordinates +/- 50 kb and overlaps
# against SNP positions. It is optional because the supplied repository does
# not contain the original geneMatrix or snp.info files.
if (nzchar(GENE_MATRIX_FILE) && nzchar(SNP_INFO_FILE)) {
  if (!file.exists(GENE_MATRIX_FILE)) stop("Gene matrix not found: ", GENE_MATRIX_FILE)
  if (!file.exists(SNP_INFO_FILE)) stop("SNP info file not found: ", SNP_INFO_FILE)
  suppressPackageStartupMessages(library(GenomicRanges))
  gene_matrix <- data.table::fread(GENE_MATRIX_FILE, data.table = FALSE)
  snp_raw <- data.table::fread(SNP_INFO_FILE, data.table = FALSE)
  if (ncol(gene_matrix) < 15L) stop("geneMatrix.tsv has fewer than 15 columns; source used columns 2,12:15.")
  genes <- gene_matrix[, c(2, 12, 13, 14, 15)]
  names(genes) <- c("gene", "chr", "start", "end", "strand")
  if (ncol(snp_raw) < 7L) stop("SNP info file has fewer than 7 columns; source selected columns 1,2,5:7.")
  snps <- snp_raw[, c(1, 2, 5, 6, 7)]
  names(snps) <- c("chr", "rsid", "pos", "A1", "A2")
  snps$chr <- sub("^chr", "", as.character(snps$chr), ignore.case = TRUE)
  snps$pos <- as.numeric(snps$pos)
  snp_ranges <- GenomicRanges::GRanges(seqnames = snps$chr,
    ranges = IRanges::IRanges(start = snps$pos, end = snps$pos))

  for (n in TOP_N) {
    top_dir <- file.path(OUTPUT_DIR, paste0("top_", n))
    top_files <- list.files(top_dir, pattern = "^top_.*\\.txt$", full.names = TRUE)
    if (!length(top_files)) next
    annotation <- data.frame(SNP = snps$rsid, stringsAsFactors = FALSE)
    for (f in top_files) {
      gl <- data.table::fread(f, data.table = FALSE)
      if (!"gene" %in% names(gl)) next
      selected <- unique(toupper(as.character(gl$gene)))
      g <- genes[toupper(as.character(genes$gene)) %in% selected, , drop = FALSE]
      g$chr <- sub("^chr", "", as.character(g$chr), ignore.case = TRUE)
      g$start <- as.numeric(g$start) - FLANK_BP
      g$end <- as.numeric(g$end) + FLANK_BP
      g <- g[is.finite(g$start) & is.finite(g$end) & g$start >= 1 & g$end >= g$start, , drop = FALSE]
      if (!nrow(g)) flag <- rep(0L, nrow(snps)) else {
        gr <- GenomicRanges::GRanges(seqnames = g$chr,
          ranges = IRanges::IRanges(start = g$start, end = g$end))
        flag <- as.integer(GenomicRanges::overlapsAny(snp_ranges, gr))
      }
      nm <- tools::file_path_sans_ext(basename(f))
      annotation[[paste0("overlap_", nm)]] <- flag
    }
    out <- file.path(OUTPUT_DIR, paste0("top_", n, "_annotation_binary.txt"))
    data.table::fwrite(annotation, out, sep = "\t", quote = FALSE)
  }
} else {
  message("Skipping SNP annotation: set both TRIAGE_GENE_MATRIX and TRIAGE_SNP_INFO to enable it.")
}

message("TRIAGE workflow outputs written to: ", normalizePath(OUTPUT_DIR))

# Original cell-cluster labels (as recorded in the source notes) are not used
# to relabel the dataset automatically; metadata already carries the chosen
# cell annotation. The source's mouse-to-human conversion attempts were not
# used in its subsequent ranking/annotation code and are intentionally omitted.
