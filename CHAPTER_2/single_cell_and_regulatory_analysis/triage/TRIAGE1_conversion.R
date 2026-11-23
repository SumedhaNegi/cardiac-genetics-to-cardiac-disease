# =============================================================================
# TRIAGE1 conversion helper
# =============================================================================

# CODE MAP
# Defines the TRIAGE1 helper: maps matrix row identifiers to a repressive-chromatin score and applies that score row-wise to expression values.
#
# Applies the supplied repressive-chromatin score to each expression-matrix
# row after mapping its gene identifier to the corresponding human/mouse gene.
# The historical method uses log(m + 1) * score when log = TRUE and m * score
# otherwise. The score is applied row-wise to every sample/column.
#
# Required reference files in `dir`:
#   repressive_hg19_cont_prop.txt  (gene identifier in V1, score in V2)
#   ensembl_id_mapping.txt         (species-specific mapping columns)
# The historical reference is labelled hg19; verify its assembly/species match
# the input expression matrix before use.

TRIAGE1 <- function(m, m_names = rownames(m), species = "Human",
                    geneid = TRUE, log = FALSE, dir,
                    repressive_file = "repressive_hg19_cont_prop.txt",
                    mapping_file = "ensembl_id_mapping.txt") {
  if (is.null(dim(m)) || length(dim(m)) != 2L) {
    stop("`m` must be a numeric matrix or two-dimensional data frame.")
  }
  m <- as.matrix(m)
  if (!is.numeric(m)) stop("`m` must contain numeric expression values.")
  if (length(m_names) != nrow(m)) {
    stop("`m_names` must contain one gene identifier per matrix row.")
  }
  if (missing(dir) || length(dir) != 1L || !dir.exists(dir)) {
    stop("Supply `dir` as an existing directory containing the two reference files.")
  }
  repressive_path <- file.path(dir, repressive_file)
  mapping_path <- file.path(dir, mapping_file)
  if (!file.exists(repressive_path)) stop("Missing TRIAGE score file: ", repressive_path)
  if (!file.exists(mapping_path)) stop("Missing gene mapping file: ", mapping_path)

  # Load references locally rather than assigning objects into .GlobalEnv.
  w <- as.data.frame(data.table::fread(repressive_path))
  g <- as.data.frame(data.table::fread(mapping_path))
  if (ncol(w) < 2L) stop("TRIAGE score file must contain at least two columns.")
  if (ncol(g) < 2L) stop("Gene mapping file must contain at least two columns.")
  if (!"V1" %in% names(w)) names(w)[1] <- "V1"

  # Preserve the source convention: Human uses 'Gene'; other species use
  # '<species> gene'. `geneid=TRUE` selects gene-symbol mapping; FALSE selects
  # stable Ensembl identifiers.
  species_column <- if (identical(species, "Human")) "Gene" else paste(species, "gene")
  identifier_column <- if (isTRUE(geneid)) "name" else "stable ID"
  mapping_column <- paste(species_column, identifier_column)
  if (!mapping_column %in% names(g)) {
    stop("Mapping file is missing expected column: ", mapping_column,
         ". Available columns: ", paste(names(g), collapse = ", "))
  }

  # Map each input identifier to the reference gene key, then retrieve its score.
  map_index <- match(as.character(m_names), as.character(g[[mapping_column]]))
  mapped_gene <- as.character(g[[2]][map_index])
  score_index <- match(mapped_gene, as.character(w[[1]]))
  scores <- suppressWarnings(as.numeric(w[[2]][score_index]))
  scores[is.na(scores)] <- 0

  # The original script multiplied each gene row by its TRIAGE score. `sweep`
  # makes that row-wise intent explicit and avoids R's column-major recycling.
  transformed <- if (isTRUE(log)) base::log(m + 1) else m
  out <- sweep(transformed, MARGIN = 1L, STATS = scores, FUN = "*")
  output_names <- rownames(m)
  if (is.null(output_names)) output_names <- as.character(m_names)
  rownames(out) <- output_names
  colnames(out) <- colnames(m)
  out
}

# Example (uncomment and replace paths/objects for a real analysis):
# expression_matrix <- as.matrix(readRDS("expression_matrix.rds"))
# triage_matrix <- TRIAGE1(expression_matrix, species = "Human",
#                          geneid = TRUE, log = FALSE, dir = "reference_data")
