# Lightweight validation of the supplied gene-list inputs.
source("R/utils.R")
files <- c(
  "HOPX_bound_genes.txt", "genes_bound_gata4.txt", "genes_bound_mef2a.txt",
  "genes_bound_nkx25.txt", "genes_bound_tbx5.txt", "genes_bound_SRF.txt"
)
for (f in files) {
  x <- read_gene_list(file.path("data/gene_lists", f))
  stopifnot(length(x) > 0L, !anyNA(x), all(nzchar(x)))
}
message("Gene-list input checks passed.")
