#!/usr/bin/env Rscript
# Visualise gene-set overlap using an UpSet plot and pairwise Jaccard heatmap.
source("R/utils.R")
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  stop("Install UpSetR with install.packages('UpSetR') before running this script.")
}
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "data/gene_lists"
out_dir <- if (length(args) >= 2) args[2] else "results/figures"
tf_files <- c(
  HOPX = "HOPX_bound_genes.txt", GATA4 = "genes_bound_gata4.txt",
  MEF2A = "genes_bound_mef2a.txt", NKX2_5 = "genes_bound_nkx25.txt",
  TBX5 = "genes_bound_tbx5.txt", SRF = "genes_bound_SRF.txt"
)
sets <- lapply(tf_files, function(f) read_gene_list(file.path(input_dir, f)))
names(sets) <- names(tf_files)
membership <- UpSetR::fromList(sets)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(out_dir, "tf_gene_set_upset.pdf"), width = 11, height = 7)
UpSetR::upset(membership, nsets = length(sets), order.by = "freq",
              sets = names(sets), keep.order = TRUE)
dev.off()

jaccard <- expand.grid(tf_1 = names(sets), tf_2 = names(sets),
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(value = length(intersect(sets[[tf_1]], sets[[tf_2]])) /
           length(union(sets[[tf_1]], sets[[tf_2]]))) |>
  ungroup()
p <- ggplot(jaccard, aes(tf_1, tf_2, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  coord_equal() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Jaccard") +
  labs(x = NULL, y = NULL, title = "Overlap among cardiac TF-associated gene sets") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "tf_gene_set_jaccard.pdf"), p, width = 7, height = 6)
