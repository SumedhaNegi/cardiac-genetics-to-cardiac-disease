#!/usr/bin/env Rscript
# Publication-oriented plots for per-SNP heritability enrichment estimates.
source("R/utils.R")
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "results/tables/sbayesrc_enrichment_all_traits.tsv"
out_dir <- if (length(args) >= 2) args[2] else "results/figures"
if (!file.exists(input_file)) stop("Run scripts/03_summarise_sbayesrc.R first.")
x <- fread(input_file, data.table = FALSE)
assert_columns(x, c("Annotation", "Enrich", "trait_file"), "enrichment table")
x <- x |>
  mutate(
    Enrich = as.numeric(Enrich),
    Annotation = factor(Annotation, levels = c("HOPX", "GATA4", "MEF2A", "NKX2_5", "TBX5", "SRF"))
  ) |>
  filter(is.finite(Enrich), !is.na(Annotation))
if (!"SD" %in% names(x)) x$SD <- NA_real_
x$SD <- as.numeric(x$SD)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- ggplot(x, aes(x = trait_file, y = Enrich, colour = Annotation, group = Annotation)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, colour = "grey45") +
  geom_point(position = position_dodge(width = 0.55), size = 1.8, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = pmax(0, Enrich - SD), ymax = Enrich + SD),
    position = position_dodge(width = 0.55), width = 0.18, na.rm = TRUE
  ) +
  labs(
    title = "Cardiac transcription-factor annotations and SNP heritability",
    subtitle = "Dashed line denotes enrichment = 1 (no enrichment relative to SNP proportion)",
    x = "GWAS summary-statistic result", y = "Per-SNP heritability enrichment",
    colour = "Annotation"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
        legend.position = "bottom")
ggsave(file.path(out_dir, "tf_heritability_enrichment.pdf"), p,
       width = 12, height = 7, limitsize = FALSE)
ggsave(file.path(out_dir, "tf_heritability_enrichment.png"), p,
       width = 12, height = 7, dpi = 320, limitsize = FALSE)
