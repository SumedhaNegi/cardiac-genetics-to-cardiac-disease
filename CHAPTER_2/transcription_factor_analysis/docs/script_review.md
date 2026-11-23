# Script selection and review notes

The archive contained multiple R and shell-script variants, including scripts
with embedded interactive console commands, hard-coded QRIS paths, repeated
code blocks, inconsistent object names, and partial exploratory analyses.

This repository replaces those fragments with modular scripts:

- `01_summarise_gene_sets.R`: gene-set sizes, membership table and Jaccard index.
- `02_build_snp_annotations.R`: SNP-to-interval overlaps using GenomicRanges.
- `03_summarise_sbayesrc.R`: collates existing SBayesRC enrichment outputs.
- `04_plot_heritability_enrichment.R`: publication-oriented enrichment plots.

## Known issues deliberately not hidden

- Some source coordinate rows store `start > end`; the new overlap routine
  normalises endpoints but does not resolve assembly or coordinate-system
  uncertainty.
- The historical archive includes both 50-kb-window and no-window coordinates.
  These are kept as distinct input files and must not be mixed.
- The historical code contains a faulty Euler/intersection construction and
  row-wise SNP interval loops. The cleaned implementation uses set operations
  and `GenomicRanges::findOverlaps`.
- Some scripts refer to annotations or input paths that are absent from the
  archive. They are not presented as executable.
- SBayesRC annotation columns with very small/zero variance caused failures in
  historical logs. Check annotation prevalence and remove or combine
  non-informative columns only under a documented, pre-specified rule.
- The source archive contains mixed-species gene identifiers and datasets.
  Species/orthology provenance must be confirmed before combining sets.
