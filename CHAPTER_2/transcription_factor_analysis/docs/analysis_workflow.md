# Analysis workflow

## Scientific question

Do genomic regions associated with cardiac transcriptional regulators contain a
disproportionate share of common-variant heritability for cardiac traits, and does
this differ between developmental stages?

The thesis compares HOPX-associated regions with GATA4, MEF2A, NKX2-5, TBX5
and SRF. It also considers adult-specific, foetal-specific and shared binding
annotations.

## Workflow modules

1. **Gene-set curation and QC** — read unique gene identifiers, report gene-set
   sizes, and calculate pairwise Jaccard similarity.
2. **Coordinate harmonisation** — validate chromosome and interval columns,
   normalise reversed endpoints, and check assembly/species compatibility.
3. **SNP annotation** — overlap reference SNP positions with genomic intervals
   using `GenomicRanges`; export a binary annotation matrix and QC summary.
4. **Heritability modelling** — run SBayesRC with compatible GWAS summary
   statistics, LD reference data and annotation matrix. This stage requires
   external inputs and an HPC-scale environment.
5. **Result collation** — combine per-trait `AnnoPerSnpHsqEnrichment` outputs.
6. **Visualisation** — plot enrichment estimates with the null reference at 1.

## Historical approach recorded in the archive

The supplied scripts describe both 1-million-SNP and approximately 7.3-million-
SNP workflows, with and without a 50-kb gene window. They include SBayesRC
`tidy`, `impute` and `sbayesrc` stages, plus a later no-window analysis. Historical
shell scripts refer to UQ HPC paths and modules and are not portable as written.

The original notes also contain attempted custom annotation-building code and
several alternative script versions. The cleaned workflow replaces the most
fragile parts with explicit input validation and interval-based overlap logic.

## Reproducibility boundary

This repository is a cleaned analysis scaffold and curated record. It is not a
one-command reproduction of the original HPC analysis because the archive does
not include all compatible GWAS summary statistics, LD reference files, complete
configuration tables, and all source annotation resources. Exact genome build,
species, coordinate convention, and SBayesRC version must be confirmed before
running the SNP-level pipeline.
