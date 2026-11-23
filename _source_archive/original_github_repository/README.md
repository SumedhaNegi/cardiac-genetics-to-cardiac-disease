# Cardiac Genetics to Cardiac Disease

Examiner-facing analysis record accompanying the thesis *Cardiac Genetics to Cardiac Disease: Statistical Genomics Analyses and Functional Characterisation of HOPX and its variants*.

## Scope of this repository version

This initial GitHub-ready build documents the Chapter 2 **genetic-correlation analysis**. It is not yet the complete thesis code repository. The folder records the analysis sequence, the scripts identified in the supplied archive, selected cleaned scripts, trait annotations, archived results and original figure outputs.

The purpose is to make the work understandable to examiners. It is not intended to recreate the original UQ computing environment or to redistribute GWAS summary statistics, institutional data, or third-party reference panels.

## Where to start

1. Read [`CHAPTER_2/genetic_correlations/docs/analysis_workflow.md`](CHAPTER_2/genetic_correlations/docs/analysis_workflow.md) for the sequence of analysis steps.
2. Read [`CHAPTER_2/genetic_correlations/docs/script_selection_and_review.md`](CHAPTER_2/genetic_correlations/docs/script_selection_and_review.md) for the selection of script versions and known caveats.
3. See [`CHAPTER_2/genetic_correlations/docs/methods_examiner_summary.md`](CHAPTER_2/genetic_correlations/docs/methods_examiner_summary.md) for an examiner-facing description of the method.
4. The cleaned plotting/clustering script is [`CHAPTER_2/genetic_correlations/scripts/plot_genetic_correlations.R`](CHAPTER_2/genetic_correlations/scripts/plot_genetic_correlations.R).

## Repository layout

```text
CHAPTER_2/genetic_correlations/
├── docs/       workflow, methods and script audit
├── figures/    archived thesis-analysis figure outputs
├── metadata/   trait annotations and source-file mapping
├── results/    archived LDSC result table and significant-count output
└── scripts/    cleaned, annotated SLURM and R scripts
```

## Important interpretation and provenance notes

- Bivariate LDSC estimates genome-wide genetic correlation; it does not identify causal variants or establish causal relationships.
- The thesis reports 1,485 pairwise correlations among quantitative cardiac traits and describes Benjamini–Hochberg FDR correction across 10,877 tested correlations in the quantitative-trait/disease analysis captions.
- The supplied `rg_results_original.tsv` contains repeated records and some repeated trait pairs with conflicting estimates. The plotting script deliberately stops rather than silently selecting or averaging conflicting results. These records need to be reconciled against the final LDSC output/logs before using this archived table to regenerate final figures.
- Ward.D2 clustering is included in the plotting workflow, using the thesis-reported dissimilarity `d = 1 - rg`. The original Ward.D2 implementation was found in the long `ldsc_corrplot.R` plotting script in the reviewed repository, not in the small standalone plotting scripts in the genetic-correlation ZIP.
- The archived reference files (`w_hm3.snplist`, European LD-score reference panels) are not redistributed here. They are third-party/reference resources and should be obtained from their original source or institutional environment.
- The SLURM scripts are sanitised templates. Replace environment variables with paths available in the computing environment in which the analysis is run. They are not represented as byte-for-byte copies of the historical jobs.

## Software

The archived SLURM scripts load the UQ LDSC module `ldsc/1.0.1--pyhdfd78af_2/module`. The plotting workflow uses R and the `corrplot` package. Record the exact R and package versions used for any future execution; the current cleaned script has not been run against the complete original HPC environment.

## Data availability and responsible sharing

Do not commit individual-level, restricted, unpublished, or collaborator-controlled data. Before making this repository public, confirm that the summary-statistic-derived result table, mapping file, and figure outputs can be shared under the relevant data-use and institutional agreements. No raw GWAS summary-statistic files or reference panels are included in this repository build.
