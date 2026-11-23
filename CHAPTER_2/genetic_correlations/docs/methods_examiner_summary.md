# Methods summary: genetic-correlation analysis

Genome-wide genetic correlations between cardiac quantitative traits and cardiovascular disease phenotypes were estimated using bivariate linkage disequilibrium score regression (LDSC) from GWAS summary statistics. Summary statistics were prepared using the LDSC `munge_sumstats.py` utility, including alignment to the HapMap3 reference SNP list. Pairwise genetic correlations were then estimated using `ldsc.py --rg`, with European-ancestry LD-score and regression-weight reference files.

The resulting estimates included genetic correlation coefficients (rg), standard errors, z-scores and p-values. Multiple testing was addressed using the Benjamini–Hochberg false-discovery-rate procedure. The thesis reports 1,485 pairwise genetic correlations among quantitative cardiac traits and describes correction across 10,877 tested correlations for the quantitative-trait/disease comparison. Trait annotations were used to organise the estimates by phenotype class and to generate correlation matrices and summaries of significant associations.

To examine relationships among traits, genetic correlations were transformed to dissimilarities using `d = 1 - rg`, and hierarchical clustering was performed using the Ward.D2 method. The resulting clusters were used to describe groups of traits with similar genome-wide genetic architectures. These clusters are descriptive and do not imply causal relationships.

**Documentation note:** This text reflects the analysis description in the thesis and the commands/scripts present in the supplied archive. The archive's result table contains repeated and conflicting records for some trait pairs; therefore, the exact final row-selection and FDR provenance should be confirmed from the final LDSC output/logs before claiming that the archived table alone reproduces the submitted figures.


## Thesis locations

The genetic-correlation workflow is described in Chapter 2, section 2.1. Figure 2.1 presents correlations among quantitative cardiac traits (thesis page 19); Figure 2.3 presents quantitative-trait/disease correlations and describes the Benjamini–Hochberg correction across 10,877 tests (thesis page 23); Figure 2.4 describes hierarchical clustering using `d = 1 - rg` and Ward.D2 (thesis page 25).
