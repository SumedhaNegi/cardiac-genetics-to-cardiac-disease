# Thesis methods overview

## Chapter 2 — From cardiac genetics to regulatory and cellular mechanisms

### Genetic correlations
The repository includes summary-statistic preparation and linkage disequilibrium score regression (LDSC) job scripts, plus R code for examining and visualising genetic-correlation results. The LDSC workflow estimates genome-wide genetic correlation between traits; it does not itself establish causality. The archived results and plotting scripts should be interpreted with the caveats in the chapter-specific documentation.

### Transcription-factor binding and regulatory annotations
The transcription-factor-analysis material includes gene-set summaries, genomic/SNP annotation preparation, SBayesRC result summaries, heritability enrichment and gene-set overlap plots. Historical DamID and enhancer-heritability scripts are retained alongside the more structured scripts. These workflows connect transcription-factor binding/gene sets and genomic annotations with genetic association or heritability evidence.

### Single-cell and regulatory-state analyses
The supplied single-cell scripts include scDRS workflows, TRIAGE conversion/discordance workflows, WSB developmental single-cell analysis, gene-expression correlation and gene identifier conversion. Some scripts mix R/Python with HPC session commands and depend on external Seurat objects, HDF5 files, gene sets and reference mappings. Species and genome-build assumptions must be checked per dataset.

## Chapter 3 — Functional characterisation of HOPX and its variants

### qPCR and gene-expression analysis
The qPCR and expression-correlation scripts record workflows for processing expression measurements and exploring relationships between genes/conditions. Input spreadsheets and experimental metadata are not all included.

### Flow cytometry
The FlowCore script documents control-based gating, including the sequence of main-cell, singlet, viable-cell and BFP-positive gates. Gate names and channels are source-specific and should be read alongside the original FCS files and controls.

### CardioExcyte / extracellular field-potential analysis
`cardioexcyte.R` reads a tab-delimited instrument export, extracts sweep/parameter/time metadata, reshapes measurements to long format and maps wells to experimental conditions. `cardioexcyte_thesis_figures.R` uses an in-memory `mydata` object to reconstruct variant and parameter labels, excludes sweeps 001–002, estimates noradrenaline addition between sweeps 006 and 007, defines baseline and post-treatment windows, calculates well-specific baselines, averages wells within variant/sweep, and creates figures and ANOVA/Dunnett tables.

The figure script's ANOVA uses the resulting sweep-level observations. Whether those observations satisfy independence assumptions depends on the experiment's repeated-measures design; this repository does not resolve that design question.

### EFPA waveform analysis
The EFPA plotting script uses mean-beat and parameter CSV files to visualise waveform features. The expected filenames and workflow assumptions are described in its source header and audit.

## Chapter 4 — Rare HOPX variant interpretation

The AlphaMissense/ClinVar heatmap script records a variant-annotation visualisation workflow. Interpretation depends on the source annotation release, transcript choice and variant nomenclature; those details should be checked against the original analysis records.

## Species and genome-build handling

The source material includes both human and mouse analyses. A shared folder structure does not constitute cross-species harmonisation. Gene symbols may differ, orthologues may not be one-to-one, and genomic coordinates are assembly-specific. Preserve species, Ensembl IDs, symbols, assembly and annotation release as metadata. Where a cross-species comparison was performed, its mapping table and version should accompany that analysis; this collection does not claim a single validated global mapping pipeline.


## Code annotation and repair status

The curated chapter copies include file-specific code maps and examiner-oriented comments. The repair record and unresolved items are listed in [`SCRIPT_AUDIT.md`](SCRIPT_AUDIT.md). The source archive is intentionally unmodified. Incomplete historical fragments and workflows requiring unavailable data/session objects remain identified as such; they have not been filled with guessed scientific steps. No code was run or tested.
