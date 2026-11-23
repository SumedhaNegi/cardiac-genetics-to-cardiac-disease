# Cardiac Genetics to Cardiac Disease — Thesis Methods Repository

This repository brings together the supplied analysis scripts and supporting material for the thesis
**Cardiac Genetics to Cardiac Disease: Statistical Genomics Analyses and Functional Characterisation of HOPX and its variants**.

## Purpose for examiners

This is a readable record of the analytical methods and the development of the thesis analyses. It is **not claimed to be a fully reproducible software pipeline** as some scripts depend on restricted datasets or UQ HPC paths/modules.

Each curated code/job file has an examiner guide and a file-specific code map at the top, with inline comments retained and expanded where the method is clear. Original source collections are retained unchanged under `_source_archive/` so that edits do not erase provenance.

## Where to start

1. [`docs/THESIS_METHODS_OVERVIEW.md`](docs/THESIS_METHODS_OVERVIEW.md) — how the code is grouped by thesis question.
2. [`docs/SCRIPT_INVENTORY.md`](docs/SCRIPT_INVENTORY.md) — complete inventory and locations.
3. [`docs/SCRIPT_AUDIT.md`](docs/SCRIPT_AUDIT.md) — known source issues and interpretation cautions.
4. Read the README in each chapter folder for workflow-specific context.

## Thesis-oriented structure

```text
CHAPTER_2/
  genetic_correlations/                  # LDSC, summary-statistic preparation and correlation figures
  transcription_factor_analysis/         # DamID/TF binding, gene sets, SNP annotations, SBayesRC and enrichment
  single_cell_and_regulatory_analysis/   # scDRS, TRIAGE, developmental single-cell and gene-ID mapping
CHAPTER_3/
  hopx_functional_characterisation/
    cardioexcyte/                         # CardioExcyte data processing and thesis EFP figures
    efpa/                                 # EFPA waveform plotting
    flow_cytometry/                       # FlowCore control-based gating and BFP analysis
    qPCR/                                 # qPCR workflow
    gene_expression_correlation/          # HOPX-related expression/correlation scripts
CHAPTER_4/
  rare_variant_interpretation/            # AlphaMissense/ClinVar visualisation
docs/
_source_archive/                          # unmodified supplied source collections
```

The chapter mapping is a practical organisation of the supplied scripts, not a claim that every script corresponds one-to-one to a thesis figure. For example the SbayesRC and enrichment scripts contribute to multiple figures throughout the thesis with different input files.


- Some historical files contain mixed R, Python, Bash, SLURM/PBS commands, terminal prompts or copied console output. The scDRS notebook has been split into language-specific components where the intent was clear; its full source remains archived. 
- `Triage_discordance_dubois.R` is the canonical reconstructed workflow. `Triage_discordance_dubois_1.R` is a short provenance pointer; the complete duplicate draft remains in `_source_archive/`.
- Historical absolute paths (for example `/QRISdata/`, `/scratch/` and user home directories) are preserved where present in originals. They identify the original computing context and will not resolve on other systems.
- No raw individual-level data, restricted GWAS summary statistics or institutional reference panels should be added to a public repository without checking data-use permissions.
