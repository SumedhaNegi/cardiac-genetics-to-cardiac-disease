# Cardiac Transcription-Factor Binding and Heritability

**Analysis scripts and workflows for the transcription-factor binding-footprint component of my PhD thesis.**

This repository documents the analysis of genomic regions associated with cardiac transcriptional regulators and their contribution to the genetic architecture of cardiac traits.

It is designed as a companion to [cardiac-genetics-to-cardiac-disease](https://github.com/SumedhaNegi/cardiac-genetics-to-cardiac-disease), using a similar examiner-facing structure: workflow documentation, methods summary, script review, curated metadata, and clearly separated archived outputs.

## Scientific scope

The analysis compares HOPX-associated genomic regions with five established cardiac transcriptional regulators:

- HOPX
- GATA4
- MEF2A
- NKX2-5
- TBX5
- SRF

The central questions are:

1. How much overlap exists among the genes and genomic regions associated with these regulators?
2. Are the associated gene sets enriched for cardiac and signalling pathways?
3. Are TF-associated genomic annotations enriched for SNP-based heritability of cardiac quantitative traits and cardiovascular disease?
4. Do adult-specific, foetal-specific and shared regulatory regions make different contributions to cardiac-trait heritability?

The thesis reports substantial overlap among the regulator-associated gene sets, enrichment for cardiac contraction and relevant signalling pathways, and consistently higher mean heritability enrichment for adult-specific than foetal-specific regulatory regions. These results are interpreted as annotation-level evidence, not proof that individual genes or intervals are causal.

## Where to start

1. Read [`docs/analysis_workflow.md`](docs/analysis_workflow.md) for the analysis sequence and inputs.
2. Read [`docs/methods_examiner_summary.md`](docs/methods_examiner_summary.md) for the thesis-aligned methods and findings.
3. Read [`docs/script_review.md`](docs/script_review.md) for the audit of the supplied scripts and known caveats.
4. Start with [`scripts/01_summarise_gene_sets.R`](scripts/01_summarise_gene_sets.R) to inspect the curated gene sets.

## Repository structure

```text
transcription-factor-analysis/
├── README.md
├── CITATION.cff
├── .gitignore
├── config/
│   └── analysis_config.yml
├── data/
│   ├── README.md
│   ├── gene_lists/
│   ├── gene_coordinates/
│   └── developmental_gene_sets/
├── docs/
│   ├── analysis_workflow.md
│   ├── methods_examiner_summary.md
│   └── script_review.md
├── R/
│   └── utils.R
├── scripts/
│   ├── 01_summarise_gene_sets.R
│   ├── 02_build_snp_annotations.R
│   ├── 03_summarise_sbayesrc.R
│   └── 04_plot_heritability_enrichment.R
├── hpc/
│   └── sbayesrc_array_template.sh
├── tests/
│   └── test_gene_sets.R
└── results/
    └── figures_archive/
```

## Software

The cleaned scripts use R and the following packages:

- `data.table`
- `dplyr`
- `stringr`
- `readr`
- `ggplot2`
- `yaml`
- `UpSetR`
- `GenomicRanges`
- `IRanges`

The historical HPC scripts used R 4.2.1 and SBayesRC in a UQ QRIS environment. Exact versions should be recorded from the environment used for any new run. The cleaned scripts have not been validated against the original HPC environment or the complete GWAS/LD inputs.

Install the CRAN dependencies:

```r
install.packages(c(
  "data.table", "dplyr", "stringr", "readr",
  "ggplot2", "yaml", "UpSetR"
))
```

Install Bioconductor dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GenomicRanges", "IRanges"))
```

## Running the cleaned workflow

Run commands from the repository root.

### 1. Summarise TF gene sets

```bash
Rscript scripts/01_summarise_gene_sets.R
```

Outputs:

- `results/tables/tf_gene_set_sizes.tsv`
- `results/tables/tf_pairwise_jaccard.tsv`
- `results/tables/tf_gene_membership.tsv`

### 2. Build SNP annotations

This step requires a SNP reference table containing chromosome, rsID and physical position, plus a confirmed genome assembly compatible with the gene-coordinate annotations.

First edit `config/analysis_config.yml` so `snp_columns` match the supplied reference table. Then run:

```bash
Rscript scripts/02_build_snp_annotations.R \
  /path/to/snp.info \
  GRCh37 \
  results/annotations
```

Replace `GRCh37` with the assembly actually used by the input files. The example is illustrative, not a claim about the supplied coordinate files.

### 3. Collate SBayesRC outputs

Place per-trait `*.AnnoPerSnpHsqEnrichment` files in `data/sbayesrc_results/` (this directory is intentionally not populated with restricted or unavailable source outputs).

```bash
Rscript scripts/03_summarise_sbayesrc.R
```

### 4. Plot heritability enrichment

```bash
Rscript scripts/04_plot_heritability_enrichment.R
```

The plot includes a reference line at enrichment = 1. Error bars are displayed only when the SBayesRC input table contains an `SD` column.

### 5. Plot gene-set overlap

```bash
Rscript scripts/05_plot_gene_set_overlap.R
```

This creates an UpSet plot and pairwise Jaccard heatmap.

### 6. Run input checks

```bash
Rscript tests/test_gene_sets.R
```

## Important interpretation and provenance notes

- **Species and genome build:** The archive contains datasets and gene identifiers from different sources, including mouse-heart occupancy datasets and human cardiac/GWAS resources. Confirm species, orthology mapping and genome assembly before combining annotations.
- **Coordinate convention:** Confirm 0-based versus 1-based coordinates and interval-end semantics for every source. The scripts normalise reversed endpoints but cannot infer the coordinate convention.
- **Window definitions:** Both expanded-window and no-window coordinate files are retained as distinct inputs. Do not mix them in a single analysis.
- **SNP reference:** The supplied historical SNP membership matrices contain rsIDs and binary membership values, not genomic coordinates. The cleaned SNP-overlap script therefore requires a compatible SNP reference file.
- **SBayesRC:** Historical logs show failures when annotation columns had zero or very low variance. Inspect annotation prevalence before model fitting and document any filtering decisions.
- **Reproducibility:** GWAS summary statistics, LD reference panels, complete HPC task lists and some source annotation resources were not included in the supplied archive. This repository is therefore a cleaned and documented analysis record, not a complete one-click reproduction of all historical results.
- **Data sharing:** Do not commit participant-level data, restricted GWAS summary statistics, collaborator-controlled data or third-party reference panels without permission.

## Archived figures

`results/figures_archive/` contains selected PDFs from the supplied working directory. These are retained as historical outputs and are not regenerated automatically by the cleaned scripts.

## Related thesis and repository

- Thesis: *Cardiac Genetics to Cardiac Disease: Statistical Genomics Analyses and Functional Characterisation of HOPX and its variants*, The University of Queensland, 2026.
- Companion repository: [cardiac-genetics-to-cardiac-disease](https://github.com/SumedhaNegi/cardiac-genetics-to-cardiac-disease).

## Citation

See [`CITATION.cff`](CITATION.cff). Add the final thesis DOI or public repository release DOI when available.

## License and release status

The scripts have been reorganised and annotated for public code review. Before making the repository public, confirm the permitted reuse of thesis-derived figures, gene annotations, and any third-party source data. A code licence should be selected once those ownership and reuse checks are complete.
