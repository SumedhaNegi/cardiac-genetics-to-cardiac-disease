# Script audit and repair record

This document describes a source-level editorial review only. **No scripts were run, syntax-checked, benchmarked or statistically tested**, as requested. The untouched supplied files remain under `_source_archive/` for provenance; curated chapter copies are the versions to read first.

## Repairs made in curated chapter copies

| File | Changes |
|---|---|
| `CHAPTER_3/.../cardioexcyte/cardioexcyte.R` | Rebuilt as a standalone R script; removed interactive `salloc`/`R` console commands; added command-line input/output paths, export-layout checks, tidy reshaping, explicit well mapping, data export and plots. The plate map remains experiment-specific and must be checked against the plate record. |
| `CHAPTER_2/.../transcription_factor_analysis/scripts/cardiac_enhancer_annpt.R` | Replaced `head()` on a path string with actual BED import; added 50 kb interval expansion, SNP coordinate import, overlap annotation and output. BED-to-GRanges coordinate conversion is explicit; genome assembly compatibility still requires confirmation. |
| `CHAPTER_2/.../triage/TRIAGE1_conversion.R` | Removed invalid placeholder code and `.GlobalEnv` side effects; added input validation and local reference loading; made gene-wise row scaling explicit with `sweep()`. The historical score/mapping file column conventions are retained and need to match the supplied reference files. |
| `CHAPTER_2/.../gene_id_mapping/ensemble to entrez id conversion.R` | Removed package installation from analysis code; corrected identifier conversion workflow, output separator and up/down output naming; made the Enrichr check optional. Assumes first worksheet column contains human Ensembl IDs. |
| `CHAPTER_2/genetic_correlations/summary_statistics/pip_snp.R` | Replaced invalid pasted shell/console fragments with native R sorting on historical PIP column 10; added file/column checks and output naming. |
| `CHAPTER_2/genetic_correlations/summary_statistics/Cleaning_sumstat.R` | Replaced mixed PBS/AWK/R fragments and malformed column selection with two explicit R formatting modes. Column positions are documented and must be checked against each source file; sample size is only overwritten when explicitly supplied. |
| `CHAPTER_4/.../alphamissense_heatmap.R` | Replaced undefined `data`/`merged` objects and interactive console fragments with a standalone plotting script. It requires V2 substitution labels and V3 scores; a ClinVar overlay is only shown if a pre-matched classification column is present. The missing ClinVar join was not invented. |
| `CHAPTER_2/genetic_correlations/historical_scripts/corr3.R` | Reconstructed the missing function header and matrix preparation around the surviving plotting body; added validation, explicit trait alignment, FDR handling and output-device selection. This is a documented reconstruction, not recovered original source. |
| `CHAPTER_2/transcription_factor_analysis/heritability/DamID_enhancer_heritability.R` | Rebuilt the BED/BIM annotation workflow; corrected the invalid SNP position indexing; added standard PLINK BIM parsing, BED coordinate conversion, chromosome-aware GenomicRanges overlap and optional WT/differential-accessibility BED inputs. |
| `CHAPTER_2/transcription_factor_analysis/heritability/Heritability_plots_M1_report.R` | Rebuilt plotting around a documented prepared-results schema; replaced undefined `case1`, `c1`, `h1` objects with a validated input table and reusable all/LV/RV/OFT plotting functions. Raw `.parSetRes` parsing remains unresolved. |
| `CHAPTER_2/single_cell_and_regulatory_analysis/triage/Triage_discordance_dubois.R` | Rebuilt the main TRIAGE workflow with configurable paths, Seurat RNA data-layer extraction, explicit mapping and weighting, sparse cell-type averaging, ranked gene outputs, top-N sets and optional SNP overlap annotation. |
| `CHAPTER_2/single_cell_and_regulatory_analysis/scdrs/` | Split clearly identifiable scDRS covariate generation and cell-ID attachment into Python scripts and provided a configurable SLURM compute-score template. Original mixed-language source is archived. |
| Mixed-language R scripts | Interactive HPC allocation/module/console-entry commands were changed to comments in curated R copies so they are not mistaken for R statements. Shell/Python commands retained as historical notes are identified as such. |

## Important historical material retained without speculative repair

- `ldsc_corrplot.R`, `ning_wallthickness.R`, `TF_damid_1.R`, `wsb_triage_heritability.R`, and other long interactive scripts still contain session fragments, pasted output, fixed institutional paths, or dependencies on objects created earlier in an R/Python session. Their visible methods are annotated, but they are not claimed to be clean standalone pipelines.
- `corr3.R` has been reconstructed as a function from the surviving plotting body. The missing original function header/setup could not be recovered; required input columns and assumptions are now explicit.
- `DamID_enhancer_heritability.R` has been rewritten as a configurable BED/BIM SNP-overlap workflow. The undefined historical objects `wt` and `c` are now optional BED inputs. Genome build compatibility remains a required scientific confirmation.
- `Heritability_plots_M1_report.R` has been rewritten to plot a prepared table with explicit columns. The original raw `.parSetRes` parser cannot be recovered because the file format/example result was not supplied. `extract_persnp_results.slurm.sh` collects raw lines only; it does not claim to parse them.
- `scDRS_multiplex.R` is now a workflow index. The source was separated into a Python covariate generator, a Python cell-ID annotation utility and a SLURM compute-score template. The original mixed-language notebook is preserved in `_source_archive`; MAGMA, Seurat visualisation and alternative batch fragments remain historical because their exact inputs/outputs are not sufficiently specified.
- `Triage_discordance_dubois.R` has been rewritten into a configurable pipeline using the source's RNA assay `data` layer, TRIAGE gene mapping/repressive scores, expression-times-score weighting, cell-type averaging, ranked gene lists and top-N subsets. Optional SNP annotation reproduces the source's gene-coordinate +/-50 kb overlap strategy. Reference table column names, geneMatrix columns, SNP-info columns and species must be confirmed against the actual source files.
- `Triage_discordance_dubois_1.R` is now a short provenance pointer to the canonical curated workflow; its full historical source remains unchanged in `_source_archive`. The available evidence did not establish that it was a separate final analysis.

## Scientific and reproducibility decisions that still need source-owner confirmation

1. **Species and gene identifiers:** Human and mouse material is present. There is no single end-to-end orthologue harmonisation workflow in the supplied source collection. Confirm species and identifier namespace at each analysis boundary.
2. **Genome assembly:** Confirm hg19/hg38 (or other build) for all BED, SNP, enhancer and gene-coordinate files before interpreting interval overlaps.
3. **Reference versions:** Record Ensembl release, gene annotation, LD reference, GWAS build, scRNA-seq object version and package versions in the thesis methods.
4. **External inputs:** Several workflows require data not included here, including GWAS summary statistics, LD panels, Seurat/H5AD objects, FCS files, qPCR exports, CardioExcyte instrument exports and ClinVar/AlphaMissense source tables.
5. **Experimental unit:** The CardioExcyte thesis-figure script averages wells and derives sweep-level summaries before statistical comparisons. Confirm the biological replicate and repeated-measures structure with the experimental design; code comments cannot establish independence.
6. **Variant matching:** The AlphaMissense source did not contain the ClinVar merge step or a documented variant-matching key. The cleaned plotting script only overlays a classification column when that match is already present in the input table.
7. **Historical summary-statistic column positions:** The two cleaning modes document the exact positional selections visible in the source. Check those positions against the original phenotype files before use.

## Annotation coverage

Every curated chapter code/job file has an examiner guide and a file-specific code map. The five priority workflows identified during review received substantive reconstruction or language separation, not only additional comments. Existing inline comments were retained. New comments explain the purpose of the major processing stages, inputs, outputs and assumptions where these can be established from the source. The `_source_archive/` copies are deliberately left untouched so the original material can be compared with the curated versions.

## What this review does not claim

- It does not claim all historical scripts are now executable.
- It does not claim missing data, objects, genome builds, orthologue maps or statistical design decisions have been reconstructed.
- It does not claim any script has been tested or any figure/result regenerated.
