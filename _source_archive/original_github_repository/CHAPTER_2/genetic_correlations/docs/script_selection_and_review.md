# Script selection and review notes

## Selected scripts in this repository

| Stage | Selected file | Basis for selection | Changes made for this repository |
|---|---|---|---|
| LDSC munging | `scripts/munge_summary_statistics.slurm.sh` | Based on the archived `munge_sumstats.sh`, which contains the six source-directory/file-pattern groups used for the broader input collection. | Recast as a parameterised template; added strict shell error handling, input checks, output-directory creation and safe glob handling. Corrected the original SLURM-array/Bash-array indexing hazard by removing the array index from the template. |
| LDSC genetic-correlation estimation | `scripts/run_ldsc_genetic_correlations.slurm.sh` | Based on the archived `corr.sh` all-files `ldsc.py --rg` approach rather than the two array scripts. | Added explicit environment variables, file-count checks, quoted comma-separated input construction and strict error handling. |
| Plotting and clustering | `scripts/plot_genetic_correlations.R` | Consolidates the most relevant sections of `Annotated_scripts/ldsc_corrplot.R` from the reviewed repository. This is the only reviewed plotting source located that explicitly contains Ward.D2 clustering code. | Removed exploratory fragments and repeated drafts; made the `d = 1 - rg` distance and `hclust(..., method = "ward.D2")` explicit; added input checks, symmetric matrix construction and duplicate-conflict checks; separates plotting from data preparation. |

## Why the other archived versions were not selected

- `array.sh` and `array2.sh`: both use SLURM array indices directly as Bash array indices. Since SLURM arrays start at 1 and Bash arrays start at 0, the indexing is offset; the scripts can skip the first file and select an empty/out-of-range item at the end. They also rely on an unquoted shell array expansion.
- `munge_sumstats.sh`: contains the six input groups but uses the SLURM task ID directly to index zero-based Bash arrays. This is an off-by-one error. The sanitised template replaces that fragile array selection with explicit per-job input variables.
- `munge.sh` and `munge_single_file_array.sh`: each is tied to one source directory and a fixed array length. They are less comprehensive than the six-group munging script and retain the same general array-indexing risk.
- `corr2 2.R`: refers to `sig_counts` without defining it in the function and contains an incorrect legend implication that square size encodes FDR. The plot's square size represents correlation magnitude; significance is represented separately.
- `corr3 2.R`: includes a barplot count calculation but still contains plotting/layout issues and does not implement the explicit thesis-reported `d = 1 - rg` Ward.D2 clustering pipeline.
- `corr.R`: is a simpler quantitative-versus-disease correlation plot, but does not include the clustering step.
- `corr2.R`, `corr3.R`, `corr4.R`, `disease.R`, `test.sh` and some RDS placeholders in the archive are empty files and were not selected.
- `debug_script.R` is an exploratory/debugging script and does not provide a reliable end-to-end plotting workflow.

## FDR handling: important caution

The reviewed `ldsc_corrplot.R` contains code that applies `p.adjust(p, method = "fdr")`, but a subsequent operation subtracts the minimum adjusted p-value from every q-value. That subtraction is not part of the Benjamini–Hochberg procedure and changes the adjusted values. I have therefore **not** copied that transformation into the cleaned plotting script. The thesis captions state that Benjamini–Hochberg FDR was used, and the archived `rg.txt` contains a `q` column; the cleaned plotting script treats that column as precomputed input.

The archived result table contains repeated rows and conflicting estimates for some canonical trait pairs. The script stops when it detects conflicting `rg`, `p` or `q` values. This is deliberate: selecting one version without the final output/logs would amount to guessing which result was used for the thesis.

## Scope of validation

The cleaned scripts have been reviewed statically for obvious shell/R logic issues. They have not been executed on Bunya, against the original full set of GWAS summary statistics, or with the institutional LD-score reference directories. The current plotting script requires a reconciled one-estimate-per-trait-pair result table before it can generate the final clustered plots.
