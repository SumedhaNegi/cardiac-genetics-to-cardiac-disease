# Multiplex scDRS workflow components

The original `scDRS_multiplex.R` was a 954-line research-session notebook mixing R, Python, Bash, SLURM directives and pasted terminal commands. It is preserved unchanged under `_source_archive/single_cell_analysis/`.

The curated components separate the parts whose intent and inputs were sufficiently clear:

- `create_scdrs_covariates.py`: detected-gene count per cell from H5AD, matching the original `.cov` preparation.
- `attach_cell_ids_to_h5ad.py`: attaches a cell annotation table to `AnnData.obs` while preserving cell order.
- `scdrs_compute_score.slurm.sh`: parameterised SLURM template for `scdrs compute-score`.

Example commands:

```bash
python create_scdrs_covariates.py --h5ad input.h5ad --output cell_covariates.tsv
python attach_cell_ids_to_h5ad.py --h5ad input.h5ad --cell-table cell_ids.tsv --output annotated.h5ad
```

Before submitting the SLURM template, export `H5AD_FILE`, `COV_FILE`, `GS_FILE`, and optionally `OUT_FOLDER`, `SPECIES`, `N_CTRL`. Select the matching species and gene-set file; the historical source includes both human and mouse contexts and multiple H5AD datasets.

The original MAGMA commands, Seurat/Harmony plots, gene-set batching and enrichment heatmaps are not presented as a single executable workflow because the original contains conflicting paths, alternative data versions, undefined objects, malformed fragments and incomplete result schemas. Consult the archived original for provenance. No analysis was run or tested during this cleanup.
