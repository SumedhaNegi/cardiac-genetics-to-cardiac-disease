# =============================================================================
# EXAMINER GUIDE
# Script: scDRS_multiplex.R
# Purpose: Index for the multiplexed scDRS analysis, historically assembled
#          from R, Python, Bash/SLURM and interactive terminal sessions.
# Provenance: The original 954-line mixed-language notebook is preserved at
# _source_archive/single_cell_analysis/scDRS_multiplex.R.
# Cleaned components: See the adjacent Python and SLURM files in this folder.
# Important: The original contained multiple alternative datasets, paths,
#            batch sizes and output folders. These are not interchangeable.
#            The split scripts expose parameters and require local review.
# =============================================================================

# This file is an R-readable workflow index, not a runnable end-to-end pipeline.
# The historical source is retained in _source_archive. Language-specific
# executable components are provided alongside this file:
#
# 1. create_scdrs_covariates.py
#    Reads an H5AD file, counts detected genes per cell and writes a covariate TSV.
# 2. attach_cell_ids_to_h5ad.py
#    Adds a cell-ID annotation table to AnnData.obs without changing cell order.
# 3. scdrs_compute_score.slurm.sh
#    SLURM template for scDRS compute-score with explicit configurable paths.
#
# The original file also contains MAGMA preparation, gene-set batching, Seurat
# visualisation, UMAP plotting, result aggregation and heatmap fragments. Their
# input schemas and exact final selection are not fully specified. They remain
# visible in the archived source and are not merged into a speculative pipeline.
