#!/usr/bin/env python3
"""Create scDRS cell-level covariates from an AnnData/H5AD file.

Reconstructs the source workflow: number of detected genes per cell, written
as a tab-separated file with cell_names and n_genes columns. This is a
cell-level QC covariate; it is not a normalised expression matrix.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def detected_gene_counts(matrix) -> np.ndarray:
    """Count nonzero genes per cell for dense or sparse AnnData matrices."""
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).sum(axis=1)).ravel().astype(np.int64)
    arr = np.asarray(matrix)
    return np.count_nonzero(arr > 0, axis=1).astype(np.int64)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Input AnnData .h5ad file")
    parser.add_argument("--output", required=True, help="Output covariate TSV")
    args = parser.parse_args()
    adata = ad.read_h5ad(args.h5ad)
    counts = detected_gene_counts(adata.X)
    cov = pd.DataFrame({"cell_names": adata.obs_names.astype(str), "n_genes": counts})
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    cov.to_csv(out, sep="\t", index=False)
    print(f"Wrote {len(cov):,} cell covariates to {out}")

if __name__ == "__main__":
    main()
