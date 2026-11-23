#!/usr/bin/env python3
"""Attach a cell annotation table to AnnData.obs using explicit cell IDs.

The historical notebook merged a CellID column into .obs. This version checks
for duplicate IDs and reports unmatched cells instead of silently changing
cell ordering or creating ambiguous duplicate rows.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import anndata as ad
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--cell-table", required=True, help="TSV/CSV containing cell ID column")
    parser.add_argument("--cell-id-column", default="CellID")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    data = ad.read_h5ad(args.h5ad)
    sep = "\t" if Path(args.cell_table).suffix.lower() in {".tsv", ".txt"} else ","
    table = pd.read_csv(args.cell_table, sep=sep)
    if args.cell_id_column not in table.columns:
        raise ValueError(f"Missing ID column {args.cell_id_column!r}; found {list(table.columns)}")
    if table[args.cell_id_column].duplicated().any():
        raise ValueError("Cell annotation table contains duplicated cell IDs")
    table = table.set_index(args.cell_id_column)
    table.index = table.index.astype(str)
    obs_index = pd.Index(data.obs_names.astype(str))
    matched = obs_index.isin(table.index)
    print(f"Matched {matched.sum():,}/{len(matched):,} cells")
    # Reindex preserves the AnnData row order; non-matching cells receive NaN.
    aligned = table.reindex(obs_index)
    aligned.index = data.obs.index
    for col in aligned.columns:
        data.obs[col] = aligned[col].to_numpy()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(out)
    print(f"Wrote annotated AnnData object to {out}")

if __name__ == "__main__":
    main()
