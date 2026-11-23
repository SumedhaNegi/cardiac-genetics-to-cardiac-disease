# Data and provenance

The archive supplied for this repository contained gene lists, gene-coordinate
tables, SNP-membership matrices, SBayesRC-related files, R/RData objects, PDFs,
and HPC shell scripts.

This public-facing build includes only selected gene-list and coordinate inputs
plus a small set of historical PDF figures. Large intermediate matrices,
workspace dumps, GWAS summary statistics, reference panels, and the nested
genetic-correlation archive are intentionally not copied into the repository.

## Included data

- `gene_lists/`: gene identifiers associated with HOPX and five comparator TFs.
- `gene_coordinates/`: source coordinate tables, including historical 50-kb
  expanded and no-window variants. Files are preserved as supplied; scripts
  normalise reversed interval endpoints before overlap calculations.
- `developmental_gene_sets/`: foetal-upregulated, foetal-downregulated,
  foetal/adult-significant and non-significant gene lists.

## Important data checks before analysis

1. Confirm species and genome assembly for every coordinate file.
2. Confirm whether coordinates are 0-based or 1-based and whether interval ends
   are inclusive. This matters for SNP overlap calculations.
3. Confirm the coordinate reference matches the GWAS SNP reference.
4. Check whether genes were assigned from mouse or human datasets and whether
   any orthology mapping was applied.
5. The archive's historical SNP membership tables contain rsIDs and binary
   membership but not chromosome/position fields. The cleaned SNP-overlap
   script therefore requires a separate compatible `snp.info`-style file.
6. Do not upload UK Biobank participant-level data, restricted summary
   statistics, controlled-access data, or third-party LD reference panels.

The repository does not claim that these inputs alone reproduce every historical
analysis result.
