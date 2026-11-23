# =============================================================================
# EXAMINER GUIDE
# Script: run_ldsc_genetic_correlations.slurm.sh
# Purpose: Supports summary-statistic preparation, LDSC genetic-correlation analysis or visualisation; see chapter workflow documentation.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# SLURM job wrapper for pairwise LDSC genetic correlations; configures Python/LDSC environment, trait pairs and output locations.
#

#!/usr/bin/env bash
#SBATCH --job-name=ldsc_rg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --output=ldsc_rg_%j.out
#SBATCH --error=ldsc_rg_%j.err

# Examiner-facing, sanitised form of the archived all-by-all LDSC command.
# LDSC --rg accepts a comma-separated list of munged summary-statistic files
# and estimates pairwise genetic correlations among them.
#
# Required environment variables:
#   SUMSTATS_DIR  directory containing *.sumstats.gz files
#   REF_LD_CHR    European LD-score reference prefix/directory
#   W_LD_CHR      regression-weight reference prefix/directory
#   OUTPUT_PREFIX output prefix for LDSC results
# Optional LDSC_MODULE as in the munging script.

set -euo pipefail
: "${SUMSTATS_DIR:?Set SUMSTATS_DIR}"
: "${REF_LD_CHR:?Set REF_LD_CHR}"
: "${W_LD_CHR:?Set W_LD_CHR}"
: "${OUTPUT_PREFIX:?Set OUTPUT_PREFIX}"
LDSC_MODULE="${LDSC_MODULE:-ldsc/1.0.1--pyhdfd78af_2/module}"
module load "$LDSC_MODULE"

shopt -s nullglob
files=( "$SUMSTATS_DIR"/*.sumstats.gz )
(( ${#files[@]} >= 2 )) || { echo "Need at least two .sumstats.gz files in $SUMSTATS_DIR" >&2; exit 1; }
[[ -d "$REF_LD_CHR" ]] || { echo "Reference LD-score path not found: $REF_LD_CHR" >&2; exit 1; }
[[ -d "$W_LD_CHR" ]] || { echo "Regression-weight path not found: $W_LD_CHR" >&2; exit 1; }

# Build the comma-separated list without word-splitting filenames.
rg_files=$(IFS=,; echo "${files[*]}")
ldsc.py \
    --rg "$rg_files" \
    --ref-ld-chr "$REF_LD_CHR" \
    --w-ld-chr "$W_LD_CHR" \
    --out "$OUTPUT_PREFIX"
