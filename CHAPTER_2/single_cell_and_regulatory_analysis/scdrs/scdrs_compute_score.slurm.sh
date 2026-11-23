#!/bin/bash
#SBATCH --job-name=scdrs_score
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=logs/scdrs_%j.out
#SBATCH --error=logs/scdrs_%j.err
#SBATCH --partition=general
#SBATCH --account=a_palpant

# Template reconstructed from the scDRS compute-score blocks in the source.
# Edit the module and paths for the selected dataset before submitting.
set -euo pipefail
mkdir -p logs

# module load python/3.10.4-gcccore-11.3.0
# source activate <environment-containing-scdrs>

H5AD_FILE="${H5AD_FILE:?Set H5AD_FILE to the input .h5ad}"
COV_FILE="${COV_FILE:?Set COV_FILE to the cell covariate TSV}"
GS_FILE="${GS_FILE:?Set GS_FILE to the scDRS .gs gene-set file}"
OUT_FOLDER="${OUT_FOLDER:-scdrs_scores}"
SPECIES="${SPECIES:-human}"
N_CTRL="${N_CTRL:-1000}"

mkdir -p "$OUT_FOLDER"
scdrs compute-score \
  --h5ad_file "$H5AD_FILE" \
  --h5ad_species "$SPECIES" \
  --cov_file "$COV_FILE" \
  --gs_file "$GS_FILE" \
  --gs_species "$SPECIES" \
  --flag_filter_data True \
  --flag_raw_count True \
  --n_ctrl "$N_CTRL" \
  --flag_return_ctrl_raw_score False \
  --flag_return_ctrl_norm_score True \
  --out_folder "$OUT_FOLDER"
