# =============================================================================
# EXAMINER GUIDE
# Script: sbayesrc_array_template.sh
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Parameterised SLURM array template for running SBayesRC across a task list; requires environment variables for summary statistics, LD, annotation, output and task configuration.
#

#!/usr/bin/env bash
#SBATCH --job-name=tf-sbayesrc
#SBATCH --array=1-29
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Generic, sanitised SLURM template. Configure paths and task list for your HPC.
set -euo pipefail

: "${SUMSTATS_DIR:?Set SUMSTATS_DIR to the directory containing formatted GWAS files}"
: "${LD_DIR:?Set LD_DIR to the compatible LD reference directory}"
: "${ANNOTATION_FILE:?Set ANNOTATION_FILE to the SNP annotation matrix}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to a writable results directory}"
: "${CONFIG_FILE:?Set CONFIG_FILE to a two-column task-index / filename table}"

mkdir -p logs "${OUTPUT_DIR}"
TASK_FILE="$(awk -v id="${SLURM_ARRAY_TASK_ID}" '$1 == id {print $2}' "${CONFIG_FILE}")"
if [[ -z "${TASK_FILE}" ]]; then
  echo "No input configured for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

module load R/4.2.1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

Rscript -e 'SBayesRC::sbayesrc(
  mafile = commandArgs(trailingOnly=TRUE)[1],
  LDdir = commandArgs(trailingOnly=TRUE)[2],
  outPrefix = commandArgs(trailingOnly=TRUE)[3],
  annot = commandArgs(trailingOnly=TRUE)[4],
  log2file = TRUE
)' "${SUMSTATS_DIR}/${TASK_FILE}" "${LD_DIR}" "${OUTPUT_DIR}/$(basename "${TASK_FILE}")_tf" "${ANNOTATION_FILE}"
