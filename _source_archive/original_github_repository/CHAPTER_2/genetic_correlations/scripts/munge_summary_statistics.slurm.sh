#!/usr/bin/env bash
#SBATCH --job-name=ldsc_munge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --output=munge_%j.out
#SBATCH --error=munge_%j.err

# Examiner-facing, sanitised version of the LDSC munging workflow.
# The original workflow processed six source directories/file patterns.
# Supply one source directory and pattern per submission. This avoids the
# original SLURM-array/Bash-index mismatch (SLURM arrays start at 1; Bash
# indexed arrays start at 0) and makes the inputs explicit in the job log.
#
# Required environment variables:
#   INPUT_DIR       directory containing source GWAS summary statistics
#   FILE_PATTERN    shell glob, e.g. '*.ma' or '*.tsv'
#   OUTPUT_DIR      destination for LDSC .sumstats.gz files
#   HM3_SNPLIST     path to LDSC HapMap3 SNP list
# Optional:
#   LDSC_MODULE     module name; defaults to the UQ LDSC module used originally

set -euo pipefail

: "${INPUT_DIR:?Set INPUT_DIR to the directory containing input summary statistics}"
: "${FILE_PATTERN:?Set FILE_PATTERN, for example '*.ma'}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR for munged summary statistics}"
: "${HM3_SNPLIST:?Set HM3_SNPLIST to the HapMap3 SNP list}"
LDSC_MODULE="${LDSC_MODULE:-ldsc/1.0.1--pyhdfd78af_2/module}"

module load "$LDSC_MODULE"
mkdir -p "$OUTPUT_DIR"
[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$HM3_SNPLIST" ]] || { echo "HapMap3 SNP list not found: $HM3_SNPLIST" >&2; exit 1; }

shopt -s nullglob
files=( "$INPUT_DIR"/$FILE_PATTERN )
(( ${#files[@]} > 0 )) || { echo "No files matched $INPUT_DIR/$FILE_PATTERN" >&2; exit 1; }

for input_file in "${files[@]}"; do
    base_name=$(basename "$input_file")
    output_prefix="$OUTPUT_DIR/$base_name"
    echo "Munging: $input_file"
    munge_sumstats.py \
        --sumstats "$input_file" \
        --out "$output_prefix" \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
done
