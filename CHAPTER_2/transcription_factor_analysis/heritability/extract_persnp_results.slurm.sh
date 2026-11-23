#!/bin/bash
#SBATCH --job-name=extract_persnp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/extract_persnp_%j.out
#SBATCH --error=logs/extract_persnp_%j.err
#SBATCH --partition=general
#SBATCH --account=a_palpant

# Collect the raw PerSnpHsq records from SBayesRC .parSetRes outputs.
# This reproduces the visible grep/find step in Heritability_plots_M1_report.R.
# It deliberately does NOT parse columns: the exact .parSetRes schema was not
# included in the supplied files and must be checked before making a TSV.
set -euo pipefail

SEARCH_DIR="${SEARCH_DIR:-.}"
OUTPUT_FILE="${OUTPUT_FILE:-ce_hopx_raw_persnp_lines.txt}"
mkdir -p "$(dirname "$OUTPUT_FILE")"

mapfile -d '' result_files < <(find "$SEARCH_DIR" -type f -name '*.parSetRes' -print0)
if ((${#result_files[@]} == 0)); then
  echo "No .parSetRes files found under: $SEARCH_DIR" >&2
  exit 1
fi

grep -H 'PerSnpHsq' "${result_files[@]}" > "$OUTPUT_FILE" || {
  echo "No PerSnpHsq records found in discovered files." >&2
  exit 1
}
echo "Collected raw PerSnpHsq records in $OUTPUT_FILE"
