#!/bin/bash
#SBATCH --job-name=enhancer_sbayes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=15:00:00
#SBATCH --output=logs/enhancer_%j.out
#SBATCH --error=logs/enhancer_%j.err
#SBATCH --partition=general
#SBATCH --account=a_palpant

# Reconstructed SBayes (GCTB) invocation from the historical source.
# Configure paths for the same genome build / LD reference / GWAS as the SNP
# annotation file. Confirm the GCTB version and command-line options locally.
set -euo pipefail

GCTB="${GCTB:?Set GCTB to the gctb executable path}"
LDM_LIST="${LDM_LIST:?Set LDM_LIST to the sparse LD matrix list}"
GWAS_SUMMARY="${GWAS_SUMMARY:?Set GWAS_SUMMARY to the formatted GWAS .ma file}"
ANNOTATION="${ANNOTATION:?Set ANNOTATION to the SNP annotation table}"
OUT_PREFIX="${OUT_PREFIX:-cardiac_enhancer}"
CHAIN_LENGTH="${CHAIN_LENGTH:-12000}"
BURN_IN="${BURN_IN:-2000}"

mkdir -p "$(dirname "$OUT_PREFIX")" logs
"$GCTB" --sbayes S \
  --mldm "$LDM_LIST" \
  --gwas-summary "$GWAS_SUMMARY" \
  --annot "$ANNOTATION" \
  --chain-length "$CHAIN_LENGTH" \
  --burn-in "$BURN_IN" \
  --out "$OUT_PREFIX" \
  > "${OUT_PREFIX}.log" 2>&1
