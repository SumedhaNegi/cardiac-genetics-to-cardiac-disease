# =============================================================================
# EXAMINER GUIDE
# Script: script.sh
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Historical HPC job script for transcription-factor analysis; inspect its module, file and output paths before adapting to a cluster.
#

#!/usr/bin/env bash
# Examiner guide: historical SLURM job template from script.sh.
# Purpose: records resource requests and commands used to launch an analysis.
# Adapt partition, account, modules, paths and executable to the target HPC system.
# The associated analysis code is documented in the relevant chapter folder.
# This job file has not been submitted or tested during repository curation.

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=150GB
#SBATCH --time=5:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=general
#SBATCH --account=a_palpant
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/imputed_sumstats

config=config.txt

cd $MY_DIR


# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/ukbEUR_HM3/"
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/top_50_annot_weighted_and_bin.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/result_50/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_50_quant_bin', annot='$annot', log2file=TRUE)"


