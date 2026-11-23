# =============================================================================
# EXAMINER GUIDE
# Script: test2.sh
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Historical second test job script; retain as provenance and verify scheduler directives and paths before reuse.
#

#!/usr/bin/env bash
# Examiner guide: historical SLURM job template from test2.sh.
# Purpose: records resource requests and commands used to launch an analysis.
# Adapt partition, account, modules, paths and executable to the target HPC system.
# The associated analysis code is documented in the relevant chapter folder.
# This job file has not been submitted or tested during repository curation.

#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=1:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=debug
#SBATCH --job-name=cluster
#SBATCH --account=a_palpant
#SBATCH --array=1-19

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/test
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/computation/dubois_development2022/supplementary_tables/annotation/table1_annotation_with_windows.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/supplementary_tables/results"
threads=4
 
export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}/${ma_file}_clusters.txt', annot='$annot', log2file=TRUE)"
