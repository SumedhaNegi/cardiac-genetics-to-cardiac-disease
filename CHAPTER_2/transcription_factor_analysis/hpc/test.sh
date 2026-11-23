# =============================================================================
# EXAMINER GUIDE
# Script: test.sh
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Historical test job script; useful as evidence of the cluster environment and commands, not a production workflow.
#

#!/usr/bin/env bash
# Examiner guide: historical SLURM job template from test.sh.
# Purpose: records resource requests and commands used to launch an analysis.
# Adapt partition, account, modules, paths and executable to the target HPC system.
# The associated analysis code is documented in the relevant chapter folder.
# This job file has not been submitted or tested during repository curation.

#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=80GB
#SBATCH --time=01:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=debug
#SBATCH --job-name=tf_card
#SBATCH --account=a_palpant

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/test/

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file="svi_bsa.ma_imp.ma"
ld_folder="/scratch/user/uqsnegi/ukbEUR_Imputed"
annot="/QRISdata/Q5149/TF_DamID/annotation/7m_nowindow.txt"
threads=4
 
export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$ld_folder', \
				 outPrefix='${ma_file}_7m_nowindow.txt', annot='$annot', log2file=TRUE)"


