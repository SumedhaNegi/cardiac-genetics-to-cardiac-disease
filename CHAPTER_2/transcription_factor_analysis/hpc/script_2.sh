# =============================================================================
# EXAMINER GUIDE
# Script: script_2.sh
# Purpose: Supports transcription-factor gene-set processing, genomic annotation, SBayesRC summaries or heritability/overlap visualisation; see the folder README.
# Inputs/outputs: Consult the script body and the chapter README for expected
# files, object names and outputs. Institutional paths and external references
# may require the original UQ environment.
# Curation: Existing analytical code is retained; this header adds navigation
# context only. No execution or testing was performed during curation.
# =============================================================================

# CODE MAP
# Historical second HPC submission/workflow variant; retained to document alternate job configuration and command sequence.
#

#!/usr/bin/env bash
# Examiner guide: historical SLURM job template from script 2.sh.
# Purpose: records resource requests and commands used to launch an analysis.
# Adapt partition, account, modules, paths and executable to the target HPC system.
# The associated analysis code is documented in the relevant chapter folder.
# This job file has not been submitted or tested during repository curation.

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=sbayesrc
#SBATCH --error=err
#SBATCH --output=out
#SBATCH --partition=general
#SBATCH --account=a_palpant

gctb=/QRISdata/Q5149/sbayesRC/gctb


gwas=/QRISdata/Q5149/sumstats/formatted/from_dalia/7m_imp_sumstats/Formatted_20015_irnt.ma_imp.ma

ldm=/scratch/user/uqsnegi/ukbEUR_Imputed
anno=/scratch/user/uqsnegi/nick_heritability/upregulated_in_replated.txt
out=Formatted_20015_irnt_nick_intercept


$gctb --ldm-eigen $ldm --gwas-summary $gwas --sbayes RC --annot $anno --thread 10 --out $out >$out.log 2>&1
