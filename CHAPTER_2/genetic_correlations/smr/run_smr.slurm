# =============================================================================
# EXAMINER GUIDE
# Script: run_smr.slurm
# Purpose: SLURM job script for an SMR (Summary-data-based Mendelian Randomization) workflow.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# SLURM job wrapper for SMR/HEIDI analysis; review the executable, reference files, summary-statistic inputs and output prefix before use.
#

#!/bin/bash

#**************************************************#
#     Run SMR  Gene -> Sz                             #
#**************************************************#


data_dir="/QRISdata/Q5149/computation/smr/"
 # create 'nickname'/variable for the directory where the data is stored

# Parameters of SMR
smr="/QRISdata/Q5149/computation/smr/smr-1.3.1-linux-x86_64/smr"
# smr="./smr_v1.3.1" # Creates variable for location of smr software
eqtl="${data_dir}/Heart_Left_Ventricle/Heart_Left_Ventricle"
  # Uses the data_dir variable to create a pathname variable for the eqtl data
gwas="/QRISdata/Q5149/sumstats/formatted/intersitital_myocardial_fibrosis.ma"

 # Creates a pathname variable for the gwas summstats
smr_result="smr_lv_interstitial_fibrosis" # This is the name we will use for the smr output files
ld_reference="${data_dir}g1000_eur/g1000_eur" # Pathname variable for the ld_reference (needed to estimate ld for SNPs in gwas summstats)

# SMR analysis with HEIDI - uses the eqtls, gwas summstats and ld reference, and writes results to the 'smr_result' variable we created
/QRISdata/Q5149/computation/smr/smr-1.3.1-linux-x86_64/smr --bfile ${ld_reference} --gwas-summary ${gwas} --beqtl-summary ${eqtl} --out ${smr_result} --diff-freq-prop 0.2 > ${smr_result}.log 2>&1

 
