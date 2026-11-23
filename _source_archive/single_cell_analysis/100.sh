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
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/top_100_annot_quant_and_binary.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/result_100/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_100_quant_bin', annot='$annot', log2file=TRUE)"
