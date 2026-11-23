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


