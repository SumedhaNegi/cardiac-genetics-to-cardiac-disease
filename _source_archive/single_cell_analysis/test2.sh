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
