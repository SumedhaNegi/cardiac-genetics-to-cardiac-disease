cd /QRISdata/Q5149/sumstats/formatted/from_dalia

ls > filenames.txt

#put filenames in excel to make a config file
#config.txt is a textfile derived from filename.txt used to order my sumstats into $ArrayTaskID where column 1 is the task number column 2 is the sumstat name and optionally you can put column 3 as the phenotype description

#example config.txt
#	ArrayTaskID	sumstat
#	1	Formatted_102_irnt.ma
#	2	Formatted_20002_1075.ma
#	3	Formatted_20002_1094.ma
#	4	Formatted_20002_1473.ma
#	5	Formatted_20003_1140860696.ma

########################################################################
#script
##########################################################################

#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=250GB
#SBATCH --time=01:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general 
#SBATCH --job-name=array_test
#SBATCH --account=a_palpant
#SBATCH --array=1-61

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/
config=config.txt



cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/sbayesRC/ukbEUR_HM3"
annot="1m_xav_chir.txt"
threads=4

export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${ma_file}_tidy.ma', log2file=TRUE)"

Rscript -e "SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$ld_folder', \
                  output='${ma_file}_imp.ma', log2file=TRUE)"


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${ma_file}_sbrc', annot='$annot', log2file=TRUE)"