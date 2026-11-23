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
