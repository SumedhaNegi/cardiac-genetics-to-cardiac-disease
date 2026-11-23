#PBS -N sbs
#PBS -A UQ-IMB-CNSG 
#PBS -o stdout 
#PBS -e stderr 
#PBS -l select=1:ncpus=2:mem=50gb -l walltime=25:0:0 

#PBS -S /bin/bash 

cd $PBS_O_WORKDIR


for file in `ls *.tsv`
do
n=`basename $file`
awk -F'\t' -v OFS="\t" '$1~/^rs/{print $1, $5, $6, $7, $11, $12, $16, $17}' $file  > $n_clean.ma

done


module load R 

MRI_svi_bsa <- read.csv("~/90days/HOPX_SBayesS/sumstats/MRI_svi_cleaned.tsv", header=TRUE, sep="\t")



colnames(MRI_svi_bsa) <- c('CHR', 'BP','SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N')



module load R 

#remove the SNPs without rs

library(dplyr)

filtered <- MRI_svi_bsa %>% 
  filter(grepl('rs', SNP))

  write.table(filtered, "~/90days/HOPX_SBayesS/sumstats/svi.ma", sep="\t", quote = FALSE , row.names = FALSE)


################ sophie code 
#PBS -N sbs
#PBS -A UQ-IMB-CNSG 
#PBS -o stdout 
#PBS -e stderr 
#PBS -l select=1:ncpus=2:mem=50gb -l walltime=25:0:0 

#PBS -S /bin/bash 

cd $PBS_O_WORKDIR

module load R 


library(dplyr)

files <- list.files()

for (file in files){
	print(file)
	x <- fread(file, header=TRUE, sep="\t")
	x_1 <- x[, c(2, 3, 1,,5, 6, 7, 11, 12, 16,18)]
	colnames(x_1) <- c("CHR", "BP",'SNP', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P', 'N')
	filtered <- x_1 %>% 
  		filter(grepl('rs', SNP))
  	newfilename <- gsub('.tsv', '.ma', file)
  	print(newfilename)
	fwrite(filtered, newfilename, sep = '\t', quote = F, row.names = F)
}











###########SURM 









##### removing the last column and putting the N in 


library(dplyr)
library(tidyverse)

files <- list.files(pattern = "\\.ma$")
for (file in files){
  print(file)
          x <- read.csv(file, header=TRUE, sep="\t")
     
        x_1 <- x[,-10]
	x_2 <- x_1 %>%
  	add_column(N = "41135")

  
        write.table(x_2, file, sep = "\t", quote = F, row.names = F)





