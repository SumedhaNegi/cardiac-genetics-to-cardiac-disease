salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:00:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l



module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)
library(tidyr)
getwd()
# /QRISdata/Q5149/TF_DamID/To_format

#extracting genes from each file since the peak annotation is based on mouse genome

#Create a list to store the first columns (gene_names) from each file

files <- list.files(pattern = "\\.txt$")

for (file in files) {
  data <- fread(file, header = TRUE, sep = "\t")
  first_column <- data[, 1]
  output_file <- paste0("genes_bound_", file)
  write.table(first_column, file = output_file, row.names = FALSE, col.names = TRUE, quote=F)
}

setwd("gene_list/")

#####

annot <- fread("/QRISdata/Q5149/tmem88_heritability/annot/geneMatrix/geneMatrix.tsv")
	annot$GeneID <- annot$gene_name
		head(annot)
		annot2 <- annot[,c("GeneID")]
		head(annot2)
		annot2$chr <- annot$hg19g0
		annot2$start <- annot$g1
		annot2$end <- annot$g2



library(dplyr)

setwd("To_format/gene_list/")

match_and_append <- function(file_name, annot2, adjustment = 50000) {
  # Read the file with gene names
  data <- fread(file_name, header = FALSE, col.names = "GeneID", sep = "\t", stringsAsFactors = FALSE)

  # Convert the "GeneID" columns to lowercase
  data$GeneID <- tolower(data$GeneID)
  annot2$GeneID <- tolower(annot2$GeneID)

  # Merge the gene names with the "annot" data based on the "GeneID" column (case-insensitive)
  merged_data <- data %>%
    left_join(annot2, by = "GeneID")

  # Adjust the start and end positions
  merged_data$start <- merged_data$start + adjustment
  merged_data$end <- merged_data$end - adjustment

    # Remove rows with NA values in any column
  merged_data <- merged_data[complete.cases(merged_data), ]

  # Count the number of matches
  num_matches <- sum(!is.na(merged_data$GeneID))

  # Save the merged data to a new file with the original file name
  new_file_name <- sub(".txt", "_location.txt", file_name)
  write.table(merged_data, new_file_name, sep = "\t", row.names = FALSE, col.names=TRUE, quote=F)

  return(num_matches)
}


# Get a list of all text files in the folder
file_names <- list.files(pattern = "\\.txt$")

# Process each of the text files and count matches
match_counts <- numeric(length(file_names))
for (i in 1:length(file_names)) {
  match_counts[i] <- match_and_append(file_names[i], annot2)
}

# Print the number of matches for each file
for (i in 1:length(file_names)) {
  cat("File:", file_names[i], "- Number of Matches:", match_counts[i], "\n")
}

save.image("updated.RData")


#Origial 																			
#File: genes_bound_HOPX.txt - Number of Matches: 16754 					
#File: genes_bound_mef2a.txt - Number of Matches: 1338 
#File: genes_bound_nkx25.txt - Number of Matches: 20574 
#File: genes_bound_SRF.txt - Number of Matches: 23807 
#File: genes_bound_tbx5.txt - Number of Matches: 55873 
#File: HOPX_bound_genes.txt - Number of Matches: 3541

#nNumber that matched with the human genes
#File: genes_bound_HOPX.txt - Number of Matches: 10562 
#File: genes_bound_mef2a.txt - Number of Matches: 868 
#File: genes_bound_nkx25.txt - Number of Matches: 12862 
#File: genes_bound_SRF.txt - Number of Matches: 15441 
#File: genes_bound_tbx5.txt - Number of Matches: 36639 
#File: HOPX_bound_genes.txt - Number of Matches: 3459 



#linux code because the first two lies are duplicated, gets rid of the second line for each file in a directory.
find . -type f -exec sed -i '2d' {} \;


#########################################################################################################

#annotating SNPs to the gene list 



snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


	dim(snp)
	#[1] 7356518      10
	snp2 <- snp[,c(1,2,5,6,7)]
	colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")



library(dplyr)
library(purrr)

getwd()
#[1] "/QRISdata/Q5149/TF_DamID/To_format/gene_list/locations"


####################################

#old way 

tbx3 <-fread ("")
tbx3$chr<- sub("^chr", "", tbx3$chr)


####################

	snp2$gata4 <- ifelse(sapply(seq_along(snp2$PhysPos), function(i) {
           inds <- gata4$start <= snp2$PhysPos[i] & gata4$end >= snp2$PhysPos[i]
           any(inds) & (snp2$Chromosome[i] == gata4$chr[which.max(inds)])
           }), "1", "0")

	table(snp2$gata4)


write.table(snp2, "/QRISdata/Q5149/TF_DamID/annotation/tbx3_snps_based_on_chipseq_location_boogaard.txt", sep="\t", quote = FALSE , col.names = TRUE, row.names = FALSE)

save.image("update.RData")


#Number of Snps

			$hopx

			      0       1 
			7058115  298403 

			$HOPX

			      0       1 
			7034353  322165 

			$srf

			      0       1 
			7025894  330624 

			$nkx25

			      0       1 
			7019671  336847 

			$tbx5

			      0       1 
			6960530  395988 

			$mef2a

			      0       1 
			7260039   96479 


			$tbx3 boogaard	
      			  0       1 
			6900022  456496 



# QC how many snps overlap ?????????/


#################### 





#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=300GB
#SBATCH --time=04:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=debug
#SBATCH --job-name=cardiac_nick
#SBATCH --account=a_palpant
#SBATCH --array=1-29

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/sbayesRC/ukbEUR_HM3"
annot=/QRISdata/Q5149/TF_DamID/annotation/annotation1.txt
threads=4
 
export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${ma_file}_tidy.ma', log2file=TRUE)"

Rscript -e "SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$ld_folder', \
                  output='${ma_file}_imp.ma', log2file=TRUE)"


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${ma_file}_sbrc_tfdamid', annot='$annot', log2file=TRUE)"

############## QPX too small with 7m snps, going to try with 1million snps and see 


#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=200GB
#SBATCH --time=03:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=debug
#SBATCH --job-name=cardiac_nick
#SBATCH --account=a_palpant
#SBATCH --array=1-29

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/sbayesRC/ukbEUR_HM3"
annot=/QRISdata/Q5149/TF_DamID/annotation/annotation1.txt
threads=4
 
export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='$ma_file', LDdir='$ld_folder', \
                  outPrefix='${ma_file}_sbrc_tfdamid', annot='$annot', log2file=TRUE)"


Annotation column Annotation column Annotation column 43 has too small XPX: 0, please remove this column from annotation data or perform re-scale.
 has too small XPX: 2 has too small XPX: 0, please remove this column from annotation data or perform re-scale.
0, please remove this column from annotation data or perform re-scale.
Annotation column 5 has too small XPX: 0, please remove this column from annotation data or perform re-scale.



/var/spool/slurmd/job6222476/slurm_script: line 31: 262115 Aborted                 (core dumped) Rscript -e "SBayesRC::sbayesrc(mafile='$ma_file', LDdir='$ld_folder',            
       outPrefix='${ma_file}_sbrc_tfdamid', annot='$annot', log2file=TRUE)"














       ########### formatting with 1m snps 

snp <- fread("/QRISdata/Q5149/sbayesRC/ukbEUR_HM3/snp.info")


	dim(snp)
	#[1] 7356518      10
	snp2 <- snp[,c(1,2,5,6,7)]
	colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")



mef2a <-fread ("")
srf$chr<- sub("^chr", "", srf$chr)


####################

	snp2$srf <- ifelse(sapply(seq_along(snp2$PhysPos), function(i) {
           inds <- srf$start <= snp2$PhysPos[i] & srf$end >= snp2$PhysPos[i]
           any(inds) & (snp2$Chromosome[i] == srf$chr[which.max(inds)])
           }), "1", "0")

	table(snp2$srf)

#mef2a   
0       1 
1137907   16,615 

#hopx3

      0       1 
1100473   54,049 

#$tbx5 
      0       1 
1085802   68,720 

#HOPX
 0       1 
1097437   57,085 

#nkx25
   0       1 
1096058   58,464 


SRF
#     0       1 
1096579   57,943


#excluding mef2a 




#######################

#find overlapping genes with a location 

library(GenomicRanges)

# Remove rows with missing values in annot2
annot2 <- annot2[complete.cases(annot2$start, annot2$end), ]

# Remove rows with missing values in locList
locList <- locList[complete.cases(locList$start, locList$end), ]

# Ensure consistent chromosome names
common_chromosomes <- intersect(annot2$chr, locList$chr)

# Filter data frames to include only common chromosomes
annot2 <- annot2[annot2$chr %in% common_chromosomes, ]
locList <- locList[locList$chr %in% common_chromosomes, ]

# Create GRanges objects for annot2 and locList
annot2_gr <- GRanges(
  seqnames = annot2$chr,
  ranges = IRanges(start = annot2$start, end = annot2$end),
  strand = "*"
)

locList_gr <- GRanges(
  seqnames = locList$chr,
  ranges = IRanges(start = locList$start, end = locList$end),
  strand = "*"
)


#QC



################################

salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=50GB --job-name=TEST --time=01:30:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l



module load r/4.2.1-foss-2021a 

		R
		library(data.table)
		library(dplyr)

		getwd()
		#[1] "/QRISdata/Q5149/TF_DamID/To_format/gene_list/locations/location"


		library(eulerr)


		# List all text files in the current working directory with a .txt extension
		file_paths <- list.files(pattern = "\\location.txt$")

		# Initialize an empty list to store vectors
		gene_id_list <- vector("list", length = length(file_paths))

		# Loop through each text file and read the GeneID names into a vector
		for (i in seq_along(file_paths)) {
		  gene_id_list[[i]] <- read.table(file_paths[i], header = TRUE, stringsAsFactors = FALSE)$GeneID
		}


		str(gene_id_list)
		#	List of 6
			 $ : chr [1:10562] "lypla1" "tcea1" "atp6v1h" "sntg1" ...
			 $ : chr [1:868] "tram2" "prim2" "dst" "dst" ...
			 $ : chr [1:12862] "lypla1" "tcea1" "atp6v1h" "pcmtd1" ...
			 $ : chr [1:15441] "lypla1" "tcea1" "tcea1" "atp6v1h" ...
			 $ : chr [1:36639] "xkr4" "sox17" "sox17" "mrpl15" ...
			 $ : chr [1:3459] "a3galt2" "aacs" "aagab" "aar2" ...


		 names(gene_id_list) <- c("HOPX","mef2a","nkx25","srf","tbx5","hopx")


		upset_plot <- upset(fromList(gene_id_list),
		      nintersects = NA, 
		      nsets = 6,
		      sets = c("HOPX","mef2a","nkx25","srf","tbx5","hopx"),
		      empty.intersections = "on",
		      order.by = "freq", 
		      decreasing = T, 
		      number.angles = 0, 
		      text.scale = 1.1, 
		      point.size = 2.8, 
		      line.size = 1)



# Use lapply to apply unique to each vector
gene_id_list_unique <- lapply(gene_id_list, unique)

# Convert the list to a matrix
gene_id_matrix <- matrix(0, nrow = length(gene_id_list_unique[[1]]), ncol = length(gene_id_list_unique))
colnames(gene_id_matrix) <- names(gene_id_list_unique)

# Fill in the matrix with 1 where the gene is present in a set
for (i in seq_along(gene_id_list_unique)) {
  gene_id_matrix[, i] <- as.integer(gene_id_list_unique[[i]] %in% gene_id_list_unique[[i]])
}

# Use the euler package to calculate intersections
gene_id_intersections <- euler(gene_id_matrix)

# Extract gene IDs for each intersection
intersection_gene_lists <- lapply(gene_id_intersections, function(intersection) {
  unique(unlist(gene_id_list_unique[intersection]))
})


##
#generating an upset plot from the SNPs overlapping 

# Load the required library
library(UpSetR)
getwd()
#[1] "/QRISdata/Q5149/TF_DamID/annotation"

data <- fread("annot.txt")
# Convert the data frame to a binary matrix
binary_matrix <- as.matrix(data[, c(-1,-2)])
binary_data <- as.data.frame(binary_matrix)
binary_data <- mutate_all(binary_data, as.numeric)

# Create the upset plot
upset_plot <- upset(binary_data, main.bar.color = "skyblue", order.by = "freq")

pdf("upset_plot.pdf", width = 10, height = 6)
print(upset_plot)
dev.off()



##############

#Doing DamID heritability without the 50kb window

match_and_append <- function(file_name, annot2, adjustment = 0) {
  # Read the file with gene names
  data <- fread(file_name, header = FALSE, col.names = "GeneID", sep = "\t", stringsAsFactors = FALSE)

  # Convert the "GeneID" columns to lowercase
  data$GeneID <- tolower(data$GeneID)
  annot2$GeneID <- tolower(annot2$GeneID)

  # Merge the gene names with the "annot" data based on the "GeneID" column (case-insensitive)
  merged_data <- data %>%
    left_join(annot2, by = "GeneID")

  # Adjust the start and end positions
  merged_data$start <- merged_data$start + adjustment
  merged_data$end <- merged_data$end - adjustment

    # Remove rows with NA values in any column
  merged_data <- merged_data[complete.cases(merged_data), ]

  # Count the number of matches
  num_matches <- sum(!is.na(merged_data$GeneID))

  # Save the merged data to a new file with the original file name
  new_file_name <- sub(".txt", "nowindow_location.txt", file_name)
  write.table(merged_data, new_file_name, sep = "\t", row.names = FALSE, col.names=TRUE, quote=F)

  return(num_matches)
}



# Get a list of all text files in the folder
file_names <- list.files(pattern = "\\.txt$")

# Process each of the text files and count matches
match_counts <- numeric(length(file_names))
for (i in 1:length(file_names)) {
  match_counts[i] <- match_and_append(file_names[i], annot2)
}

# Print the number of matches for each file
for (i in 1:length(file_names)) {
  cat("File:", file_names[i], "- Number of Matches:", match_counts[i], "\n")
}


#qc


		#locations

		File: genes_bound_HOPX.txt - Number of Matches: 10562 
		File: genes_bound_mef2a.txt - Number of Matches: 868 
		File: genes_bound_nkx25.txt - Number of Matches: 12862 
		File: genes_bound_SRF.txt - Number of Matches: 15441 
		File: genes_bound_tbx5.txt - Number of Matches: 36639 
		File: HOPX_bound_genes.txt - Number of Matches: 3459 


# annotating snps to genes without window - using 7m 
 snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


        dim(snp)
        #[1] 7356518      10
        snp2 <- snp[,c(1,2,5,6,7)]
        colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")




tbx5 <-fread ("genes_bound_tbx5nowindow_location.txt")
tbx5$chr<- sub("^chr", "", tbx5$chr)



	snp2$tbx5 <- ifelse(sapply(seq_along(snp2$PhysPos), function(i) {
           inds <- tbx5$start <= snp2$PhysPos[i] & tbx5$end >= snp2$PhysPos[i]
           any(inds) & (snp2$Chromosome[i] == tbx5$chr[which.max(inds)])
           }), "1", "0")

	table(snp2$tbx5)


	#nkx25
	      0       1 
		6876618  479900 

	#tbx5
	      0       1 
		6842356  514162

	 #mef2a
	     0       1 
	 7187109  169409

	 #hopx
	   0       1 
	  6927523  428995
	 
	 #srf
	     0       1 
	6881900  474618 
 
 		#gata4 
	       0       1 
		6889106  467412
save.image("update2.RData")




#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=80GB
#SBATCH --time=03:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general 
#SBATCH --job-name=tf_card
#SBATCH --account=a_palpant
#SBATCH --array=1-61

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/tf_damid_7m
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/scratch/user/uqsnegi/ukbEUR_Imputed"
annot="/QRISdata/Q5149/TF_DamID/annotation/7m_nowindow.txt"
threads=4
 
export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$ld_folder', \
				 outPrefix='${ma_file}_7m_nowindow.txt', annot='$annot', log2file=TRUE)"






##### Plotting 

file_names <- list.files(pattern = ".AnnoPerSnpHsqEnrichment", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Function to extract parameter from file name
extract_parameter <- function(file_name) {
  sub("^\\./(.+?)\\.ma.*$", "\\1", file_name)
}

# Read and combine the data from each file
for (file in file_names) {
  df <- read.table(file, header = TRUE, sep = "\t")  # Adjust sep based on your actual delimiter
  df$Parameter <- extract_parameter(file)
  combined_data <- rbind(combined_data, df)
}

#replace code names with ukb names 
parameter <- fread ("/QRISdata/Q5149/sumstats/formatted/from_dalia/ukb31063_h2_all.02Oct2019.tsv")
parameter <- parameter[,c(1,8)]
colnames(parameter) <- c("Parameter", "description")
combined_data$Parameter <- gsub ("Formatted_", "",combined_data$Parameter)

  x <- combined_data %>%
  mutate(Parameter = trimws(Parameter)) %>%
  left_join(parameter %>% mutate(Parameter = trimws(Parameter)), by = "Parameter") %>%
  mutate(Full_Name = ifelse(is.na(description), Parameter, description)) %>%
  select(-Parameter, -description) 



combined_data <- x

# Create a bar chart using ggplot2
desired_order <- c("hopx", "gata4", "srf", "nkx25", "tbx5")

mean_values <- combined_data %>%
  group_by(Annotation) %>%
  summarise(mean_enrich = mean(Enrich))



# Convert the Annotation column to a factor with the desired order
combined_data$Annotation <- factor(combined_data$Annotation, levels = desired_order)


ggplot(combined_data, aes(x = Parameter, y = Enrich, fill = Annotation)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Enrich - sd(SD), ymax = Enrich + sd(SD)),
      position = position_dodge(width = 0.8),  # Adjust width to match geom_bar
    width = 0.25
  ) +
  labs(
    title = "Per Snp Enrichment",
    x = "Full_Name",
    y = "Enrichment Value"
  ) +
  scale_fill_brewer(palette = "Set3") +  # You can choose a different palette
   geom_hline(data = mean_data, aes(yintercept = Enrich, color = Annotation),
             linetype = "dashed", size = 1) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))




#without window for cardiac traits 

#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=80GB
#SBATCH --time=01:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general 
#SBATCH --job-name=tf_mri
#SBATCH --account=a_palpant
#SBATCH --array=1-19
MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/test
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/scratch/user/uqsnegi/ukbEUR_Imputed"
annot="/QRISdata/Q5149/TF_DamID/annotation/7m_nowindow.txt"
threads=4
 
export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='$ma_file', LDdir='$ld_folder', \
                  outPrefix='${ma_file}_nowindow', annot='$annot', log2file=TRUE)"






############
#Unique SNPS
# Assuming your dataframe is named 'df'
module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)
library(tidyr)
getwd()

setwd("/QRISdata/Q5149/TF_DamID/annotation")

x <- fread ("7m_nowindow.txt")
x <- x[,-1]
x <- x[,-1]

 dt <- x
 dt_processed <- copy(dt)  # Create a copy of the original data.table

# Check if the data.table has rows
if (nrow(dt) > 0) {
  # Iterate over each row
  for (i in 1:nrow(dt)) {
    # Check if the row has more than one column with the value 1
    if (sum(dt[i, .SD, .SDcols = -1] == 1) > 1) {
      # If yes, set all values in that row to 0
      dt_processed[i, names(dt)[-1] := 0]
    }
  }
}















############################


visualization of no window results 

# using facet wrap to get the results 



setwd("/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/test/results")


# Assuming your files are named in a specific way, adjust the pattern accordingly
file_names <- list.files(pattern = ".AnnoPerSnpHsqEnrichment", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Function to extract parameter from file name
extract_parameter <- function(file_name) {
  sub("^\\./(.+?)\\.ma.*$", "\\1", file_name)
}

# Read and combine the data from each file
for (file in file_names) {
  df <- read.table(file, header = TRUE, sep = "\t")  # Adjust sep based on your actual delimiter
  df$Parameter <- extract_parameter(file)
  combined_data <- rbind(combined_data, df)

}



desired_order <- c("HOPX", "Gata4", "SRF", "nkx25", "tbx5", "MEF2A")


combined_data$Annotation <- factor(
  combined_data$Annotation,
  levels = desired_order
)

combined_data <- na.omit(combined_data)

table(combined_data$Annotation)

     HOPX     Gata4       SRF     nkx25      tbx5     MEF2A Intercept 
       54        54        54        54        54        36        18 

# Create a bar chart using ggplot2
mean_values <- aggregate(Enrich ~ Annotation, data = combined_data, mean, na.rm = TRUE)


ggplot(combined_data, aes(x = Parameter, y = Enrich, fill = Annotation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = Enrich - SD, ymax = Enrich + SD),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  labs(
    title = "Per Snp Enrichment of heritability",
    x = "Parameter",
    y = "Enrichment Value"
  ) +
  scale_fill_manual(
    values = c("HOPX" = "black", "Gata4" = "lightcoral", "SRF" = "lightblue", "nkx25" = "cornflowerblue", "tbx5" = "darkgreen", "MEF2A" = "orange"),
    breaks = desired_order
  ) +
  scale_color_manual(
    values = c("HOPX" = "black", "Gata4" = "lightcoral", "SRF" = "lightblue", "nkx25" = "cornflowerblue", "tbx5" = "darkgreen", "MEF2A" = "orange"),
    breaks = desired_order
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
  guides(fill = guide_legend(title = "Annotation")) 



ggsave("Full_barplot_TF_binding.pdf", width = 30, height = 25)






###########################################################################

# dalia GWAS sumstats 1million 
/QRISdata/Q5149/sumstats/formatted/from_dalia/1m_imputed_sumstats

# file_list
/QRISdata/Q5149/sumstats/formatted/from_dalia/ukb31063_h2_all.02Oct2019.tsv

# non imputed gwas sumstats rv lv 
[uqsnegi@bunya2 tf_damid_7m_nowindow]$ pwd
/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow


# imputed 7m 
/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/imputed_sumstats





#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=150GB
#SBATCH --time=3:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=young_old
#SBATCH --job-name=general
#SBATCH --account=a_palpant
#SBATCH --array=1-61

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/1m_imputed_sumstats
config=config.txt

cd $MY_DIR


# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/ukbEUR_HM3/"
annot="/QRISdata/Q5149/TF_DamID/M2/annot/fetal_vs_adult_human_vs_common.txt"
out_prefix="/QRISdata/Q5149/TF_DamID/M2/annot/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_fetal_adult_common', annot='$annot', log2file=TRUE)"



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
annot="/QRISdata/Q5149/TF_DamID/M2/annot/fetal_vs_adult_human_vs_common.txt"
out_prefix="/QRISdata/Q5149/TF_DamID/M2/annot/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_fetal_adult_common', annot='$annot', log2file=TRUE)"









# venn diagram schematic ?


# Install and load the VennDiagram package if you haven't already
if (!require("VennDiagram")) {
  install.packages("VennDiagram")
}



# Load library
library(VennDiagram)
 Akerberg_et_al <- c("GATA4", "NKX25", "MEF2A", "MEF2C", "SRF", "TBX5", "TEAD1")
He_et_al <- c("GATA4", "MEF2A", "NKX25", "SRF", "TBX5", "P300")
Schlesinger_et_al <- c("GATA4", "MEF2A", "NKX25", "SRF")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "PuBu")# Chart
venn.diagram(
  x = list(
    "Akerberg et al" = Akerberg_et_al,
    "He et al" = He_et_al,
    "Schlesinger et al" = Schlesinger_et_al
  ),        
  category.names = c("Arkerberg et al" , "He et al " , "Schlesinger et al"),
   
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)


venn <- venn.diagram(
  x = list(
    "Akerberg et al" = Akerberg_et_al,
    "He et al" = He_et_al,
    "Schlesinger et al" = Schlesinger_et_al
  ), filename=Venn,   
  category.names = c("Akerberg et al", "He et al", "Schlesinger et al"),
        
  # Circles
  fill = myCol
)










############

/QRISdata/Q5149/TF_DamID/M2/annot


#general fetal vs adult for TFs de genes 

# Load required library
library(biomaRt)

# Define function to convert mouse gene names to human gene names
convertMouseGeneList <- function(musGenes) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "dec2021.archive.ensembl.org")

  # Query Ensembl using the main server
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = musGenes, mart = mouse, attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows = TRUE)

  # Initialize a vector of NA
  human_symbols <- rep(NA, length(musGenes))

  # Identify indices where conversion is successful
  valid_indices <- nrow(genesV2) > 0

  # Assign human gene symbols to the valid indices
  human_symbols[valid_indices] <- unique(genesV2[valid_indices, 2])

  # Return the result
  return(human_symbols)
}

# List of txt files
files <- list.files(pattern = "\\.txt$")

# Initialize an empty list to store results
results_list <- list()

# Iterate over each file
for (file in files) {
  # Read the file
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)

  # Extract mouse gene names from the first column
  musGenes <- data[[1]]

  # Convert mouse gene names to human gene names
  human_symbols <- convertMouseGeneList(musGenes)

  # Get the base file name
  base_name <- tools::file_path_sans_ext(basename(file))

  # Assign the result to the list with base file name as vector name
  results_list[[base_name]] <- human_symbols
}



str(results_list)
List of 4
 $ cm_downregulated_in_fetal_ms: chr [1:418] "ACSS2" "KCNB1" "JAM2" "SIX4" ...
 $ cm_foetal_vs_adult_sig_ms   : chr [1:652] "ACSS2" "SIX4" "FLRT3" "IL13RA1" ...
 $ cm_upregulated_in_fetal_ms  : chr [1:233] "GINS1" "RPL3" "PGRMC1" "KLHL13" ...
 $ nonsig_ms                   : chr [1:416] "GPHB5" "ICOSLG" "TMEM33" "C21orf62" ...


results_list <- lapply(results_list, unique)
str(results_list)
  # Store the result in the merged_results_list

annot <- fread("/QRISdata/Q5149/tmem88_heritability/annot/geneMatrix/geneMatrix.tsv")
	annot$GeneID <- annot$gene_name
		head(annot)
		annot2 <- annot[,c("GeneID")]
		head(annot2)
		annot2$chr <- annot$hg19g0
		annot2$start <- annot$g1
		annot2$end <- annot$g2




# Create a new list to store the adjusted start, end, and chromosome values




 # Iterate over each file in results_list
for (file_name in names(results_list)) {
  # Get gene IDs from results_list
  gene_ids <- results_list[[file_name]]
  
  # Find matching rows in annot2
  matching_rows <- annot2[GeneID %in% gene_ids]
  
  # Append annot2 columns to results_list
  if (nrow(matching_rows) > 0) {
    # Find the unique matching rows
    unique_matching_rows <- matching_rows[!duplicated(matching_rows$GeneID), ]
    
    # Get the indices of matching gene IDs in results_list
    indices <- match(unique_matching_rows$GeneID, gene_ids)
    
    # Create or update chr, start, and end columns in results_list
    results_list[[file_name]] <- data.frame(
      GeneID = gene_ids,
      chr = NA_character_,
      start = NA_integer_,
      end = NA_integer_
    )
    results_list[[file_name]]$chr[indices] <- unique_matching_rows$chr
    results_list[[file_name]]$start[indices] <- unique_matching_rows$start
    results_list[[file_name]]$end[indices] <- unique_matching_rows$end
  }
}




























List of 4
 $ cm_downregulated_in_fetal_ms:List of 2
  ..$ merged_data:'data.frame':	397 obs. of  5 variables:
  .. ..$ GeneID  : chr [1:397] "glrx5" "corin" "phyh" "fundc1" ...
  .. ..$ GeneID.y: chr [1:397] "GLRX5" "CORIN" "PHYH" "FUNDC1" ...
  .. ..$ chr     : chr [1:397] "chr14" "chr4" "chr10" "chrX" ...
  .. ..$ start   : num [1:397] 96049840 47646016 13369796 44432885 34907837 ...
  .. ..$ end     : num [1:397] 95961061 47790123 13294412 44352105 34814023 ...
  ..$ num_matches: int 397
 $ cm_foetal_vs_adult_sig_ms   :List of 2
  ..$ merged_data:'data.frame':	609 obs. of  5 variables:
  .. ..$ GeneID  : chr [1:609] "dnajc28" "mex3c" "cenpc" "glrx5" ...
  .. ..$ GeneID.y: chr [1:609] "DNAJC28" "MEX3C" "CENPC" "GLRX5" ...
  .. ..$ chr     : chr [1:609] "chr21" "chr18" "chr4" "chr14" ...
  .. ..$ start   : num [1:609] 34907837 48750920 68384480 96049840 47959931 ...
  .. ..$ end     : num [1:609] 34814023 48694703 68361221 95961061 48050834 ...
  ..$ num_matches: int 609
 $ cm_upregulated_in_fetal_ms  :List of 2
  ..$ merged_data:'data.frame':	211 obs. of  5 variables:
  .. ..$ GeneID  : chr [1:211] "tceal8" "gins1" "rpl3" "pgrmc1" ...
  .. ..$ GeneID.y: chr [1:211] "TCEAL8" "GINS1" "RPL3" "PGRMC1" ...
  .. ..$ chr     : chr [1:211] "chrX" "chr20" "chr22" "chrX" ...
  .. ..$ start   : num [1:211] 1.03e+08 2.54e+07 3.98e+07 1.18e+08 3.41e+07 ...
  .. ..$ end     : num [1:211] 1.02e+08 2.54e+07 3.97e+07 1.18e+08 3.41e+07 ...
  ..$ num_matches: int 211
 $ nonsig_ms                   :List of 2
  ..$ merged_data:'data.frame':	334 obs. of  5 variables:
  .. ..$ GeneID  : chr [1:334] "gphb5" "icoslg" "tmem33" "c21orf62" ...
  .. ..$ GeneID.y: chr [1:334] "GPHB5" "ICOSLG" "TMEM33" "C21orf62" ...
  .. ..$ chr     : chr [1:334] "chr14" "chr21" "chr4" "chr21" ...
  .. ..$ start   : num [1:334] 63829549 45686897 41987146 34212984 2009403 ...
  .. ..$ end     : num [1:334] 63735649 45610849 41912820 34136053 1924931 ...
  ..$ num_matches: int 334

> str(result_list) #mouse 
List of 13
 $ Adult_Gata4 : chr [1:19712] "Secisbp2l" "4933406G16Rik" "Srd5a3" "Agap2" ...
 $ Adult_Mef2a : chr [1:11103] "Btrc" "Mir3091" "4930594O21Rik" "9230112D13Rik" ...
 $ Adult_Nkx2-5: chr [1:20600] "Arhgef1" "Erbb4" "Asap2" "Gadd45b" ...
 $ Adult_Srf   : chr [1:12485] "Dlc1" "Myh11" "Srsf5" "Trio" ...
 $ Adult_Tbx5  : chr [1:6541] "Sgms1" "Asph" "Ift80" "1110058D11Rik" ...
 $ Adult_Tead  : chr [1:25442] "Gphn" "Noct" "Repin1" "Misp" ...
 $ Fetal_Gata4 : chr [1:44184] "Fbxo33" "Dok5" "Grb10" "Igsf6" ...
 $ Fetal_Mef2a : chr [1:15933] "Nedd4" "Napg" "Il20ra" "Galnt5" ...
 $ Fetal_Mef2c : chr [1:30178] "Arhgef16" "Kirrel" "Angptl2" "Col16a1" ...
 $ Fetal_Nkx2-5: chr [1:25864] "Gprin1" "Sobp" "Wdr7" "Celf2" ...
 $ Fetal_Srf   : chr [1:18521] "Prr30" "Cic" "Napb" "Scube3" ...
 $ Fetal_Tbx5  : chr [1:30789] "Prkag2" "Fabp3" "Gm9054" "Anp32b" ...
 $ Fetal_Tead  : chr [1:45806] "Mir599" "Itpr1" "4933411G06Rik" "Gpc4" ...

saveRDS(result_list,"list_of_fetal_vs_adult_binding_genes_TFs_humans.RDS")








   snp <- fread("/QRISdata/Q5149/ukbEUR_Imputedsnp.info")

                   
        dim(snp)
        #[1] 7356518      10
        snp2 <- snp[,c(1,2,5,6,7)]
        colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")

# Iterate over each vector in results_list
for (vector_name in names(results_list)) {
  # Filter out rows where the chr value is "chrX" and where any value is NA
  results_list[[vector_name]] <- results_list[[vector_name]][complete.cases(results_list[[vector_name]]), ]
}

# Iterate over each vector in results_list
for (vector_name in names(results_list)) {
  # Filter out rows where the chr value is "chrX"
  results_list[[vector_name]] <- results_list[[vector_name]][results_list[[vector_name]]$chr != "chrX", ]
}
for (vector_name in names(results_list)) {
  # Filter out rows where the chr value is "chrX"
  results_list[[vector_name]]$chr <- results_list[[vector_name]][results_list[[vector_name]]$chr != "chrX", ]
}

  results_list[[vector_name]]$chr <- gsub("^chr", "", results_list[[vector_name]]$chr)


        # Load required libraries
library(GenomicRanges)
library(data.table)

# Convert 'Chromosome' column in snp2 to character format
snp2$Chromosome <- as.character(snp2$Chromosome)

# Create GenomicRanges object for snp2
snp2_ranges <- GRanges(seqnames = snp2$Chromosome, ranges = IRanges(start = snp2$PhysPos, end = snp2$PhysPos))

# Iterate over each file in merged_results_list
for (file_name in names(merged_results_list)) {
  # Extract merged data from merged_results_list
  merged_data <- merged_results_list[[file_name]]$merged_data
  
  # Convert 'seqnames' column in merged data to character format
  merged_data$seqnames <- as.character(merged_data$seqnames)
  
  # Create GenomicRanges object for the merged data
  merged_ranges <- GRanges(seqnames = merged_data$seqnames, ranges = IRanges(start = merged_data$start, end = merged_data$end))
  
  # Create column name prefix for overlap columns
  overlap_col_prefix <- gsub("_", "", file_name)
  
  # Check for overlaps and create new columns in snp2
  overlap_col <- paste0("overlap_", overlap_col_prefix)
  
  # Remove rows with NAs in 'start' or 'end' columns
  merged_data <- na.omit(merged_data)
  
  # Create GenomicRanges object for the current merged data
  snp2[[overlap_col]] <- as.integer(overlapsAny(snp2_ranges, merged_ranges))
}

# List of overlap columns
overlap_columns <- grep("^overlap_", names(snp2), value = TRUE)

# Create a data.table to store the results
result_dt <- data.table(column_name = character(), table_result = list())

# Iterate over each overlap column
for (col in overlap_columns) {
  # Apply the table function to the column
  table_result <- table(snp2[[col]], useNA = "ifany")
  
  # Append results to result_dt
  result_dt <- rbind(result_dt, data.table(column_name = col, table_result = list(table_result)))
}

# Print result_dt


# if error remove 

> result_dt
                        column_name table_result
                             <char>       <list>
1: overlap_cmdownregulatedinfetalms 7254710,....
2:     overlap_cmfoetalvsadultsigms 7215204,....
3:   overlap_cmupregulatedinfetalms 7317094,....
4:                 overlap_nonsigms 7311521,....




print(result_dt)







            column_name table_result
                 <char>       <list>
 1:  overlap_AdultGata4 5202962,....
 2:  overlap_AdultMef2a 5878358,....
 3: overlap_AdultNkx2-5 5500617,....
 4:    overlap_AdultSrf 5815965,....
 5:   overlap_AdultTbx5 6331555,....
 6:   overlap_AdultTead 5261705,....
 7:  overlap_FetalGata4 4869803,....
 8:  overlap_FetalMef2a 5423876,....
 9:  overlap_FetalMef2c 5056457,....
10: overlap_FetalNkx2-5 5165097,....
11:    overlap_FetalSrf 5279177,....
12:   overlap_FetalTbx5 4980593,....
13:   overlap_FetalTead 4795143,....









salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)




R
#gene ontology of overlap

"/QRISdata/Q5149/TF_DamID/To_format/gene_list"



library(data.table)

file_names <- list.files()
# Split the input string into a vector of filenames

file_names <- unlist(strsplit(file_names, "\\s+"))

filtered_files <- grep("^genes_bound(?!.*nowindow_location)", file_names, value = TRUE, perl = TRUE)

filtered_files
	[1] "genes_bound_gata4.txt"  "genes_bound_gata4.txt~" "genes_bound_mef2a.txt" 
	[4] "genes_bound_nkx25.txt"  "genes_bound_SRF.txt"    "genes_bound_tbx5.txt"  




# Initialize an empty list to store vectors
gene_vectors <- list()

# Read the first column of each file and store in the list
for (file in filtered_files) {
  # Check if the file exists
  if (file.exists(file)) {
    # Read the first column
    data <- fread(file, select = 1, header = FALSE)
    # Store the column as a vector in the list
    gene_vectors[[file]] <- data[[1]]
  } else {
    warning(paste("File not found:", file))
  }
}

# Print the list of vectors
print(gene_vectors)

> str(gene_vectors)
List of 7
 $ genes_bound_gata4.txt : chr [1:16753] "AK149000" "Lypla1" "Tcea1" "Atp6v1h" ...
 $ genes_bound_gata4.txt~: chr [1:16753] "AK149000" "Lypla1" "Tcea1" "Atp6v1h" ...
 $ genes_bound_mef2a.txt : chr [1:1337] "Tram2" "Prim2" "1700001G17Rik" "Dst" ...
 $ genes_bound_nkx25.txt : chr [1:20573] "mKIAA1889" "mKIAA1889" "AK149000" "AK149000" ...
 $ genes_bound_SRF.txt   : chr [1:23806] "AK149000" "Lypla1" "Tcea1" "Tcea1" ...
 $ genes_bound_tbx5.txt  : chr [1:55872] "Xkr4" "Sox17" "Sox17" "Mrpl15" ...
 $ genes_bound_hopx.txt  : chr [1:3459] "a3galt2" "aacs" "aagab" "aar2" ...


 vector_to_remove <- "genes_bound_gata4.txt~"
if (vector_to_remove %in% names(gene_vectors)) {
  gene_vectors[[vector_to_remove]] <- NULL
}







# Load required library
library(biomaRt)

# Define function to convert mouse gene names to human gene names
convertMouseGeneList <- function(musGenes) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "dec2021.archive.ensembl.org")

  # Query Ensembl using the main server
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = musGenes, mart = mouse, attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows = TRUE)

  # Initialize a vector of NA
  human_symbols <- rep(NA, length(musGenes))

  # Identify indices where conversion is successful
  valid_indices <- nrow(genesV2) > 0

  # Assign human gene symbols to the valid indices
  human_symbols[valid_indices] <- unique(genesV2[valid_indices, 2])

  # Return the result
  return(human_symbols)
}

gene_vectors <- lapply(gene_vectors, unique)
x

for (file in names(gene_vectors)) {
  gene_vectors[[file]] <- convertMouseGeneList(gene_vectors[[file]])
}




> str(gene_vectors)
List of 6
 $ genes_bound_gata4.txt: chr [1:4543] "MAFIP" "RELL1" "UCHL1" "RRP1" ...
 $ genes_bound_mef2a.txt: chr [1:730] "ZCCHC2" "AKAP6" "NCOA3" "PTPN21" ...
 $ genes_bound_nkx25.txt: chr [1:5026] "MAFIP" "PSPC1" "SLC37A1" "CDK5RAP1" ...
 $ genes_bound_SRF.txt  : chr [1:6914] "MAFIP" "NRIP1" "SEZ6L" "MCAT" ...
 $ genes_bound_tbx5.txt : chr [1:9743] "MAFIP" "DAZ4" "DAZ3" "DAZ1" ...
 $ genes_bound_hopx.txt : chr [1:3225] "JAKMIP1" "TMEM184B" "RGS6" "CNMD" ...




		upset_plot <- upset(fromList(gene_vectors),
		      sets = names(gene_vectors),
		      empty.intersections = "on",
		      order.by = "freq", 
		      decreasing = T, 
		      number.angles = 0, 
		      text.scale = 1.1, 
		      point.size = 2.8, 
		      line.size = 1)





set_names <- c("hopx", "mef2a", "gata4", "nkx25", "SRF", "tbx5")

# Initialize an empty list to store the row names for each set
set_row_names <- list()

# Iterate over each set name
for (set_name in set_names) {
  # Get the intersection members for the current set
  z <- get_intersect_members(x, paste0("genes_bound_", set_name, ".txt"))
  
  # Extract row names for rows with at least one value of 1
  row_names <- rownames(z)[rowSums(z) > 0]
  
  # Store the row names in the results list
  set_row_names[[set_name]] <- row_names
}



 $ hopx : chr [1:833] "CNMD" "N6AMT1" "ELFN2" "TFAP2C" ...
 $ mef2a: chr "TYRP1"
 $ gata4: chr [1:125] "SNRPB2" "LRRC9" "CUBN" "DEFB118" ...
 $ nkx25: chr [1:210] "TRDV3" "GABRB1" "PIGN" "NKX3-2" ...
 $ SRF  : chr [1:203] "GPR33" "RBBP9" "TGDS" "UFM1" ...
 $ tbx5 : chr [1:1636] "DAZ4" "DAZ3" "DAZ1" "DAZ2" ...


# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Define a function to convert gene symbols to Entrez IDs
convertToEntrez <- function(gene_symbols) {
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Return the Entrez IDs
  return(entrez_ids$ENTREZID)
}

# Apply the function to each element of the list
entrez_ids_list <- lapply(set_r, convertToEntrez)
# Compare enrichment results between sets
compare_result_kegg <- compareCluster(enrich_results_kegg, fun = "enrichKEGG")

# Visualize the results
dotplot(compare_result_kegg)
simple <- simplify(
  compare_result_kegg,
  cutoff = 0.1,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)


library(clusterProfiler)


compare_cluster_results_up <- compareCluster(geneCluster = gl, 
                                             fun = "enrichGO",
                                             OrgDb = org.Hs.eg.db,
                                             ont = "ALL", 
                                             pvalueCutoff = 1.0,
                                             qvalueCutoff = 1.0)




simple <- simplify(
  compare_cluster_results_up,
  cutoff = 0.1,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)



dotplot(simple)





dotplot(compare_cluster_results_up)

saveRDS(compare_cluster_results_up, "GO_enrichment_for_uniquetfboundgenes_results.RDS")


compareCluster(geneCluster=gene_list_up, 
             ont ="ALL", 
             pvalueCutoff = 1.0, 
            OrgDb = org.Hs.eg.db,

             pAdjustMethod = "none")






# Perform GO enrichment analysis for downregulated genes
compare_cluster_results_down <- compareCluster(geneCluster = gene_list_down, 
                                               fun = "enrichGO",
                                               OrgDb = org.Hs.eg.db,
                                               ont = "BP", 
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.05)

# Generate dot plots
dotplot(compare_cluster_results_up) + ggtitle("GO Enrichment Analysis for Upregulated Genes")
dotplot(compare_cluster_results_down) + ggtitle("GO Enrichment Analysis for Downregulated Genes")
l















#### overlapping genes ( all tf's except mef2c, all tf's except hopx, all tf's)


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=30GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)




> # Define the factors you're interested in
factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt")

# Apply the function to the data frame x
common_genes <- get_intersect_members(x, factors)

> head(common_genes)
        genes_bound_gata4.txt genes_bound_mef2a.txt genes_bound_nkx25.txt
RALY                        1                     0                     1
CERK                        1                     0                     1
HRH4                        1                     1                     1
TNRC6B                      1                     0                     1
SLC23A2                     1                     0                     1
WDR19                       1                     0                     1
        genes_bound_SRF.txt genes_bound_tbx5.txt genes_bound_hopx.txt
RALY                      1                    1                    1
CERK                      1                    1                    1
HRH4                      1                    1                    1
TNRC6B                    1                    1                    1
SLC23A2                   1                    1                    1
WDR19                     1                    1                    1
> dim(common_genes)
[1] 816   6
> factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt", "genes_bound_mef2a.txt")
> # Define the factors you're interested in
factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt")

# Apply the function to the data frame x
common_genes <- get_intersect_members(x, factors)
> factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt", "genes_bound_mef2a.txt")

common_all <- get_intersect_members(x, factors)
Warning message:
Using `across()` in `filter()` was deprecated in dplyr 1.0.8.
ℹ Please use `if_any()` or `if_all()` instead.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 
> dim(common_all)
[1] 159   6
> dim(common_genes)
[1] 816   6
> history()

  result <- x %>%
    filter(across(all_of(factors), ~ . == 1)) %>%
    column_to_rownames(var = "gene")
  
  return(result)
}

# Define the factors you're interested in
factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt")

# Apply the function to the data frame x
common_genes <- get_intersect_members(x, factors)

# Print the resulting data frame
print(common_genes)

head(common_genes)
dim(common_genes)
factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt", "genes_bound_mef2a.txt")
factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt", "genes_bound_mef2a.txt")

common_all <- get_intersect_members(x, factors)

dim(common_all)
dim(common_genes)
history()






q


factors <- c("genes_bound_SRF.txt", "genes_bound_tbx5.txt", "genes_bound_gata4.txt", "genes_bound_hopx.txt", "genes_bound_nkx25.txt", "genes_bound_mef2a.txt")

common_all <- get_intersect_members(x, factors)
#159   6


> dim(common_all)
[1] 159   6
> dim(common_genes)
[1] 816   6






setwd("gene_list/")

#####

annot <- fread("/QRISdata/Q5149/tmem88_heritability/annot/geneMatrix/geneMatrix.tsv")
	annot$GeneID <- annot$gene_name
		head(annot)
		annot2 <- annot[,c("GeneID")]
		head(annot2)
		annot2$chr <- annot$hg19g0
		annot2$start <- annot$g1
		annot2$end <- annot$g2



library(dplyr)


match_and_append <- function(file_name, annot2, adjustment = 50000) {
  # Read the file with gene names
  data <- fread(file_name, header = FALSE, col.names = "GeneID", sep = "\t", stringsAsFactors = FALSE)

  # Convert the "GeneID" columns to lowercase
  data$GeneID <- tolower(data$GeneID)
  annot2$GeneID <- tolower(annot2$GeneID)

  # Merge the gene names with the "annot" data based on the "GeneID" column (case-insensitive)
  merged_data <- data %>%
    left_join(annot2, by = "GeneID")

  # Adjust the start and end positions
  merged_data$start <- merged_data$start + adjustment
  merged_data$end <- merged_data$end - adjustment

    # Remove rows with NA values in any column
  merged_data <- merged_data[complete.cases(merged_data), ]

  # Count the number of matches
  num_matches <- sum(!is.na(merged_data$GeneID))

  # Save the merged data to a new file with the original file name
  new_file_name <- sub(".txt", "_location.txt", file_name)
  write.table(merged_data, new_file_name, sep = "\t", row.names = FALSE, col.names=TRUE, quote=F)

  return(num_matches)
}









snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


	dim(snp)
	#[1] 7356518      10
	snp2 <- snp[,c(1,2,5,6,7)]
	colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")



library(data.table)
library(GenomicRanges)


# Assuming snp2 is a data.table containing information about SNPs with columns Chromosome and PhysPos
# Convert 'Chromosome' column to character format
snp2$Chromosome <- as.character(snp2$Chromosome)

# Create GenomicRanges object for snp2
snp2_ranges <- GRanges(seqnames = snp2$Chromosome, ranges = IRanges(start = snp2$PhysPos, end = snp2$PhysPos))


# Convert 'chr' column in each dataframe to character format and create GenomicRanges object
for (df_name in names(annotated_data)) {
  df <- annotated_data[[df_name]]
  df$chr <- as.character(df$chr)
  
  # Create GenomicRanges object for the dataframe
  df_ranges <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  # Create column name prefix for overlap columns
  overlap_col_prefix <- gsub("\\.txt$", "", df_name)
  
  # Check for overlaps and create new columns in snp2
  overlap_col <- paste0("overlap_", overlap_col_prefix)
  
  # Remove rows with NAs in 'start' or 'end' columns
  df <- na.omit(df)
  
  # Create GenomicRanges object for the current dataframe
  snp2[[overlap_col]] <- as.integer(overlapsAny(snp2_ranges, df_ranges))
}


print(head(snp2))



#################

 table(snp2$all)

      0       1 
7239696  116822 








x <- snp2[,c(2,53,6:52)]
setnames(x, "rsid", "SNP")


 fwrite(x, file = "top_250_annotation_quantitative.txt", append = FALSE, quote =FALSE,sep = "\t", col.names=TRUE)



setwd("/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/")

save.image("update.RData")

## jaccard index 

salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

	R
	library(data.table)
	library(dplyr)


	 gene_sets <- readRDS("tf_bound_genes.RDS")
	# Function to calculate Jaccard Index
	jaccard_index <- function(list1, list2) {
	  intersection <- length(intersect(list1, list2))
	  union <- length(union(list1, list2))
	  return(intersection / union)
	}

	# Calculate Jaccard Index for each pair of gene sets
	combinations <- combn(names(gene_sets), 2, simplify = FALSE)
	jaccard_results <- sapply(combinations, function(pair) {
	  list1 <- gene_sets[[pair[1]]]
	  list2 <- gene_sets[[pair[2]]]
	  index <- jaccard_index(list1, list2)
	  return(index)
	})

	# Convert results to a matrix
	jaccard_matrix <- matrix(0, nrow = length(gene_sets), ncol = length(gene_sets),
	                         dimnames = list(names(gene_sets), names(gene_sets)))

	for (i in 1:length(combinations)) {
	  pair <- combinations[[i]]
	  jaccard_matrix[pair[1], pair[2]] <- jaccard_results[i]
	  jaccard_matrix[pair[2], pair[1]] <- jaccard_results[i]
	}

	# Load the necessary library
	library(ggplot2)
	library(reshape2)

	# Melt the matrix for ggplot2
	jaccard_melted <- melt(jaccard_matrix)

	# Plot the heatmap
	x <-ggplot(data = jaccard_melted, aes(Var1, Var2, fill = value)) +
	  geom_tile(color = "white") +
	  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
	                       midpoint = 0.5, limit = c(0, 1), space = "Lab",
	                       name = "Jaccard Index") +
	  theme_minimal() +
	  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
	  coord_fixed() +
	  labs(title = "Heatmap of Jaccard Indices",
	       x = "Gene Set 1",
	       y = "Gene Set 2")


#######



######### adult vs fetal annotation and heritability ##############



salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)
getwd()
#/QRISdata/Q5149/TF_DamID/M2/biotin



x <- readRDS("humanGenesList.RDS")





# Define function to annotate each transcription factor list in x and adjust start/end values
annotate_transcription_factors_adjusted <- function(x, annot2) {
  # Iterate through each transcription factor list
  for (tf_name in names(x)) {
    # Get gene IDs from x[[tf_name]]
    gene_ids <- x[[tf_name]]
    
    # Find matching rows in annot2
    matching_rows <- annot2[annot2$GeneID %in% gene_ids, ]
    
    # Append annot2 columns to x[[tf_name]]
    if (nrow(matching_rows) > 0) {
      # Find the unique matching rows
      unique_matching_rows <- matching_rows[!duplicated(matching_rows$GeneID), ]
      
      # Get the indices of matching gene IDs in x[[tf_name]]
      indices <- match(unique_matching_rows$GeneID, gene_ids)
      
      # Create or update chr, start, and end columns in x[[tf_name]]
      x[[tf_name]] <- data.frame(
        GeneID = gene_ids,
        chr = NA_character_,
        start = NA_integer_,
        end = NA_integer_
      )
      x[[tf_name]]$chr[indices] <- unique_matching_rows$chr
      x[[tf_name]]$start[indices] <- unique_matching_rows$start - 50000  # Adjust start
      x[[tf_name]]$end[indices] <- unique_matching_rows$end + 50000      # Adjust end
    }
  }
  
  return(x)
}

# Apply the adjusted annotation function to your x list
x <- annotate_transcription_factors_adjusted(x, annot2)









overlap_columns <- grep("^overlap_", names(snp2), value = TRUE)

# Create a data.table to store the results
result_dt <- data.table(column_name = character(), table_result = list())

# Iterate over each overlap column
for (col in overlap_columns) {
  # Apply the table function to the column
  table_result <- table(snp2[[col]], useNA = "ifany")
  
  # Append results to result_dt
  result_dt <- rbind(result_dt, data.table(column_name = col, table_result = list(table_result)))
}


> result_dt
              column_name table_result
                   <char>       <list>
 1:   overlap_Adult_Gata4 6743786,....
 2:   overlap_Adult_Mef2a 6905410,....
 3:  overlap_Adult_Nkx2_5 6843457,....
 4:     overlap_Adult_Srf 6881556,....
 5:    overlap_Adult_Tbx5 7258811,....
 6:    overlap_Adult_Tead 7059586,....
 7:   overlap_Fetal_Gata4 6177327,....
 8:   overlap_Fetal_Mef2a 6110006,....
 9:   overlap_Fetal_Mef2c 3673071,....
10:  overlap_Fetal_Nkx2_5 6180521,....
11:     overlap_Fetal_Srf 6048487,....
12:    overlap_Fetal_Tbx5 4858212,....
13:    overlap_Fetal_Tead 6150930,....
14:  overlap_Common_Gata4 4165560,....
15:  overlap_Common_Mef2a 5195348,....
16:    overlap_Common_Srf 4940712,....
17:   overlap_Common_Tbx5 5675517,....
18:   overlap_Common_Tead 4063778,....
19: overlap_Common_Nkx2_5 4601758,....


7356518      25



