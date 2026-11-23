

cd /QRISdata/Q4879/WSB_Seurat_integrated
salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=0:45:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l


R
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


x <- readRDS("common_orig_log_FPKM2L_integrated.RDS")
library(Seurat)
dim(x)
str(x)




setwd("/QRISdata/Q5149/computation/dubois_development2022")





#I want to go rna and then the data slot for normalized data, the counts slot is raw data. 

subset <- x@assays$RNA$data
head(subset)
dim(subset)

#[1] 13656  4279



dir="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score"

	TRIAGE1 <- function(m, m_names = NULL, species = NULL, geneid = NULL, 
	log = NULL, dir = NULL){
			if (length(grep('^w$', ls(.GlobalEnv))) == 0){
				w <<- as.data.frame(data.table::fread(paste0(dir, 'repressive_hg19_cont_prop.txt')))
		        g <<- as.data.frame(data.table::fread(paste0(dir, 'ensembl_id_mapping.txt')))
			}
			if (length(m_names) == 0){m_names <- row.names(m)}

			#check species of input data
			if (species == "Human") {
				s.col <- "Gene"
			} else{
				s.col <- paste(species, "gene")
			}

		    #check ensembl id or gene names of input data gene 
			if (geneid == TRUE){
		    g.col <- "name"
		    } else{
		    g.col <- "stable ID"
		    }
		    
		    #get the colnames should be used to select column from reference data "ensembl_id_mapping.txt"
		    select.col <- paste(s.col, g.col)

		    colfn <- function(x){return(g[which(x == g[, select.col]), ][1,2])}
		    hgname <- sapply(m_names, colfn)

			mappable <- which(hgname %in% w$V1)
			l <- rep(0, nrow(m))

			RTSfn <- function(gene){return(as.numeric(w[which(w$V1 == gene),][1,2]))}
			l[mappable] <- sapply(hgname[mappable], RTSfn)
		    
		    #check if input normalised data is log transformed or not
		    if (log == TRUE){
		    	out <- log(m+1)*l

		    } else {out <- m*l}

			# row.names(out) <- m_names
			return(out)
		}



df1 <- subset 
df2 <- TRIAGE1(df1, species = "Mouse", geneid = FALSE, log = FALSE, dir = "/QRISdata/Q5149/TRIAGE1/")

# if geneid is ensemble change geneid parameter to false
#df2 <- TRIAGE1(df1, species = "Mouse", geneid = FALSE, log = FALSE, dir = "/QRISdata/Q5149/TRIAGE1/")

metadata <- x@meta.data
# Load necessary libraries
library(Matrix)

# Assuming df2 is your sparse matrix and metadata is your data frame
# First, transpose the sparse matrix to have cells as rows and genes as columns
df2_transposed <- Matrix::t(df2)

# Convert the transposed matrix to a data frame
df2_transposed_df <- as.data.frame(as(df2_transposed, "matrix"))

# Add the time point and cell type columns from metadata
df2_transposed_df$time_point <- metadata$time
df2_transposed_df$cell_type <- metadata$cell

# Use the aggregate function to calculate the mean for each gene, grouped by time point and cell type
df2_transposed_df$identifier <- paste(df2_transposed_df$time_point, df2_transposed_df$cell_type, sep = "_")

mean_expression_by_identifier <- aggregate(df2_transposed_df[numeric_columns], 
                                            by = list(identifier = df2_transposed_df$identifier), 
                                            FUN = mean, na.rm = TRUE)

# Filter metadata to include only cells from wild-type tissue
metadata_wt <- metadata[metadata$type == "WT", ]

# Filter the mean expression data to include only cells from wild-type tissue
mean_expression_by_identifier_wt <- mean_expression_by_identifier[mean_expression_by_identifier$identifier %in% paste(metadata_wt$time, metadata_wt$cell, sep = "_"), ]

 mean_exp <- as.data.frame(t (mean_expression_by_identifier_wt))




		# Assuming mean_exp is your data frame
gene_names <- rownames(mean_exp)

head(gene_names)
colnames(mean_exp)
mean_exp[1,]
colnames(mean_exp) <- mean_exp[1,]
colnames(mean_exp)
dim(mean_exp)
mean_exp[1:10,1:10]
#mean_exp <- mean_exp[-1,]




#redundant
	#colnames(mean_exp) <- mean_exp[1,]
	expression_values <- as.data.frame(mean_exp)


	expression_values <- expression_values[-1,-1]



	#convert mouse genes to human genes
	musGenes <- rownames(expression_values)


				# Basic function to convert mouse to human gene names
				convertMouseGeneList <- function(x){

						require("biomaRt")
						human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "dec2021.archive.ensembl.org")
						mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "dec2021.archive.ensembl.org")

					
				  # Query Ensembl using the main server
genesV2 <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                  values = musGenes, mart = mouse, attributesL = c("hgnc_symbol"), 
                  martL = human, uniqueRows = TRUE)

				  # Initialize a vector of NA
				  human_symbols <- rep(NA, length(x))

				  # Identify indices where conversion is successful
				  valid_indices <- nrow(genesV2) > 0

				  # Assign human gene symbols to the valid indices
				  human_symbols[valid_indices] <- unique(genesV2[valid_indices, 2])

				  # Print the result
				  head(human_symbols)

				  return(human_symbols)
				}



				expression_values$genes <- rownames(expression_values)

				# Remove duplicates based on gene IDs
expression_values <- expression_values[!duplicated(expression_values$genes), ]

# Match mouse gene IDs with human gene symbols
mapping <- setNames(genesV2$HGNC.symbol, genesV2$Gene.stable.ID)

# Replace mouse gene IDs in the "genes" column with human gene symbols
expression_values$genes <- mapping[expression_values$genes]

# Check if there are any NA values in the "genes" column
missing_symbols <- sum(is.na(expression_values$genes))
if (missing_symbols > 0) {
  # Remove rows with missing human gene symbols
  expression_values <- expression_values[!is.na(expression_values$genes), ]
}

# Assign the "genes" column as row names
rownames(expression_values) <- expression_values$genes

# Remove the "genes" column
expression_values$genes <- NULL


#removing genes column if it still exists
mean_exp <- expression_values[,-48]


  gene_names <- rownames(mean_exp)

#convert columns to numeric
# Create a list to store the ordered genes and expression values for each cell type
ordered_genes_list <- lapply(1:ncol(mean_exp), function(i) {
  cell_type <- colnames(mean_exp)[i]

  # Check the class of the selected column
  if (!is.numeric(mean_exp[, i])) {
    # If not numeric, attempt to convert to numeric
    mean_exp[, i] <- as.numeric(as.character(mean_exp[, i]))
  }
  
  # Check if conversion was successful
  if (is.numeric(mean_exp[, i])) {
    # Sort the indices based on expression values
    sorted_indices <- order(-mean_exp[, i])
    genes <- gene_names[sorted_indices]
    expression_vals <- mean_exp[sorted_indices, i]
  
    return(data.frame(cell_type = cell_type, gene = genes, expression_value = expression_vals))
  } else {
    # Return NA if the column couldn't be converted to numeric
    return(NA)
  }
})

# Remove any elements with NA (indicating conversion failure)
ordered_genes_list <- ordered_genes_list[sapply(ordered_genes_list, function(x) !is.null(x$cell_type))]

setwd("/QRISdata/Q5149/computation/dubois_development2022")



# Save each data frame as a separate txt file in the new folder

# Create a new folder in the current directory
wsb_cell_type_time <- "ordered_genes_files"
dir.create(wsb_cell_type_time, showWarnings = FALSE)

# Save each data frame as a separate txt file in the new folder
for (i in seq_along(ordered_genes_list)) {
  cell_type <- ordered_genes_list[[i]]$cell_type[1]
  file_path <- file.path(wsb_cell_type_time, paste0("ordered_genes_", cell_type, ".txt"))
  write.table(ordered_genes_list[[i]], file = file_path, sep = "\t", row.names = FALSE)
  cat("Data frame for", cell_type, "saved as", file_path, "\n")
}


setwd("wsb_ordered_genes_files/")



#@@@@@@@@@@@@@@@@@@@@


files <- list.files()

# Set the topN values
topN_values <- c( 25, 50, 250, 500)

# Loop through topN values
for (topN in topN_values) {
  # Process each file for the current topN value
  for (file in files) {
    process_file(file, topN)
  }
}

# Function to process each file
process_file <- function(file, topN) {
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(FALSE)
  }

  file_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(file_data) == 0) {
    warning(paste("File is empty:", file))
    return(FALSE)
  }

  top_rows <- file_data[1:min(nrow(file_data), topN), ]

  # Extract the base name of the file without extension
  base_name <- tools::file_path_sans_ext(basename(file))

  # Construct the destination file path with a subfolder
  dest_folder <- file.path("ordered_genes_files", paste0("top_", topN))
  dest_path <- file.path(dest_folder, paste0("top_", topN, "_", gsub(" ", "_", base_name), ".txt"))

  # Create the destination folder if it doesn't exist
  dir.create(dest_folder, recursive = TRUE, showWarnings = FALSE)

  # Print information for debugging
  cat("Processing:", file, "TopN:", topN, "\n")
  cat("Destination folder:", dest_folder, "\n")
  cat("Destination path:", dest_path, "\n")

  # Write the top rows to the destination file
  write.table(top_rows, dest_path, sep = "\t", row.names = FALSE)
  
  if (file.exists(dest_path)) {
    cat("File processed:", dest_path, "\n")
    return(TRUE)
  } else {
    warning(paste("Failed to write file:", dest_path))
    return(FALSE)
  }
}


#################################

setwd("ordered_genes_files/top_50")



# Get list of text files in the directory
file_list <- list.files(pattern = "\\.txt$", full.names = TRUE)

# Load the geneMatrix data
mtrix <- fread("/QRISdata/Q5149/tmem88_heritability/annot/geneMatrix/geneMatrix.tsv")
m <- mtrix[, c(2, 12, 13, 14, 15)]
colnames(m) <- c("gene", "chr", "start", "stop", "strand")
annotated_data <- list()


# Initialize the list of annotated data
annotated_data <- list()

# Loop through each unique file
for (file in unique(file_list)) {
  # Read the gene list from the file
  genelist <- fread(file)

  # Rename columns
  colnames(genelist) <- c("cell_type", "gene", "triage_discordance_score")

  # Convert gene names to uppercase
  genelist$gene <- toupper(genelist$gene)

  # Convert gene names in m to uppercase
  m$gene <- toupper(m$gene)

  # Merge gene list with m
  genelist1 <- merge(x = genelist, y = m, by = "gene", all.x = TRUE)

  # Modify start and stop positions
  genelist1$start <- genelist1$start - 50000
  genelist1$stop <- genelist1$stop + 50000

  # Remove the './' prefix from file name
  file <- sub("^./", "", file)

  print(colnames(genelist1))

  # Count the number of rows with at least one NA value
  num_na_rows <- sum(rowSums(is.na(genelist1)) > 0)
  cat("Number of rows with at least one NA value:", num_na_rows, "\n")

  # Remove rows with NA values
  genelist1 <- genelist1[complete.cases(genelist1), ]

  # Prepare data for annotation
  x <- genelist1[, 3:6]
  x$chr <- sub("chr", "", x$chr)

  # Write annotated data to a new file
  annotated_data[[file]] <- x
}
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
  df_ranges <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$stop))
  
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


#(47 different annotations)



>table(snp2$overlap_top_50_ordered_genes_p7_cardiomyocyte)

      0       1 
7313816   42702 




snp2$Intercept <- "1"

###################### NEEDS TO BE EDITED each time 
x <- snp2[,c(2,53,6:52)]
setnames(x, "rsid", "SNP")


 fwrite(x, file = "top_250_annotation_quantitative.txt", append = FALSE, quote =FALSE,sep = "\t", col.names=TRUE)



setwd("/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/")

save.image("update.RData")



read.table("filename", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

######## BINARY annotation(?)

  file
#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=250GB
#SBATCH --time=7:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=mri
#SBATCH --account=a_palpant
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/imputed_sumstats/
config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M/"
annot="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_top50', annot='$annot', log2file=TRUE)"







# trouble shooting 

> setwd="/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/imputed_sumstats/"
> setwd(setwd)
> library(SBayesRC)
Loading required package: data.table
data.table 1.15.0 using 1 threads (see ?getDTthreads).  Latest news: r-datatable.com
> mafile = "svi_bsa.ma_imp.ma"
> x <- fread(mafile)


#--------- if you get illegal coredumped it might be a memory error






















############## IMPORTANT DIRECTORIES 


* MRI 7m impute files 

/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/imputed_sumstats



*MRI files before tidy 

/QRISdata/Q4879/cardiac_gwas_sumstats










#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=180GB
#SBATCH --time=12:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general 
#SBATCH --job-name=cardiac
#SBATCH --account=a_palpant
#SBATCH --array=1

MY_DIR=/QRISdata/Q4879/cardiac_gwas_sumstats/
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/scratch/user/uqsnegi/ukbEUR_Imputed"
annot="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results"
threads=4
 
export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${ma_file}_tidy.ma', log2file=TRUE)"

Rscript -e "SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$ld_folder', \
                  output='${ma_file}_imp.ma', log2file=TRUE)"

Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"



#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=180GB
#SBATCH --time=01:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=debug 
#SBATCH --job-name=cardiac
#SBATCH --account=a_palpant
#SBATCH --array=-19

MY_DIR=/QRISdata/Q4879/cardiac_gwas_sumstats/
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M/"
annot="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}${ma_file}_tidy.ma', log2file=TRUE)"

Rscript -e "SBayesRC::impute(mafile='${out_prefix}${ma_file}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}imp_${ma_file}', log2file=TRUE)"

Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}imp_${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}${ma_file}_sbrc', annot='$annot', log2file=TRUE)"








#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=1:00:00
#SBATCH --job-name=sbayesrc
#SBATCH --error=err
#SBATCH --output=out
#SBATCH --partition=debug
#SBATCH --account=a_palpant
#SBATCH --array=1-19

MY_DIR=/QRISdata/Q4879/cardiac_gwas_sumstats/
config=config.txt

cd $MY_DIR

gctb=/scratch/user/uqsnegi/gctb/gctb_2.5.1_Linux/gctb

gwas=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ldm=/scratch/user/uqsnegi/ukbEUR_Imputed
anno=/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt
out=/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results/xd

$gctb --ldm-eigen $ldm --gwas-summary $gwas --sbayes RC --annot $anno --thread 10 --out $out$gwas >$out$gwas.log 2>&1















############# ACTUALLY WORKS ######





#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --job-name=sbayesrc
#SBATCH --error=err
#SBATCH --output=out
#SBATCH --partition=general
#SBATCH --account=a_palpant
#SBATCH --array=1-19

MY_DIR=/QRISdata/Q4879/cardiac_gwas_sumstats/
config=config.txt

cd $MY_DIR

gctb=/scratch/user/uqsnegi/gctb/gctb_2.5.1_Linux/gctb

gwas=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ldm=/scratch/user/uqsnegi/ukbEUR_Imputed
anno=/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt
out=/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results/new_

$gctb --ldm-eigen $ldm --gwas-summary $gwas --sbayes RC --annot $anno --thread 10 --out $out$gwas >$out$gwas.log 2>&1




















######################### CELL TYPE SPECIFICITY VISUALIZATION ###########################

salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=50GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

 

R
library(data.table)

setwd("/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/")




# Load required libraries
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)

# List all text files in the directory
file_list <- list.files(pattern = "*.txt")

# Initialize an empty data frame
gene_presence_df <- data.frame()

# Iterate through each file in the directory
for (file in file_list) {
  # Read the genes from the file
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract relevant columns (cell type and gene)
  cell_type <- unique(data$cell_type)
  genes <- data$gene
  
  # Create a data frame with genes as rows and file name as a column
  gene_df <- data.frame(Gene = genes, Presence = 1, File = cell_type)
  
  # Append to the main data frame
  gene_presence_df <- bind_rows(gene_presence_df, gene_df)
}

# Spread the data frame to create a binary matrix
gene_presence_matrix <- gene_presence_df %>%
  spread(File, Presence, fill = 0)

# Set the gene names as rownames and remove the Gene column
rownames(gene_presence_matrix) <- gene_presence_matrix$Gene
gene_presence_matrix$Gene <- NULL

# Convert the data frame to a matrix
gene_presence_matrix <- as.matrix(gene_presence_matrix)

# Plot the heatmap with adjusted font size
ht <- pheatmap(gene_presence_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("white", "blue"))(50),
         main = "Gene Presence Across Text Files",
         fontsize_row = 6, # Adjust this value to change the font size of row names
         fontsize_col = 8, # Adjust this value to change the font size of column names
         cellwidth = 10,   # Adjust this value to change the width of cells
         cellheight = 10)  # Adjust this value to change the height of cells

pdf("results/cell_type_specificity_heatmap.pdf", width = 105/25.4, height = 297/25.4) # Width and height in inches
ht
dev.off()





upset(gene_df,
	nsets = 47,
	sets = c("E10.5_endothelial", "E10.5_epicardial", "E10.5_mesenchymal_cell",
                "E11.5_cardiomyocyte", "E11.5_endothelial", "E11.5_fibroblast",
                "E14.5_cardiomyocyte", "E14.5_endothelial", "E14.5_fibroblast",
                "E18.5_cardiomyocyte", "E18.5_endothelial", "E18.5_fibroblast",
                "E6.5_epiblast_cells", "E6.5_others", "E6.75_endothelial",
                "E6.75_epiblast_cells", "E6.75_notochord_or_DEF_endoderm",
                "E6.75_others", "E7.25_anterior_SHF_or_pharyngeal_mesoderm",
                "E7.25_cardiomyocytes_or_mesenchymal", "E7.25_endothelial",
                "E7.25_epiblast_cells", "E7.25_notochord_or_DEF_endoderm",
                "E7.25_others", "E7.5_anterior_SHF_or_pharyngeal_mesoderm",
                "E7.5_cardiomyocytes_or_mesenchymal", "E7.5_endothelial",
                "E7.5_epiblast_cells", "E7.5_notochord_or_DEF_endoderm",
                "E7.5_others", "E7.5_posterior_SHF_or_somatic_mesoderm",
                "E8.5_cardiomyocyte", "E9.5_cardiomyocyte", "E9.5_endothelial",
                "E9.5_mesenchymal_cell", "p0_cardiomyocyte", "p0_endothelial",
                "p0_fibroblast", "p21_cardiomyocyte", "p21_endothelial",
                "p21_fibroblast", "p3_cardiomyocyte", "p3_endothelial",
                "p3_fibroblast", "p7_cardiomyocyte", "p7_endothelial",
                "p7_fibroblast"),
	order.by = "freq",
    sets.bar.color = "skyblue",
    matrix.color = "lightgray",
    text.scale = 1.5,
    point.size = 3,
    empty.intersections = "off",
      )








#################


# Load necessary libraries
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
library(biomaRt)
library(data.table)
library(parallel)

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(mouse_genes) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  # Query Ensembl using the main server
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = mouse_genes, mart = mouse, attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows = TRUE)
  
  # Create a named vector for fast lookup
  mouse_to_human <- setNames(genesV2$HGNC.symbol, genesV2$MGI.symbol)
  
  # Convert the mouse genes to human genes
  human_genes <- mouse_to_human[mouse_genes]
  
  # Replace NA with the original mouse gene name if no match found
  human_genes[is.na(human_genes)] <- mouse_genes[is.na(human_genes)]
  
  return(human_genes)
}

# Function to process each file
processFile <- function(file_path) {
  # Read the file
  data <- fread(file_path, header = TRUE, sep = "\t")
  
  # Original mouse gene names
  original_genes <- data$gene
  
  # Convert mouse gene names to human gene names
  human_genes <- convertMouseGeneList(original_genes)
  
  # Identify successful conversions
  successful <- original_genes != human_genes
  successful_indices <- which(successful)
  
  # Limit to first 50 successful conversions
  if (length(successful_indices) > 100) {
    selected_indices <- successful_indices[1:100]
  } else {
    selected_indices <- successful_indices
  }
  
  # Subset the data to include only the first 50 successful conversions
  data_subset <- data[selected_indices, ]
  data_subset$gene <- human_genes[selected_indices]
  
  # Create new file name with "human_" prefix
  file_dir <- dirname(file_path)
  file_name <- basename(file_path)
  new_file_name <- paste0("human100_", file_name)
  new_file_path <- file.path(file_dir, new_file_name)
  
  # Write the updated data to the new file
  fwrite(data_subset, new_file_path, sep = "\t", quote = FALSE)
  
  # Print conversion summary
  cat(sprintf("File: %s - Successfully converted %d out of %d genes (%.2f%%)\n", 
              new_file_name, length(selected_indices), length(original_genes), 
              (length(selected_indices) / length(original_genes)) * 100))
}

# Directory containing the files

# List all txt files in the directory
files <- list.files(input_directory, pattern = "*.txt", full.names = TRUE)

# Parallel processing of files
num_cores <- detectCores() - 1
mclapply(files, processFile, mc.cores = num_cores)

#################




# Load necessary libraries
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
library(biomaRt)
library(data.table)
library(parallel)

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(mouse_genes) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  # Query Ensembl using the main server
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = mouse_genes, mart = mouse, attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows = TRUE)
  
  # Create a named vector for fast lookup
  mouse_to_human <- setNames(genesV2.symbol, genesV2.symbol)
  
  # Convert the mouse genes to human genes
  human_genes <- mouse_to_human[mouse_genes]
  
  # Replace NA with the original mouse gene name if no match found
  human_genes[is.na(human_genes)] <- mouse_genes[is.na(human_genes)]
  
  return(human_genes)
}

# Function to process each file
processFile <- function(file_path) {
  # Read the file
  data <- fread(file_path, header = TRUE, sep = "\t")
  
  # Original mouse gene names
  original_genes <- data
  
  # Convert mouse gene names to human gene names
  human_genes <- convertMouseGeneList(original_genes)
  
  # Identify successful conversions
  successful <- original_genes != human_genes
  successful_indices <- which(successful)
  
  # Limit to first 50 successful conversions
  if (length(successful_indices) > 50) {
    selected_indices <- successful_indices[1:50]
  } else {
    selected_indices <- successful_indices
  }
  
  # Subset the data to include only the first 50 successful conversions
  data_subset <- data[selected_indices, ]
  data_subset <- human_genes[selected_indices]
  
  # Create new file name with "human_" prefix
  file_dir <- dirname(file_path)
  file_name <- basename(file_path)
  new_file_name <- paste0("human50_", file_name)
  new_file_path <- file.path(file_dir, new_file_name)
  
  # Write the updated data to the new file
  fwrite(data_subset, new_file_path, sep = "\t", quote = FALSE)
  
  # Print conversion summary
  cat(sprintf("File: %s - Successfully converted %d out of %d genes (%.2f%%)\n", 
              new_file_name, length(selected_indices), length(original_genes), 
              (length(selected_indices) / length(original_genes)) * 100))
}

# Directory containing the files

# List all txt files in the directory
files <- list.files(input_directory, pattern = "^ordered.*\\.txt$"*.txt", full.names = TRUE)

# Parallel processing of files
num_cores <- detectCores() - 1
mclapply(files, processFile, mc.cores = num_cores)"









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
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/1m_imputed_sumstats
config=config.txt

cd $MY_DIR


# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/ukbEUR_HM3/"
annot="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/annot/top_50_annot.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/wsb_ordered_genes_files/ordered_genes_files/top_50/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_dalia', annot='$annot', log2file=TRUE)"






####### JACCARD INDEX SIMILARIY 

gene_sets <- readRDS("upset_gene_sets.rds")




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)



# Define the list of gene sets
gene_sets <- list(
  genes_bound_gata4 = c("MAFIP", "RELL1", "UCHL1", "RRP1"),
  genes_bound_mef2a = c("ZCCHC2", "AKAP6", "NCOA3", "PTPN21"),
  genes_bound_nkx25 = c("MAFIP", "PSPC1", "SLC37A1", "CDK5RAP1"),
  genes_bound_SRF = c("MAFIP", "NRIP1", "SEZ6L", "MCAT"),
  genes_bound_tbx5 = c("MAFIP", "DAZ4", "DAZ3", "DAZ1"),
  genes_bound_hopx = c("JAKMIP1", "TMEM184B", "RGS6", "CNMD")
)

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
ggplot(data = jaccard_melted, aes(Var1, Var2, fill = value)) +
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
