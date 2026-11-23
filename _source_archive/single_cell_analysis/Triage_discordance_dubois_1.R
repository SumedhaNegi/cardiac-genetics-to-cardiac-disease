#Discordance gene transformation 

cd /QRISdata/Q5149/computation/dubois_development2022

cd triage_discordance_score/

salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=00:20:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l


R
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


setwd("/QRISdata/Q5149/computation/dubois_development2022")


x <- readRDS("GSE205950_TP_merge_cardiac.rds")
library(Seurat)
dim(x)
str(x)


#I want to go rna and then the data slot for normalized data, the counts slot is raw data. 

subset <- x@assays$RNA$data
head(subset)
dim(subset)



# x@assays$integrated$scale.data used for WSB

#21738 12832

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
df2 <- TRIAGE1(df1, species = "Mouse", geneid = TRUE, log = FALSE, dir = "/QRISdata/Q5149/TRIAGE1/")

# if geneid is ensemble change geneid parameter to false
#df2 <- TRIAGE1(df1, species = "Mouse", geneid = FALSE, log = FALSE, dir = "/QRISdata/Q5149/TRIAGE1/")


dim(df2)

#21738 12832

#save df2 as the output

#finding the mean expresiion value for each cell type and sorting them in a decreasing order 

		metadata <- x@meta.data

		# Assuming df2 is your sparse matrix and metadata is your data frame
		# First, transpose the sparse matrix to have cells as rows and genes as columns

		library(Matrix)

		# Transpose the sparse matrix
		df2_transposed <- Matrix::t(df2)

		# Add a column to the transposed matrix with cell types from the metadata
		df2_transposed_df <- as.data.frame(as(df2_transposed, "matrix"))

		# Add the cell type column from metadata
		df2_transposed_df$cell_type <- metadata$cell.type.anno.man

		# Use the aggregate function to calculate the mean for each cell type
		mean_expression_by_cell_type <- aggregate(. ~ cell_type, data = df2_transposed_df, mean)

		# View the result
		head(mean_expression_by_cell_type)


		mean_exp <- as.data.frame(t(mean_expression_by_cell_type))


		# Assuming mean_exp is your data frame
gene_names <- rownames(mean_exp)

head(gene_names)
colnames(mean_exp)
mean_exp[1,]
colnames(mean_exp) <- mean_exp[1,]
colnames(mean_exp)
dim(mean_exp)
mean_exp[1:10,1:10]
mean_exp <- mean_exp[-1,]




#redundant
	colnames(mean_exp) <- mean_exp[1,]
	expression_values <- as.data.frame(mean_exp)

	colnames(expression_values) <- expression_values[1,]
	colnames(expression_values)
	rownames(expression_values)


	expression_values <- expression_values[-1,-1]

	exprs_val <- apply(expression_values, 2, as.numeric)


	#convert mouse genes to human genes
	musGenes <- rownames(expression_values)


				# Basic function to convert mouse to human gene names
				convertMouseGeneList <- function(x){

						require("biomaRt")
						human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "dec2021.archive.ensembl.org")
						mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "dec2021.archive.ensembl.org")

					
				  # Query Ensembl using the main server
				  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)

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
						


				human_expression_values <- convertMouseGeneList(musGenes)
				length(human_expression_values) 

				#21738

				length(musGenes)
				#21738

#took too long
	# Define the function to convert mouse genes to human genes
	convert_mouse_to_human <- function(gene_list) { 
	  output = c()
	  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

	  for(gene in gene_list) {
	    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
	    if( !identical(class_key, integer(0)) ) {
	      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
	      for(human_gene in human_genes) {
	        output = rbind(c(gene, human_gene), output)
	      }
	    }
	  }
	  return (output)
	}

	# Define the function to convert mouse genes to human genes in batches
	convert_mouse_to_human_batches <- function(gene_list, batch_size) {
	  # Initialize an empty list to store converted genes
	  converted_genes_list <- list()
	  
	  # Determine the number of batches
	  num_batches <- ceiling(length(gene_list) / batch_size)
	  
	  # Split the gene list into batches
	  gene_batches <- split(gene_list, ceiling(seq_along(gene_list) / batch_size))
	  
	  # Iterate over each batch
	  for (i in seq_len(num_batches)) {
	    # Get the current batch
	    current_batch <- gene_batches[[i]]
	    
	    # Convert the current batch of genes
	    converted_genes <- lapply(current_batch, convert_mouse_to_human)
	    
	    # Combine the converted genes into a single data frame
	    converted_genes_df <- do.call(rbind, converted_genes)
	    
	    # Store the converted genes in the list
	    converted_genes_list[[i]] <- converted_genes_df
	    
	    # Print a message indicating that the batch has been converted
	    cat("Batch", i, "has been converted.\n")
	  }
	  
	  # Return the list of converted genes
	  return(converted_genes_list)
	}

	# Apply the function to all row names of mean_exp
	converted_genes <- convert_mouse_to_human_batches(rownames(mean_exp), batch_size)



#quick fix

> x <- toupper(rownames(mean_exp))
> rownames(mean_exp) <- x
> head(rownames(mean_exp))






# Create a list to store the ordered genes and expression values for each cell type
ordered_genes_list <- lapply(1:ncol(mean_exp), function(i) {
  cell_type <- colnames(mean_exp)[i]
  sorted_indices <- order(-mean_exp[, i])
  genes <- gene_names[sorted_indices]
  expression_vals <- mean_exp[sorted_indices, i]
  
  return(data.frame(cell_type = cell_type, gene = genes, expression_value = expression_vals))
})

# View the first element of the resulting list
print(ordered_genes_list[[1]])

for (i in seq_along(ordered_genes_list)) {
  cat("Cell Type:", ordered_genes_list[[i]]$cell_type[1], "\n")
  print(head(ordered_genes_list[[i]], 10))
  cat("\n")
}





		#############saving each cell type as a different txt file 


# Save each data frame as a separate txt file in the new folder

		# Create a new folder in the current directory
		mean_ds_by_cell_type <- "ordered_genes_files"
		dir.create(mean_ds_by_cell_type, showWarnings = FALSE)

		# Function to clean cell type names for file paths
		clean_cell_type <- function(cell_type) {
		  cleaned_name <- gsub(" ", "_", cell_type)
		  cleaned_name <- gsub("/", "_", cleaned_name)
		  return(cleaned_name)
		}

		# Save each data frame as a separate txt file in the new folder
		for (i in seq_along(ordered_genes_list)) {
		  cell_type <- ordered_genes_list[[i]]$cell_type[1]
		  cleaned_cell_type <- clean_cell_type(cell_type)
		  file_path <- file.path(mean_ds_by_cell_type, paste0("ordered_genes_", cleaned_cell_type, ".txt"))
		  write.table(ordered_genes_list[[i]], file = file_path, sep = "\t", row.names = FALSE)
		  cat("Data frame for", cell_type, "saved as", file_path, "\n")
		}



save.image("update.RData")

#Cell clusters 
				Sublist 1 : 1 posterior second heart field 
				Sublist 2 : 2 RV 
				Sublist 3 : 3 hematoendothelial progenitor 
				Sublist 4 : 4 atrial CM 
				Sublist 5 : 5 aSHF myocardial cells/OFT 
				Sublist 6 : 6 lateral plate mesoderm 
				Sublist 7 : 7 atrioventricular canal 
				Sublist 8 : 8 differentiating aSHF 
				Sublist 9 : 9 endothelial cells 
				Sublist 10 : 10 endothelial cells 
				Sublist 11 : 11 paraxial mesoderm 
				Sublist 12 : 12 sinus venosus 
				Sublist 13 : 13 endothelial cells 
				Sublist 14 : 14 differentiating pSHF 
				Sublist 15 : 15 cardiac committed mesodermal progenitor 
				Sublist 16 : 16 endothelial cells 
				Sublist 17 : 17 LV 
				Sublist 18 : 18 neural crest 
				Sublist 19 : 19 anterior second heart field 
				Sublist 20 : 20 paraxial mesoderm 
				Sublist 21 : 21 posterior second heart field 
				Sublist 22 : 22 anterior second heart field 
				Sublist 23 : 23 endocardial cells 
				Sublist 24 : 24 hematoendothelial progenitor 
				Sublist 25 : 25 epicardium 
				Sublist 26 : 26 lateral plate mesoderm 





#################################################################################################################################################
###############################################################################################################################################

# extracting top 20 , 50 , 100 , 250 and 500 genes from each cell cluster 

setwd("/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/")
# List all files in the current directory with a .txt extension

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



##################




setwd("ordered_genes_files/top_20")

# let's annotate top 20, top 100, top 500 and run heritability.

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


       table_result			
 1: 7351128,   5390
 2: 7351919,   4599
 3: 7350732,   5786
 4: 7351293,   5225
 5: 7350936,   5582
 6: 7347897,   8621
 7: 7350027,   6491
 8: 7351962,   4556
 9: 7348863,   7655
10: 7351058,   5460
11: 7352523,   3995
12: 7352357,   4161
13: 7350360,   6158
14: 7352131,   4387
15: 7351337,   5181
16: 7352557,   3961
17: 7350853,   5665
18: 7349980,   6538
19: 7352154,   4364
20: 7352386,   4132
21: 7350243,   6275
22: 7351341,   5177
23: 7351714,   4804
24: 7352614,   3904
25: 7351374,   5144
26: 7353251,   3267
27: 7349508,   7010

 [1] "Chromosome"                                                             
 [2] "rsid"                                                                   
 [3] "PhysPos"                                                                
 [4] "A1"                                                                     
 [5] "A2"                                                                     
 [6] "overlap_top_20_ordered_genes_0_early_myocardial_progenitors"            
 [7] "overlap_top_20_ordered_genes_1_posterior_second_heart_field"            
 [8] "overlap_top_20_ordered_genes_10_endothelial_cells"                      
 [9] "overlap_top_20_ordered_genes_11_paraxial_mesoderm"                      
[10] "overlap_top_20_ordered_genes_12_sinus_venosus"                          
[11] "overlap_top_20_ordered_genes_13_endothelial_cells"                      
[12] "overlap_top_20_ordered_genes_14_differentiating_pSHF"                   
[13] "overlap_top_20_ordered_genes_15_cardiac_committed_mesodermal_progenitor"
[14] "overlap_top_20_ordered_genes_16_endothelial_cells"                      
[15] "overlap_top_20_ordered_genes_17_LV"                                     
[16] "overlap_top_20_ordered_genes_18_neural_crest"                           
[17] "overlap_top_20_ordered_genes_19_anterior_second_heart_field"            
[18] "overlap_top_20_ordered_genes_2_RV"                                      
[19] "overlap_top_20_ordered_genes_20_paraxial_mesoderm"                      
[20] "overlap_top_20_ordered_genes_21_posterior_second_heart_field"           
[21] "overlap_top_20_ordered_genes_22_anterior_second_heart_field"            
[22] "overlap_top_20_ordered_genes_23_endocardial_cells"                      
[23] "overlap_top_20_ordered_genes_24_hematoendothelial_progenitor"           
[24] "overlap_top_20_ordered_genes_25_epicardium"                             
[25] "overlap_top_20_ordered_genes_26_lateral_plate_mesoderm"                 
[26] "overlap_top_20_ordered_genes_3_hematoendothelial_progenitor"            
[27] "overlap_top_20_ordered_genes_4_atrial_CM"                               
[28] "overlap_top_20_ordered_genes_5_aSHF_myocardial_cells_OFT"               
[29] "overlap_top_20_ordered_genes_6_lateral_plate_mesoderm"                  
[30] "overlap_top_20_ordered_genes_7_atrioventricular_canal"                  
[31] "overlap_top_20_ordered_genes_8_differentiating_aSHF"                    
[32] "overlap_top_20_ordered_genes_9_endothelial_cells" 






x <- snp2[,c(2,33,6:32)]
setnames(x, "rsid", "SNP")


 fwrite(x, file = "top_20_annotation_binary.txt", append = FALSE, quote =FALSE,
                                  sep = "\t", col.names=TRUE)



# top 20 


#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=top50
#SBATCH --account=a_palpant
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/
config=config.txt

cd $MY_DIR


module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_20/top_20_annotation_binary.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/ordered_genes_files/top_20/results/"
threads=4

export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}${ma_file}_top50', annot='$annot', log2file=TRUE)"








 # top 250 

        table_result
 1: 7288187,  68331
 2: 7289135,  67383
 3: 7285727,  70791
 4: 7283340,  73178
 5: 7281952,  74566
 6: 7279362,  77156
 7: 7288785,  67733
 8: 7291255,  65263
 9: 7275820,  80698
10: 7281113,  75405
11: 7283478,  73040
12: 7278791,  77727
13: 7281288,  75230
14: 7281498,  75020
15: 7291594,  64924
16: 7288398,  68120
17: 7276424,  80094
18: 7281850,  74668
19: 7282630,  73888
20: 7286654,  69864
21: 7274725,  81793
22: 7281405,  75113
23: 7282278,  74240
24: 7287673,  68845
25: 7278292,  78226
26: 7291520,  64998
27: 7282868,  73650
       table_result


 [1] "Chromosome"                                                              
 [2] "rsid"                                                                    
 [3] "PhysPos"                                                                 
 [4] "A1"                                                                      
 [5] "A2"                                                                      
 [6] "overlap_top_250_ordered_genes_0_early_myocardial_progenitors"            
 [7] "overlap_top_250_ordered_genes_1_posterior_second_heart_field"            
 [8] "overlap_top_250_ordered_genes_10_endothelial_cells"                      
 [9] "overlap_top_250_ordered_genes_11_paraxial_mesoderm"                      
[10] "overlap_top_250_ordered_genes_12_sinus_venosus"                          
[11] "overlap_top_250_ordered_genes_13_endothelial_cells"                      
[12] "overlap_top_250_ordered_genes_14_differentiating_pSHF"                   
[13] "overlap_top_250_ordered_genes_15_cardiac_committed_mesodermal_progenitor"
[14] "overlap_top_250_ordered_genes_16_endothelial_cells"                      
[15] "overlap_top_250_ordered_genes_17_LV"                                     
[16] "overlap_top_250_ordered_genes_18_neural_crest"                           
[17] "overlap_top_250_ordered_genes_19_anterior_second_heart_field"            
[18] "overlap_top_250_ordered_genes_2_RV"                                      
[19] "overlap_top_250_ordered_genes_20_paraxial_mesoderm"                      
[20] "overlap_top_250_ordered_genes_21_posterior_second_heart_field"           
[21] "overlap_top_250_ordered_genes_22_anterior_second_heart_field"            
[22] "overlap_top_250_ordered_genes_23_endocardial_cells"                      
[23] "overlap_top_250_ordered_genes_24_hematoendothelial_progenitor"           
[24] "overlap_top_250_ordered_genes_25_epicardium"                             
[25] "overlap_top_250_ordered_genes_26_lateral_plate_mesoderm"                 
[26] "overlap_top_250_ordered_genes_3_hematoendothelial_progenitor"            
[27] "overlap_top_250_ordered_genes_4_atrial_CM"                               
[28] "overlap_top_250_ordered_genes_5_aSHF_myocardial_cells_OFT"               
[29] "overlap_top_250_ordered_genes_6_lateral_plate_mesoderm"                  
[30] "overlap_top_250_ordered_genes_7_atrioventricular_canal"                  
[31] "overlap_top_250_ordered_genes_8_differentiating_aSHF"                    
[32] "overlap_top_250_ordered_genes_9_endothelial_cells"    



snp2$Intercept <- "1"
x <- snp2[,c(2,33,6:32)]
setnames(x, "rsid", "SNP")


 fwrite(x, file = "top_250_annotation_binary.txt", append = FALSE, quote =FALSE,
                                  sep = "\t", col.names=TRUE)



#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=top250
#SBATCH --account=a_palpant
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/
config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_250/annot/top_250_annotation_binary.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_20/results/"
threads=4

export OMP_NUM_THREADS=$threads

# Corrected the closing parenthesis for the Rscript command
Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}${ma_file}_top250', annot='$annot', log2file=TRUE)"




########## 


quantitative annotation 



snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


	dim(snp)
	#[1] 7356518      10
	snp2 <- snp[,c(1,2,5,6,7)]
	colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")




file_list <- list.files(pattern = "\\.txt$")
              
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

  # Write annotated data to a new fileX
  annotated_data[[file]] <- x
}



library(data.table)
library(GenomicRanges)

# Assuming snp2 is a data.table containing information about SNPs with columns Chromosome and PhysPos
# Convert 'Chromosome' column to character format
snp2$Chromosome <- as.character(snp2$Chromosome)

# Create GenomicRanges object for snp2
snp2_ranges <- GRanges(seqnames = snp2$Chromosome, ranges = IRanges(start = snp2$PhysPos, end = snp2$PhysPos))



# Loop through annotated dataframes and annotate SNPs in snp2
for (df_name in names(annotated_data)) {
  df <- annotated_data[[df_name]]
  df$chr <- as.character(df$chr)
  
  # Create GenomicRanges object for the current dataframe
  df_ranges <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$stop))
  
  # Check for overlaps between SNP positions and genomic ranges
  overlaps <- findOverlaps(snp2_ranges, df_ranges)
  
  # Extract scores associated with overlapping regions
  overlapping_scores <- df$triage_discordance_score[subjectHits(overlaps)]
  
  # Create column name for overlap
  overlap_col <- paste0("quantitative_", df_name)
  
  # Initialize the 'overlap_' column with 0
  snp2[[overlap_col]] <- 0
  
  # Assign scores to the 'overlap_' columns in snp2 where there are overlaps
  snp2[queryHits(overlaps), (overlap_col) := ifelse(overlapping_scores == 0, 0, overlapping_scores)]
}






head(snp2)
result_dt <- data.table(column_name = character(), num_na = integer(), num_zero = integer(), num_nonzero = integer())

# Iterate over each overlap column
for (col in overlap_columns) {
  # Count the number of NA's in the column
  num_na <- sum(is.na(snp2[[col]]))
  
  # Count the number of 0 values in the column
  num_zero <- sum(snp2[[col]] == 0, na.rm = TRUE)
  
  # Count the number of non-zero values in the column
  num_nonzero <- sum(snp2[[col]] != 0, na.rm = TRUE)
  
  # Append results to result_dt
  result_dt <- rbind(result_dt, data.table(column_name = col, num_na = num_na, num_zero = num_zero, num_nonzero = num_nonzero))
}




    num_na num_zero num_nonzero
     <int>    <int>       <int>
 1:      0  7288187       68331
 2:      0  7289135       67383
 3:      0  7285727       70791
 4:      0  7283340       73178
 5:      0  7281952       74566
 6:      0  7279362       77156
 7:      0  7288785       67733
 8:      0  7291255       65263
 9:      0  7275820       80698
10:      0  7281113       75405
11:      0  7283478       73040
12:      0  7278791       77727
13:      0  7281288       75230
14:      0  7281498       75020
15:      0  7291594       64924
16:      0  7288398       68120
17:      0  7276424       80094
18:      0  7281850       74668
19:      0  7282630       73888
20:      0  7286654       69864
21:      0  7274725       81793
22:      0  7281405       75113
23:      0  7282278       74240
24:      0  7287673       68845
25:      0  7278292       78226
26:      0  7291520       64998
27:      0  7282868       73650
    num_na num_zero num_nonzero





snp2$Intercept <- "1"


x <- snp2[,c(2,60,33:59)]
setnames(x, "rsid", "SNP")


 fwrite(x, file = "top_250_annotation_quantitative.txt", append = FALSE, quote =FALSE,
                                  sep = "\t", col.names=TRUE)




#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=top250
#SBATCH --account=a_palpant
#SBATCH --array=1-18

MY_DIR=/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/
config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_250/annot/top_250_annotation_quantitative.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_20/results/"
threads=4

export OMP_NUM_THREADS=$threads

# Corrected the closing parenthesis for the Rscript command
Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}${ma_file}_top250_quant', annot='$annot', log2file=TRUE)"





Plot top 20 of discordance score vs top 250 vs rna_velocity 

Plot top 250 binary vs top 250 quant 



# Doing RNA velocity first 

R

setwd("3")


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)




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


# Manual replacement mapping
replacement_map <- c(
  "overlap_V2"  = "0_early_myocardial_progenitors",
  "overlap_V3"  = "1_posterior_second_heart_field",
  "overlap_V4"  = "2_RV",
  "overlap_V5"  = "3_hematoendothelial_progenitor",
  "overlap_V6"  = "4_atrial_CM",
  "overlap_V7"  = "5_aSHF_myocardial_cells_OFT",
  "overlap_V8"  = "6_lateral_plate_mesoderm",
  "overlap_V9"  = "7_atrioventricular_canal",
  "overlap_V10" = "8_differentiating_aSHF",
  "overlap_V11" = "9_endothelial_cells",
  "overlap_V12" = "10_endothelial_cells",
  "overlap_V13" = "11_paraxial_mesoderm",
  "overlap_V14" = "12_sinus_venosus",
  "overlap_V15" = "13_endothelial_cells",
  "overlap_V16" = "14_differentiating_pSHF",
  "overlap_V17" = "15_cardiac_committed_mesodermal_progenitor",
  "overlap_V18" = "16_endothelial_cells",
  "overlap_V19" = "17_LV",
  "overlap_V20" = "18_neural_crest",
  "overlap_V21" = "19_anterior_second_heart_field",
  "overlap_V22" = "20_paraxial_mesoderm",
  "overlap_V23" = "21_posterior_second_heart_field",
  "overlap_V24" = "22_anterior_second_heart_field",
  "overlap_V25" = "23_endocardial_cells",
  "overlap_V26" = "24_hematoendothelial_progenitor",
  "overlap_V27" = "25_epicardium",
  "overlap_V28" = "26_lateral_plate_mesoderm"
)

# Perform the replacement
combined_data$Annotation <- replacement_map[as.character(combined_data$Annotation)]

# Verify the result
table(combined_data$Annotation)








library(ggplot2)

# Your ggplot code
p <- ggplot(combined_data, aes(x = Annotation, y = Enrich, color = Annotation)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = Enrich - sd(SD), ymax = Enrich + sd(SD)),
    position = position_dodge(width = 0.5),
    width = 0.25
  ) +
  labs(
    title = "Per Snp Enrichment of Heritability",
    y = "Per SNP Enrichment Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"  # Move the legend to the bottom
  ) +
  guides(fill = guide_legend(title = "Annotation")) +
  facet_wrap(~Parameter, scales = "free_y", ncol = 4, nrow = 5)

# Save the ggplot as a PDF in the "visualization" subdirectory
ggsave("visualization/your_plot.pdf", plot = p, width = 20, height = 30, units = "cm")





#heatmap


library(ggplot2)

# Create a heatmap
p <- ggplot(combined_data, aes(x = Annotation, y = Parameter, fill = Enrich)) +
  geom_tile() +
  labs(
    title = "Per Snp Enrichment of Heritability",
    x = "Annotation",
    y = "Parameter",
    fill = "Per SNP Enrichment Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom"  # Move the legend to the bottom
  ) +
  scale_fill_gradient2(low = "white", high = "red")

p


ggsave("visualization/heatmap.pdf", plot = p, width = 70, height = 40, units = "cm")



# Remove prefixes and suffixes from parameter names
combined_data$Parameter <- gsub("invnorm_|indexed_|cleaned|_bsa", "", combined_data$Parameter)

combined_data$Parameter <- gsub("_$", "", combined_data$Parameter)
# Verify the updated parameter names
unique(combined_data$Parameter)



 saveRDS(combined_data,"combine_data.RDS")











# Load the reshape2 package
library(reshape2)

# Create the matrix without aggregation, specifying the aggregation function as mean
enrichment_matrix <- dcast(combined_data, Annotation ~ Parameter, value.var = "Enrich", fun.aggregate = mean)

# Convert Annotation column to character if needed
enrichment_matrix$Annotation <- as.character(enrichment_matrix$Annotation)

# Check the structure of the matrix
str(enrichment_matrix)





data <- enrichment_matrix 
rownames(data) <- data[,1]
data <- data[,-1]


 # Load required library
library(pheatmap)

# Define a custom function to format breaks with decimal points 0 or 5


# Generate breaks
breaks <- seq(min(data), max(data), length.out = 5)
labels <- sprintf("%.1f", breaks)

pheatmap(data,
         cluster_rows = TRUE,    # Cluster rows (cell types)
         cluster_cols = TRUE,   # Do not cluster columns (parameters)
         show_rownames = TRUE,   # Show row names (cell type names)
         show_colnames = TRUE,   # Show column names (parameter names)
         main = "Heatmap of Enrichment by Parameter and Annotation",
         fontsize = 8,           # Adjust font size if needed
         angle_col = 90,         # Rotate x-axis labels by 90 degrees
         legend_breaks = breaks,
         legend_labels = labels)



saveRDS(data, "heatmap_data.RDS")










# top_20 discordance score anno 

"/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/divided_files/top_20/results/top_20_Anno"


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


combined_data$Annotation <- sub("overlap_top_20_ordered_genes_", "", combined_data$Annotation)

# Check the updated combined_data
combined_data


combined_data <- combined_data %>%
  filter(Annotation != "Intercept")




combined_data <- combined_data %>%
  mutate(Annotation = factor(Annotation, levels = unique(Annotation[order(as.numeric(gsub("_[^_]+$", "", Annotation)))])))



p <- ggplot(combined_data, aes(x = Annotation, y = Enrich, color = Annotation)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = Enrich - SD, ymax = Enrich + SD),
    position = position_dodge(width = 0.5),
    width = 0.25
  ) +
  labs(
    title = "Per SNP Enrichment of Heritability",
    y = "Per SNP Enrichment Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"  # Move the legend to the bottom
  ) +
  guides(fill = guide_legend(title = "Annotation")) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 4, nrow = 5)

ggsave("top_20_genes.pdf", plot = p, width = 20, height = 30, units = "cm")





########### p value from enrichment data 

x <- readRDS("combine_data.RDS")

# Load necessary libraries if not already loaded
library(dplyr)

# Full dataframe x
# Assuming x is already loaded with the necessary data

# Calculate z-scores
x$z_score <- (x$Enrich - 1) / x$SD

# Calculate one-tailed p-values
x$one_tail_pvalue <- pnorm(x$z_score, lower.tail = FALSE)

# Calculate two-tailed p-values
x$two_tail_pvalue <- 2 * pnorm(abs(x$z_score))

# Print the dataframe with z-scores and p-values
print(x)



x_sorted <- x[order(x$two_tail_pvalue), ]

saveRDS(x, "combine_data.RDS")




x <- readRDS("combined_data.RDS")
combined_data <- x

# heatmap 

# Load required library



# Load the reshape2 package
library(reshape2)

# Create the matrix without aggregation, specifying the aggregation function as mean
pvalue_matrix <- dcast(combined_data, Annotation ~ Parameter, value.var = "one_tail_pvalue", fun.aggregate = mean)

# Convert Annotation column to character if needed
pvalue_matrix$Annotation <- as.character(pvalue_matrix$Annotation)

# Check the structure of the matrix
str(pvalue_matrix)

asterisk_data <- pvalue_matrix 
rownames(asterisk_data) <- asterisk_data[,1]
asterisk_data <- asterisk_datcd a[,-1]





# Generate breaks
breaks <- seq(min(data), max(data), length.out = 5)
labels <- sprintf("%.1f", breaks)

pheatmap(df,
            # Cluster rows (cell types)
       # Do not cluster columns (parameters)
         show_rownames = TRUE,   # Show row names (cell type names)
         show_colnames = TRUE,   # Show column names (parameter names)
         main = "Heatmap of Enrichment by Parameter and Annotation",
         fontsize = 8,           # Adjust font size if needed
         angle_col = 90,         # Rotate x-axis labels by 90 degrees
      )



saveRDS(data, "heatmap_data.RDS")


data_matrix <- dcast(combined_data, Annotation ~ File, value.var = "Enrich", fun.aggregate = mean)

pvalue_matrix <- dcast(combined_data, Annotation ~ File, value.var = "one_tail_pvalue", fun.aggregate = mean)


library(circlize)

col_fun = colorRamp2(c(0,), c("blue", "red"))

library(ComplexHeatmap)
Heatmap(name= "mat", col = col_fun, cell_fun = function(j, i, x, y, w, h, fill) {
    if(mat2[i, j] < 0.001) {
        grid.text("***", x, y)
    } else if(mat2[i, j] < 0.01) {
        grid.text("**", x, y)
    }else if(mat2[i, j] < 0.05) {
        grid.text("*", x, y)
    ....
}

library(ComplexHeatmap)

Heatmap(mat, name = "mat", 
    cell_fun = function(j, i, x, y, w, h, fill) {
        if (mat2[i, j] < 0.001) {
            grid.text("***", x, y)
        } else if (mat2[i, j] < 0.01) {
            grid.text("**", x, y)
        } else if (mat2[i, j] < 0.05) {
            grid.text("*", x, y)
        } else {
            grid.text("", x, y)  # Add any additional conditions here
        }
    }
)






"/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/top_100_annot_quant_and_binary.txt"




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
annot="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/top_50_annot_weighted_and_bin.txt"
out_prefix="/QRISdata/Q5149/computation/dubois_development2022/triage_discordance_score/ordered_genes_files/result_50/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_50_quant_bin', annot='$annot', log2file=TRUE)"






<svg id="Layer_3" data-name="Layer 3" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 17.35 7.84"><defs><style>.cls-1,.cls-2{fill:#e5d8bd;stroke:#231f20;stroke-miterlimit:10;}.cls-1{stroke-width:0.93px;}.cls-2{stroke-width:0.77px;}</style></defs><line class="cls-1" x1="8.93" y1="0.17" x2="8.93" y2="7.84"/><line class="cls-2" y1="0.39" x2="17.35" y2="0.39"/></svg>





normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply the normalization function to each data frame in the list
normalized_data <- lapply(annotated_data, function(df) {
  df$triage_discordance_score <- normalize(df$triage_discordance_score)
  return(df)
})









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
annot="/QRISdata/Q5149/TF_DamID/To_format/gene_list/all_common_tfboundgene_snps.txt"
out_prefix="/QRISdata/Q5149/TF_DamID/M2/annot/results/"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', outPrefix='${out_prefix}${ma_file}_allcommon', annot='$annot', log2file=TRUE)"


