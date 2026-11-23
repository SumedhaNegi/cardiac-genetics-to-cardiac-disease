 salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=200GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l



#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=500GB
#SBATCH --time=5:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=download
#SBATCH --account=a_palpant



# Set the directory where you want to download the files
output_dir="/QRISdata/Q5149/sumstats/raw_files/"

urls=(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000291/GCST90000291_GRCh37.tsv"
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000292/GCST90000292_GRCh37.tsv"
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000293/GCST90000293_GRCh37.tsv"
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000294/GCST90000294_GRCh37.tsv"
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000295/GCST90000295_GRCh37.tsv"
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000296/GCST90000296_GRCh37.tsv"
)

# Load the necessary module (if required)
# module load wget

# Download files in parallel using GNU Parallel
parallel --jobs 4 wget -P "$output_dir" ::: "${urls[@]}"

echo "Download complete!"


#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=500GB
#SBATCH --time=5:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=download
#SBATCH --account=a_palpant



# Set the directory where you want to download the files
output_dir="/QRISdata/Q5149/sumstats/raw_files/"

urls=(

"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278508/GCST90278508.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278509/GCST90278509.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278510/GCST90278510.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278511/GCST90278511.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278512/GCST90278512.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278513/GCST90278513.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278514/GCST90278514.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278515/GCST90278515.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278516/GCST90278516.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278517/GCST90278517.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278518/GCST90278518.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278519/GCST90278519.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90296001-GCST90297000/GCST90296096/GCST90296096.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90296001-GCST90297000/GCST90296097/GCST90296097.tsv"
"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90296001-GCST90297000/GCST90296098/GCST90296098.tsv"


parallel --jobs 4 wget -P "$output_dir" ::: "${urls[@]}"

echo "Download complete!"



# Load the necessary library
library(data.table)


# Load the necessary library
library(data.table)

# Get a list of all TSV files in the current working directory
file_list <- list.files(pattern = "\\.tsv$", full.names = TRUE)

# Function to process each TSV file
process_file <- function(file_path) {
    # Read the original summary statistics file
    df <- fread(file_path)
    
    # Reformat the data
  # Reformat the data
    df_formatted <- df[, .(SNP = MarkerName,
                           A1 = Allele1,
                           A2 = Allele2,
                           freq = Freq1,
                           b = Effect,
                           se = StdErr,
                           p = p_value,
                           N = 19260)] 

output_file <- paste0("formatted_", tools::file_path_sans_ext(basename(file_path)), ".ma")
    fwrite(df_formatted, output_file, sep = "\t", quote = FALSE)
}

       
  # Use the provided value for N
    # Write to a new file with the same name but prefixed with "formatted_"
    output_file <- paste0("formatted_", basename(file_path))
    fwrite(df_formatted, output_file, sep = "\t", quote = FALSE)
}







# Process each TSV file in the current working directory
for (file_path in file_list) {
    process_file(file_path)
}

print("Conversion complete!")



process_file <- function(file_path) {
    # Read the original summary statistics file
    df <- fread(file_path)
    
    # Reformat the data
 



#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=fractal
#SBATCH --account=a_palpant
#SBATCH --array=9

MY_DIR=/QRISdata/Q5149/sumstats/raw_files/fractal_structure_heart_Meyer_2020/formatted_sumstats

config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/TF_DamID/annotation/annot.txt"
out_prefix="/QRISdata/Q5149/sumstats/raw_files/fractal_structure_heart_Meyer_2020/results"
threads=4

export OMP_NUM_THREADS=$threads


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}${ma_file}_TF_', annot='$annot', log2file=TRUE)"












#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=fractal
#SBATCH --account=a_palpant
#SBATCH --array=1-9

MY_DIR=/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023


config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/TF_DamID/annotation/annot.txt"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${ma_file}_tidy.ma', log2file=TRUE)"

Rscript -e "SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$ld_folder', \
                  output='${ma_file}_imp.ma', log2file=TRUE)"


Rscript -e "SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$ld_folder', \
                outPrefix='${ma_file}_result', log2file=TRUE)"  











# HOPX bound and regulated genes annotation 



salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=500GB --job-name=TEST --time=00:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l
module load r/4.2.1-foss-2021a

# let's annotate top 20, top 100, top 500 and run heritability.

# Get list of text files in the directory
file_list <- list.files(pattern = "\\.txt$", full.names = TRUE)

# Load the geneMatrix data
mtrix <- fread("/QRISdata/Q5149/tmem88_heritability/annot/geneMatrix/geneMatrix.tsv")
m <- mtrix[, c(1, 12, 13, 14, 15)]
colnames(m) <- c("gene_id", "chr", "start", "stop", "strand")
annotated_data <- list()


# Initialize the list of annotated data
annotated_data <- list()
# Define a function to process each file
process_file <- function(file_path, file_name) {
  # Read the gene list from the file
  genelist <- fread(file_path)

  # Rename the first column to "gene"
  names(genelist)[1] <- "gene"

  # Convert gene names to uppercase
  genelist$gene <- toupper(genelist$gene)

  # Convert gene names in m to uppercase
  m$gene_id <- toupper(m$gene_id)

  # Merge gene list with m
  genelist1 <- merge(x = genelist, y = m, by.x = "gene", by.y = "gene_id", all.x = TRUE)

  # Modify start and stop positionshea
  genelist1$start <- genelist1$start - 50000
  genelist1$stop <- genelist1$stop + 50000

  # Remove the './' prefix from file name
  file_name <- sub("^./", "", file_name)

  # Print the first few rows of the annotated data
  print(head(genelist1))

  # Count the number of rows with at least one NA value
  num_na_rows <- sum(rowSums(is.na(genelist1)) > 0)
  cat("Number of rows with at least one NA value:", num_na_rows, "\n")

  # Remove rows with NA values
  genelist1 <- genelist1[complete.cases(genelist1), ]

  # Prepare data for annotation
  x <- genelist1  # Keep all columns from genelist1

  # Return annotated data with file name
  return(list(file_name = file_name, annotated_data = x))
}

# Process each file using lapply
annotated_data_list <- lapply(file_list, function(file) {
  file_name <- basename(file)
  process_file(file, file_name)
})

# Convert the list to a named list using file names
annotated_data <- setNames(annotated_data_list, sapply(annotated_data_list, "[[", "file_name"))



# Define a function to remove "chr" from the chr column
remove_chr_prefix <- function(data) {
  data$chr <- gsub("^chr", "", data$chr)
  return(data)
}

# Use lapply to apply the function to each data.table in the annotated_data list
annotated_data_clean <- lapply(annotated_data, function(entry) {
  entry$annotated_data <- remove_chr_prefix(entry$annotated_data)
  return(entry)
})




snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


	dim(snp)
	#[1] 7356518      10
	snp2 <- snp[,c(1,2,5,6,7)]
	colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")

# Convert 'chr' column in each dataframe to character format and create GenomicRanges object
for (df_name in names(annotated_data_clean)) {
  # Access the annotated_data element within the list
  df <- annotated_data_clean[[df_name]]$annotated_data
  
  # Convert 'chr' column to character format
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

result_dt <- data.table(column_name = character(), table_result = list())

# Iterate over each overlap column
for (col in overlap_columns) {
  # Apply the table function to the column
  table_result <- table(snp2[[col]], useNA = "ifany")
  
  # Append results to result_dt
  result_dt <- rbind(result_dt, data.table(column_name = col, table_result = list(table_result)))
}

result_dt



# ADD INTERCEPT
x <- snp2[,c(2,7:9)]

x
setnames(x, "rsid", "SNP")


 fwrite(x, file = "/QRISdata/Q5149/computation/heritability/hopx_annot/hopxdamid_annot.txt", append = FALSE, quote =FALSE,
                                  sep = "\t", col.names=TRUE)



/QRISdata/Q5149/computation/heritability/hopx_annot








#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=110GB
#SBATCH --time=4:30:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=lv_wall
#SBATCH --account=a_palpant
#SBATCH --array=1-15

MY_DIR=/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/imputed_sumstats


config=config.txt

cd $MY_DIR

module load r/4.2.1-foss-2021a

# Extract the ma_file_name for the current $SLURM_ARRAY_TASK_ID
ma_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ld_folder="/QRISdata/Q5149/snps_7M"
annot="/QRISdata/Q5149/computation/heritability/hopx_annot/hopxdamid_annot.txt"
threads=4

export OMP_NUM_THREADS=$threads

Rscript -e "SBayesRC::sbayesrc(mafile='$ma_file', LDdir='$ld_folder', \
                outPrefix='${ma_file}_hopx', annot='$annot', log2file=TRUE)"  




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
asterisk_data <- asterisk_data[,-1]





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

col_fun = colorRamp2(c(0,8), c("blue", "red"))


Heatmap(name= "mat", col = col_fun, cell_fun = function(j, i, x, y, w, h, fill) {
    if(mat2[i, j] < 0.001) {
        grid.text("***", x, y)
    } else if(mat2[i, j] < 0.01) {
        grid.text("**", x, y)
    }else if(mat2[i, j] < 0.05) {
        grid.text("*", x, y)
    ....
}


Heatmap(mat, name = "Heatmap", col=col_fun,
    cell_fun = function(j, i, x, y, w, h, fill) {
        if (mat2[i, j] < 0.001) {
            grid.text("***", x, y)
        } else if (mat2[i, j] < 0.01) {
            grid.text("**", x, y)
        } else if (mat2[i, j] < 0.05) {
            grid.text("*", x, y)
        } else {
            grid.text("ns", x, y)  # "ns" for p-values above 0.05
        }







#mapping the GWAS TRAIT t
#"/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/imputed_sumstats/AnnoPerSnpHsq"

# Load the dplyr library
library(dplyr)

# Directory containing the files
files_directory <- "path/to/your/files/directory"

# Get a list of filenames from the directory
file_list <- list.files()

# Extract accession IDs from filenames
accession_ids <- sub(".*_(GCST\\d+)_.*", "\\1", file_list)

# Extract accession IDs from filenames

# Dictionary mapping accession IDs to reported traits
3    "90278508" = "Left ventricular end systole inferoseptal wall thickness",
    "90278509" = "Left ventricular end systole inferior wall thickness",
    "90278510" = "Left ventricular end systole inferolateral wall thickness",
    "90278511" = "Left ventricular end systole anterolateral wall thickness",
    "90278512" = "Left ventricular end systole anterior wall thickness",
    "90278513" = "Left ventricular end systole anterospetal wall thickness",
    "90278514" = "Left ventricular end diastole inferoseptal wall thickness",
    "90278515" = "Left ventricular end diastole inferior wall thickness",
    "90278516" = "Left ventricular end diastole inferolateral wall thickness",
    "90278517" = "Left ventricular end diastole anterolateral wall thickness",
    "90278518" = "Left ventricular end diastole anterior wall thickness",
    "90278519" = "Left ventricular end diastole anterospetal wall thickness",
    "90296096" = "Hypertrophic cardiomyopathy",
    "90296097" = "Dilated cardiomyopathy",
    "90296098" = "Pulmonary hypertension"
    # Add more mappings as needed
)

# Get reported trait names using accession IDs
reported_traits <- sapply(accession_ids, function(id) accession_to_trait[[id]])

# Create a data frame with reported trait names and accession IDs
file_data <- data.frame(
    filename = file_list,
    accessionId = accession_ids,
    reportedTrait = reported_traits,
    stringsAsFactors = FALSE
)

# Print the file data
print(file_data)



#HEATMAP 

 modified_files <- list.files(pattern = "^formatted_", full.names = TRUE)

file_labels <- gsub("^formatted_(.*?)\\.tsv.*", "\\1", basename(modified_files))


data_list <- list()

# Read each file and store it in the list
for (file in file_labels) {
    # Read the file
    data <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
    # Store the data frame in the list
    data_list[[file]] <- data
}



transposed_data <- t(combined_data)

# Create the heatmap using pheatmap
pheatmap(transposed_data, cluster_rows = FALSE, cluster_cols = FALSE)




# Load necessary libraries
library(dplyr)
library(reshape2)

# Directory containing the files
files_directory <- "path/to/your/files/directory"

# Get list of files
file_list <- list.files(pattern = "^formatted_", full.names = TRUE)

# Read data from all files and combine into a single dataframe
combined_data <- lapply(file_list, read.table, header = TRUE, sep = "\t")
combined_data <- bind_rows(combined_data, .id = "File")

# Calculate z-scores
combined_data$z_score <- (combined_data$Enrich - 1) / combined_data$SD

# Calculate one-tailed p-values
combined_data$one_tail_pvalue <- pnorm(combined_data$z_score, lower.tail = FALSE)

# Calculate two-tailed p-values
combined_data$two_tail_pvalue <- 2 * pnorm(abs(combined_data$z_score))

# Reshape data to create a matrix of p-values
pvalue_matrix <- dcast(combined_data, Annotation ~ Parameter, value.var = "two_tail_pvalue", fun.aggregate = mean)

# Convert Annotation column to character if needed
pvalue_matrix$Annotation <- as.character(pvalue_matrix$Annotation)

# Save the p-value matrix as an RDS file
saveRDS(pvalue_matrix, "pvalue_matrix.RDS")







# Khurshid mass
       df_formatted <- x[, .(SNP = SNP,
                           A1 = ALLELE0,
                           A2 = ALLELE1,
                           freq = A0FREQ,
                           b = BETA,
                           se = SE,
                           p = P_BOLT_LMM,
                           N = 32239 )]






