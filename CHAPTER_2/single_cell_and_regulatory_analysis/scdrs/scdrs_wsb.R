# =============================================================================
# EXAMINER GUIDE
# Script: scdrs_wsb.R
# Purpose: Runs or prepares scDRS analysis on the WSB single-cell dataset; source includes HPC module/session commands.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# WSB single-cell workflow: loads the Seurat object, converts it for Python/scDRS, maps mouse genes to human orthologues and prepares trait-score visualisation.
#

# HISTORICAL SOURCE CAVEAT: This file contains interactive HPC allocation, HPC module command.
# These commands/output are retained for provenance; see docs/SCRIPT_AUDIT.md.

########## scDRS WSB 

#load hd5fr module on terminal before r package if you need to use seuratdisk 
# HPC / shell command from the interactive source session (not R code): module load hdf5/1.14.0-gompi-2023a

# HPC / shell command from the interactive source session (not R code): salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=60GB --job-name=TEST --time=01:00:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

# The original session entered the interactive R console at this point.

library(Seurat)
library(SeuratData)
library(SeuratDisk)


getwd()
#[1] "/QRISdata/Q4879/WSB_Seurat_integrated"

wsb_seurat <- readRDS("/QRISdata/Q4879/WSB_Seurat_integrated/common_orig_log_FPKM2L_integrated.RDS")

class(wsb_seurat)
#[1] "Seurat"
attr(,"package")
#[1] "Seurat"


wsb_seurat<- UpdateSeuratObject(wsb_seurat)

wsb_seurat
	An object of class Seurat 
	43711 features across 4279 samples within 2 assays 
	Active assay: integrated (13656 features, 13656 variable features)
	 2 layers present: data, scale.data
	 1 other assay present: RNA
	 2 dimensional reductions calculated: pca, umap






SaveH5Seurat(wsb_seurat, filename = "/QRISdata/Q5149/computation/WSB_dataset/wsb.h5Seurat")

Convert("/QRISdata/Q5149/computation/WSB_dataset/wsb.h5Seurat", dest = "h5ad")


# /QRISdata/Q5149/computation/WSB_dataset/ WSB_gene_number.cov generated in jupyter WSB_testing 


# HPC / shell command from the interactive source session (not R code): module load python/3.10.4-gcccore-11.3.0


# HPC / shell command from the interactive source session (not R code): scdrs compute-score \
    --h5ad-file /QRISdata/Q5149/computation/WSB_dataset/wsb.h5ad\
    --h5ad-species mouse\
    --gs-file nicktraits.gs\
    --gs-species human\
    --out-folder res\
    --cov-file /QRISdata/Q5149/computation/WSB_dataset/WSB_gene_number.cov\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True




###### scaled the dataset to have no negatives (consider running harmony next time)


# Reference URL recorded in the original workflow: https://portals.broadinstitute.org/harmony/articles/quickstart.html







adata <- anndata::read_h5ad("/scratch/user/uqsnegi/scDRS/scdrs/data/wsb_scaled.h5ad")



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

mouse_genes<- rownames(wsb_seurat@assays$RNA)


# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(mouse_genes) {
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  # Query Ensembl using the main server
  genesV2 <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
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
  if (length(successful_indices) > 100) {
    selected_indices <- successful_indices[1:100]
  } else {
    selected_indices <- successful_indices
  }
  
  # Subset the data to include only the first 50 successful conversions
  data_subset <- data[selected_indices, ]
  data_subset <- human_genes[selected_indices]
  
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

