############ generating a .cov file 


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:00:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l
#python3

import scanpy as sc
import pandas as pd

# Load the .h5ad file
adata = sc.read_h5ad("/scratch/user/uqsnegi/scDRS/scdrs/data/edited_gene_names.h5ad")

# Extract cell names
cell_names = adata.obs_names

# Calculate the number of genes for each cell
n_genes = (adata.X > 0).sum(axis=1).A1  # Convert to 1-dimensional array

# Create a DataFrame
cov_df = pd.DataFrame({
    'cell_names': cell_names,
    'n_genes': n_genes
})

# Save the DataFrame as a .tsv file
cov_df.to_csv('/scratch/project/palpant_scratch/SNegi/edit_gene_names.cov', sep='\t', index=False)



################


#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=60GB
#SBATCH --time=06:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=test
#SBATCH --account=a_palpant



module load python/3.10.4-gcccore-11.3.0


H5AD_FILE=/scratch/user/uqsnegi/scDRS/scdrs/data/edited_gene_names.h5ad
COV_FILE=/scratch/project/palpant_scratch/SNegi/edit_gene_names.cov

# zscore
GS_FILE=/scratch/project/palpant_scratch/SNegi/hopx_scdrs/filtered_traits.gs
OUT_FOLDER=/scratch/project/palpant_scratch/SNegi/hopx_scdrs/individual



scdrs compute-score \
    --h5ad_file $H5AD_FILE\
    --h5ad_species human\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter_data True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER




#array script 
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=100GB
#SBATCH --time=05:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=test
#SBATCH --account=a_palpant
#SBATCH --array=1-75
# Load necessary 
module load python/3.10.4-gcccore-11.3.0

# Path to your trait_info.tsv file
TRAIT_INFO_FILE="/scratch/project/palpant_scratch/SNegi/hopx_scdrs/traitlist_cleaned.tsv"                                                                                                                                                                       
# Read the trait identifier from trait_info.tsv based on the SLURM_ARRAY_TASK_ID
TRAIT_IDENTIFIER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TRAIT_INFO_FILE}")

# Construct the full path to the score file
SCORE_FILE="/scratch/project/palpant_scratch/SNegi/hopx_scdrs/individual/${TRAIT_IDENTIFIER}.full_score.gz"

scdrs perform-downstream \
--h5ad-file /scratch/user/uqsnegi/scDRS/scdrs/data/edited_gene_names.h5ad \
--score-file "${SCORE_FILE}" \
--out-folder /scratch/project/palpant_scratch/SNegi/hopx_scdrs/group/ \
--group-analysis new_clusters \
--corr-analysis None \
--gene-analysis True \
--flag-filter-data True \
--flag-raw-count True



grep "HOPX" *scdrs_gene




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=100GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)
library(Seurat)
getwd()
[1] "/QRISdata/Q4879/scRNAseq_Data"



harmony <- readRDS('harmony_pc20_norm_seu_240619.RDS')
harmony[["RNA"]] <- JoinLayers(harmony[["RNA"]])
DefaultAssay(harmony) <- 'RNA'

library(Seurat)

# Set the cell identities to a specific resolution
# Here we choose the resolution 0.6 as an example
Idents(harmony) <- harmony@meta.data$SCT_snn_res.0.6

# Visualize the UMAP with legends and labels
DimPlot(harmony, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 12))






library(Seurat)

# Set the cell identities to the 'sample' column
Idents(harmony) <- harmony@meta.data$sample

# Visualize the UMAP with legends and labels by sample
DimPlot(harmony, reduction = "umap", label = TRUE, group.by = "sample", pt.size = 0.5) + 
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 12))







library(Seurat)

# Visualize the expression of the gene "HOPX" on the UMAP
FeaturePlot(harmony, features = "HOPX--ENSG00000171476", reduction = "umap", pt.size = 0.5) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 12))






library(Seurat)
library(ggplot2)

# Define the gene of interest
gene_of_interest <- "A# Extract the feature names
features <- rownames(harmony[["RNA"]])

# Remove the "--ENSG000001XXXX" part
clean_features <- sub("--ENSG[0-9]+", "", features)

# Update the feature names in the Seurat object
rownames(harmony[["RNA"]]) <- clean_features

# Verify the changes
head(rownames(harmony[["RNA"]]))
"

# Generate the violin plot for the gene of interest without dots
p <- VlnPlot(harmony, features = gene_of_interest, pt.size = 0, combine = FALSE) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Print the plot
print(p)

# Save the plot to a file
ggsave("ACTC1_violin_plot_no_dots.png", plot = p, width = 8, height = 6, dpi = 300)
library(Seurat)
library(ggplot2)

# Define the gene of interest
gene_of_interest <- "ACTC1--ENSG00000159251"

# Generate the violin plot for the gene of interest without dots
p <- VlnPlot(harmony, features = gene_of_interest, pt.size = 0, combine = FALSE) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Print the plot
print(p)

# Save the plot to a file
ggsave("ACTC1_violin_plot_no_dots.png", plot = p, width = 8, height = 6, dpi = 300)






#module load hdf5/1.14.0-gompi-2023a

# Install and load SeuratDisk package
install.packages("SeuratDisk")
library(Seurat)
library(SeuratDisk)

# Load your Seurat object
harmony <- readRDS("harmony_seurat_object.rds")

# Save Seurat object as .h5seurat file
SaveH5Seurat(harmony, filename = "harmony_seurat_object.h5seurat")

# Convert .h5seurat file to .h5ad file
Convert("harmony_seurat_object.h5seurat", dest = "h5ad")












salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=150GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)








############################## MAGMA #####################



./magma \
    --bfile g1000_eur \
    --pval /QRISdata/Q5149/sumstats/raw_files/francis_hcm_dcm/formatted_meanLVwallthickness_metal.ma use='SNP,p' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/meanLVwallthickness_metal





/QRISdata/Q5149/computation/magma
/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/*ma

    mkdir -p out/step2
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${trait}.pval use='SNP,p' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/${trait}














pwd 


/scratch/project/palpant_scratch/Multiplex/dalia_outputs/Label_Transfer/Multiplex_Dif_Data_Full_110924.h5ad 




scdrs compute-score \
    --h5ad-file /scratch/project/palpant_scratch/Multiplex/dalia_outputs/Label_Transfer/Multiplex_Dif_Data_Full_110924.h5ad \
    --h5ad-species human\
    --gs-file nicktraits.gs\
    --gs-species human\
    --out-folder res\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True








#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=400GB
#SBATCH --time=06:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=test
#SBATCH --account=a_palpant



module load python/3.10.4-gcccore-11.3.0


H5AD_FILE=/scratch/project/palpant_scratch/Multiplex/dalia_outputs/Label_Transfer/Multiplex_Dif_Data_Full_110924.h5ad 


# zscore
GS_FILE=/scratch/project/palpant_scratch/SNegi/hopx_scdrs/filtered_traits.gs
OUT_FOLDER=/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test



scdrs compute-score \
    --h5ad_file $H5AD_FILE\
    --h5ad_species human\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter_data True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER












######### make cov file ~ 
######## adjust for cell type proportions??





salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=300GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)



#SBATCH --array=1-36      # Set the job array size (adjust this based on the number of files)


#!/bin/bash
#SBATCH --job-name=magma_array
#SBATCH --output=logs/magma_%A_%a.out  # Save the output log
#SBATCH --error=logs/magma_%A_%a.err   # Save the error log
#SBATCH --ntasks=4
#SBATCH --time=01:00:00  # Adjust time limit as needed
#SBATCH --mem=8G         # Adjust memory based on your needs
#SBATCH --account=a_palpant
#SBATCH --partition=general


# Directories
sumstats_dir="/scratch/project/palpant_scratch/SNegi/sumstats"
magma_dir="/QRISdata/Q5149/computation/magma"
output_dir="out/step2"

# Create output directory
mkdir -p ${output_dir}

# Get list of all *.tsv_tidy.ma files
files=(${sumstats_dir}/*.ma)

# Get the file corresponding to the current array task
file=${files[$SLURM_ARRAY_TASK_ID]}

# Extract the trait name from the filename
trait=$(basename ${file}.ma)

# Run MAGMA
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${file} use='SNP,p' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out ${output_dir}/${trait}






salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)

setwd("/QRISdata/Q5149/computation/magma/out/step2/")

library(dplyr)

zscore_files <- list.files(pattern = "out")

# Initialize an empty list to store the data frames
df_list <- list()

# Loop through each file
for (file in zscore_files) {
  # Read the file
  df <- read.table(file, header = TRUE, sep = "")
  
  # Select only the first and eighth qcolumns
  selected_columns <- df %>% select(GENE, col_8th = colnames(df)[8])
  
  # Rename the 8th column based on the file name (or another unique identifier)
  colnames(selected_columns)[2] <- paste0("zscore_", sub(".txt", "", file))  # Remove file extension for column name
  
  # Add the data frame to the list
  df_list[[file]] <- selected_columns
}

# Merge all data frames by the first column
merged_df <- Reduce(function(x, y) merge(x, y, by = "GENE", all = TRUE), df_list)


  write.table(merged_df, file = "120924_merged_zscore_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)





module load python/3.10.4-gcccore-11.3.0

    # Select top 1,000 genes and use z-score weights
scdrs munge-gs \
    --out-file 120924_cardiac_traits.gs \
    --zscore-file /QRISdata/Q5149/computation/magma/out/step2/120924_cleaned_zscore_data.tsv\
    --weight zscore \
    --n-max 1000



import pandas as pd
import numpy as np
import os

# Path to the original .gs file
gs_file = 'magma_10kb_top1000_zscore.74_traits.rv1.gs'

# Load the data
data = pd.read_csv(gs_file, sep='\t')

# Number of batches
num_batches = 74

# Split the data into batches
batches = np.array_split(data, num_batches)

# Save each batch to a separate file
for i, batch in enumerate(batches):
    batch_file = os.path.join(f'magma_scdrs_paper_batch{i}.gs')
    batch.to_csv(batch_file, sep='\t', index=False)
    print(f'Batch {i} saved to {batch_file}')







###############################################################################

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=80GB
#SBATCH --time=02:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=test
#SBATCH --account=a_palpant
#SBATCH --array=0-52



module load python/3.10.4-gcccore-11.3.0

BATCH_NUM=$SLURM_ARRAY_TASK_ID


H5AD_FILE=/scratch/project/palpant_scratch/Multiplex/dalia_outputs/Label_Transfer/Multiplex_Dif_Data_Full_110924.h5ad 
COV_FILE=/scratch/project/palpant_scratch/SNegi/multiplex_study/multiplex_120924.cov 


# zscore
GS_FILE=/scratch/project_mnt/S0138/SNegi/gs_files/magma_scdrs_paper_batch$BATCH_NUM.gs
OUT_FOLDER=/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test



scdrs compute-score \
    --h5ad_file $H5AD_FILE\
    --h5ad_species human\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter_data True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER



salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=200GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l






import anndata as ad
import pandas as pd

# Load the H5AD file
adata = ad.read_h5ad('/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test/multiplex_120924_umap.h5ad')
# Load the cell identities from the text file
cell_ids = pd.read_csv('/scratch/project_mnt/S0138/Multiplex/Cardiac_Cell_IDs_scGPT_Cluster_Overlap_130924.txt', 
                       delim_whitespace=True, header=None, skiprows=1, names=['x', 'CellID'])

# Check if the cell IDs match with the H5AD object's obs dataframe
print(adata.obs.head())
print(cell_ids.head())

# Merge the cell identity information into the adata.obs dataframe
adata.obs = adata.obs.merge(cell_ids[['CellID']], left_index=True, right_on='CellID', how='left')

# Save the annotated H5AD file
adata.write('/scratch/project/palpant_scratch/Multiplex/overlapping_annotated.h5ad')







/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test/multiplex_120924_umap.h5ad



############### group #############


#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=25GB
#SBATCH --time=01:00:00
#SBATCH -o group_out_%x_%j.txt
#SBATCH -e group_error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=test
#SBATCH --account=a_palpant
#SBATCH --array=1-52


module load python/3.10.4-gcccore-11.3.0

score_dir="/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test/full_scores"
H5AD_FILE="/scratch/project/palpant_scratch/Multiplex/scDRS/umap_annotation_and_connectivities.h5ad"
score_files=($(ls "$score_dir"/*.full_score.gz))
score_file="${score_files[$SLURM_ARRAY_TASK_ID]}"
file_out_folder="/scratch/project/palpant_scratch/Multiplex/scDRS/group"
mkdir -p "$file_out_folder"

echo "Processing file: $score_file"

# Perform downstream analysis
scdrs perform-downstream \
    --h5ad-file "$H5AD_FILE" \
    --score-file "$score_file" \
    --out-folder "$file_out_folder" \
    --group-analysis Cardiomyocytes \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False

echo "Finished processing file: $score_file"




df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=dict_df_score["SCZ"],
    group_cols=["level1class"],
)["level1class"]

display(df_stats.style.set_caption("Group-level statistics for SCZ"))



# Set your working directory to the subdirectory containing the files
directory = "/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test"  # Change this to your directory path
file_pattern = "*.full_score.gz"  # File pattern to match full_score.gz files

adata = sc.read_h5ad("multiplex_120924_umap.h5ad")

# Loop through all files that match the pattern in the directory

import pandas as pd
import scanpy as sc
import os
import glob
import subprocess

# Set your working directory to the subdirectory containing the files
directory = "/scratch/project_mnt/S0138/Multiplex/scDRS/individual_test"  # Change this to your directory path
file_pattern = "*.full_score.gz"  # File pattern to match full_score.gz files

# Loop through all files that match the pattern in the directory
for file_path in glob.glob(os.path.join(directory, file_pattern)):
# Extract the trait name from the file name (without the full_score.gz extension)
trait = os.path.basename(file_path).replace(".full_score.gz", "")

try:
# Read the file into a DataFrame
score_df = pd.read_csv(file_path, sep="\t", index_col=0)

# Check if the DataFrame is empty
if score_df.empty:
print(f"Warning: File {file_path} is empty.")
continue

# Add the norm_score to adata.obs if it exists
if "norm_score" in score_df.columns:
adata.obs[trait] = score_df["norm_score"]
else:
print(f"Warning: 'norm_score' column missing in file {file_path}")
continue  # Skip this file if the column is missing

# Set figure parameters for visualization
sc.set_figure_params(figsize=[2.5, 2.5], dpi=150)

# Define the path for saving the plot
plot_path = os.path.join(directory, f"_{trait}_umap.png")

# Plot UMAP with the specific trait
sc.pl.umap(
adata,
color=trait,
color_map="RdBu_r",
vmin=-5,
vmax=5,
s=20,
save=f"_{trait}_umap.png"  # Save each plot as a PNG with the trait name
)
        
        # Open the plot with the default image viewer (which can be XQuartz on macOS)
    






adata.write('umap_annotation_and_connectivities.h5ad')

















library(ggplot2)
library(ggrastr)
library(dplyr)

# Load the UMAP coordinates and metadata
umap_data <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load the trait scores
trait_data <- read.table("zscore_LVESV.genes.out.score", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

 trait_data$cell_barcode <- trait_data$X
# Merge on cell_barcode
merged_data <- umap_data %>%
  left_join(trait_data %>% select(cell_barcode, norm_score), by = "cell_barcode")

# Create the UMAP plot
umap_plot <- ggplot(merged_data, aes(x = umap_1, y = umap_2, color = norm_score)) +
  geom_point(size = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "UMAP Plot Colored by Norm Score", x = "UMAP 1", y = "UMAP 2", color = "Norm Score")
# Save the rasterized plot
ggsave(filename = "umap_plot.pdf", plot = umap_plot, device = "pdf")





########### single cell level plotting


library(ggplot2)
library(ggrastr)
library(dplyr)
library(tools)  



trait_file <- "/scratch/project/palpant_scratch/Multiplex/scDRS/individual_test/full_scores/zscore_LVESV.genes.out.score"


trait_data <- read.table(trait_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract the base name of the file (without extension)
base_name <- file_path_sans_ext(basename(trait_file))

# Function to calculate convex hulls for each cluster group
# Create the UMAP plot
umap_plot <- ggplot(merged_data, aes(x = umap_1, y = umap_2, color = norm_score)) +
  geom_point(size = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                        limits = c(-5, 5), breaks = c(-5, 0, 5), labels = c("-5", "0", "+5")) +
  theme_minimal() +
  labs(title = paste("sCDRS -", base_name), 
       x = "UMAP 1", y = "UMAP 2", color = "Norm Score") +
  # Add labels for each cluster at the median position
  geom_text(data = merged_data %>%
              group_by(SCT_snn_res.0.2) %>%
              summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)),
            aes(label = SCT_snn_res.0.2), size = 3, vjust = -1, hjust = 0.5, color = "black")

# Rasterize the points layer
rasterized_plot <- rasterize(umap_plot, layers = 'Point', dpi = 300)

# Save the rasterized plot
ggsave(filename = paste0("umap_plot_with_cluster_annotations_", base_name, ".pdf"), 
       plot = rasterized_plot, device = "pdf")









##################################################################


group plotting








R
library(data.table)
library(dplyr)

# Define the directory containing your files

file_directory <- getwd()
#[1] "/scratch/project_mnt/S0138/Multiplex/scDRS/group

# List all files in the directory
file_list <- list.files(file_directory, pattern = "^zscore.*\\.Cardiomyocytes$", full.names = TRUE)


# Initialize an empty list to store data from each file
data_list <- list()

# Loop through each file and read the data
for (file in file_list) {
  data <- fread(file, header = TRUE)
  trait_name <- sub(".*/|\\.scdrs_group\\.cell_ontology_class$", "", file)
  data$trait <- trait_name
  data_list[[trait_name]] <- data
}

# Combine all data into a single data frame
combined_data <- rbindlist(data_list)

# Apply FDR correction within each trait group for assoc_mcp and hetero_mcp
combined_data[, `:=`(assoc_mcp_fdr_zscore = p.adjust(assoc_mcp, method = "fdr"),
                     hetero_mcp_fdr_zscore = p.adjust(hetero_mcp, method = "fdr")), 
              by = trait]



# Perform the join


# Remove the suffix from the trait column
combined_data[, trait := sub("\\Cardiomyocytes", "", trait)]

# Create the new column for proportion of significant cells
combined_data[, prop_sig_cells_zscore := n_fdr_0.2 / n_cell]


# Check the updated data frame with the suffix removed and new column added
head(combined_data)



library(ggplot2)


# Add a border color and dot column based on conditions
combined_data <- combined_data %>%
  mutate(
    border_color = ifelse(assoc_mcp_fdr_zscore <= 0.10, "black", "transparent"),
    dot = ifelse(hetero_mcp_fdr_zscore <= 0.10, "black", NA)
  )

# Plot the heatmap
ggplot(combined_data, aes(x = trait, y = group, fill = prop_sig_cells_zscore)) +
  geom_tile(color = combined_data$border_color) +  # Add border color
  geom_point(aes(color = dot), size = 2, na.rm = TRUE) +  # Add dots for hetero_mcp_fdr
  scale_fill_gradient(low = "blue", high = "red") +  # Adjust color scale for fill
  scale_color_identity() +  # Use the color values directly
  theme_minimal() +
  labs(x = "Trait", y = "Group", fill = "Proportion of Significant Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels






library(ggplot2)
library(ggforce)

# Prepare the data
filtered_data <- sorted_combined_data[grepl("cardio", group)]
filtered_data <- filtered_data %>%
  mutate(
    fill_color = ifelse(assoc_mcp <= 0.05, "red", "white"),  # Set fill color to red where assoc_mcp is <= 0.05
    dot = ifelse(hetero_mcp <= 0.05, "black", NA)  # Set dot color based on hetero_mcp
  )

# Plot the heatmap
ggplot(filtered_data, aes(x = trait, y = group)) +
  geom_tile(aes(fill = fill_color), color = "black") +  # Set fill color and border color
  geom_point(aes(color = dot), size = 2, na.rm = TRUE) +  # Add dots for hetero_mcp
  scale_fill_identity() +  # Use the fill values directly
  scale_color_identity() +  # Use the color values directly for dots
  theme_minimal() +
  labs(x = "Trait", y = "Group", fill = "Proportion of Significant Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels


