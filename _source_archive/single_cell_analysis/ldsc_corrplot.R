# LDSC 

######## directories 

#RV, outflow 
/QRISdata/Q5149/sumstats/Right_ventricle/ 





salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)




module load ldsc/1.0.1--pyhdfd78af_2/module

 /QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/formatted_GCST90278509.tsv

 /QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/formatted_GCST90278510.tsv






module load ldsc/1.0.1--pyhdfd78af_2/module

munge_sumstats.py \
--sumstats LV_mass.ma \
--out LV_mass \
--merge-alleles w_hm3.snplist \ 
--chunksize 500000





######### consider the need for merge-alleles 


	Reading list of SNPs for allele merge from w_hm3.snplist
	Read 1217311 SNPs for allele merge.
	Reading sumstats from formatted_GCST90278509.tsv into memory 500000 SNPs at a time.
	................. done
	Read 8456866 SNPs from --sumstats file.
	Removed 7346330 SNPs not in --merge-alleles.
	Removed 0 SNPs with missing values.
	Removed 0 SNPs with INFO <= 0.9.
	Removed 0 SNPs with MAF <= 0.01.
	Removed 0 SNPs with out-of-bounds p-values.
	Removed 2 variants that were not SNPs or were strand-ambiguous.
	1110534 SNPs remain.
	Removed 0 SNPs with duplicated rs numbers (1110534 SNPs remain).
	Removed 0 SNPs with N < 27989.3333333 (1110534 SNPs remain).
	Median value of b was -4.59517e-06, which seems sensible.
	Removed 2 SNPs whose alleles did not match --merge-alleles (1110532 SNPs remain).
	Writing summary statistics for 1217311 SNPs (1110532 with nonmissing beta) to /scratch/user/uqsnegi/ldsc/test/GCST90278509.sumstats.gz.

	Metadata:
	Mean chi^2 = 1.132
	Lambda GC = 1.097
	Max chi^2 = 62.02
	19 Genome-wide significant SNPs (some may have been removed by filtering).

	Conversion finished at Sat May 25 13:56:09 2024
	Total time elapsed: 30.41s





ldsc.py --rg GCST90278510.sumstats.gz,GCST90278509.sumstats.gz \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out corr


ldsc.py --rg invnorm_RVSV.ma.sumstats.gz,invnorm_lv_sv.ma.sumstats.gz,longcovid_vs_covid.ma.sumstats.gz \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out corr








correlation = 0.82







try with kidney size ----------> 

/QRISdata/Q5149/tmem88_heritability/traits/kidney_size.ma



 munge_sumstats.py \
--sumstats kidney_size.ma \
--out /scratch/user/uqsnegi/ldsc/test/kidney \
--merge-alleles w_hm3.snplist \
--N 38000 \
--chunksize 500000






ldsc.py --rg kidney.sumstats.gz,GCST90278509.sumstats.gz \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out corr_kid




correlation ---- 
0.2632











###################### downloading summary statistics ##########


------------> munge sumstats on these summary statistics 


#garnier dcm 

/QRISdata/Q5149/sumstats/raw_files/garnier_dcm/garner_dcm_formatted.ma 

==done

#fractal structure of the heart (trabeculae)

/QRISdata/Q5149/sumstats/raw_files/fractal_structure_heart_Meyer_2020/formatted_sumstats/*imp.ma
== done

# Ning lW wall thickness (DCM/HCM/pulomary)
/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/*ma

== half 

# hermes and right ventricular and left ventricular and outflow tract 
/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/*ma



/QRISdata/Q5149/sumstats/Right_ventricle/*tsv
	file_list <- list.files(pattern="tsv$", full.names=TRUE)
		for (file in file_list) {
			  # Read the file
			  data <- fread(file)
			  
			  # Rename the column P_BOLT_LMM to P
			  if ("P_BOLT_LMM" %in% names(data)) {
			    setnames(data, "P_BOLT_LMM", "P")
			  }
			  
			  # Append the column "N" with the value 32860
			  data[, N := 36041]
			  
			  # Write the modified data back to the file, overwriting the original
			  fwrite(data, file)
			}


#== P p bolt lmm 

# could not find A1 and A2 columns (consult cleaning_sumstat.R)

#liver and kidney and lung and covid (random?)
/QRISdata/Q5149/tmem88_heritability/traits/*cleaned.ma

		setwd("/QRISdata/Q5149/tmem88_heritability/traits")
		file_list <- list.files(pattern = "cleaned\\.ma$", full.names = TRUE)

					# Loop through each file
					for (file in file_list) {
					  # Read the file
					  data <- fread(file)
					  
					  # Append the column "N" with the value 32860
					  data[, N := 32860]
					  
					  # Write the modified data back to the file, overwriting the original
					  fwrite(data, file) # check code make sure the written file is tab delimited 
					}
/QRISdata/Q5149/computation/heritability/long_covid/formatted_summary_statistics/longcovid_vs_covid.ma





# v's traits (?)

for FILE in 
*tsv; do
 # Check if the file exists
   if [ -f "$FILE" ]; then
     # Convert the file to tab-delimited and overwrite the original file
     sed 's/,/\t/g' "$FILE" > temp && mv temp "$FILE"
   fi
done

 


####################interactive job 

#!/bin/bash

module load ldsc/1.0.1--pyhdfd78af_2/module


OUTPUT_DIR="/scratch/user/uqsnegi/ldsc/sumstat_correlations/"

# Define an array with directory paths and corresponding file patterns
declare -A FILE_PATTERNS=(
    ["/QRISdata/Q5149/sumstats/raw_files/garnier_dcm/"]="garner_dcm_formatted.ma"
    ["/QRISdata/Q5149/sumstats/raw_files/fractal_structure_heart_Meyer_2020/formatted_sumstats/"]="*imp.ma"
    ["/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/"]="*ma"
    ["/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/"]="*ma"
    ["/QRISdata/Q5149/sumstats/Right_ventricle/"]="*tsv"
    ["/QRISdata/Q5149/tmem88_heritability/traits/"]="*cleaned.ma"
    ["/QRISdata/Q5149/computation/heritability/long_covid/formatted_summary_statistics/"]="longcovid_vs_covid.ma"
    ["/QRISdata/Q5149/tmem88_heritability/traits/"]="*P.ma"
)

# Iterate over each directory and process the corresponding file pattern
for DIRECTORY in "${!FILE_PATTERNS[@]}"; do
    # Change directory to the current directory
    cd "$DIRECTORY" || { echo "Failed to change directory to $DIRECTORY"; exit 1; }

    # Extract the directory name
    DIR_NAME=$(basename "$DIRECTORY")

    # Get the corresponding file pattern for the current directory
    FILE_PATTERN="${FILE_PATTERNS[$DIRECTORY]}"

    # Iterate over files matching the file pattern
    for FILENAME in $FILE_PATTERN; do
        # Check if the current file actually exists
        if [ -f "$FILENAME" ]; then
            # Extract the filename without the path
            FILENAME=$(basename "$FILENAME")
            echo "Processing file: $FILENAME"
            # Run munge_sumstats.py for the current file
            munge_sumstats.py \
            --sumstats "$DIRECTORY""$FILENAME" \
            --out "$OUTPUT_DIR""$FILENAME" \
            --merge-alleles /scratch/user/uqsnegi/ldsc/w_hm3.snplist \
            --chunksize 500000
        fi
    done
done


###### munging francis_hcm_dcm_strain 
salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

#!/bin/bash

module load ldsc/1.0.1--pyhdfd78af_2/module


OUTPUT_DIR="/scratch/user/uqsnegi/ldsc/sumstat_correlations/"

# Define an array with directory paths and corresponding file patterns
declare -A FILE_PATTERNS=(
    ["/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/"]="*.ma"
)

# Iterate over each directory and process the corresponding file pattern
for DIRECTORY in "${!FILE_PATTERNS[@]}"; do
    # Change directory to the current directory
    cd "$DIRECTORY" || { echo "Failed to change directory to $DIRECTORY"; exit 1; }

    # Extract the directory name
    DIR_NAME=$(basename "$DIRECTORY")

    # Get the corresponding file pattern for the current directory
    FILE_PATTERN="${FILE_PATTERNS[$DIRECTORY]}"

    # Iterate over files matching the file pattern
    for FILENAME in $FILE_PATTERN; do
        # Check if the current file actually exists
        if [ -f "$FILENAME" ]; then
            # Extract the filename without the path
            FILENAME=$(basename "$FILENAME")
            echo "Processing file: $FILENAME"
            # Run munge_sumstats.py for the current file
            munge_sumstats.py \
            --sumstats "$DIRECTORY""$FILENAME" \
            --out "$OUTPUT_DIR""$FILENAME" \
            --merge-alleles /scratch/user/uqsnegi/ldsc/w_hm3.snplist \
            --chunksize 500000
        fi
    done
done





######### script where each array is a directory

#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH --mem=50GB
#SBATCH --time=01:00:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=ldsc
#SBATCH --account=a_palpant
#SBATCH --array=1-7
# Load the required module
module load ldsc/1.0.1--pyhdfd78af_2/module

OUTPUT_DIR="/scratch/user/uqsnegi/ldsc/sumstat_correlations/"

# Define arrays with directory paths and corresponding file patterns
DIR_ARRAY=(
    "/QRISdata/Q5149/sumstats/raw_files/LV_wall_thickness_Ning_2023/"
    "/QRISdata/Q5149/sumstats/formatted/from_dalia/cardiac/tf_damid_7m_nowindow/"
    "/QRISdata/Q5149/sumstats/Right_ventricle/"
    "/QRISdata/Q5149/tmem88_heritability/traits/"
    "/QRISdata/Q5149/computation/heritability/long_covid/formatted_summary_statistics/"
    "/QRISdata/Q5149/tmem88_heritability/traits/"
)

PATTERN_ARRAY=(
    "*ma"
    "*ma"
    "*tsv"
    "*cleaned.ma"
    "longcovid_vs_covid.ma"
    "*P.txt"
)

# Get the directory and file pattern for the current SLURM_ARRAY_TASK_ID
DIRECTORY="${DIR_ARRAY[$SLURM_ARRAY_TASK_ID]}"
FILE_PATTERN="${PATTERN_ARRAY[$SLURM_ARRAY_TASK_ID]}"

# Change directory to the current directory
cd "$DIRECTORY" || { echo "Failed to change directory to $DIRECTORY"; exit 1; }

# Iterate over files matching the file pattern
for FILENAME in $FILE_PATTERN; do
    # Check if the current file actually exists
    if [ -f "$FILENAME" ]; then
        # Extract the filename without the path
        FILENAME=$(basename "$FILENAME")
        echo "Processing file: $FILENAME"
        # Run munge_sumstats.py for the current file
        munge_sumstats.py \
        --sumstats "$DIRECTORY$FILENAME" \
        --out "$OUTPUT_DIR$FILENAME" \
        --merge-alleles /scratch/user/uqsnegi/ldsc/w_hm3.snplist \
        --chunksize 500000
    fi
done










###### Didn't work for the tsv files so cleaning script here first 

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10GB
#SBATCH --time=00:20:00
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --partition=general
#SBATCH --job-name=cleaning
#SBATCH --account=a_palpant
#SBATCH --array=1-34

# Load R module
R

# Define function to run R script
run_r_script() {
    Rscript my_script.R "$1"
}

# Run R script for each file in the array
run_r_script "${SLURM_ARRAY_TASK_ID}"





# my_script.R

library(data.table)
library(dplyr)

files <- list.files(pattern="tsv$", full.names=TRUE)
file_index <- commandArgs(trailingOnly = TRUE)

file <- files[as.numeric(file_index)]

print(file)
x <- fread(file, header=TRUE, sep="\t")
x_1 <- x[, c(2, 3, 1, 5, 6, 7, 11, 12, 16, 18)]
colnames(x_1) <- c("CHR", "BP", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")
filtered <- x_1 %>% 
    filter(grepl('rs', SNP))
newfilename <- gsub('.tsv', '.ma', file)
print(newfilename)
fwrite(filtered, newfilename, sep = '\t', quote = FALSE, row.names = FALSE)






salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=30GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l


R
library(data.table)
library(dplyr)


###################


#munge sumstats for tsv file redo 
# now each file in the directory is 1 array task 


#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=20GB
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --job-name=munge
#SBATCH --account=a_palpant
#SBATCH --array=1-83
# Load the required module
module load ldsc/1.0.1--pyhdfd78af_2/module


# Set the directory where files are located
DIRECTORY="/QRISdata/Q5149/sumstats/Right_ventricle/"

# Get the list of .ma files in the directory
FILES=$(ls "$DIRECTORY"*.ma)

# Get the file for the current SLURM_ARRAY_TASK_ID
FILE=$(echo "$FILES" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Check if the file actually exists
if [ -f "$FILE" ]; then
    # Change directory to the directory of the file
    FILE_DIR=$(dirname "$FILE")
    cd "$FILE_DIR" || { echo "Failed to change directory to $FILE_DIR"; exit 1; }
    
    # Extract the filename without the path
    FILENAME=$(basename "$FILE")
    echo "Processing file: $FILENAME"
    
    # Run munge_sumstats.py for the current file
    munge_sumstats.py \
    --sumstats "$FILENAME" \
    --out "${FILENAME}" \
    --merge-alleles /scratch/user/uqsnegi/ldsc/w_hm3.snplist \
    --chunksize 500000
else
    echo "File not found for task ID: $SLURM_ARRAY_TASK_ID"
fi



    


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=20GB
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --job-name=munge
#SBATCH --account=a_palpant
#SBATCH --array=1-50
# Load the required module
module load ldsc/1.0.1--pyhdfd78af_2/module

OUTPUT_DIR="/scratch/user/uqsnegi/ldsc/sumstat_correlations/"

# Set the directory where files are located
DIRECTORY="/QRISdata/Q5149/sumstats/formatted/from_dalia/1m_imputed_sumstats/"

# Get the list of .ma files in the directory
FILES=$(ls "$DIRECTORY"*.ma)

# Get the file for the current SLURM_ARRAY_TASK_ID
FILE=$(echo "$FILES" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Check if the file actually exists
if [ -f "$FILE" ]; then
    # Change directory to the directory of the file
    FILE_DIR=$(dirname "$FILE")
    cd "$FILE_DIR" || { echo "Failed to change directory to $FILE_DIR"; exit 1; }
    
    # Extract the filename without the path
    FILENAME=$(basename "$FILE")
    echo "Processing file: $FILENAME"
    
    # Run munge_sumstats.py for the current file
    munge_sumstats.py \
    --sumstats "$FILENAME" \
    --out "${OUTPUT_DIR}${FILENAME}" \
    --merge-alleles /scratch/user/uqsnegi/ldsc/w_hm3.snplist \
    --chunksize 500000
else
    echo "File not found for task ID: $SLURM_ARRAY_TASK_ID"
fi






salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)


#####################missing files 
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=20GB
#SBATCH --time=00:15:00
#SBATCH --partition=general
#SBATCH --job-name=corr
#SBATCH --account=a_palpant
#SBATCH --array=1-15
#SBATCH -o out_%x_%j.txt  # Standard output will be written to out_corr_<job_id>.txt
#SBATCH -e error_%x_%j.txt  # Standard error will be written to error_corr_<job_id>.txt

# Load necessary modules (if needed)
module load ldsc/1.0.1--pyhdfd78af_2/module

# Define the directory containing the input .gz files
input_dir="/QRISdata/Q5149/sumstats/all_sumstats_celltype/sumstat_correlations"

# Read the missing files into an array
mapfile -t files < ${input_dir}/missing_files.txt

# Get the current array index (SLURM_ARRAY_TASK_ID) and the corresponding file
current_file=${files[$SLURM_ARRAY_TASK_ID-1]}

# Generate the --rg argument
rg_files=$(echo ${current_file} ${files[@]} | tr ' ' ',' | sed 's/,,/,/g' | sed 's/,$//')

# Extract the filename without path and extension to use as the output name
current_file_base=$(basename ${current_file} .sumstats.gz)

# Run ldsc.py command with the constructed --rg argument
ldsc.py --rg "${rg_files}" \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out "corr_${current_file_base}"
########################


#########

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=50GB
#SBATCH --time=00:10:00
#SBATCH --partition=debug
#SBATCH --job-name=corr
#SBATCH --account=a_palpant
#SBATCH --array=1-147

# Load necessary modules (if needed)
module load ldsc/1.0.1--pyhdfd78af_2/module

# Define the directory containing the input .gz files
input_dir="/scratch/user/uqsnegi/ldsc/sumstat_correlations"

# Get list of .gz files in the input directory
files=(${input_dir}/*.gz)

# Get the total number of files
num_files=${#files[@]}

# Get the current array index (SLURM_ARRAY_TASK_ID) and the corresponding first file
first_file=${files[$SLURM_ARRAY_TASK_ID]}

# Generate the --rg argument
rg_files=$(echo ${first_file} ${files[@]} | tr ' ' ',' | sed 's/,,/,/g' | sed 's/,$//')

# Extract the filename without path and extension to use as the output name
first_file_base=$(basename ${first_file} .sumstats.gz)

# Run ldsc.py command with the constructed --rg argument
ldsc.py --rg "${rg_files}" \
        --ref-ld-chr eur_w_ld_chr/ \
        --w-ld-chr eur_w_ld_chr/ \
        --out "corr_${first_file_base}"





salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)

# saving the .log files to a folder names correlation results 

R
library(data.table)

# Define the dir]ectory containing the log files
log_files <- list.files(pattern = "*.log", full.names = TRUE)


# Initialize an empty list to store the summaries
summaries_list <- list()

# Function to extract summary from a log file
extract_summary <- function(file_path) {
  # Read the log file
  lines <- readLines(file_path)
  
  # Find the line numbers for the summary table
  start_line <- grep("Summary of Genetic Correlation Results", lines) + 2
  end_line <- grep("^\\s*$", lines[start_line:length(lines)])[1] + start_line - 1
  
  # Extract the summary lines
  summary_lines <- lines[start_line:end_line]
  
  # Create a data table from the summary lines
  summary_dt <- fread(text = paste(summary_lines, collapse = "\n"), fill = TRUE, header = FALSE)
  
  # Set the column names
  setnames(summary_dt, c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"))
  
  # Remove the specified string from p1 and p2 columns
  summary_dt[, p1 := gsub("/scratch/user/uqsnegi/ldsc/sumstat_correlations/", "", p1)]
  summary_dt[, p2 := gsub("/scratch/user/uqsnegi/ldsc/sumstat_correlations/", "", p2)]
  
  # Ensure p1 is always the "smaller" value and p2 the "larger" to remove order relevance
  summary_dt[, c("p1", "p2") := .(pmin(p1, p2), pmax(p1, p2))]
  
  return(summary_dt)
}



# Loop over all log files and extract summaries
for (log_file in log_files) {
  summary_dt <- extract_summary(log_file)
  summaries_list <- append(summaries_list, list(summary_dt))
}

# Combine all summaries into a single data table
combined_summaries <- rbindlist(summaries_list)
dim(combined_summaries)

final_summary <- rbindlist(summaries_list, fill = TRUE)


# length(unique(combined_summaries$p1))
	[1] 84
	 length(unique(combined_summaries$p2))
	[1] 83
#[1] 7055   12
# Remove duplicate rows based on p1 and p2
library(dplyr)
library(stringr)

# Assuming combined_summaries is your data frame
combined_summaries <- combined_summaries %>%
  mutate(pair = if_else(p1 < p2, 
                        paste(p1, p2, sep = "_"), 
                        paste(p2, p1, sep = "_")))

# View the updated data frame
head(combined_summaries)

combined_summaries_unique <- combined_summaries %>%
  distinct(pair, .keep_all = TRUE)



x <- combined_summaries_unique[,c(1,2,3,4,5,6)]

x <- x[!grepl("cleaned", x$p1) & !grepl("cleaned", x$p2), ]

x <- x %>%
  mutate(q = p.adjust(p, method = "fdr"))

# Set the smallest adjusted p-value to 0
min_q <- min(x$q)
x <- x %>%
  mutate(q = q - min_q)


> dim(combined_summaries_unique)
[1] 10893    13


 

original order


data <- data.frame(Filename = filenames, Trait_Name = trait_names)


merged_data <- merge(x, data, by.x = "p1", by.y = "Filename", all.x = TRUE)
names(merged_data)[names(merged_data) == 'Trait_Name'] <- 'p1_trait'


merged_data2 <- merge(merged_data, data, by.x = "p2", by.y = "Filename", all.x = TRUE)
names(merged_data2)[names(merged_data2) == 'Trait_Name'] <- 'p2_trait'

#p1_category p1  p2_category p2  rg  se  z   p   q

p1_category p1  p2_category p2  rg  se  z   p   q
Anthropometric  Height  Anthropometric  BMI -0.0587 0.024   -2.4421 0.014602    0.136797684
Anthropometric  Height  Metabolic   TC  -0.0778 0.0344  -2.2634 0.023611    0.18869

x <- merged_data2
saveRDS(x,"rg.RDS")





CATEGORY    TRAIT   COLOR
Anthropometric  Height  #4C72AE
Anthropometric  BMI #4C72AE
Metabolic   TC  #E4812F
Metabolic   HDL-C   #E4812




# Create a named vector mapping each trait to its corresponding category
# Create the mapping
# Load necessary library
library(dplyr)

# Given category mapping
# Create a named vector mapping variables to categories
category_map <- c(
  LVEDV_indexed = "LV_MRI",
  LVESV_indexed = "LV_MRI",
  MFD1 = "Fractal",
  MFD2 = "Fractal",
  MFD3 = "Fractal",
  MFD4 = "Fractal",
  MFD5 = "Fractal",
  MFD6 = "Fractal",
  MFD7 = "Fractal",
  MFD8 = "Fractal",
  MFD9 = "Fractal",
  LVESIFWT = "Wall thickness",
  LVESIWT = "Wall thickness",
  LVESILWT = "Wall thickness",
  LVESALWT = "Wall thickness",
  LVESAWT = "Wall thickness",
  LVESASWT = "Wall thickness",
  LVEDISWT = "Wall thickness",
  DCM = "CVD",
  Heart_failure = "CVD",
  LV_SV_indexed = "LV_MRI",
  LV_SV = "LV_MRI",
  LVEF = "LV_MRI",
  max_ascending_aorta_diam_indexed = "OF_MRI",
  max_ascending_aorta_diam = "OF_MRI",
  max_lv_indexed = "LV_MRI",
  max_lv = "LV_MRI",
  min_lv_indexed = "LV_MRI",
  min_lv = "LV_MRI",
  pa_ao_diastole = "OF_MRI",
  pa_ao = "OF_MRI",
  rafac = "RA_MRI",
  ramax_area_indexed = "RA_MRI",
  ramax_area = "RA_MRI",
  ramin_area_indexed = "RA_MRI",
  ramin_area = "RA_MRI",
  RVEDV_indexed = "RV_MRI",
  RVEDV_LVEDV_ratio = "V_MRI",
  RVEDV = "RV_MRI",
  RVEF_LVEF_ratio = "V_MRI",
  RVEF = "RV_MRI",
  RVESV_indexed = "RV_MRI",
  RVESV_LVESV_ratio = "V_MRI",
  RVESV = "RV_MRI",
  RVSV_indexed = "RV_MRI",
  RVSV_LVSV_ratio = "RV_MRI",
  RVSV = "RV_MRI",
  short_axis_cm_pa_diastole_indexed = "OF_MRI",
  short_axis_cm_pa_diastole = "OF_MRI",
  short_axis_cm_pa_indexed = "OF_MRI",
  short_axis_cm_pa_strain = "OF_MRI",
  short_axis_cm_pa = "OF_MRI",
  short_axis_cm_root_indexed = "OF_MRI",
  long_covid_vs_covid = "Other",
  LV_SV_indexed2 = "LV_MRI",
  Cardiac_arrhythmia_COPD_comorbidities = "CVD",
  Coronary_atherosclerosis = "CVD",
  Diagnoses_main_ICD10_I20_Angina_pectoris = "CVD",
  Diagnoses_main_ICD10_I21_Acute_myocardial_infarction = "CVD",
  Diagnoses_main_ICD10_I25_Chronic_ischaemic_heart_disease = "CVD",
  Diagnoses_main_ICD10_I48_Atrial_fibrillation_and_flutter = "CVD",
  Diastolic_blood_pressure_automated_reading = "BP_HR",
  Diastolic_blood_pressure_manual_reading = "BP_HR",
  ECG_heart_rate = "BP_HR",
  ECG_load = "BP_HR",
  Hip_circumference = "PM",
  Ischaemic_heart_disease_wide_definition = "CVD",
  Major_coronary_heart_disease_event = "CVD",
  Maximum_heart_rate_during_fitness_test = "CVD",
  Non_cancer_illness_code_self_reported_high_cholesterol = "Other",
  Myocardial_infarction = "CVD",
  Systolic_blood_pressure_automated_reading = "BP_HR",
  Systolic_blood_pressure_manual_reading = "BP_HR",
  Waist_circumference = "PM",
  Weight = "PM",
  Whole_body_fat_mass = "PM"
)

# Display the category mapping vector
print(category_map)

# Convert the mapping to a data frame
df <- data.frame(
  TRAIT = names(category_map),
  CATEGORY = unname(category_map)
)

# Assign arbitrary colors to each unique category
unique_categories <- unique(df$CATEGORY)
colors <- c("#4C72AE", "#E4812F", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33") # Add more colors if needed


#disease 
#unique_categories <- c("LV_MRI", "Fractal", "Wall thickness", "CVD", "OF_MRI", "RA_MRI", "RV_MRI", "V_MRI", "Other", "BP_HR", "PM")
#colors <- c("#4C72AE", "#E4812F", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33","#A65628", "#F781BF", "#FFEDA0")




# Ensure there are enough colors for the categories
if (length(unique_categories) > length(colors)) {
  stop("Not enough colors for the number of unique categories.")
}

# Create a data frame with categories and their assigned colors
category_colors <- data.frame(
  CATEGORY = unique_categories,
  COLOR = colors[1:length(unique_categories)]
)

# Merge the data frames to assign colors to each trait
df <- merge(df, category_colors, by = "CATEGORY")

# Print the resulting data frame



# Create two new columns p1_category and p2_category in x
x <- x %>%
  left_join(df, by = c("p1_trait" = "TRAIT")) %>%
  rename(p1_category = CATEGORY) %>%
  left_join(df, by = c("p2_trait" = "TRAIT")) %>%
  rename(p2_category = CATEGORY)

#na removal
x <- x[complete.cases(x), ]


 x <- x[,c(10,8,11,9,3,4,5,6,7)]

 names(x) <- c("p1_category", "p1", "p2_category", "p2", "rg","se","z","p","q")

setorder(x, p1_category, p2_category)
  dim(x)
[1] 1485    9

saveRDS(x, "formattedRG.RDS")


fwrite(x,"formatted_RG.txt", sep="\t")
##### sort by category!!!!! 
#fwrite - formatted_RG.txt
# get chatgpt to help with category?
rg.txt 

p1_category	p1	p2_category	p2	rg	se	z	p	q
Anthropometric	Height	Anthropometric	BMI	-0.0587	0.024	-2.4421	0.014602	0.136797684
Anthropometric	Height	Metabolic	TC	-0.0778	0.0344	-2.2634	0.023611	0.188695257
Anthropometric	Height	Metabolic	HDL-C	-0.0045	0.0364	-0.123	0.90208	0.971590187









category_colors <- data.frame(
  CATEGORY = c("Disease", "LV volume", "Fractal slice", "LV wall thickness", "RV volume", "Ventricular", "Contraction", "OF", "RA"),
  COLOR = c("#4C72AE", "#E4812F", "#8C61B6", "#FFA500", "#008000", "#FFFF00", "#FF4500", "#1E90FF", "#FF1493")
)

# Merge category_colors with df to add color codes
df_with_colors <- df %>%
  left_join(category_colors, by = "CATEGORY")

# Print the resulting data frame
print(df_with_colors)
saveRDS(df_with_colors,"traitlist.RDS")
fwrite(df_with_colors,"traitlist.txt", sep="\t")



traitlist.txt 

	CATEGORY	TRAIT	COLOR
	Anthropometric	Height	#4C72AE
	Anthropometric	BMI	#4C72AE
	Metabolic	TC	#E4812F
	Metabolic	HDL-C	#E4812F
	Metabolic	LDL-C	#E4812F
	Metabolic	TG	#E4812F



######## script  corr.R

library(corrplot)
library(reshape2)

# quantitative traits

TRAIT_CATEGORY1 = c("LV volume", "Fractal slice", "LV wall thickness", "RV volume", "Ventricular", "Contraction", "OF", "RA"),


# diseases
TRAIT_CATEGORY2 = c("Disease")

################################################################################
args = commandArgs(trailingOnly=T)
if (identical(args, character(0))) {
  args = paste0("./input_example/", c("input_rg.txt", "traitlist.txt"))
}
rg_fname = args[1]
traitlist_fname = args[2]


################################################################################
# load data
rg = read.table(rg_fname, T, sep='\t', as.is = T)
trait_all = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', fileEncoding='utf-8', comment.char="")

traitlist1 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY1)
traitlist2 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY2)

# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)

################################################################################

corrplot_nsquare = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  x2 = dcast(rg, p1 ~ p2, value.var = "rg")
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1

  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  x2 = dcast(rg, p1 ~ p2, value.var = "q")
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1 = match(trait1, rownames(mat2))
  trait2 = match(trait2, colnames(mat2))
  mat2 = mat2[trait1, trait2]
  qmat2 = qmat2[trait1, trait2]

  corrplot(mat2, method = "psquare", order = order,
           p.mat = qmat2, sig.level = 0.05, sig = "pch",
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")
  # return(list(mat = mat2, qmat = qmat2))
}
x11()

pdf("./output/corrplot_rg_disease.pdf", width = 16, height = 9, family = "Helvetica")
  corrplot_nsquare(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()

png("./output/corrplot_rg_disease.png", width = 16, height = 9, units = 'in', res = 300, family = "Helvetica")
  corrplot_nsquare(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()


warnings()



## corr2.R
library(corrplot)
library(reshape2)

# quantitative traits

TRAIT_CATEGORY1 = c("BP", "Heart rate","contractility", "MFD","CONtractility",  "LVWT", "LV", "RV", "RV_LV_ratio", "OFT","Comorbidity")



# diseases
TRAIT_CATEGORY2 = c("BP", "Heart rate", "contractility", "CONtractility", "MFD",  "LVWT", "LV", "RV", "RV_LV_ratio", "OFT","Comorbidity")



################################################################################
args = commandArgs(trailingOnly=T)
if (identical(args, character(0))) {
  args = paste0("./input_example/", c("input_rg.txt", "traitlist.txt"))
}
rg_fname = args[1]
traitlist_fname = args[2]
g

################################################################################
# load data
rg = read.table(rg_fname, T, sep='\t', as.is = T)
trait_all = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', fileEncoding='utf-8', comment.char="")

traitlist1 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY1)
traitlist2 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY2)

# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)
corrplot_nsquare_with_barplot = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  
  # Prepare correlation matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "rg",fun.aggregate=mean)
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  # Prepare FDR matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate=mean)
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1_idx = match(trait1, rownames(mat2))
  trait2_idx = match(trait2, colnames(mat2))
  mat2 = mat2[trait1_idx, trait2_idx]
  qmat2 = qmat2[trait1_idx, trait2_idx]

  # Calculate the number of significant correlations per trait
  sig_counts = apply(qmat2, 1, function(x) sum(x < 0.05))

  # Set up the layout to include the barplot
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

  # Plot correlation matrix with significance
  corrplot(mat2, method = "square", order = order,
           p.mat = qmat2, sig.level = 0.05, insig = "label_sig", 
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")

  # Add legend for square size
  par(new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  legend("bottomleft", legend = c("FDR < 0.05", "FDR < 0.1", "FDR >= 0.1"), 
         pch = 22, pt.cex = c(2, 1.5, 1), col = "black", 
         title = "Square Size", cex = 1.2, bty = "n")

  # Plot the barplot
  par(mar = c(5, 4, 0, 2))  # Adjust margins to fit the barplot
  barplot(sig_counts, main = "Significant Correlations per Trait",
          xlab = "Trait", ylab = "Number of Significant Correlations",
          las = 2, col = "skyblue")
}

# Example usage


png("corrplot2_rg_disease_with_barplot.png", width = 40, height = 20, units = 'in', res = 300, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()



Rscript corr2.R formatted_RG.txt traitlist.txt 




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)




# corr3.R
library(corrplot)
library(reshape2)

TRAIT_CATEGORY1 = c("Fractal","LV_MRI","OF_MRI","RV_MRI","V_MRI","Wall thickness")


# diseases
TRAIT_CATEGORY2 = c("BP_HR","CVD","Other","PM")

################################################################################
args = commandArgs(trailingOnly=T)
if (identical(args, character(0))) {
  args = paste0("./input_example/", c("input_rg.txt", "traitlist.txt"))
}
rg_fname = args[1]
traitlist_fname = args[2]


################################################################################
# load data
rg = read.table(rg_fname, T, sep='\t', as.is = T)
trait_all = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', fileEncoding='utf-8', comment.char="")

traitlist1 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY1)
traitlist2 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY2)

# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)
corrplot_nsquare_with_barplot = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  
  # Prepare correlation matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "rg",fun.aggregate=mean)
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  # Prepare FDR matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate=mean)
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1_idx = match(trait1, rownames(mat2))
  trait2_idx = match(trait2, colnames(mat2))
  mat2 = mat2[trait1_idx, trait2_idx]
  qmat2 = qmat2[trait1_idx, trait2_idx]

  # Calculate the number of significant correlations per trait
  sig_counts = apply(qmat2, 1, function(x) sum(x < 0.05))

  # Set up the layout to include the barplot
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

  # Plot correlation matrix with significance
  corrplot(mat2, method = "square", order = order,
           p.mat = qmat2, sig.level = 0.05, insig = "label_sig", 
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")

  # Add legend for square size
  par(new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  legend("bottomleft", legend = c("FDR < 0.05", "FDR < 0.1", "FDR >= 0.1"), 
         pch = 22, pt.cex = c(2, 1.5, 1), col = "black", 
         title = "Square Size", cex = 1.2, bty = "n")

  # Plot the barplot
  par(mar = c(5, 4, 0, 2))  # Adjust margins to fit the barplot
  barplot(sig_counts, main = "Significant Correlations per Trait",
          xlab = "Trait", ylab = "Number of Significant Correlations",
          las = 2, col = "skyblue")
}

# Example usage
x11()
pdf("corrplot2_rg_disease_with_barplot.pdf", width = 16, height = 12, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()

png("corrplot2_rg_disease_with_barplot.png", width = 16, height = 12, units = 'in', res = 300, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=200GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)


hclust_corrplot_nsquare_with_barplot = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  
  # Prepare correlation matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "rg",fun.aggregate=mean)
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  # Prepare FDR matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate=mean)
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1_idx = match(trait1, rownames(mat2))
  trait2_idx = match(trait2, colnames(mat2))
  mat2 = mat2[trait1_idx, trait2_idx]
  qmat2 = qmat2[trait1_idx, trait2_idx]

  # Calculate the number of significant correlations per trait
  sig_counts = apply(qmat2, 1, function(x) sum(x < 0.05))

  # Set up the layout to include the barplot
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

  # Plot correlation matrix with significance
  corrplot(mat2, method = "square", order = 'hclust',
           p.mat = qmat2, sig.level = 0.05, insig = "label_sig", 
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")

  # Add legend for square size
  par(new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  legend("bottomleft", legend = c("FDR < 0.05", "FDR < 0.1", "FDR >= 0.1"), 
         pch = 22, pt.cex = c(2, 1.5, 1), col = "black", 
         title = "Square Size", cex = 1.2, bty = "n")

  # Plot the barplot
  par(mar = c(5, 4, 0, 2))  # Adjust margins to fit the barplot
  barplot(sig_counts, main = "Significant Correlations per Trait",
          xlab = "Trait", ylab = "Number of Significant Correlations",
          las = 2, col = "skyblue")
}

(order.hc2 = corrMatOrder(mat2, order = 'hclust', hclust.method = 'ward.D2'))


# same as: corrplot(M, order = 'hclust', hclust.method = 'ward.D2', addrect = 2)
corrplot(M.hc2)
corrRect.hclust(M.hc2, k = 2, method = 'ward.D2')


########### amazing 

corrplot(M.hc2, method = "square", order = 'hclust', hclust.method = 'ward.D2', addrect = 11, type = 'lower',                                                      
         p.mat = q.hc2, sig.level = 0.05, insig = "label_sig", 
         pch = "*", pch.cex = 1.5, full_col = FALSE,
         na.label = "square", na.label.col = "grey30",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45)  # Adjust tl.srt for text angle
##########
sig_counts_ordered <- sig_counts[order.hc2]

# Open a PDF device for saving plots
pdf("barplot_significant_correlations.pdf")

# Plot the barplot with ordered sig_counts
barplot(sig_counts_ordered, main = "Significant Correlations per Trait",
        xlab = "Trait", ylab = "Number of Significant Correlations",
        las = 2, col = "skyblue")

# Close the PDF device
dev.off()

########################









salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)

bash script to format a file n

#!/bin/bash
#!/bin/bash

# Function to clean filenames
clean_filename() {
    local filename="$1"
    # Remove spaces, colons, and commas from the beginning and end
    filename=$(echo "$filename" | sed 's/^[ ,:]*//; s/[ ,:]*$//')
    # Replace spaces, colons, and commas in the middle with underscores
    filename=$(echo "$filename" | sed 's/[ ,:]/_/g')
    # Replace multiple underscores with a single underscore
    filename=$(echo "$filename" | sed 's/_\+/_/g')
    echo "$filename"
}

# Directory to process
directory="/QRISdata/Q5149/sumstats/formatted/from_dalia/1m_imputed_sumstats"

# Iterate through each file in the directory
for filepath in "$directory"/*; do
    filename=$(basename "$filepath")
    new_filename=$(clean_filename "$filename")
    new_filepath="$directory/$new_filename"
    
    # Rename the file if the new filename is different
    if [[ "$new_filename" != "$filename" ]]; then
        mv "$filepath" "$new_filepath"
        echo "Renamed: \"$filename\" to \"$new_filename\""
    else
        echo "No change needed for: \"$filename\""
    fi
done








######################### 


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=200GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
library(data.table)
library(dplyr)





############ total number of analyzed 
3175 pairwise correlations for binary traits and disease traits 






## corr2.R
library(corrplot)
library(reshape2)

# quantitative traits

TRAIT_CATEGORY1 = c("BP", "Heart rate","contractility", "MFD","CONtractility",  "LVWT", "LV", "RV", "RV_LV_ratio", "OFT","Comorbidity")



# diseases
TRAIT_CATEGORY2 = c("disease")



################################################################################
################################################################################
# load data
rg = read.table(rg_fname, T, sep='\t', as.is = T)
trait_all = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', fileEncoding='utf-8', comment.char="")

traitlist1 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY1)
traitlist2 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY2)


hclust_corrplot_nsquare_with_barplot = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  


# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)

#########################
# Prepare correlation matrix
  x2 <- reshape2::dcast(rg, p1 ~ p2, value.var = "rg", fun.aggregate = mean)
  mat2 <- as.matrix(x2[, -1])
  rownames(mat2) <- x2$p1
  mat2[mat2 > 1] <- 1
  mat2[mat2 < -1] <- -1
  mat2[is.na(mat2)] <- 0
  
  # Prepare FDR matrix
  x2 <- reshape2::dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate = mean)
  qmat2 <- as.matrix(x2[, -1])
  rownames(qmat2) <- x2$p1
  qmat2[is.na(qmat2)] <- 1
  
  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) <- 1
    diag(qmat2) <- -1
  }
  
  # Match trait indices based on formatted_trait

  # Calculate the number of significant correlations per trait
  sig_counts <- apply(qmat2, 1, function(x) sum(x < 0.05))
  
(order.hc2 = corrMatOrder(mat2, order = 'hclust', hclust.method = 'ward.D2'))


M.hc2 = mat2[order.hc2, order.hc2]
q.hc2 = qmat2[order.hc2, order.hc2]


########### amazing 
pdf("final_quantiative_corrplot.pdf")

corrplot(M.hc2, method = "square", order = 'hclust', hclust.method = 'ward.D2', addrect = 11, type = 'lower',                                                      
         p.mat = q.hc2, sig.level = 0.05, insig = "label_sig", 
         pch = "*", pch.cex = 1.5, full_col = FALSE,
         na.label = "square", na.label.col = "grey30",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45)  # Adjust tl.srt for text angle

dev.off()
##########

sig_counts_ordered <- sig_counts[order.hc2]

# Open a PDF device for saving plots
pdf("barplot_significant_correlations.pdf")

# Plot the barplot with ordered sig_counts
barplot(sig_counts_ordered, main = "Significant Correlations per Trait",
        xlab = "Trait", ylab = "Number of Significant Correlations",
        las = 2, col = "violet")

# Close the PDF device
dev.off()






# Example usage:
# Assuming traitlist1$TRAIT and traitlist2$TRAIT are vectors of original trait names









library(corrplot)
library(reshape2)

TRAIT_CATEGORY1 = c("BP", "Heart rate","contractility", "MFD","CONtractility",  "LVWT", "LV", "RV", "RV_LV_ratio", "OFT","Comorbidity")


# diseases
TRAIT_CATEGORY2 = c("Disease")

################################################################################

################################################################################
# load data

traitlist1 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY1)
traitlist2 = subset(trait_all, CATEGORY %in% TRAIT_CATEGORY2)

# duplicate lines
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
tmp$p1_category = rg$p2_category
tmp$p2_category = rg$p1_category
rg = rbind(rg, tmp)
corrplot_nsquare_with_barplot = function(rg, trait1, trait2, traits_use = NULL, order = "original", landscape=FALSE) {
  if (landscape & length(trait1) > length(trait2)) {
    tmp = trait1
    trait1 = trait2
    trait2 = tmp
  }
  if (!is.null(traits_use)) {
    rg = subset(rg, p1 %in% traits_use & p2 %in% traits_use)
  }
  
  # Prepare correlation matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "rg",fun.aggregate=mean)
  mat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(mat2) = x2$p1
  mat2[mat2 > 1] = 1
  mat2[mat2 < -1] = -1
  mat2[is.na(mat2)] = 0

  # Prepare FDR matrix
  x2 = dcast(rg, p1 ~ p2, value.var = "q", fun.aggregate=mean)
  qmat2 = as.matrix(x2[, 2:ncol(x2)])
  rownames(qmat2) = x2$p1
  qmat2[is.na(qmat2)] = 1

  if (nrow(mat2) == ncol(mat2)) {
    diag(mat2) = 1
    diag(qmat2) = -1
  }

  trait1_idx = match(trait1, rownames(mat2))
  trait2_idx = match(trait2, colnames(mat2))
  mat2 = mat2[trait1_idx, trait2_idx]
  qmat2 = qmat2[trait1_idx, trait2_idx]

  # Calculate the number of significant correlations per trait
  sig_counts = apply(qmat2, 1, function(x) sum(x < 0.05))

  # Set up the layout to include the barplot
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))

  # Plot correlation matrix with significance
  corrplot(mat2, method = "square", order = order,
           p.mat = qmat2, sig.level = 0.05, insig = "label_sig", 
           pch = "*", pch.cex = 1.5, full_col=FALSE,
           na.label = "square", na.label.col = "grey30")

  # Add legend for square size
  par(new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  legend("bottomleft", legend = c("FDR < 0.05", "FDR < 0.1", "FDR >= 0.1"), 
         pch = 22, pt.cex = c(2, 1.5, 1), col = "black", 
         title = "Square Size", cex = 1.2, bty = "n")

  # Plot the barplot
  par(mar = c(5, 4, 0, 2))  # Adjust margins to fit the barplot
  barplot(sig_counts, main = "Significant Correlations per Trait",
          xlab = "Trait", ylab = "Number of Significant Correlations",
          las = 2, col = "skyblue")

}
# Example usage
x11()
pdf("corrplot2_rg_disease_with_barplot.pdf", width = 16, height = 12, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()

png("corrplot2_rg_disease_with_barplot.png", width = 16, height = 12, units = 'in', res = 300, family = "Helvetica")
  corrplot_nsquare_with_barplot(rg, traitlist1$TRAIT, traitlist2$TRAIT, landscape=TRUE)
dev.off()


}




salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=20GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

R
library(data.table)
library(dplyr)





