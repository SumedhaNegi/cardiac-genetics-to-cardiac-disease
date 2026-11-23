# =============================================================================
# EXAMINER GUIDE
# Script: correlation_hopx_paper.R
# Purpose: Historical workflow for expression/correlation analysis associated with the HOPX paper; includes session and directory commands.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# Performs the HOPX-focused gene-expression correlation analysis used for paper/thesis comparisons.
#

# HISTORICAL SOURCE CAVEAT: This file contains HPC module command, PBS scheduler command.
# These commands/output are retained for provenance; see docs/SCRIPT_AUDIT.md.


#BSW data set
##########################

ssh -XYC delta2.imb.uq.edu.au

cd 90days/analysis/scRNA_correlation/

qsub -XI -l mem=120gb,walltime=4:00:00


##
# HPC / shell command from the interactive source session (not R code): module load R/4.1.2
export R_LIBS=/scratch/90days/uqsnegi/RLIBS

# The original session entered the interactive R console at this point.
library(Seurat)
library(biomaRt)
library(plyr)

dir = '/scratch/90days/uqsnegi/analysis/scRNA_correlation/'


#read functions for GeneCorrelationAnalysis.R
mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

genes <- c('MTF1','HOPX','BCL6B','GATA3','FOSL1','KLF4','TBX3','TEAD4','PRDM6','FOXJ3','SP4','PATZ1', "TTN", "ATXN2", "ELN", "BAG3" )


bsw <- readRDS('/afm01/UQ/Q4879/WSB_Seurat_integrated/common_orig_log_FPKM2L_integrated.RDS')


library(ggrepel)
library(ggrastr)

# 210125 plot genes ranked by correlation 
# mark position of HOPX, draw hashed line at 10% most up and most down correlated genes then save as pdf
library(ggrepel)
library(ggrastr)
genes <- tools::toTitleCase(tolower(c('MTF1','HOPX','BCL6B','GATA3','FOSL1','KLF4','TBX3','TEAD4','PRDM6','FOXJ3','SP4','PATZ1', "TTN", "ATXN2", "ELN", "BAG3")))

for (gene in genes){pltord(gene, df = bsw_df)}

# 210127 blanpain data, save the ordered expression/correlation plots for correlation with MYL2 for MYL2, MYH7 and ACTB
genes <- c('Myl2', 'Myh7', 'Actb')
myl2_corr <- runCorr(bsw_df, 'Myl2')

pdf('outs/210127_NEW_Myl2CorrsBSW.pdf', width = 21)
print(pltgenes(genes, myl2_corr, rast = TRUE))
dev.off(

for(gene in genes){pltord(gene = 'HOPX', df = bsw_df, plotgene = gene, , saveprefix = paste0(dir, 'outs/bsw_runx3'))}


######## this code works :) 


genex <- c('^LRP5_', '^TCF7_', '^GSK3B_')
for(genx in genex){pltord(gene = 'HOPX', plotgene = genx, outt = dim_ot, save = TRUE, saveprefix = paste0(dir, 'outs/210304_dimmeler_'))}


# Genes not present in data 																	#
#ZNF384	
#"Onecut3
#Znf524
#"Zn770"
#Zsc22
#Andr""Zn770" 
#"Zn264
#"Znf384" 

#################

Doing for claytons d30 dataset 

cf30 <- readRDS('/afm01/UQ/Q4879/Friedman2018_scRNAseq_MaikaProcessed/day30_separate_final.RDS')
cf_df <- data.frame(t(cf30$Exprs))[which(cf30$Clusters == '2'),]

dir="/home/uqsnegi/90days/analysis/scRNA_correlation/"


genes <- (c('MTF1','HOPX','GATA3','FOSL1','KLF4','TBX3','TEAD4','PRDM6','FOXJ3','SP4','PATZ1', "TTN", "ATXN2", "ELN", "BAG3"))

for(gene in genes){pltord(gene = 'HOPX', df = cf_df, plotgene = gene, , saveprefix = paste0(dir, 'outs/CF30_'))}


#########################################################################

# Doing gene ontology of top 10% most correlated and bottom 10$ most correlated 

df <- read.csv("210120_NEW_BSW_Hopx_ENSMUSG00000059325_corr.txt", header=T, sep =
+ "\t")


> dim(df)
[1] 13656     2

TOP 100 GENES
df2 <- df[1:100,]

bottom 100 genes
df3 <- df[13557:13656,]

#only keep ensembl ids (remove everything before underscore)
df3$gene <- gsub(".*_","",df3$gene)

write.table(df3, "", sep="\t", quote = FALSE , col.names = TRUE, row.names = FALSE)


