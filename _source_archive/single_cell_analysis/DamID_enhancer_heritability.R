


hopx_enriched_enhancer <- read.csv("hopx_wt_ratio1_average_peaks2_enhancers.bed", header= F ,sep="\t")

dim (hopx_enriched_enhance
[1] 174   4

 cardiac_enhancers <- read.csv("cardiac_enhancer_hg38.bed", header=F , sep="\t")
 dim (cardiac_enhancers)
[1] 45887     5


snps <- read.csv("ukbEURu_v3_hm3_maf01.bim", header=F, sep="\t")
dim (snps)
[1] 1130559   6    

snps2 <- snps[, c(1,2,4)]

dim(snps2)
[1] 1130559       3


library(tidyverse)

#making bed file from chr to number only 
cardiac_enhancers$V1 <- gsub('[chr]', '', cardiac_enhancers$V1)


####if the snps overlap with a cardiac enhancer location put 1  or else put zero 

snps2$overlapping_with_cardiac_enhancers <- ifelse(sapply(seq_along(snps2$V4), function(i) {
           inds <- cardiac_enhancers$V2 <= snps2$V4[i] & cardiac_enhancers$V3 >= snps2$V4[i]
           any(inds) & (snps2$V1[i] == cardiac_enhancers$V1[which.max(inds)])
           }), "1", "0")

snps2$overlapping_with_cardiac_enhancers

table(snps2$overlapping_with_cardiac_enhancers)

#      0       1 
#1092341   29862 
#2.6% of snps lie in cardiac enhancers 


##############################
#hopx_enriched_enhancer
###############################


making bed file from chr to number only 
hopx_enriched_enhancer$V1 <- gsub('[chr]', '', hopx_enriched_enhancer$V1)


####if the snps overlap with a cardiac enhancer location put 1  or else put zero 

snps2$overlapping_with_hopxenriched <- ifelse(sapply(seq_along(snps2$V4), function(i) {
           inds <- hopx_enriched_enhancer$V2 <= snps2$V4[i] & hopx_enriched_enhancer$V3 >= snps2$V4[i]
           any(inds) & (snps2$V1[i] == hopx_enriched_enhancer$V1[which.max(inds)])
           }), "1", "0")

snps2$overlapping_with_hopxenriched

table(snps2$overlapping_with_hopxenriched)

#      0       1 
#1130426     133 

####################
#WT
##################

wt$V1 <- gsub('[chr]', '', wt$V1)


snps2$overlapping_with_wtenriched <- ifelse(sapply(seq_along(snps2$V4), function(i) {
           inds <- wt$V2 <= snps2$V4[i] & wt$V3 >= snps2$V4[i]
           any(inds) & (snps2$V1[i] == wt$V1[which.max(inds)])
           }), "1", "0")



#      0       1 
#1130526      33 


######################################################
#combined 
#########################################

c$V1 <- gsub('[chr]', '', c$V1)


snps2$overlapping_with_differentialaccessability <- ifelse(sapply(seq_along(snps2$V4), function(i) {
           inds <- c$V2 <= snps2$V4[i] & c$V3 >= snps2$V4[i]
           any(inds) & (snps2$V1[i] == c$V1[which.max(inds)])
           }), "1", "0")

table(snps2$overlapping_with_differentialaccessability)
#      0       1 
#1130393     166 



########################################################################################################
###########################################################################################################
#making annotation files
cardiac_enhancer_snps <- snps2[,c(2,4)]

cardiac_enhancer_snps$overlapping_with_cardiac_enhancers <- as.numeric(cardiac_enhancer_snps$overlapping_with_cardiac_enhancers)
 cardiac_enhancer_snps$opposite <- (1 - cardiac_enhancer_snps$overlapping_with_cardiac_enhancers)


na.omit(cardiac_enhancer_snps)
write.table(cardiac_enhancer_snps, "hm3_cardiac_enhancer_snps.txt", sep="\t", quote = FALSE , col.names = TRUE, row.names = FALSE)


awk '{print $3, $2}'  hm3_cardiac_enhancer_snps.txt | sort | uniq




meta <- read.csv("all_enhancer_annotations_to_hm3_snps.txt", header=T, sep="\t")


differentially_accessible_enhancer_snps <- meta[,c(2,7)]
differentially_accessible_enhancer_snps$overlapping_with_differentialaccessability <- as.numeric(differentially_accessible_enhancer_snps$overlapping_with_differentialaccessability)
 differentially_accessible_enhancer_snps$opposite <- (1 - differentially_accessible_enhancer_snps$overlapping_with_differentialaccessability )


 write.table(differentially_accessible_enhancer_snps,"differentially_accyessible_enhancer_snps_genomewide.txt", sep="\t", quote = FALSE , col.names = TRUE, row.names = FALSE)

####################################################################################

#annotation file for differentially accesible enhancer snps among snps in cardiac enhancer regions ONLY
# DA_snps_in_CE_regions
#####################################################################################

meta2 <- meta[,c(2,4,7)]

#code below didn't work only made a new column
##protein_coding_snps_hopx <- Overlapping_protein_coding_snps %>%  mutate_all(opposite = if_else (Overlapping_protein_coding_snps$V1 %in% hopx$V4, "0","1" )










##########################################
##########SBAYESS#######################
#######################################
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=50gb
#SBATCH --time=15:00:00
#SBATCH -o 
#SBATCH -e 
#SBATCH --partition=general
#SBATCH --job-name=enhancer
#SBATCH --account=a_palpant





MY_DIR=/QRISdata/Q4879/Cardiac_enhancer_heritability_SN/

cd $MY_DIR


cd $SLURM_SUBMIT_DIR
gctb="/gctb"

$gctb --sbayes S --mldm LDM/ukbEURu_hm3_sparse_mldm_list.txt --gwas-summary sumstats/esv_bsa.ma --annot hm3_cardiac_enhancer_snps.txt --chain-length 12000 --burn-in 2000 --out cardiac_enhancer > cardiac_enhancer.log 2>&1




QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr10_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr11_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr12_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr13_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr14_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr15_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr16_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr17_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr18_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr19_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr1_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr20_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr21_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr22_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr2_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr3_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr4_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr5_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr6_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr7_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr8_v3_50k.ldm.sparse
QRISdata/Q4789/cardiac_enhancer_heritability_SN/LDM/ukbEURu_hm3_chr9_v3_50k.ldm.sparse



 cardiac_enhancer_snps$opposite <- (1 - cardiac_enhancer_snps$overlapping_with_cardiac_enhancers)











##########################################################
#SNPS within HOPX bound peaks table S2
####################################################
meta$overlapping_with_hopx_bound_peaks <- ifelse(sapply(seq_along(meta$V4), function(i) {
           inds <-h$V2 <= meta$V4[i] & h$V3 >= meta$V4[i]
           any(inds) & (meta$V1[i] == h$V1[which.max(inds)])
           }), "1", "0")

table(meta$overlapping_with_hopx_bound_peaks)

> h2$overlapping_with_hopx_bound_peaks<- as.numeric(h2$overlapping_with_hopx_bound_peaks)
> h2$opposite <- (1- h2$overlapping_with_hopx_bound_peaks)
> head(h2)


write.table(h2, "snps_within_hopx_bound_peaks.txt", header=T, sep="\t", quote=F ,col.names = TRUE, row.names = FALSE)


#####################################################################
#SNPS within hOPX bound peaks that are also cardiac enhancers 
#####################################################################

>hopx_bound_peaks_and_ce <- meta[(meta$overlapping_with_cardiac_enhancers == '1' & meta$overlapping_with_hopx_bound_peaks == '1'),]
dim(hopx_bound_peaks_and_ce)
[1] 372   8


write.table(h2, "snps_within_hopx_bound_peaks_overlapping_ce.txt", header=T, sep="\t", quote=F ,col.names = TRUE, row.names = FALSE)
