# =============================================================================
# EXAMINER GUIDE
# Script: sumedha_qpcr_pipelie.R
# Purpose: Processes qPCR measurements; inspect the code for the specific normalisation and comparison steps.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# Imports qPCR setup, melt-curve and results exports; derives sample groups and HPRT-normalised expression, then explores replicate distributions and group comparisons.
#

# HISTORICAL SOURCE CAVEAT: This file contains interactive HPC allocation, HPC module command.
# These commands/output are retained for provenance; see docs/SCRIPT_AUDIT.md.

Sumedha version 


# HPC / shell command from the interactive source session (not R code): salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=100GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

# HPC / shell command from the interactive source session (not R code): module load r/4.2.1-foss-2021a

# The original session entered the interactive R console at this point.
library(data.table)
library(dplyr)
library(tidyr)

setwd("/QRISdata/Q5149/qpcr/exported_files_to_format")


'%!in%' <- function(x,y)!('%in%'(x,y))
fcn <- function(df){setnames(df, gsub(' ', '_', colnames(df)))} # remove spaces from colnames


files <- list.files(pattern = '.*Setup.*')
meta <- list()
for (i in 1:min(6, length(files))) {
  data <- fcn(fread(files[i], skip = '[Sample Setup]'))
 meta[[i]] <- data
}
meta <- rbindlist(meta)


files <- list.files(pattern = '.*MeltCurve.*')
melt <- list()
for (i in 1:min(6, length(files))) {
  data <- fcn(fread(files[i]))
 melt[[i]] <- data
}
melt <- rbindlist(melt)


files <- list.files(pattern = '.*Results.*')
resu <- list()
for (i in 1:min(6, length(files))) {
  data <- fcn(fread(files[i], skip = '[Results]'))
 resu[[i]] <- data
}
resu <- rbindlist(resu, fill=TRUE)

resu$Target_Name <- gsub("HOPX WITH HA|HOPX HA TEST", "HOPX_HA", resu$Target_Name)
genes <- unique(resu$Target_Name)[which(unique(resu$Target_Name) %!in% c('HPRT', 'HPRT1', 'Hprt', 'Hprt1'))]

# [1] "MYL2"         "MYL7"         "MYH7"         "MYH6"         "TNNI3"       
# [6] "TNNI1"        "HOPX"         "HMGCR"        "CCNB1"        "SRF"         
#[11] "NPPB"         "NFATC"        "SERCA2"       "HOPX HA TEST" "HOPX WITH HA"

meta$Biogroup_Name <- gsub('NA_', '', ifelse(meta$Sample_Name == '', NA, paste0(meta$Biogroup_Name, '_', ifelse(grepl('K', meta$Sample_Name), 'KD', 'CTRL'))))
meta <- meta[is.na(Biogroup_Name) == FALSE] 


meta$Biogroup_Name <- gsub('NA_', '', ifelse(meta$Sample_Name == '', NA, gsub('.*(HOPX|D12Y|W58|KD|WT).*', '\\1', meta$Sample_Name)))
meta$Biogroup_Name <- sub("HOPX", "WT", meta$Biogroup_Name)





# HPC / shell command from the interactive source session (not R code): salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=100GB --job-name=TEST --time=01:00:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

# HPC / shell command from the interactive source session (not R code): module load r/4.2.1-foss-2021a

# The original session entered the interactive R console at this point.
library(data.table)
library(dplyr)
library(tidyr)
library (stringr)

load("formatted_qpcr.RData")

resu$br <- str_extract(resu$Sample_Name, "BR\\d")
resu$Sample_Name <- sub("BR\\d", "", resu$Sample_Name)

meta$br<- str_extract(meta$Sample_Name, "BR\\d")
meta$Sample_Name <- sub("BR\\d", "", meta$Sample_Name)


#get rid of all non alphabetical columns
meta$condition <- gsub("[A-Za-z]", "", meta$Sample_Name)
meta$condition <- gsub("\\b(?:12|58)\\b", "", meta$condition)


resu$Biogroup_Name <- gsub('NA_', '', ifelse(resu$Sample_Name == '', NA, gsub('.*(HOPX|D12Y|W58|KD|WT).*', '\\1', resu$Sample_Name)))
resu$Biogroup_Name <- sub("HOPX", "WT", resu$Biogroup_Name)

resu$condition <- gsub("[A-Za-z]", "", resu$Sample_Name)
resu$condition <- gsub("\\b(?:12|58)\\b", "", resu$condition)



resu$Target_Name <- gsub("HOPX_HA", "HA", resu$Target_Name)
meta$Target_Name <- gsub("HOPX WITH HA|HOPX HA TEST", "HA", meta$Target_Name)

genes <- unique(resu$Target_Name)[which(unique(resu$Target_Name) %!in% c('HPRT', 'HPRT1', 'Hprt', 'Hprt1'))]


# RESULTS
# for each sample, normalise its gene expression by hprt expression: norm = 2**(HPRT_CtMean - Gene_ctMean)

res <- unique(merge(meta, resu,
                     by = c("Well", "condition", "Biogroup_Name", "Target_Name", "br"),
                     suffixes = c(".meta", ".result")))
# Both source tables contain Sample_Name; retain the result-table label where
# available, falling back to the setup metadata label.
if (all(c("Sample_Name.meta", "Sample_Name.result") %in% names(res))) {
  res$Sample_Name <- ifelse(!is.na(res$Sample_Name.result) & nzchar(res$Sample_Name.result),
                            res$Sample_Name.result, res$Sample_Name.meta)
} else if (!"Sample_Name" %in% names(res)) {
  stop("Merged qPCR table has no Sample_Name field; inspect the source exports.")
}
res <- as.data.table(res)
columns_to_keep <- c("Biogroup_Name", "Target_Name", "Ct_Mean", "condition", "br", "Sample_Name")

res <- subset(res, select = columns_to_keep)
res$Name <- paste(res$Biogroup_Name, res$Sample_Name, sep = '_')



hres <- res[Target_Name == 'HPRT', .(Biogroup_Name=Biogroup_Name,condition=condition,br=br, HPRT_Mean = Ct_Mean)]
hres <- unique(hres)
hres <- hres %>%
  group_by(Biogroup_Name, condition, br) %>%
  summarise(HPRT_Mean = mean(HPRT_Mean, na.rm = TRUE)) %>%
  ungroup()
res <- merge(res[res$Target_Name != "HPRT"], hres,
             by = c("Biogroup_Name", "condition", "br"))
res <- as.data.table(res)

res$Norm <- 2**(res$HPRT_Mean - res$Ct_Mean)

# add run identifier column so that they're easy to rbind

# save formatted output
write.table(res,'br1_br7_formatted.txt', sep = '\t', quote = F, col.names = T, row.names = F)
save.image("start_plots.RData")

#############################################################################################

#PLOTTING
res$name <- ifelse(grepl('1:1000', res$condition), 'high', 'low')

res$name <- paste(res$Biogroup_Name,res$name, sep = " ")
res <- as.data.frame(res)

library(scales)
library(patchwork)
library(ggpubr)

# boxplot labelling all samples to see outliers
rawresults <- function(gn, df){ 
	df <- df[df$Target_Name == gn, , drop = FALSE]
	p2 <- ggplot(df, aes(x= name, y= Norm, colour = br)) +
			geom_boxplot(alpha = 0.5, width = 0.1, colour = 'black') + 
			geom_text(aes(label = br),size = 3) +  
			ggtitle(gn) + 
			        labs(x = NULL, y = "fold change")+
			theme_classic() + 
			theme(legend.position = 'none', 
					axis.text = element_text(size = 8),
					axis.title = element_text(size = 10))			        
	return(p2)
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
wrap_plots(lapply(genes, rawresults, df = res))
ggsave("results/all_raw_results.pdf", width=25)





#remove anomaly (?)

res_filter <- res

res_filter$edit <- paste(res_filter$br,res_filter$name, sep = "_")


res_filter <- res_filter[!(res_filter$edit %in% c('BR6_KD high','BR7_WT low','BR6_WT low', 'BR5_WT high', 'BR5_WT low', 'BR5_D12Y high', 'BR5_D12Y low')), , drop = FALSE]
wrap_plots(lapply(genes, rawresults, df = res_filter))
ggsave("results/removed_major_anomalies.pdf", width=20)


# Compare high- and low-concentration conditions within the same biological
# group, biological replicate and target gene. The source code referenced
# res_h/res_l without creating them; derive them explicitly from res_filter.
res_high <- res_filter[grepl("high", res_filter$name, ignore.case = TRUE), , drop = FALSE]
res_low  <- res_filter[grepl("low",  res_filter$name, ignore.case = TRUE), , drop = FALSE]
res_high$condition <- "high"
res_low$condition <- "low"

# Construct a stable pairing key and retain only pairs present in both groups.
res_high$smoosh <- paste(res_high$Biogroup_Name, res_high$br,
                         res_high$Target_Name, sep = "_")
res_low$smoosh <- paste(res_low$Biogroup_Name, res_low$br,
                        res_low$Target_Name, sep = "_")

div <- merge(res_high, res_low, by = "smoosh", suffixes = c(".high", ".low"))
div <- div[!duplicated(div$smoosh), , drop = FALSE]
div$Div <- div$Norm.high / div$Norm.low


	boxsig <- function(gene, rdf = div){ 
	ggboxplot(rdf[rdf$Target_Name.high == gene, , drop = FALSE], x = 'Biogroup_Name.low', y = 'Div', add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}


wrap_plots(lapply(genes, boxsig, rdf = div))




res_high <- res_high %>%
  distinct(Biogroup_Name, br, Target_Name, condition, .keep_all = TRUE)




	boxsig <- function(gene, rdf = res_high){ 
	ggboxplot(rdf[rdf$Target_Name == gene, , drop = FALSE], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}
##########################


boxsig <- function(gene, rdf = res_filter){ 
	ggboxplot(rdf[rdf$Target_Name == gene, , drop = FALSE], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}
wrap_plots(lapply(genes, boxsig, rdf = res_filter))
ggsave("results/comparing_biogroups.pdf", width=20, height=20)





boxsig <- function(gene, rdf = res_filter){ 
	ggboxplot(rdf[rdf$Target_Name == gene, , drop = FALSE], x = 'name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}

wrap_plots(lapply(genes, boxsig, rdf = res_filter))


high <- data.frame()
low <- data.frame()

# Copy rows based on the Description column
for (i in 1:nrow(res_filter)) {
  if (grepl("high", res_filter$name[i], ignore.case = TRUE)) {
    high <- rbind(high, res_filter[i, ])
  } else if (grepl("low", res_filter$name[i], ignore.case = TRUE)) {
    low <- rbind(low, res_filter[i, ])
  }
}

boxsig <- function(gene, rdf = high){ 
	ggboxplot(rdf[rdf$Target_Name == gene, , drop = FALSE], x = 'name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT high') + ggtitle(gene)
}

wrap_plots(lapply(genes, boxsig, rdf = high))
ggsave("results/comparing_high.pdf", width=20, height=20)


##################################################

#doing a graph withut BR6 AND BR7
res_filter2 <- res_filter
res_filter2 <- res_filter[!(res_filter$br %in% c('BR6', 'BR7')), , drop = FALSE]

boxsig <- function(gene, rdf = res_filter2){ 
	ggboxplot(rdf[rdf$Target_Name == gene, , drop = FALSE], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}

wrap_plots(lapply(genes, boxsig, rdf = res_filter2))
ggsave("results/no_br6_br7.pdf", width =20, height=20)

