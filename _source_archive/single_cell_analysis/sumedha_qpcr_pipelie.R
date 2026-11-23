Sumedha version 


salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=100GB --job-name=TEST --time=01:55:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
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





salloc --x11 --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=100GB --job-name=TEST --time=01:00:00 --partition=general srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load r/4.2.1-foss-2021a 

R
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

res <- unique(merge(meta,resu,by = c("Well","condition","Biogroup_Name", "Target_Name", "br")))
columns_to_keep <- c('Biogroup_Name','Target_Name','Ct_Mean', "condition", "br")

res <- subset(res, select = columns_to_keep)
res$Name <- paste(res$Biogroup_Name, res$Sample_Name, sep = '_')



hres <- res[Target_Name == 'HPRT', .(Biogroup_Name=Biogroup_Name,condition=condition,br=br, HPRT_Mean = Ct_Mean)]
hres <- unique(hres)
hres <- hres %>%
  group_by(Biogroup_Name, condition, br) %>%
  summarise(HPRT_Mean = mean(HPRT_Mean, na.rm = TRUE)) %>%
  ungroup()
res <- merge(res[res$Target_Name != 'HPRT'], hres, by = c("Biogroup_Name", "condition", "br"))

res$Norm <- 2**(res$HPRT_Mean - res$Ct_Mean)

# add run identifier column so that they're easy to rbind

# save formatted output
write.table(res,'br1_br7_formatted.txt', sep = '\t', quote = F, col.names = T, row.names = F)
save.image("start_plots.RData")

#############################################################################################

#PLOTTING
res$name <- ifelse(grepl('1:1000', res$condition), 'high', 'low')

res$name <- paste(res$Biogroup_Name,res$name, sep = " ")

library(scales)
library(patchwork)
library(ggpubr)

# boxplot labelling all samples to see outliers
rawresults <- function(gn, df){ 
	df <- df[Target_Name == gn]
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

wrap_plots(lapply(genes, rawresults, df = res))
ggsave("results/all_raw_results.pdf", width=25)





#remove anomaly (?)

res_filter <- res

res_filter$edit <- paste(res_filter$br,res_filter$name, sep = "_")


res_filter <- res_filter[edit %!in% c('BR6_KD high','BR7_WT low','BR6_WT low', 'BR5_WT high', 'BR5_WT low', 'BR5_D12Y high', 'BR5_D12Y low')]
wrap_plots(lapply(genes, rawresults, df = res_filter))
ggsave("results/removed_major_anomalies.pdf", width=20)


#comparimng only high
res_high$condition <- str_extract(res_high$name, "high|low")

#divided dataframe into 2 

res_h
res_l 
res_l$smoosh <- paste0(res_l$Biogroup_Name, "_", res_l$br, "_", res_l$Target_Name)

res_h$smoosh <- paste0(res_h$Biogroup_Name, "_", res_h$br, "_", res_h$Target_Name)





res_h <- res_h %>%  mutate(Div = if (res_h$smoosh %in% res_l$smoosh, (res_h$Norm/res_l$Norm))
)                      

# Display the result_df



res_h <- res_h %>%
  inner_join(res_l, by = "smoosh") %>%
  mutate(Div = Norm.x / Norm.y) %>%
  select(-Norm.x, -Norm.y)


div <- div %>%
  distinct(smoosh, .keep_all = TRUE)



	boxsig <- function(gene, rdf = div){ 
	ggboxplot(rdf[Target_Name.x == gene], x = 'Biogroup_Name.y', y = 'Div', add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}


wrap_plots(lapply(genes, boxsig, rdf = div)




res_h$Norm <- 2**(res_h$HPRT_Mean.y - res_h$Ct_Mean.y)

res_high <- res_high %>%
  distinct(Norm, .keep_all = TRUE)




	boxsig <- function(gene, rdf = res_high){ 
	ggboxplot(rdf[Target_Name == gene], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}
##########################


boxsig <- function(gene, rdf = res_filter){ 
	ggboxplot(rdf[Target_Name == gene], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}
wrap_plots(lapply(genes, boxsig, rdf = res_filter))
ggsave("results/comparing_biogroups.pdf", width=20, height=20)





boxsig <- function(gene, rdf = res_filter){ 
	ggboxplot(rdf[Target_Name == gene], x = 'name', y = 'Norm', colour='br' ,add = c('jitter')) + 
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
	ggboxplot(rdf[Target_Name == gene], x = 'name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT high') + ggtitle(gene)
}

wrap_plots(lapply(genes, boxsig, rdf = high))
ggsave("results/comparing_high.pdf", width=20, height=20)


##################################################

#doing a graph withut BR6 AND BR7
res_filter2 <- res_filter
res_filter2 <- res_filter[br %!in% c('BR6', 'BR7')]

boxsig <- function(gene, rdf = res_filter2){ 
	ggboxplot(rdf[Target_Name == gene], x = 'Biogroup_Name', y = 'Norm', colour='br' ,add = c('jitter')) + 
	stat_compare_means(label = 'p.signif', ref.group = 'WT') + ggtitle(gene)
}

wrap_plots(lapply(genes, boxsig, rdf = res_filter2))
ggsave("results/no_br6_br7.pdf", width =20, height=20)

