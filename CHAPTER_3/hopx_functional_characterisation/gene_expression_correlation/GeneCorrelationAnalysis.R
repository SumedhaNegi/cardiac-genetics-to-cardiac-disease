# =============================================================================
# EXAMINER GUIDE
# Script: GeneCorrelationAnalysis.R
# Purpose: Creates ordered gene-expression correlation plots; source comments refer to cardiomyocyte datasets.
# Provenance: Adapted as an annotated copy from the supplied source archive.
# Inputs/outputs: Confirm exact filenames, genome build, species and output paths
# in the original analysis environment before interpreting or re-running.
# Reproducibility: This is a methodological record, not a guaranteed standalone
# pipeline. UQ HPC paths, modules and external datasets may need local equivalents.
# Review note: No analysis was executed or statistically re-analysed during curation.
# =============================================================================

# CODE MAP
# Analyses correlations among gene-expression measurements and generates gene-level comparison plots/statistics.
#

# ordered gene expression correlation plots 
# want to do this with CSC and Seidman cardios data (before sept4)

library(ggplot2)
library(patchwork)
library(ggrastr)
library(ggrepel)

# you need to first read in the data and have it be a data.table with cells = rows and genes = columns. genes need to be in column names

runCorr <- function(df, findGene){
	gene <- grep(findGene, colnames(df), value = TRUE)[1]
	print(paste0('Gene being used: ', gene))

	df <- df[order(-df[[gene]]),]
	df$ordcells <- 1:nrow(df)

	corrs <- unlist(lapply(colnames(df), function(gn){cor(df[,c(gene, gn)], method = 'pearson')[2]}))
	names(corrs) <- colnames(df)
	outs <- list(df, corrs)
	names(outs) <- c('df', 'corrs')
	return(outs)
}

plt <- function(gene, corr, do.grep = TRUE, suffix = ''){
	if(do.grep){
		gene <- grep(paste0(gene, suffix), names(corr$corrs), value = TRUE)[1]
		print(paste0('Gene being used: ', gene))
	}
	ggplot(corr$df, aes(ordcells, !!sym(gene))) + 
	ggtitle(round(corr$corrs[[gene]], 4)) + 
	ylab(gsub('_.*', '', gene))
}

pltn <- function(corr, n = 16, low = FALSE, rast = FALSE){
	if (low == FALSE){
		topn <- head(names(corr$corrs[order(-corr$corrs)]), n)} else {
		topn <- head(names(corr$corrs[order(corr$corrs)][-1]),n)}

	plts <- lapply(topn, plt, corr = corr, do.grep = FALSE)
	p <- wrap_plots(plts) & theme_minimal() & theme(axis.text.x = element_blank())
	if (rast == TRUE){
		p <- p & geom_point_rast(size = 1)
	} else {p <- p & geom_point(size = 0.1)}
	return(p)
}

pltgenes <- function(genelist, corr, do.grep = TRUE, rast = FALSE, suffix = ''){
	if (do.grep){
		usegenes <- grep(paste(paste0(genelist, suffix), collapse = '|'), names(corr$corrs), value = TRUE)
	} else {usegenes <- genelist}
	
	p <- wrap_plots(lapply(usegenes, plt, corr = corr, do.grep = FALSE)) & theme_minimal() & theme(axis.text.x = element_blank()) 
	if (rast == TRUE){
		p <- p & geom_point_rast(size = 1)
	} else {p <- p & geom_point(size = 0.1)}

	return(p)
}


plthist <- function(corr, bins = 100){
	ggplot(data.frame('correlation' = corr$corrs), aes(x = correlation)) + geom_histogram(bins = bins) + ggtitle(paste0('Distribution of correlations in ', bins, ' bins')) + theme_classic()
}


mkot <- function(corr, saveprefix = NULL){
	outtable <- data.frame('gene' = names(corr$corrs), 'correlation' = corr$corrs)[1:(length(corr$corrs)-1),]
	outtable <- outtable[order(-outtable$correlation),]

	if (length(saveprefix) != 0){
		write.table(outtable, paste0(saveprefix, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
	}
	return(outtable)
}


saveOuts <- function(prefix, corr = NULL, df = NULL, gene = NULL, plots = TRUE, n = 25, h = 1750, w = 1750, r = 200, bins = 100){

	if (length(corr) == 0){corr <- runCorr(df, gene)}
	write.table(mkot(corr), paste0(prefix, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

	if (plots == TRUE){
		jpeg(paste0(prefix, '_top', n, '.jpeg'), height = h, width = w, res = r)
		print(pltn(corr, n = n))
		dev.off()

		jpeg(paste0(prefix, '_low', n, '.jpeg'), height = h, width = w, res = r)
		print(pltn(corr, n = n, low = TRUE))
		dev.off()

		jpeg(paste0(prefix, '_corrdist.jpeg'), height = 750, width = 1500, r = 200)
		suppressWarnings(print(plthist(corr, bins = bins)))
		dev.off()
	}
	
}


saveMany <- function(genes, df, saveprefix, collapsepatt = '_|', dir_ = dir, plots = FALSE){
	for (gene in grep(paste(genes, collapse = collapsepatt), colnames(df), value = TRUE)){
		# thiscorr <- runCorr(df, gene)
		saveOuts(paste0(dir_, saveprefix, '_', gene, '_corr'), df = df, gene = gene, plots = plots)
	}
}


# labelled ordered correlation plots (pltord)
pltord <- function(gene, plotgene = 'HOPX', outt = NULL, df = NULL, saveprefix = NULL){
	if (length(outt) == 0){
		corr <- runCorr(df, gene)
		outt <- mkot(corr)}
	outt$gene <- factor(outt$gene, levels = outt$gene)
	if (length(which(is.na(outt$correlation))) > 0){outt <- outt[-(which(is.na(outt$correlation))),]}
	pltgene <- grep(plotgene, row.names(outt), value = TRUE)[1]
	print(paste0('Gene being plotted: ', pltgene))

	p <- ggplot(outt, aes(gene, correlation, colour = ifelse(gene == pltgene, 'plotgene', 'not'))) +
			geom_point_rast(size = ifelse(outt$gene == pltgene, 3, 1)) + 
			geom_vline(xintercept = 0.1*nrow(outt), linetype = 'dashed') +
			geom_vline(xintercept = 0.9*nrow(outt), linetype = 'dashed') +
			geom_label_repel(data = outt[outt$gene == pltgene,], aes(gene, correlation, label = gsub('_.*|--.*', '', gene), size = 0.1), box.padding = 1, segment.colour = NA) +
			coord_cartesian(clip = 'off') +
			scale_colour_manual(values = c('black', 'red')) + 
			theme_classic() +
			theme(legend.position = 'none',
				  axis.text.x = element_blank(),
				  axis.ticks.x = element_blank(),
				  axis.line = element_line(linewidth = 0.5)) + 
			expand_limits(x = c(-150,nrow(outt)+150), y = 0) +
			ggtitle(gene)

	if (length(saveprefix) != 0){
		pdf(paste0(saveprefix,'_', gsub('\\^|-.*|_.*', '', pltgene), '_corrRank.pdf'))
		print(p)
		dev.off()
	}
	return(p)
}
#pltord('Hopx', 'Hopx', df = bsw_df)

# # example use
# day30 <- readRDS('/scratch/90days/uqsshen2/ANALYSIS/day30_separate_final.RDS')
# df30 <- data.frame(t(day30$Exprs))
# corr30 <- runCorr(df30, 'HOPX')
# head(as.data.frame(corr30$corrs))

# pltn(corr30, n = 16)
# pltn(corr30, low = TRUE, n = 16)

# saveOuts('/scratch/90days/uqsshen/ANALYSIS/corrplots', corr30)

