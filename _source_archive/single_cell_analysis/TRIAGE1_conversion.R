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

df1 <- your data matrix 
df2 <- TRIAGE1(df1, species = "Human", geneid = TRUE, log = FALSE, dir = dir1)

