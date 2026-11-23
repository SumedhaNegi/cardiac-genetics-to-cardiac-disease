


library(readxl)
upregulated <- read_xlsx(path = "Desktop/HOPX_bound_genes_TFs.xlsx", sheet = "HOPX bound upregulated genes", col_names = FALSE)

downregulated <- read_xlsx(path = "Desktop/HOPX_bound_genes_TFs.xlsx", sheet = "HOPX bound downregulated genes", col_names = FALSE)


# convert from ensemle to entrez id's 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

install("ChIPpeakAnno")
library(ChIPpeakAnno)


ensemblIDs = upregulated$...1

library(org.Hs.eg.db)
entrezIDs_upregulated <- (as.data.frame(convert2EntrezID(IDs=ensemblIDs, orgAnn="org.Hs.eg.db",ID_type="ensembl_gene_id")))



ensemblIDs = downregulated$...1

entrezIDs_downregulated <- (as.data.frame(convert2EntrezID(IDs=ensemblIDs, orgAnn="org.Hs.eg.db",ID_type="ensembl_gene_id")))
                             
write.table(entrezIDs_downregulated, file = "Desktop/dry lab results/hopx_upregulated_entrezid.txt", quote = FALSE, sep = "/t", col.names = FALSE, row.names = FALSE)
                             






library(devtools)
install_github("wjawaid/enrichR")

library(enrichR)
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)