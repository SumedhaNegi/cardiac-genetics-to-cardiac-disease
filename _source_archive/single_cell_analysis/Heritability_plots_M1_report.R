data visualization of HERITABILITY 

qsub -XI -l mem=120gb,walltime=4:00:00
##
module load R/4.1.2
export R_LIBS=/scratch/90days/uqsnegi/RLIBS

Compiling all the results 

grep 'PerSnpHsq' $(find . -name "*.parSetRes") > ce_hopx.results.txt

R
#con1 

con1 <- read.csv("/con1/con1_persnp.txt", sep="", header = F)


c <- case1[,c(1,3,4,5)]
c2 <- c[!(c$V3 =="opposite"),]


colnames(c2) <- c('trait','annot', 'Per_SNP_enrichment', 'error')

c2$trait<-gsub("|./invnorm_|.ma.parSetRes:|./|_cleaned","",as.character(c2$trait))

 
 library(ggplot2)
 med <- median(c2$V4)

 pdf("proteincoding_persnp.pdf")
 ggplot(c2, aes(Per_SNP_enrichment, trait)) + geom_point() + geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error)) + geom_vline(xintercept = med) 
dev.off()



####################################

#dividing the data 

hopx_med <- median(c1$Per_SNP_enrichment)
genome_med <- median(c2$Per_SNP_enrichment)

 ggplot(h1, aes(Per_SNP_enrichment, trait, colour=annot)) + geom_point() + 
 geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error, width=0.5)) +
 geom_vline(xintercept = hopx_med, colour="purple", alpha=0.3) + 
 geom_vline(xintercept=genome_med, colour="grey", alpha=0.5)+
 scale_color_manual(values = c("purple","grey"))+ 
 theme_classic()
 

 ggsave("enhancer_all_traits.pdf", width = 8, height = 8, dpi = 1000)
   
lv <- c("bsa", "lv_indexed", "lvef", "svi", "ratio")

rv <- c("RVEDV_indexed", "ratio", "RVEF", "RVESV_indexed", "RV_SV_indexed")

oft <- c("aa_diam_indexed", "pa_ao")


lv_only <-subset(h1,grepl(paste0(lv, collapse = "|"), trait))

rv_only <-subset(h1,grepl(paste0(rv, collapse = "|"), trait))

oft_only <-subset(h1,grepl(paste0(oft, collapse = "|"), trait))

############################################################


######### DO a position switch next time 


#plotting individual lv

genome <- lv_only[grepl("^overlapping_with_cardiac", lv_only$annot), ]
prot_med <- median(genome$Per_SNP_enrichment)
hopx <- lv_only[grepl("^overlapping_with_hopx", lv_only$annot), ]
hopx_med <- median(hopx$Per_SNP_enrichment)


 ggplot(lv_only, aes(Per_SNP_enrichment, trait, colour=annot)) + geom_point() + 
 geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error, width=0.3)) +
 geom_vline(xintercept = hopx_med, colour="purple", alpha=0.3) + 
 geom_vline(xintercept=prot_med, colour="grey", alpha=0.5)+
 scale_color_manual(values = c("purple","grey"))+ 
 theme_classic()
ggsave("lv_traits_enhancer.pdf", width = 8, height = 6, dpi = 1000)


genome <- rv_only[grepl("^overlapping_with_cardiac", rv_only$annot), ]
prot_med <- median(genome$Per_SNP_enrichment)
hopx <- rv_only[grepl("^overlapping_with_hopx", rv_only$annot), ]
hopx_med <- median(hopx$Per_SNP_enrichment)


  ggplot(rv_only, aes(Per_SNP_enrichment, trait, colour=annot)) + geom_point() + 
 geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error, width=0.3)) +
 geom_vline(xintercept = hopx_med, colour="purple", alpha=0.3) + 
 geom_vline(xintercept=prot_med, colour="grey", alpha=0.5)+
 scale_color_manual(values = c("purple","grey"))+ 
 theme_classic()

ggsave("rv_traits_enhancer.pdf", width = 8, height = 4.5, dpi = 1000)




	
genome <- oft_only[grepl("^overlapping_with_cardiac", oft_only$annot), ]
prot_med <- median(genome$Per_SNP_enrichment)
hopx <- oft_only[grepl("^overlapping_with_hopx", oft_only$annot), ]
hopx_med <- median(hopx$Per_SNP_enrichment)
  ggplot(oft_only, aes(Per_SNP_enrichment, trait, colour=annot)) + geom_point() + 
 geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error, width=0.1)) +
 geom_vline(xintercept = hopx_med, colour="purple", alpha=0.3) + 
 geom_vline(xintercept=prot_med, colour="grey", alpha=0.5)+
 scale_color_manual(values = c("purple","grey"))+ 
 theme_classic()
ggsave("oft_traits_enhancer.pdf", width = 8, height = 3, dpi = 1000)


#############################################################################



genome <- h1[grepl("^overlapping_with_cardiac", h1$annot), ]
prot_med <- median(genome$Per_SNP_enrichment)
hopx <- h1[grepl("^overlapping_with_hopx", h1$annot), ]
hopx_med <- median(hopx$Per_SNP_enrichment)


 ggplot(h1, aes(Per_SNP_enrichment, trait, colour=annot)) + geom_point() + 
 geom_errorbar(aes(xmin = Per_SNP_enrichment - error, xmax = Per_SNP_enrichment + error, width=0.5)) +
 geom_vline(xintercept = hopx_med, colour="purple", alpha=0.3) + 
 geom_vline(xintercept=genome_med, colour="grey", alpha=0.5)+
 scale_color_manual(values = c("purple","grey"))+ 
 theme_classic()

ggsave("enhancer_all_traits_enhancer.pdf", width = 8, height = 6, dpi = 1000)











