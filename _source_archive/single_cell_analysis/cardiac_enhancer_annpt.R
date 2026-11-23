#Cardiac_enhancer with 50kb window


x <- head("/QRISdata/Q5149/computation/heritability/hopx_annot/hg19_cardiac_enhancer (1).bed")

x$V2 <- x$V2 - 50000
x$V3 <- x$V3 + 50000


snp <- fread("/scratch/user/uqsnegi/ukbEUR_Imputed/snp.info")


        dim(snp)
        #[1] 7356518      10
        snp2 <- snp[,c(1,2,5,6,7)]
        colnames(snp2) <- c ("Chromosome", "rsid", "PhysPos", "A1", "A2")



library(data.table)
library(GenomicRanges)

snp2$Chromosome <- as.character(snp2$Chromosome)
snp2_ranges <- GRanges(seqnames = snp2$Chromosome, ranges = IRanges(start = snp2$PhysPos, end = snp2$PhysPos))

df <- x
 df$V1 <- as.character(df$V1)
 colnames(df) <- c("chr", "start","stop")
 df$chr <- sub("^chr", "", df$chr)

 df_ranges <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$stop))

 snp2$enhancer <- as.integer(overlapsAny(snp2_ranges, df_ranges))


 	