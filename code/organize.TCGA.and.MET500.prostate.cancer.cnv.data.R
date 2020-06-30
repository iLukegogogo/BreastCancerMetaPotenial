load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
require(CNTools)



################# organize hg19 gene coordinates information #####################
hg19.RefSeq.gene.coordinates <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hg19.RefSeq.gene.coordinates.txt", stringsAsFactors=FALSE)
flag                         <- grepl(x=hg19.RefSeq.gene.coordinates$chrom,pattern = '_')
hg19.RefSeq.gene.coordinates <- hg19.RefSeq.gene.coordinates[!flag,] #Let us remove the genes which are on the un-assembled contigs!
tmp                          <- ddply(hg19.RefSeq.gene.coordinates,.(name2),function(x) data.frame(chrom=x$chrom[1], start=min(x$txStart),end=max(x$txEnd),geneid=x$name2[1],genename=x$name2[1]))
tmp$name2                    <- NULL
hg19.gene.info               <- tmp
hg19.gene.info$chrom         <- gsub(x=hg19.gene.info$chrom,pattern='chr',replacement = '')
hg19.gene.info$chrom         <- as.character(hg19.gene.info$chrom)
hg19.gene.info$geneid        <- as.character(hg19.gene.info$geneid)
hg19.gene.info$genename      <- as.character(hg19.gene.info$genename)

################# organize TCGA PRAD cnv data #####################
data               <- read.table("client-side/Data/gdac.broadinstitute.org_PRAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/PRAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", stringsAsFactors=FALSE,header = TRUE)
data               <- data[,c('Sample','Chromosome','Start','End','Segment_Mean')]
colnames(data)     <- c('ID','chrom','loc.start','loc.end','seg.mean')

get.TCGA.sample.id <- function(x) {
  tmp <- strsplit(x=x,split='-') %>% unlist 
  s   <-paste(tmp[1:4],collapse = '-',sep='')
  s   <- gsub(x = s,pattern = 'A$',replacement = '')
  s   <- gsub(x = s,pattern = 'B$',replacement = '')
  s
}
data$ID                               <- sapply(data$ID,get.TCGA.sample.id)
TCGA.seg                              <- CNSeg(data)
TCGA.prostate.cancer.gene.cnv           <- getRS(object = TCGA.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
TCGA.prostate.cancer.gene.cnv           <- TCGA.prostate.cancer.gene.cnv[,5:ncol(TCGA.prostate.cancer.gene.cnv)]
rownames(TCGA.prostate.cancer.gene.cnv) <- TCGA.prostate.cancer.gene.cnv$genename
TCGA.prostate.cancer.gene.cnv$genename  <- NULL
TCGA.prostate.cancer.cnv.matrix           <- as.matrix(TCGA.prostate.cancer.gene.cnv)




################# organize MET500 cnv data #####################

get.id <- function(x){
  l <- strsplit(x = x,split='\\.')  %>% unlist 
  l[1]
}
flag                                    <- grepl(x=MET500.sample.meta$cancer.type,pattern='Prostate Adenocarcinoma') 
prostate.cancer.sample.MET500.id          <- MET500.sample.meta$MET500.id[flag] %>% as.character %>% unique
MET500.cnv.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.cnv.data                         <- MET500.cnv.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.cnv.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.cnv.data$ID                      <- sapply(MET500.cnv.data$ID,get.id)
MET500.prostate.cancer.cnv.data           <- MET500.cnv.data[MET500.cnv.data$ID %in% prostate.cancer.sample.MET500.id,]
MET500.prostate.cancer.seg                <- CNSeg(MET500.prostate.cancer.cnv.data)
MET500.prostate.cancer.gene.cnv           <- getRS(object = MET500.prostate.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
MET500.prostate.cancer.gene.cnv           <- MET500.prostate.cancer.gene.cnv[,5:ncol(MET500.prostate.cancer.gene.cnv)]
rownames(MET500.prostate.cancer.gene.cnv) <- MET500.prostate.cancer.gene.cnv$genename
MET500.prostate.cancer.gene.cnv$genename  <- NULL
MET500.prostate.cancer.cnv.matrix         <- as.matrix(MET500.prostate.cancer.gene.cnv)

########################################################
save(file='client-side/output/organize.TCGA.and.MET500.prostate.cancer.cnv.data.R.output/organize.TCGA.and.MET500.prostate.cancer.cnv.data.RData',list=c('TCGA.prostate.cancer.cnv.matrix','MET500.prostate.cancer.cnv.matrix','hg19.gene.info'))
