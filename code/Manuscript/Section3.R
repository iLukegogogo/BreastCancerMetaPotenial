require(ggplot2)
require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')


################################################################################   
#----------------------------- Figure 3 ---------------------------------------#
################################################################################  

##############################################
# Fig 3a: expressin pattern of upregulated tumor-intrinsic DE genes (LuminalB)
##############################################
median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
get.gene.id <- function(x) {
  strsplit(x=x,split='\\.') %>% unlist %>% head(1)
}
rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)


draw.color.bar <- function(value.range, color.range){
  f <- colorRamp2(value.range,color.range)
  x <- sapply(value.range, f)
  plot(c(0, 40*length(x)), c(0, 100), type = "n")
  i <- 0:(length(x)-1)
  rect(i*40, 0, (i+1)*40, 100, col=x)
}
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/color.bar.pdf',width = 20,height=20)
draw.color.bar(value.range = seq(from=-5,to=5,by=0.5),
               color.range = colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
)
dev.off()

lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x

lumb.gene  <- intersect(c(lumb.up.gene),  rownames(median.tpm.matrix))
her2.gene  <- intersect(c(her2.up.gene),  rownames(median.tpm.matrix))
basal.gene <- intersect(c(basal.up.gene), rownames(median.tpm.matrix))


tmp           <- apply((median.tpm.matrix[lumb.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/lumb.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()




##############################################
# Fig 3b: comparision GO annotation results before/after DEBoost
##############################################
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


perform.GO.enrichment <- function(DEBoost.gene,DESeq2.gene){
    DEBoost.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DEBoost.gene)
    DEBoost.GO.BP               <- enrichGO(gene=DEBoost.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    DEBoost.GO.BP               <- DEBoost.GO.BP[DEBoost.GO.BP$p.adjust < 0.05 & DEBoost.GO.BP$Count >=5,]             
  
    DESeq2.gene.annotation.df   <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DESeq2.gene)
    DESeq2.GO.BP                <- enrichGO(gene=DESeq2.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    DESeq2.GO.BP                <- DESeq2.GO.BP[DESeq2.GO.BP$p.adjust < 0.05 & DESeq2.GO.BP$Count >=5,]             

    GO.gene.symbol <- sapply(DESeq2.GO.BP$geneID , function(x) strsplit(x,split='/') %>% unlist) %>% unlist %>% unique
    r              <- 1- intersect(GO.gene.symbol,DEBoost.gene.annotation.df$SYMBOL) %>% length / length(GO.gene.symbol)
    list(DEBoost.GO.BP=DEBoost.GO.BP,DESeq2.GO.BP=DESeq2.GO.BP,uDE.gene.ratio = r)
}
  
# Basal-like, upregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.basal$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.basal$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.basal)[flag]
basal.up.rs                 <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)

# Basal-like, dnregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.basal$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.basal$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.basal)[flag]
basal.dn.rs                 <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)


# Lumb, upregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.lumb$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.lumb)[flag]
lumb.up.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)

# Lumb, dnregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.lumb$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.lumb)[flag]
lumb.dn.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)


# Her2, upregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.her2$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.her2$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.her2)[flag]
her2.up.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)

# Her2, dnregulated gene
DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv", stringsAsFactors=FALSE)$x
flag                        <- de.res.metastasis.liver.vs.breast.her2$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.her2$padj < 0.05
DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.her2)[flag]
her2.dn.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)














################################################################################   
#----------------------------- Figure S3 ---------------------------------------#
################################################################################  

##############################################
# Fig S3a, S3b: expressin pattern of upregulated tumor-intrinsic DE genes (LuminalB)
##############################################
tmp           <- apply((median.tpm.matrix[her2.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file= '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/her2.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()

tmp           <- apply((median.tpm.matrix[basal.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/basal.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()













