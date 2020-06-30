require(clusterProfiler)
library(org.Hs.eg.db)
require(dplyr)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')


DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, COAD.DE.rs, PRAD.DE.rs, NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2', 'BRCA.LumB', 'COAD','PRAD', 'NET.PAAD', 'NET.SI')


up.gene                  <- lapply(DE.rs.list,function(x) x$tumor.intrinsic.DE.gene.rs$up.gene) %>% unlist %>% unique
up.gene.matrix           <- matrix(0,nrow = length(up.gene),ncol = length(DE.rs.list))
rownames(up.gene.matrix) <- up.gene
colnames(up.gene.matrix) <- names(DE.rs.list)
for(i in 1:length(DE.rs.list) )  {
    gene   <- DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$up.gene
    cancer <- names(DE.rs.list)[i]
    up.gene.matrix[gene,cancer] <- 1
}
up.gene.freq <- apply(up.gene.matrix,1,sum)
up.gene.freq <- sort(up.gene.freq,decreasing = TRUE)


dn.gene                  <- lapply(DE.rs.list,function(x) x$tumor.intrinsic.DE.gene.rs$dn.gene) %>% unlist %>% unique
dn.gene.matrix           <- matrix(0,nrow = length(dn.gene),ncol = length(DE.rs.list))
rownames(dn.gene.matrix) <- dn.gene
colnames(dn.gene.matrix) <- names(DE.rs.list)
for(i in 1:length(DE.rs.list) )  {
  gene   <- DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$dn.gene
  cancer <- names(DE.rs.list)[i]
  dn.gene.matrix[gene,cancer] <- 1
}
dn.gene.freq <- apply(dn.gene.matrix,1,sum)
dn.gene.freq <- sort(dn.gene.freq,decreasing = TRUE)











# ########### Analyze TF ################
# 
# require(dplyr)
# require(stringr)
# 
# file.name <- "client-side/Data/JASPAR.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# get.TF.name <- function(x) {
#   if(grepl(x = x,pattern='>')){
#     l <- strsplit(x=x,split='\t')  %>% unlist
#     toupper(l[2])
#   }else{
#     NA
#   }
# }
# 
# TF.list <- sapply(line,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# JASPAR.TF.list <- TF.list
# close(conn)
# 
# 
# file.name <- "client-side/Data/CISTROME.factor.line.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# tmp       <- str_extract_all(pattern="id=\"[:alnum:]+\"",string=line[1],simplify = TRUE)
# get.TF.name <- function(x) {
#   x <- str_remove_all(string = x,pattern="\"")
#   x <- str_remove_all(string = x,pattern="id=")
#   x
# }
# TF.list <- sapply(tmp,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# CISTROME.TF.list <- TF.list
# close(conn)
# 
# TF.list <- c(CISTROME.TF.list,JASPAR.TF.list) %>% unique
# 
# 
# c.up.gene        <- names(up.gene.freq)[up.gene.freq >=4]
# c.up.gene.symbol <- mapping.df[c.up.gene,'symbol']
# c.up.TF          <- intersect(c.up.gene.symbol,TF.list)
# 
# 
# 
# 
# c.dn.gene        <- names(dn.gene.freq)[dn.gene.freq >=4]
# c.dn.gene.symbol <- mapping.df[c.dn.gene,'symbol']
# c.dn.TF          <- intersect(c.dn.gene.symbol,TF.list)
# 
# 
# 
# 
# 
# 
# cosine_sim <- function(a, b) crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))
# 
# 
# 
# 
# 
# 
# 
# PRRX1 <- "ENSG00000116132" # now, focus on PRRX1
# EGR1  <- "ENSG00000120738" # now, focus on EGR1
# SPARCL1 <- "ENSG00000152583"
# 
# 
# 
# PRRX1.expr   <- PRI.log2.fpkm.matrix[PRRX1,] 
# q            <- quantile(PRRX1.expr)
# h.sample     <- names(PRRX1.expr)[PRRX1.expr >= q['75%']]
# l.sample     <- names(PRRX1.expr)[PRRX1.expr <= q['25%']]
# 
# 
# rs <- perform.DE.analysis.between.TRE.and.CON(CON.log2.fpkm.matrix = PRI.log2.fpkm.matrix[,h.sample],
#                                               TRE.log2.fpkm.matrix = PRI.log2.fpkm.matrix[,l.sample],
#                                               CON.log2.read.count.matrix = PRI.log2.read.count.matrix[,h.sample],
#                                               TRE.log2.read.count.matrix = PRI.log2.read.count.matrix[,l.sample]
#                                               )
# 
# up.gene <- rownames(rs)[rs$log2FoldChange >  1 & rs$padj < 0.05]
# dn.gene <- rownames(rs)[rs$log2FoldChange < -1 & rs$padj < 0.05]
# 
# 
# w.test.p.value <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
#     wilcox.test(PRI.log2.fpkm.matrix[g,h.sample],PRI.log2.fpkm.matrix[g,l.sample])$p.value  
# }
# names(w.test.p.value) <- c(up.gene,dn.gene)
# 
# de.gene <- names(w.test.p.value)[w.test.p.value < 0.05]
# 
# # m1 <- apply(PRI.log2.fpkm.matrix[de.gene,h.sample],1,median)
# # m2 <- apply(PRI.log2.fpkm.matrix[de.gene,l.sample],1,median)
# # 
# # max.m <- apply(rbind(m1,m2),2,max)
# # de.gene <- names(max.m)[max.m > log2(11)]
# 
# 
# up.gene <- intersect(up.gene,de.gene)
# dn.gene <- intersect(dn.gene,de.gene)
# 
# 
# 
# g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　dn.gene)
# 
# write.csv(x=g.df,quote = FALSE,file='~/Desktop/g.df.csv')
# 
# 
# # cor.with.SPARCL1 <- cor(PRI.log2.fpkm.matrix[c(dn.gene),] %>% t,method='spearman')[SPARCL1,]
# # 
# # #cor.with.SPARCL1 <- cor(PRI.log2.fpkm.matrix[c(NET.SI.DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene),] %>% t,method='spearman')[SPARCL1,]
# # 
# # cor.with.SPARCL1[c.dn.gene] %>% boxplot
# # 
# # 
# # hh <- names(cor.with.SPARCL1)[cor.with.SPARCL1 >= 0.5]
# 
# 
# g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　hh)
# 
# 
# plot(PRI.log2.fpkm.matrix[c(PRRX1,dn.gene[1]),] %>% t)
# 
# y <- PRI.log2.fpkm.matrix[PRRX1,]
# 
# dn.gene <- setdiff(dn.gene,PRRX1)
# r.value <- foreach(g = dn.gene,.combine='c') %do% {
#     df <- data.frame(x=PRI.log2.fpkm.matrix[g,],y)  
#     fit.rs <- lm(df,formula = y ~ x)
#     mad(fit.rs$residuals)
# }
# 
# names(r.value) <- dn.gene



#####################
# median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
# get.gene.id <- function(x) {
#   strsplit(x=x,split='\\.') %>% unlist %>% head(1)
# }
# rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)
# 
# 
# g <- intersect(rownames(median.tpm.matrix),BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),BRCA.LumB.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),COAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),PRAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))


perform.GO.analysis <- function(DE.rs) {
    up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,mapping.df$ensembl_gene_id))
    GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]

    dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,mapping.df$ensembl_gene_id))
    GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]
    
    list(GO.rs.up = GO.rs.1.up,GO.rs.dn = GO.rs.1.dn)
}




rs <- perform.GO.analysis(BRCA.LumB.DE.rs)

rs <- perform.GO.analysis(BRCA.Basal.DE.rs)

rs <- perform.GO.analysis(BRCA.Her2.DE.rs)

rs <- perform.GO.analysis(PRAD.DE.rs)

rs <- perform.GO.analysis(COAD.DE.rs)


rs <- perform.GO.analysis(NET.PAAD.DE.rs)
rs <- perform.GO.analysis(NET.SI.DE.rs)



g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　c.dn.gene)
g.df <- g.df[complete.cases(g.df),]
GO.rs             <- enrichGO(gene=g.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]


