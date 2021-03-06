#Aim: DE analysis between primary and metastatic breast cancer (liver) in subtype-specific manner

require(DESeq2)
require(dplyr)
library (VennDiagram)
require(gplots)
require(foreach)
source('client-side/code/util.R')

protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


################################################################################################################
#Perform DE analysis between normal tissues, liver vs breast
################################################################################################################
load('server-side//RData//Breast - Mammary Tissue.RData')
GTex.breast.log2.read.count.matrix <- log2.read.count.matrix
GTex.breast.log2.fpkm.matrix       <- log2.fpkm.matrix
breast.expressed.gene              <- get.expressed.gene(GTex.breast.log2.fpkm.matrix)


load('server-side//RData//Liver.RData')
GTex.liver.log2.read.count.matrix <- log2.read.count.matrix
GTex.liver.log2.fpkm.matrix       <- log2.fpkm.matrix
liver.expressed.gene              <- get.expressed.gene(GTex.liver.log2.fpkm.matrix)


expressed.gene <- intersect(c(breast.expressed.gene,liver.expressed.gene) %>% unique,protein.coding.gene.id)
expr.matrix    <- cbind(GTex.breast.log2.read.count.matrix[expressed.gene,],GTex.liver.log2.read.count.matrix[expressed.gene,])
expr.matrix    <- 2^expr.matrix - 1
df             <- data.frame(condition=c(rep(x='breast',times=GTex.breast.log2.read.count.matrix %>% ncol), rep(x='liver',times=GTex.liver.log2.read.count.matrix %>% ncol)))
df$condition   <- factor(df$condition,levels = c('breast','liver'))
dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                       colData = df,
                                       design = ~ condition)

dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','liver','breast')) %>% as.data.frame
res <- res[order(res$pvalue),]
res <- res[complete.cases(res),]
cut.off                        <- 0.05
de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange >  1  & res$padj < cut.off]
de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1  & res$padj < cut.off]
de.res.liver.vs.breast         <- res



#######################################################################################################################################################################################
###### Prepare the data 
#######################################################################################################################################################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')

TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 

#######################################################################################################################################################################################
###### Function to perform DE analysis
#######################################################################################################################################################################################
perform.DE.analysis <- function(){
    g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
    g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
    expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
  
    expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
    expr.matrix    <- 2^expr.matrix - 1
  
    purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
    purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
  
    df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
    df$purity      <- c(purity.MET500,purity.TCGA)
    df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
  
    dds           <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                            colData = df,
                                            design = ~ condition + purity)
    dds <- DESeq(dds)
    res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
    res <- res[order(res$pvalue),]
    res <- res[complete.cases(res),]     
    res
}

fpkm.cut.off <- log2(2+1)


#######################################################################################################################################################################################
###### DE analysis, Basal subtype
#######################################################################################################################################################################################
MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
res            <- perform.DE.analysis()

##### remove DE genes identified from normal tissue DE analysis ############3
de.res.metastasis.liver.vs.breast.basal.up.gene <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
de.res.metastasis.liver.vs.breast.basal.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
basal.up.gene                                   <- setdiff(de.res.metastasis.liver.vs.breast.basal.up.gene,de.res.liver.vs.breast.up.gene)
basal.dn.gene                                   <- setdiff(de.res.metastasis.liver.vs.breast.basal.dn.gene,de.res.liver.vs.breast.dn.gene)
de.res.metastasis.liver.vs.breast.basal         <- res

######  wilcoxon rank test to assign p-value ###########
purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$purity      <- c(purity.MET500,purity.TCGA)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
wilcox.test.p.value.vec <- foreach(g=c(basal.up.gene,basal.dn.gene),.combine='c') %do% {
    df$expr           <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
    df$expr           <- df$expr + 0.0001
    loess.fit         <- loess(data = df,formula=expr~purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    adjusted.expr     <- loess.fit$residuals
    wilcox.test(adjusted.expr[f1],adjusted.expr[f2])$p.value
}
names(wilcox.test.p.value.vec) <- c(basal.up.gene,basal.dn.gene)
basal.wilcox.test.p.value.vec  <- wilcox.test.p.value.vec
# deseq.test.p.value.vec         <- de.res.metastasis.liver.vs.breast.basal[c(basal.up.gene,basal.dn.gene),'pvalue']
# names(deseq.test.p.value.vec)  <- c(basal.up.gene,basal.dn.gene)

######## remove genes whose DE appears to be due to outlier
basal.up.gene.p.value <- wilcox.test.p.value.vec[basal.up.gene]
basal.dn.gene.p.value <- wilcox.test.p.value.vec[basal.dn.gene]
basal.up.gene         <- names(basal.up.gene.p.value)[basal.up.gene.p.value < 0.05]
basal.dn.gene         <- names(basal.dn.gene.p.value)[basal.dn.gene.p.value < 0.05]


######## I only want genes whose expr is high in MET500 samples
m.value               <- apply(MET500.log2.fpkm.matrix[basal.dn.gene,MET500.sample],1,median)
basal.dn.gene         <- basal.dn.gene[m.value > fpkm.cut.off]
m.value               <- apply(MET500.log2.fpkm.matrix[basal.up.gene,MET500.sample],1,median)
basal.up.gene         <- basal.up.gene[m.value > fpkm.cut.off]

write.csv(x=basal.up.gene,file = 'client-side/output/DE.breast.cancer.R.output/basal.up.csv',quote=FALSE)
write.csv(x=basal.dn.gene,file = 'client-side/output/DE.breast.cancer.R.output/basal.dn.csv',quote=FALSE)


#LPA receptor 1 mediates LPA-induced ovarian cancer metastasis: an in vitro and in vivo study  LPAR2 de
#https://www.biorxiv.org/content/biorxiv/early/2015/08/04/023911.full.pdf  HCFC1 and MYC
#######################################################################################################################################################################################
###### DE analysis, Her2 subtype
#######################################################################################################################################################################################
MET500.sample <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample   <- pure.TCGA.breast.cancer.polyA.Her2.sample
res           <- perform.DE.analysis()


de.res.metastasis.liver.vs.breast.her2.up.gene  <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
de.res.metastasis.liver.vs.breast.her2.dn.gene  <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
her2.up.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.her2.up.gene,de.res.liver.vs.breast.up.gene)
her2.dn.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.her2.dn.gene,de.res.liver.vs.breast.dn.gene)
de.res.metastasis.liver.vs.breast.her2          <- res


purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$purity      <- c(purity.MET500,purity.TCGA)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
wilcox.test.p.value.vec <- foreach(g=c(her2.up.gene,her2.dn.gene),.combine='c') %do% {
    df$expr           <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
    df$expr           <- df$expr + 0.0001
    loess.fit         <- loess(data = df,formula=expr~purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    adjusted.expr     <- loess.fit$residuals
    wilcox.test(adjusted.expr[f1],adjusted.expr[f2])$p.value
}
names(wilcox.test.p.value.vec) <- c(her2.up.gene,her2.dn.gene)
her2.wilcox.test.p.value.vec   <- wilcox.test.p.value.vec

# deseq.test.p.value.vec         <- de.res.metastasis.liver.vs.breast.her2[c(her2.up.gene,her2.dn.gene),'pvalue']
# names(deseq.test.p.value.vec)  <- c(her2.up.gene,her2.dn.gene)


her2.up.gene.p.value <- wilcox.test.p.value.vec[her2.up.gene]
her2.dn.gene.p.value <- wilcox.test.p.value.vec[her2.dn.gene]
her2.up.gene         <- names(her2.up.gene.p.value)[her2.up.gene.p.value < 0.05]
her2.dn.gene         <- names(her2.dn.gene.p.value)[her2.dn.gene.p.value < 0.05]

m.value              <- apply(MET500.log2.fpkm.matrix[her2.dn.gene,MET500.sample],1,median)
her2.dn.gene         <- her2.dn.gene[m.value > fpkm.cut.off]
m.value              <- apply(MET500.log2.fpkm.matrix[her2.up.gene,MET500.sample],1,median)
her2.up.gene         <- her2.up.gene[m.value > fpkm.cut.off]

write.csv(x=her2.up.gene,file = 'client-side/output/DE.breast.cancer.R.output/her2.up.csv',quote=FALSE)
write.csv(x=her2.dn.gene,file = 'client-side/output/DE.breast.cancer.R.output/her2.dn.csv',quote=FALSE)


#######################################################################################################################################################################################
######  DE analysis, Luminal A subtype
#######################################################################################################################################################################################
MET500.sample <- c(MET500.breast.cancer.polyA.LumA.sample)
MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample   <- pure.TCGA.breast.cancer.polyA.LumA.sample
res           <- perform.DE.analysis()


de.res.metastasis.liver.vs.breast.luma.up.gene  <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
de.res.metastasis.liver.vs.breast.luma.dn.gene  <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
luma.up.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.luma.up.gene,de.res.liver.vs.breast.up.gene)
luma.dn.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.luma.dn.gene,de.res.liver.vs.breast.dn.gene)
de.res.metastasis.liver.vs.breast.luma          <- res


purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$purity      <- c(purity.MET500,purity.TCGA)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
wilcox.test.p.value.vec <- foreach(g=c(luma.up.gene,luma.dn.gene),.combine='c') %do% {
    df$expr           <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
    df$expr           <- df$expr + 0.0001
    loess.fit         <- loess(data = df,formula=expr~purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    adjusted.expr     <- loess.fit$residuals
    wilcox.test(adjusted.expr[f1],adjusted.expr[f2])$p.value
}
names(wilcox.test.p.value.vec) <- c(luma.up.gene,luma.dn.gene)
luma.wilcox.test.p.value.vec   <- wilcox.test.p.value.vec

# deseq.test.p.value.vec         <- de.res.metastasis.liver.vs.breast.luma[c(luma.up.gene,luma.dn.gene),'pvalue']
# names(deseq.test.p.value.vec)  <- c(luma.up.gene,luma.dn.gene)


luma.up.gene.p.value <- wilcox.test.p.value.vec[luma.up.gene]
luma.dn.gene.p.value <- wilcox.test.p.value.vec[luma.dn.gene]
luma.up.gene         <- names(luma.up.gene.p.value)[luma.up.gene.p.value < 0.05]
luma.dn.gene         <- names(luma.dn.gene.p.value)[luma.dn.gene.p.value < 0.05]

m.value              <- apply(MET500.log2.fpkm.matrix[luma.dn.gene,MET500.sample],1,median)
luma.dn.gene         <- luma.dn.gene[m.value > fpkm.cut.off]
m.value              <- apply(MET500.log2.fpkm.matrix[luma.up.gene,MET500.sample],1,median)
luma.up.gene         <- luma.up.gene[m.value > fpkm.cut.off]

write.csv(x=luma.up.gene,file = 'client-side/output/DE.breast.cancer.R.output/luma.up.csv',quote=FALSE)
write.csv(x=luma.dn.gene,file = 'client-side/output/DE.breast.cancer.R.output/luma.dn.csv',quote=FALSE)


#######################################################################################################################################################################################
#DE analysis, Luminal B subtype
#######################################################################################################################################################################################
MET500.sample <- c(MET500.breast.cancer.polyA.LumB.sample)
MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample   <- pure.TCGA.breast.cancer.polyA.LumB.sample
res           <- perform.DE.analysis()


de.res.metastasis.liver.vs.breast.lumb.up.gene  <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
de.res.metastasis.liver.vs.breast.lumb.dn.gene  <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
lumb.up.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.lumb.up.gene,de.res.liver.vs.breast.up.gene)
lumb.dn.gene                                    <- setdiff(de.res.metastasis.liver.vs.breast.lumb.dn.gene,de.res.liver.vs.breast.dn.gene)
de.res.metastasis.liver.vs.breast.lumb          <- res


purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$purity      <- c(purity.MET500,purity.TCGA)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
wilcox.test.p.value.vec <- foreach(g=c(lumb.up.gene,lumb.dn.gene),.combine='c') %do% {
    df$expr           <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
    df$expr           <- df$expr + 0.0001
    loess.fit         <- loess(data = df,formula=expr~purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
    adjusted.expr     <- loess.fit$residuals
    wilcox.test(adjusted.expr[f1],adjusted.expr[f2])$p.value
}
names(wilcox.test.p.value.vec) <- c(lumb.up.gene,lumb.dn.gene)
lumb.wilcox.test.p.value.vec   <- wilcox.test.p.value.vec

# deseq.test.p.value.vec         <- de.res.metastasis.liver.vs.breast.lumb[c(lumb.up.gene,lumb.dn.gene),'pvalue']
# names(deseq.test.p.value.vec)  <- c(lumb.up.gene,lumb.dn.gene)


lumb.up.gene.p.value <- wilcox.test.p.value.vec[lumb.up.gene]
lumb.dn.gene.p.value <- wilcox.test.p.value.vec[lumb.dn.gene]
lumb.up.gene         <- names(lumb.up.gene.p.value)[lumb.up.gene.p.value < 0.05]
lumb.dn.gene         <- names(lumb.dn.gene.p.value)[lumb.dn.gene.p.value < 0.05]

m.value              <- apply(MET500.log2.fpkm.matrix[lumb.dn.gene,MET500.sample],1,median)
lumb.dn.gene         <- lumb.dn.gene[m.value > fpkm.cut.off]
m.value              <- apply(MET500.log2.fpkm.matrix[lumb.up.gene,MET500.sample],1,median)
lumb.up.gene         <- lumb.up.gene[m.value > fpkm.cut.off]

write.csv(x=lumb.up.gene,file = 'client-side/output/DE.breast.cancer.R.output/lumb.up.csv',quote=FALSE)
write.csv(x=lumb.dn.gene,file = 'client-side/output/DE.breast.cancer.R.output/lumb.dn.csv',quote=FALSE)



#Lumb: cxcl1 DE  Ref:Tumor cell-intrinsic factors underlie heterogeneity of immune cell infiltration and response to immunotherapy
#In vivo screening identifies GATAD2B as a metastasis driver in KRAS-driven lung cancer




####################

save(file='client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData',list=c('de.res.liver.vs.breast','de.res.metastasis.liver.vs.breast.basal','de.res.metastasis.liver.vs.breast.lumb','de.res.metastasis.liver.vs.breast.luma','de.res.metastasis.liver.vs.breast.her2'))


##################### Trash code ####################################

# ####################################################################################################################################
# ############# compute correlation values with cell line, as estimate of tumor purity
# load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
# CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
# CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
# tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
# tmp.rank                    <- apply(tmp,2,rank)
# rank.mean                   <- apply(tmp.rank,1,mean)
# rank.sd                     <- apply(tmp.rank,1,sd)
# plot(x=rank.mean,y=rank.sd)
# lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
# CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
# 
# 
# 
# 
# 
# res <- de.res.metastasis.liver.vs.breast.her2
# de.res.metastasis.liver.vs.breast.her2.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < cut.off]
# de.res.metastasis.liver.vs.breast.her2.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
# 
# 
# res <- de.res.metastasis.liver.vs.breast.luma
# de.res.metastasis.liver.vs.breast.luma.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < cut.off]
# de.res.metastasis.liver.vs.breast.luma.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
# 
# res <- de.res.metastasis.liver.vs.breast.lumb
# de.res.metastasis.liver.vs.breast.lumb.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < cut.off]
# de.res.metastasis.liver.vs.breast.lumb.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
# 
# 
# 
# 
# require(gplots)
# venn(list(normal.up.gene=de.res.liver.vs.breast.up.gene,
#           de.res.metastasis.liver.vs.breast.lumb.up.gene=de.res.metastasis.liver.vs.breast.lumb.up.gene,
#           de.res.metastasis.liver.vs.breast.luma.up.gene=de.res.metastasis.liver.vs.breast.luma.up.gene,
#           de.res.metastasis.liver.vs.breast.basal.up.gene=de.res.metastasis.liver.vs.breast.basal.up.gene,
#           de.res.metastasis.liver.vs.breast.her2.up.gene=de.res.metastasis.liver.vs.breast.her2.up.gene
#           )
#     )
# 
# 
# 
# 
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.luma.up.gene,de.res.liver.vs.breast.up.gene),file = 'client-side/output/DE.breast.cancer.R.output/luma.up.csv',quote=FALSE)
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.lumb.up.gene,de.res.liver.vs.breast.up.gene),file = 'client-side/output/DE.breast.cancer.R.output/lumb.up.csv',quote=FALSE)
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.her2.up.gene,de.res.liver.vs.breast.up.gene),file = 'client-side/output/DE.breast.cancer.R.output/her2.up.csv',quote=FALSE)
# 
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.luma.dn.gene,de.res.liver.vs.breast.dn.gene),file = 'client-side/output/DE.breast.cancer.R.output/luma.dn.csv',quote=FALSE)
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.lumb.dn.gene,de.res.liver.vs.breast.dn.gene),file = 'client-side/output/DE.breast.cancer.R.output/lumb.dn.csv',quote=FALSE)
# write.csv(x=setdiff(de.res.metastasis.liver.vs.breast.her2.dn.gene,de.res.liver.vs.breast.dn.gene),file = 'client-side/output/DE.breast.cancer.R.output/her2.dn.csv',quote=FALSE)
# 
# 
# 
# 
# 
# 
# 
# s <- intersect(de.res.metastasis.liver.vs.breast.lumb.up.gene,de.res.metastasis.liver.vs.breast.luma.up.gene)
# s <- intersect(s,de.res.metastasis.liver.vs.breast.her2.up.gene) # in this step, pay attentation to ADAT1
# s <- intersect(s,de.res.metastasis.liver.vs.breast.basal.up.gene)
# 
# s <- setdiff(s,de.res.liver.vs.breast.up.gene) 
# View(s) # pay attentation to ADAT3 and HES4
# #ref: The T Box Transcription Factor TBX2 Promotes Epithelial-Mesenchymal Transition and Invasion of Normal and Malignant Breast Epithelial Cells
# #ref: ID1 promotes breast cancer metastasis by S100A9 regulation
# #ref: Smad6 determines BMP-regulated invasive behaviour of breast cancer cells in a zebrafish xenograft model
# #ref: KLF15 promotes the proliferation and metastasis of lung adenocarcinoma cells and has potential as a cancer prognostic marker
# #ref: Interplay between Notch1 and Notch3 promotes EMT and tumor initiation in squamous cell carcinoma
# #ref: NOTCH3 expression is linked to breast cancer seeding and distant metastasis (for Basal)
# #ref: Overexpression of NOTCH-regulated Ankyrin Repeat Protein is associated with papillary thyroid carcinoma progression
# #ref: Metastatic effect of LY-6K gene in breast cancer cells (LY6K,lumb subtype)
# #ref: Hes4: A potential prognostic biomarker for newly diagnosed patients with high-grade osteosarcoma.
# #ref: Overexpression of NOTCH-regulated Ankyrin Repeat Protein is associated with papillary thyroid carcinoma progression
# #ref: MiR-126-3p suppresses tumor metastasis and angiogenesis of hepatocellular carcinoma by targeting LRP6 and PIK3R2
# #ref: Bromodomain 4 activation predicts breast cancer survival
# #ref: BRD4 Regulates Breast Cancer Dissemination through Jagged1/Notch1 Signaling
# #ref: DOT1L cooperates with the c-Myc-p300 complex to epigenetically derepress CDH1 transcription factors in breast cancer progression (basal subtype)
# #ref: Inhibition of histone H3K79 methylation selectively inhibits proliferation, self-renewal and metastatic potential of breast cancer.
# #ref:Cyclin-dependent kinase 11p110 (CDK11p110) is crucial for human breast cancer cell proliferation and growth
# #GSE43837:brain metastasis
# #CXCL3 is a potential target for breast cancer metastasis.
# 
# s <- intersect(de.res.metastasis.liver.vs.breast.lumb.dn.gene,de.res.metastasis.liver.vs.breast.luma.dn.gene)
# s <- intersect(s,de.res.metastasis.liver.vs.breast.her2.dn.gene)
# s <- intersect(s,de.res.metastasis.liver.vs.breast.basal.dn.gene)
# 
# setdiff(s,de.res.liver.vs.breast.dn.gene) %>% View
# 
# 
# 
# intersect(de.res.metastasis.liver.vs.breast.her2.up.gene,de.res.liver.vs.breast.dn.gene) %>% View
# intersect(de.res.metastasis.liver.vs.breast.her2.dn.gene,de.res.liver.vs.breast.up.gene) %>% View
# 
# 
# intersect(de.res.metastasis.liver.vs.breast.basal.up.gene,de.res.liver.vs.breast.dn.gene) %>% View
# intersect(de.res.metastasis.liver.vs.breast.basal.dn.gene,de.res.liver.vs.breast.up.gene) %>% View
# 
# 
# intersect(de.res.metastasis.liver.vs.breast.luma.up.gene,de.res.liver.vs.breast.dn.gene) %>% View
# intersect(de.res.metastasis.liver.vs.breast.luma.dn.gene,de.res.liver.vs.breast.up.gene) %>% View
# 
# 
# 
# 
# 
# require(gplots)
# venn(list())
# # 
# 
# 
# # TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
# # TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
# # 
# # 
# # MET500.breast.cancer.polyA.sample.biopsy.site   <-  MET500.sample.meta[MET500.breast.cancer.polyA.sample,'biopsy.site']
# # MET500.breast.cancer.polyA.liver.sample         <-  MET500.breast.cancer.polyA.sample[MET500.breast.cancer.polyA.sample.biopsy.site == 'LIVER']
# # MET500.breast.cancer.polyA.lymph.sample         <-  MET500.breast.cancer.polyA.sample[MET500.breast.cancer.polyA.sample.biopsy.site == 'LYMPH_NODE']
# # MET500.breast.cancer.polyA.non.Basal.sample     <-  c(MET500.breast.cancer.polyA.Her2.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.Normal.sample)
# # 
# # TCGA.breast.cancer.polyA.non.Basal.sample       <- c(TCGA.breast.cancer.polyA.Her2.sample,TCGA.breast.cancer.polyA.LumA.sample,TCGA.breast.cancer.polyA.LumB.sample,TCGA.breast.cancer.polyA.Normal.sample)
# # 
# # 
# # # Perform DE analysis, LIVER metastasis vs primary cancer, non Basal-like subtypes 
# # MET500.idx <- intersect(MET500.breast.cancer.polyA.liver.sample,MET500.breast.cancer.polyA.non.Basal.sample)
# # TCGA.idx   <- TCGA.breast.cancer.polyA.non.Basal.sample[1:70]
# # 
# # TCGA.breast.cancer.read.count.matrix   <- TCGA.breast.cancer.log2.read.count.matrix[,TCGA.idx] 
# # MET500.breast.cancer.read.count.matrix <- MET500.log2.read.count.matrix[,MET500.idx]  
# # TCGA.breast.cancer.read.count.matrix   <- TCGA.breast.cancer.read.count.matrix[rownames(TCGA.breast.cancer.read.count.matrix) %in% protein.coding.gene.id,]
# # MET500.breast.cancer.read.count.matrix <- MET500.breast.cancer.read.count.matrix[rownames(MET500.breast.cancer.read.count.matrix) %in% protein.coding.gene.id,]
# # 
# # 
# # 
# # 
# # expr.matrix  <- cbind(MET500.breast.cancer.read.count.matrix,TCGA.breast.cancer.read.count.matrix)
# # expr.matrix  <- 2^expr.matrix - 1
# # expr.matrix  <- expr.matrix[rownames(expr.matrix) %in% protein.coding.gene.id,]
# # flag         <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
# # expr.matrix  <- expr.matrix[flag == ncol(expr.matrix),]
# # df           <- data.frame(condition=c(rep(x='MET500',times=length(MET500.idx)), rep(x='TCGA',times=length(TCGA.idx))))
# # df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
# # dds          <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
# #                                        colData = df,
# #                                        design = ~ condition)
# # 
# # dds <- DESeq(dds)
# # res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
# # res <- res[order(res$pvalue),]
# # res <- res[complete.cases(res),]
# # res.liver.metastasis.vs.primary.tumor.res <- res
# # 
# # #liver.metastasis.vs.primary.tumor.non.Basal.up.de.gene   <- rownames(res)[res$padj < 0.0001 & res$log2FoldChange >  1.5]
# # #liver.metastasis.vs.primary.tumor.non.Basal.down.de.gene <- rownames(res)[res$padj < 0.0001 & res$log2FoldChange < -1.5]
# # 
# # c.gene       <- intersect(rownames(res.liver.metastasis.vs.primary.tumor.res),rownames(res.liver.vs.breast))
# # dd           <- data.frame(res.liver.vs.breast = res.liver.vs.breast[c.gene,2], res.liver.metastasis.vs.primary.tumor.res=res.liver.metastasis.vs.primary.tumor.res[c.gene,2])
# # rownames(dd) <- c.gene
# # View(dd[dd$res.liver.vs.breast > 1.0 & dd$res.liver.metastasis.vs.primary.tumor.res < -1.0,])
# # View(dd[dd$res.liver.vs.breast < -1.0 & dd$res.liver.metastasis.vs.primary.tumor.res > 1.0,])
# # 
# # 
# # #############################################################################################################################################
# # 
# # # Perform DE analysis, LIVER metastasis vs primary cancer,  Basal-like subtypes 
# # MET500.idx <- intersect(MET500.breast.cancer.polyA.liver.sample,MET500.breast.cancer.polyA.Basal.sample)
# # TCGA.idx   <- TCGA.breast.cancer.polyA.Basal.sample
# # 
# # TCGA.breast.cancer.read.count.matrix   <- TCGA.breast.cancer.log2.read.count.matrix[,TCGA.idx] 
# # MET500.breast.cancer.read.count.matrix <- MET500.log2.read.count.matrix[,MET500.idx]  
# # TCGA.breast.cancer.read.count.matrix   <- TCGA.breast.cancer.read.count.matrix[rownames(TCGA.breast.cancer.read.count.matrix) %in% protein.coding.gene.id,]
# # MET500.breast.cancer.read.count.matrix <- MET500.breast.cancer.read.count.matrix[rownames(MET500.breast.cancer.read.count.matrix) %in% protein.coding.gene.id,]
# # 
# # 
# # expr.matrix  <- cbind(MET500.breast.cancer.read.count.matrix,TCGA.breast.cancer.read.count.matrix)
# # expr.matrix  <- 2^expr.matrix - 1
# # expr.matrix  <- expr.matrix[rownames(expr.matrix) %in% protein.coding.gene.id,]
# # flag         <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
# # expr.matrix  <- expr.matrix[flag == ncol(expr.matrix),]
# # df           <- data.frame(condition=c(rep(x='MET500',times=length(MET500.idx)), rep(x='TCGA',times=length(TCGA.idx))))
# # df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
# # dds          <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
# #                                        colData = df,
# #                                        design = ~ condition)
# # 
# # dds <- DESeq(dds)
# # res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
# # res <- res[order(res$pvalue),]
# # res <- res[complete.cases(res),]
# # 
# # 
# # res.liver.metastasis.vs.primary.tumor.basal <- res
# # 
# # 
# # c.gene       <- intersect(rownames(res.liver.metastasis.vs.primary.tumor.basal),rownames(res.liver.vs.breast))
# # dd           <- data.frame(res.liver.vs.breast = res.liver.vs.breast[c.gene,2], res.liver.metastasis.vs.primary.tumor.basal=res.liver.metastasis.vs.primary.tumor.basal[c.gene,2])
# # rownames(dd) <- c.gene
# # View(dd[dd$res.liver.vs.breast > 1.5 & dd$res.liver.metastasis.vs.primary.tumor.basal < -1.5,])
# # View(dd[dd$res.liver.vs.breast < -1.5 & dd$res.liver.metastasis.vs.primary.tumor.basal > 1.5,])
# # 
# # plot(x=dd$res.liver.vs.breast,y=dd$res.liver.metastasis.vs.primary.tumor.basal,xlim=c(-10,10),ylim=c(-10,10))
# # 
# # #liver.metastasis.vs.primary.tumor.Basal.up.de.gene   <- rownames(res)[res$padj < 0.0001 & res$log2FoldChange >  1.5]
# # #liver.metastasis.vs.primary.tumor.Basal.down.de.gene <- rownames(res)[res$padj < 0.0001 & res$log2FoldChange < -1.5]
# # 
# # 
# # 
# # #################
# # require(gplots)
# # venn(list(normal.up.de.gene     = liver.vs.breast.up.de.gene,  non.basal.up.de.gene   = liver.metastasis.vs.primary.tumor.non.Basal.up.de.gene,  basal.up.de.gene  =  liver.metastasis.vs.primary.tumor.Basal.up.de.gene))
# # venn(list(normal.down.de.gene   = liver.vs.breast.down.de.gene,non.basal.down.de.gene = liver.metastasis.vs.primary.tumor.non.Basal.down.de.gene,basal.down.de.gene = liver.metastasis.vs.primary.tumor.Basal.down.de.gene))
# # 
# # ########
# # 
# # o.gene <- intersect(liver.metastasis.vs.primary.tumor.Basal.down.de.gene,liver.metastasis.vs.primary.tumor.non.Basal.down.de.gene)
# # o.gene <- setdiff(o.gene,liver.vs.breast.down.de.gene)
# 
# #save(file='client-side/output/DE.breast.cancer.R.output/DE.RData',list=c('liver.vs.breast.up.de.gene','liver.vs.breast.down.de.gene','liver.metastasis.vs.primary.tumor.Basal.down.de.gene','liver.metastasis.vs.primary.tumor.Basal.up.de.gene','liver.metastasis.vs.primary.tumor.non.Basal.up.de.gene','liver.metastasis.vs.primary.tumor.non.Basal.down.de.gene'))
