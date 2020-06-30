require(DESeq2)
require(dplyr)
library (VennDiagram)
require(gplots)
require(foreach)
source('client-side/code/util.R')

protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


################################################################################################################
#Load reference tissue data
################################################################################################################
load('server-side//RData//Liver.RData')
female.sample                     <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
Ref.liver.log2.read.count.matrix  <- log2.read.count.matrix[,female.sample]
Ref.liver.log2.fpkm.matrix        <- log2.fpkm.matrix[,female.sample]
liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)





#######################################################################################################################################################################################
###### Prepare the data 
#######################################################################################################################################################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')

TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 

#######################################################################################################################################################################################
###### Function to perform DE analysis between tumor samples
#######################################################################################################################################################################################
perform.DE.analysis.between.metastatic.and.primary.cancer <- function(){
  g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
  g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
  expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
  
  expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
  expr.matrix    <- 2^expr.matrix - 1
  
  #purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
  #purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
  
  df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
  df$purity      <- c(purity.MET500,purity.TCGA)
  df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
  
  dds           <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                          colData = df,
                                          design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
  res <- res[order(res$pvalue),]
  res <- res[complete.cases(res),]     
  res
}



fpkm.cut.off <- log2(10+1)
cut.off      <- 0.05

#######################################################################################################################################################################################
###### DE analysis, LumB subtype
#######################################################################################################################################################################################
MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample

#TM7SF3        <- 'ENSG00000205089'
CSF3R          <- 'ENSG00000119535'
CSF3R.expr     <- MET500.log2.fpkm.matrix[CSF3R,MET500.sample] %>% sort(decreasing = TRUE)
max.sample     <- names(CSF3R.expr)[1]
MET500.sample  <- setdiff(MET500.sample,max.sample)
deseq2.res     <- perform.DE.analysis.between.metastatic.and.primary.cancer()
deseq2.res[CSF3R,]

save(file='client-side/output/re.perform.DE.for.CSF3R.R.output/re.perform.DE.for.CSF3R.RData',list=c('deseq2.res'))




############## run edgeR ###########
require(edgeR)
MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample


g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)

expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
expr.matrix    <- 2^expr.matrix - 1

group <- c(rep(x = 'MET500',length(MET500.sample)), rep(x = 'TCGA',length(TCGA.sample)))
group <- factor(group,levels = c('TCGA','MET500'))

y <- DGEList(counts=expr.matrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
rs.table <- topTags(lrt,p.value = 1,n= nrow(expr.matrix))$table
rs.table[CSF3R,]



########## exclude the outlier then run edgeR ######
MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample
CSF3R.expr     <- MET500.log2.fpkm.matrix[CSF3R,MET500.sample] %>% sort(decreasing = TRUE)
max.sample     <- names(CSF3R.expr)[1]
MET500.sample  <- setdiff(MET500.sample,max.sample)

g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)

expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
expr.matrix    <- 2^expr.matrix - 1

group <- c(rep(x = 'MET500',length(MET500.sample)), rep(x = 'TCGA',length(TCGA.sample)))
group <- factor(group,levels = c('TCGA','MET500'))

y <- DGEList(counts=expr.matrix,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
rs.table <- topTags(lrt,p.value = 1,n= nrow(expr.matrix))$table
rs.table[CSF3R,]
