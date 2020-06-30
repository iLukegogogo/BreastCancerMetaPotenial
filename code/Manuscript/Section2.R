require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')

##############################################################################
#---------------- Figure 2 --------------------------------------------------#
##############################################################################

#######################################################
# Fig 2a: motivation to use piecewise linear regression to remove liver speicfic genes (LuminalB subtype)
######################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.lumb),rownames(de.res.metastasis.liver.vs.breast.lumb))
x             <- de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.log2FC.MP.and.normal.pdf',width = 20,height=20)


#######################################################
# Fig 2b: motivation to use piecewise linear regression to remove liver speicfic genes (GSE58708 dataset)
######################################################
load('client-side/output/validation.of.confounding.R.output/validation.of.confounding.RData')

c.gene        <- intersect(rownames(SRP043470.de.res.liver.vs.breast),rownames(SRP043470.de.res.metastasis.liver.vs.breast))
x             <- SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange']
y             <- SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange'],y=SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/SRP043470.log2FC.MP.and.normal.pdf',width = 20,height=20)


#######################################################
# Fig 2c: DEBoost pipeline
######################################################


##############################################################################
#---------------- Figure S2 --------------------------------------------------#
##############################################################################

#######################################################
# Fig S2a: motivation to use piecewise linear regression to remove liver speicfic genes (Her2 subtype)
######################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.her2),rownames(de.res.metastasis.liver.vs.breast.her2))
x             <- de.res.liver.vs.breast.her2[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.her2[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.log2FC.MP.and.normal.pdf',width = 20,height=20)


#######################################################
# Fig S2b: motivation to use piecewise linear regression to remove liver speicfic genes (Basal-like subtype)
######################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.basal),rownames(de.res.metastasis.liver.vs.breast.basal))
x             <- de.res.liver.vs.breast.basal[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.basal[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/basal.log2FC.MP.and.normal.pdf',width = 20,height=20)







####################################################### Trash ####################################################### 
#######################################################
# FigS2a: Bar plot to show subtype-sepcificity of DE genes
######################################################
# up.gene <- c(lumb.up.gene,her2.up.gene,basal.up.gene)
# dn.gene <- c(lumb.dn.gene,her2.dn.gene,basal.dn.gene)
# 
# up.gene.freq.df <- table(up.gene) %>% as.data.frame
# dn.gene.freq.df <- table(dn.gene) %>% as.data.frame
# up.gene.freq.df$color <- 'up'
# dn.gene.freq.df$color <- 'dn'
# colnames(up.gene.freq.df)[1] <- 'gene'
# colnames(dn.gene.freq.df)[1] <- 'gene'
# 
# ggplot(rbind(up.gene.freq.df,dn.gene.freq.df)) + 
#   geom_bar( aes( x=factor(Freq),fill=color),position= 'dodge') + 
#   ggplot.style + 
#   xlab('Number of subtypes') + 
#   scale_fill_manual(values=c('up'='red','dn'='blue'))
# 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/DE.gene.subtype.specificity.pdf',width = 20,height=20)

# subtype.sample          <- intersect(pure.TCGA.breast.cancer.polyA.Her2.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
# her2.eg.up.gene.KM.plot <- draw.KM.plot(her2.eg.up.gene)
# her2.eg.dn.gene.KM.plot <- draw.KM.plot(her2.eg.dn.gene)
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.eg.up.gene.KM.plot.pdf',width=20,height=15)
# print(her2.eg.up.gene.KM.plot[[1]])
# dev.off()
# her2.survival.rs.up[her2.eg.up.gene,]
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.eg.dn.gene.KM.plot.pdf',width=20,height=15)
# print(her2.eg.dn.gene.KM.plot[[1]])
# dev.off()
# her2.survival.rs.dn[her2.eg.dn.gene,]
# 
# 
# 
# subtype.sample          <- intersect(pure.TCGA.breast.cancer.polyA.LumB.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
# lumb.eg.up.gene.KM.plot <- draw.KM.plot(lumb.eg.up.gene)
# lumb.eg.dn.gene.KM.plot <- draw.KM.plot(lumb.eg.dn.gene)
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.eg.up.gene.KM.plot.pdf',width=20,height=15)
# print(lumb.eg.up.gene.KM.plot[[1]])
# dev.off()
# lumb.survival.rs.up[lumb.eg.up.gene,]
# 
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.eg.dn.gene.KM.plot.pdf',width=20,height=15)
# print(lumb.eg.dn.gene.KM.plot[[1]])
# dev.off()
# lumb.survival.rs.dn[lumb.eg.dn.gene,]

# 
# 
# 
# 
# 
# 
# 
# tmp <- basal.survival.rs.up[basal.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# basal.eg.up.gene <- 'ENSG00000103253'
# 
# 
# tmp <- her2.survival.rs.up[her2.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# her2.eg.up.gene <- 'ENSG00000156453'
# 
# 
# tmp <- lumb.survival.rs.up[lumb.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# lumb.eg.up.gene <- 'ENSG00000132613'
# 
# 
# tmp <- basal.survival.rs.dn[basal.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# basal.eg.dn.gene <- 'ENSG00000160307'
# 
# 
# tmp <- her2.survival.rs.dn[her2.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# her2.eg.dn.gene <- 'ENSG00000182326'
# 
# 
# tmp <- lumb.survival.rs.dn[lumb.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# lumb.eg.dn.gene <- 'ENSG00000188001'

