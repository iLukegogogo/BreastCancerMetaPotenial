
########### analysis of GDSC IC50##############

c.cell.line      <- c('HCC1143','HDQ-P1','MDA-MB-436','HCC1395','HCC1937')
c.cell.line      <- intersect(c.cell.line,colnames(GDSC.IC50.matrix))

pan.cell.line.score        <- apply(GDSC.IC50.matrix,1,function(x) x[is.na(x) == FALSE] %>% median)
amp.cell.line.score        <- apply(GDSC.IC50.matrix[,c.cell.line],1,function(x) x[is.na(x) == FALSE] %>% median)

plot(x=pan.cell.line.score,y=amp.cell.line.score,xlim=c(-8,8),ylim=c(-8,8))
lines(c(-8,8),c(-8,8))
lm.mod <- lm(amp.cell.line.score ~ pan.cell.line.score)
r      <- lm.mod$residuals
GDSC.IC50.result.df <- data.frame(pan.cell.line.score=pan.cell.line.score[names(r)],amp.cell.line.score=amp.cell.line.score[names(r)],
                             residual = r
)






clean.name <- function(x) {
    x <- str_remove_all(x,'_')  
    toupper(x)
}
colnames(GDSC.AUC.matrix) <- sapply(colnames(GDSC.AUC.matrix),clean.name)

########### analysis of GDSC AUC##############
NOTCH3.amp.cell.line <- str_remove_all(NOTCH3.amp.cell.line,'_BREAST')
NOTCH3.oth.cell.line <- str_remove_all(NOTCH3.oth.cell.line,'_BREAST')

c.cell.line           <- intersect(NOTCH3.amp.cell.line,colnames(GDSC.AUC.matrix))
amp.cell.line.score   <- apply(GDSC.AUC.matrix[,c.cell.line],1,function(x) x[is.na(x) == FALSE] %>% median)

c.cell.line           <- intersect(NOTCH3.oth.cell.line,colnames(GDSC.AUC.matrix))
pan.cell.line.score        <- apply(GDSC.AUC.matrix[,c.cell.line],1,function(x) x[is.na(x) == FALSE] %>% median)

plot(x=pan.cell.line.score,y=amp.cell.line.score,xlim=c(-8,8),ylim=c(-8,8))
lines(c(-8,8),c(-8,8))
lm.mod <- lm(amp.cell.line.score ~ pan.cell.line.score)
r      <- lm.mod$residuals
GDSC.AUC.result.df <- data.frame(pan.cell.line.score=pan.cell.line.score[names(r)],amp.cell.line.score=amp.cell.line.score[names(r)],
                                  residual = r
)



########### analysis of CTRP AUC##############

c.cell.line      <- c('HCC1143','HDQP1')
c.cell.line      <- intersect(c.cell.line,colnames(CTRP.AUC.matrix))

pan.cell.line.score        <- apply(CTRP.AUC.matrix,1,function(x) x[is.na(x) == FALSE] %>% median)
amp.cell.line.score        <- apply(CTRP.AUC.matrix[,c.cell.line],1,function(x) x[is.na(x) == FALSE] %>% median)

plot(x=pan.cell.line.score,y=amp.cell.line.score,xlim=c(0,18),ylim=c(0,18))
lines(c(-8,28),c(-8,28))
lm.mod <- lm(amp.cell.line.score ~ pan.cell.line.score)
r      <- lm.mod$residuals
CTRP.AUC.result.df <- data.frame(pan.cell.line.score=pan.cell.line.score[names(r)],amp.cell.line.score=amp.cell.line.score[names(r)],
                                 residual = r
)





########

screen.matrix           <- primary.drug.screen.matrix
breast.id               <- sample.info$DepMap_ID[sample.info$lineage == 'breast' & sample.info$lineage_subtype == 'TNBC'] %>% as.character()
breast.id               <- intersect(breast.id,rownames(screen.matrix))
screen.matrix           <- screen.matrix[breast.id,]
rownames(screen.matrix) <- sample.info[breast.id,'CCLE.Name']
breast.cancer.primary.drug.screen.matrix    <- screen.matrix


screen.matrix           <- secondary.drug.screen.matrix
breast.id               <- sample.info$DepMap_ID[sample.info$lineage == 'breast' & sample.info$lineage_subtype == 'TNBC'] %>% as.character()
breast.id               <- intersect(breast.id,rownames(screen.matrix))
screen.matrix           <- screen.matrix[breast.id,]
rownames(screen.matrix) <- sample.info[breast.id,'CCLE.Name']
breast.cancer.secondary.drug.screen.matrix    <- screen.matrix




require(foreach)
compute.correlation <- function(expr.vec, screen.matrix) {
  cell.line     <- intersect(names(expr.vec),rownames(screen.matrix)) 
  cor.vec <- foreach(j = 1:ncol(screen.matrix),.combine='c') %do% {
    tmp <-     screen.matrix[cell.line,j]
    if( sum(is.na(tmp) == FALSE) >= 5){
      cor(expr.vec[cell.line],screen.matrix[cell.line,j],method='spearman',use='pairwise.complete.obs')
    }else{
      NA      
    }
  }
  names(cor.vec) <- colnames(screen.matrix)
  cor.vec
}

compute.z.score <- function(x) {
  m  <- median(x)
  sd <- 1.4826 * mad(x)
  (x - m) /sd
}

cor.vec         <- compute.correlation(NOTCH3.HES4.score,breast.cancer.crispr.screen.matrix)
cor.vec         <- sort(cor.vec)
z.score         <- compute.z.score(cor.vec)
gene            <- names(z.score)[z.score < -1.5]
crispr.screen.notch3.hes4.cor.vec <- cor.vec
#targetable.gene <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/targetable_GeneNames_chembl25.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
#intersect(gene,targetable.gene)



cor.vec.primary       <- compute.correlation(NOTCH3.HES4.score,breast.cancer.primary.drug.screen.matrix)
cor.vec.secondary     <- compute.correlation(NOTCH3.HES4.score,breast.cancer.secondary.drug.screen.matrix)
cor.vec               <- c(cor.vec.primary,cor.vec.secondary)
cor.vec               <- cor.vec[is.na(cor.vec) == FALSE]
drug.screen.notch3.hes4.cor.vec <- cor.vec
z.score               <- compute.z.score(cor.vec)
df1                   <- data.frame(pert.id = primary.drug.screen.meta %>% rownames,   name=primary.drug.screen.meta$name,stringsAsFactors = FALSE)
df2                   <- data.frame(pert.id = secondary.drug.screen.meta %>% rownames, name=secondary.drug.screen.meta$name,stringsAsFactors = FALSE)
screen.meta           <- rbind(df1,df2)
screen.meta           <- screen.meta[complete.cases(screen.meta),]

df <- foreach(drug.name = screen.meta$name %>% unique,.combine='rbind') %do% {
  pert.id <- screen.meta$pert.id[screen.meta$name == drug.name]
  data.frame(drug.name=drug.name,mean=median(cor.vec[pert.id]), sd = mad(cor.vec[pert.id]),cnt=pert.id %>% length)
  #data.frame(drug.name=drug.name,mean=median(z.score[pert.id]), sd = mad(z.score[pert.id]),cnt=pert.id %>% length)
}

cp.with.replicates.df <- df[df$cnt >=2,]



save(file = 'client-side/output/pharmacological.interaction.R.output/pharmacological.interaction.RData',
     list=c('NOTCH3.HES4.score','breast.cancer.crispr.screen.matrix','breast.cancer.primary.drug.screen.matrix','breast.cancer.secondary.drug.screen.matrix',
            'crispr.screen.notch3.hes4.cor.vec','drug.screen.notch3.hes4.cor.vec','cp.with.replicates.df','primary.drug.screen.meta','secondary.drug.screen.meta')
)

#cp.without.replicates.df <- df[df$cnt == 1,]





# cor.vec  <- sort(cor.vec)
# z.score  <- compute.z.score(cor.vec)
# primary.pert.id <- names(z.score)[z.score < -1.5]
# primary.cp.name <- primary.drug.screen.meta[primary.pert.id,'name']
# 
# 
# cor.vec  <- compute.correlation(NOTCH3.HES4.score,breast.cancer.secondary.drug.screen.matrix)
# cor.vec  <- sort(cor.vec)
# z.score  <- compute.z.score(cor.vec)
# secondary.pert.id <- names(z.score)[z.score < -1.5]
# secondary.cp.name <- secondary.drug.screen.meta[secondary.pert.id,'name']
# 
# 
# # NOTCH3.crispr.association.vec <- compute.correlation(NOTCH3.expr,breast.cancer.crispr.screen.matrix)
# # HES4.crispr.association.vec   <- compute.correlation(HES4.expr,breast.cancer.crispr.screen.matrix)
# # 
# # plot(NOTCH3.crispr.association.vec,HES4.crispr.association.vec)
# # 
# # gene1 <- names(NOTCH3.crispr.association.vec)[NOTCH3.crispr.association.vec < -0.5]
# # gene2 <- names(HES4.crispr.association.vec)[HES4.crispr.association.vec     < -0.5]
# # gene <- intersect(gene1,gene2)




# cor.matrix            <- cor(screen.matrix[breast.cell.line,],BRACA.log2.fpkm.matrix[HES4,breast.cell.line],method='spearman',use = 'complete.obs')
# cor.vec.HES4          <- c(cor.matrix)
# names(cor.vec.HES4)   <- rownames(cor.matrix)
# cor.vec.HES4          <- sort(cor.vec.HES4)
# 
# cor.matrix            <- cor(screen.matrix[breast.cell.line,],BRACA.log2.fpkm.matrix[NOTCH3,breast.cell.line],method='spearman',use = 'complete.obs')
# cor.vec.NOTCH3        <- c(cor.matrix)
# names(cor.vec.NOTCH3) <- rownames(cor.matrix)
# cor.vec.NOTCH3        <- sort(cor.vec.NOTCH3)
# 
# gene1 <- names(cor.vec.HES4)[cor.vec.HES4     < -0.4]
# gene2 <- names(cor.vec.NOTCH3)[cor.vec.NOTCH3 < -0.4]
# gene <- intersect(gene1,gene2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# RNAi.dep.score <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/DepMap/D2_combined_gene_dep_scores.csv", stringsAsFactors=FALSE)
# breast.cell.line <- colnames(BRACA.log2.fpkm.matrix)
# flag             <- colnames(RNAi.dep.score) %in% breast.cell.line
# screen.matrix    <- as.matrix(RNAi.dep.score[,flag])
# rownames(screen.matrix) <- (RNAi.dep.score$X)
# screen.matrix           <- screen.matrix %>% t
# breast.cell.line        <- intersect(breast.cell.line,rownames(screen.matrix))
# 
# cor.matrix            <- cor(screen.matrix[breast.cell.line,],BRACA.log2.fpkm.matrix[HES4,breast.cell.line],method='spearman',use = 'complete.obs')
# cor.vec.HES4          <- c(cor.matrix)
# names(cor.vec.HES4)   <- rownames(cor.matrix)
# cor.vec.HES4          <- sort(cor.vec.HES4)
# 
# cor.matrix            <- cor(screen.matrix[breast.cell.line,],BRACA.log2.fpkm.matrix[NOTCH3,breast.cell.line],method='spearman',use = 'complete.obs')
# cor.vec.NOTCH3        <- c(cor.matrix)
# names(cor.vec.NOTCH3) <- rownames(cor.matrix)
# cor.vec.NOTCH3        <- sort(cor.vec.NOTCH3)
# 
# gene1 <- names(cor.vec.HES4)[cor.vec.HES4     < -0.4]
# gene2 <- names(cor.vec.NOTCH3)[cor.vec.NOTCH3 < -0.4]
# gene <- intersect(gene1,gene2)




############################ Trash code ############################ 
# breast.id <- sample.info$DepMap_ID[sample.info$lineage == 'breast']
# breast.id <- intersect(breast.id,rownames(Achilles.gene.effect.matrix))
# 
# dep.score <- apply(Achilles.gene.effect.matrix[breast.id,],2,function(x) median(x[is.na(x) == FALSE]))
# mad.score <- apply(Achilles.gene.effect.matrix[breast.id,],2,function(x) mad(x[is.na(x) == FALSE]))
# 
# library(mixtools)
# x       <- dep.score[is.na(dep.score) == FALSE]
# mixmdl  <- normalmixEM(x)
# plot(mixmdl,which=2,breaks=30)
# lower   <- mixmdl$mu[2] - 3 * mixmdl$sigma[2]
# 


# ######### Check those pan-cell-line essential genes were still expressed ##########
# 
# 
# basal.up.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.annotation.df.csv",stringsAsFactors = FALSE)$SYMBOL
# basal.up.gene.dep.score <- dep.score[basal.up.gene]
# basal.up.gene.dep.score <- basal.up.gene.dep.score[is.na(basal.up.gene.dep.score) == FALSE]
# basal.up.gene.dep.score[basal.up.gene.dep.score < lower]
# 
# 
# 
# basal.dn.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv",stringsAsFactors = FALSE)$SYMBOL
# basal.dn.gene.dep.score <- dep.score[basal.dn.gene]
# basal.dn.gene.dep.score <- basal.dn.gene.dep.score[is.na(basal.dn.gene.dep.score) == FALSE]
# basal.dn.gene.dep.score[basal.dn.gene.dep.score < lower]
# 
# 
# 
# lumb.up.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.annotation.df.csv",stringsAsFactors = FALSE)$SYMBOL
# lumb.up.gene.dep.score <- dep.score[lumb.up.gene]
# lumb.up.gene.dep.score <- lumb.up.gene.dep.score[is.na(lumb.up.gene.dep.score) == FALSE]
# lumb.up.gene.dep.score[lumb.up.gene.dep.score < lower]
# 
# 
# 
# her2.up.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.annotation.df.csv",stringsAsFactors = FALSE)$SYMBOL
# her2.up.gene.dep.score <- dep.score[her2.up.gene]
# her2.up.gene.dep.score <- her2.up.gene.dep.score[is.na(her2.up.gene.dep.score) == FALSE]
# her2.up.gene.dep.score[her2.up.gene.dep.score < lower]


# basal.up.gene.dep.score <- dep.score[basal.up.gene.annotation.df$SYMBOL %>% as.character()]  %>% sort
# lumb.up.gene.dep.score  <- dep.score[lumb.up.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# her2.up.gene.dep.score  <- dep.score[her2.up.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# boxplot(basal.up.gene.dep.score,lumb.up.gene.dep.score,her2.up.gene.dep.score)
# 
# 
# basal.dn.gene.dep.score <- dep.score[basal.dn.gene.annotation.df$SYMBOL %>% as.character()]  %>% sort
# lumb.dn.gene.dep.score  <- dep.score[lumb.dn.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# her2.dn.gene.dep.score  <- dep.score[her2.dn.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# boxplot(basal.dn.gene.dep.score,lumb.dn.gene.dep.score,her2.dn.gene.dep.score)
# 
# boxplot(basal.up.gene.dep.score,basal.dn.gene.dep.score)
# boxplot(her2.up.gene.dep.score,her2.dn.gene.dep.score)
# boxplot(lumb.up.gene.dep.score,lumb.dn.gene.dep.score)
#data.matrix <- drug.screen.matrix
