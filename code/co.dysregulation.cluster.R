require(clusterProfiler)
library(org.Hs.eg.db)
require(dplyr)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output//DE.prostate.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')





DE.rs             <- PRAD.DE.rs
up.gene.symbol    <- mapping.df[DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,'symbol']
up.gene.symbol    <- up.gene.symbol[is.na(up.gene.symbol) == FALSE]
idx               <- match(up.gene.symbol,hg19.gene.info$genename)
up.gene.coor.df   <- hg19.gene.info[idx,]     
up.gene.coor.df   <- up.gene.coor.df[complete.cases(up.gene.coor.df),]

bed.region        <- paste("chr",up.gene.coor.df$chrom,':',up.gene.coor.df$start,"-",up.gene.coor.df$end,sep='')
names(bed.region) <- up.gene.coor.df$genename

bed.sorted.region <- bedr.sort.region(bed.region)
merged.region     <- bedr.merge.region(x=bed.sorted.region,distance=50000,check.chr = TRUE,check.sort = TRUE)

co.up.clusters    <- setdiff(merged.region,bed.sorted.region)
region.size       <- bedr:::size.region(co.up.clusters)
df                <- data.frame(co.up.clusters=co.up.clusters,region.size=region.size)
df$co.up.clusters <- as.character(df$co.up.clusters)

df$gene.name <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
    f <- bed.sorted.region %in.region% x
    s <- bed.sorted.region[f]
    paste(names(bed.region)[match(s,bed.region)],collapse =":")
}

df$gene.number <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
  f <- bed.sorted.region %in.region% x
  sum(f)
}

df <- df[order(df$gene.number,decreasing = TRUE),]
View(df)

