library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)
allowWGCNAThreads()

# where are kalliso files?
fp<-"~/Documents/FROG/samples/"
# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files
samples<-dir(fp) # where the directory 'samples'
# contains the kallisto output directories - 1 per sample.
files <- file.path(samples, "abundance.h5")
setwd(fp)
names(files) <- paste0(samples)
all(file.exists(files)) # need to navigate into samples directory for this to work!!
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.csv("~/Documents/FROG/metadata.csv")
setwd("~/Documents/FROG/")
save(txi.kallisto.tsv, colData, file="txi_and_colData.RData")
colData$Timepoint<-as.factor(colData$Timepoint)
# ==================================================================================


# if already have a txi object, load it with the metadata (colData)
load("txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same

# ==================================================================================
# read count stats chunk starts here

expressed_genes<-txi.kallisto.tsv$abundance
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),]
expressed_genes<-expressed_genes[,c(19, 1:18)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:19)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Sample")
sub<-all_wheat_genes[,c(8, 2, 3, 4)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>=0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>=0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>=0.5, 1,0)
rep_wise$Sum<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3

expressed<-rep_wise[(rep_wise$Sum >=2),]

for(i in unique(expressed$Factor)){
  data<-expressed[(expressed$Factor==i),]
  factor<-paste(i)
  write.csv(data, file=paste("~/Documents/FROG/", factor, ".csv", sep=""))
  assign(factor, data)
}
write.csv(expressed_genes, file="~/Documents/FROG/all_gene_counts.csv")
write.csv(tpm, file="~/Documents/FROG/all_gene_tpm.csv", row.names=T)


# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Genotype"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Genotype, shape=Treatment), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black"))

vst_counts<-assay(vsd)

