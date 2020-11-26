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

setwd("~/Documents/FROG/")
# ==================================================================================
# if already have a txi object, load it with the metadata (colData)
colData<-read.csv("~/Documents/FROG/metadata.csv")
load("~/Documents/FROG/txi.RData")
colData$Timepoint<-as.factor(colData$Timepoint)
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

# N1O vs N1W comparison is called "N1"
# M2O vs M2W comparison is called "M2"
# F2O vs F2W comparison is called "F2"
# F2O vs M2O comparison is called "FO"
# F2W vs M2W comparison is called "W2"

N1<-rbind(N1O,N1W)
M2<-rbind(M2O,M2W)
F2<-rbind(F2O,F2W)
FO<-rbind(F2O,M2O)
W2<-rbind(F2W,M2W)

N1<-N1[!(duplicated(N1$GeneID)),]
M2<-M2[!(duplicated(M2$GeneID)),]
F2<-F2[!(duplicated(F2$GeneID)),]
FO<-FO[!(duplicated(FO$GeneID)),]
W2<-W2[!(duplicated(W2$GeneID)),]

N1$Comparison<-"N1"
M2$Comparison<-"M2"
F2$Comparison<-"F2"
FO$Comparison<-"FO"
W2$Comparison<-"W2"

N1<-N1[,c(2, 10)]
M2<-M2[,c(2, 10)]
F2<-F2[,c(2, 10)]
FO<-FO[,c(2, 10)]
W2<-W2[,c(2, 10)]

all_filtered_lists<-rbind(N1,
                          M2,
                          F2,
                          FO,
                          W2)

write.csv(all_filtered_lists, file="~/Documents/FROG/all_lists_filtered.csv", row.names = F)


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


write.csv(pcaData, file="~/Documents/FROG/pcadata.csv")
vst_counts<-assay(vsd)

