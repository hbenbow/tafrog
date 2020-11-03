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
library(tximport)

setwd("~/Documents/FROG/")
load("txi_and_colData.RData")
colData<-read.csv("~/Documents/FROG/metadata.csv")

# check that order of samples in metadata and txi object are the same
order<-as.data.frame(colnames(txi.kallisto.tsv$abundance))
colnames(order)<-"Sample"
meta<-merge(order, colData, by="Sample", sort=F)
colData<-meta
colData$Rep<-as.factor(colData$Rep)
colData$Timepoint<-as.factor(colData$Timepoint)
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor )
dds<-DESeq(dds)


R1<-as.data.frame(results(dds, contrast=c("Factor","N1O","N1W")))
R2<-as.data.frame(results(dds, contrast=c("Factor","M2O","M2W")))
R3<-as.data.frame(results(dds, contrast=c("Factor","F2O","F2W")))
R4<-as.data.frame(results(dds, contrast=c("Factor","F2O","M2O")))
R5<-as.data.frame(results(dds, contrast=c("Factor","F2W","M2W")))

R1$Factor<-"N1O vs N1W"
R2$Factor<-"M2O vs M2W"
R3$Factor<-"F2O vs F2W"
R4$Factor<-"F2O vs M2O"
R5$Factor<-"F2W vs M2W"

R1$Gene<-row.names(R1)
R2$Gene<-row.names(R2)
R3$Gene<-row.names(R3)
R4$Gene<-row.names(R4)
R5$Gene<-row.names(R5)

results<-rbind(R1, R2, R3, R4, R5)
results<-na.omit(results)
significant<-results[(results$padj < 0.05),]

write.csv(results, file="~/Documents/FROG/all_DE_results.csv")
write.csv(significant, file="~/Documents/FROG/Significant_DEGs.csv")
