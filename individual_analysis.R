# de testing for individual datasets that were filtered before

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

setwd("data/")

# Import kallisto files with txi 
# ==================================================================================
colData<-read.csv("~/Documents/FROG/metadata.csv")
load("~/Documents/FROG/txi.RData")
colData$Timepoint<-as.factor(colData$Timepoint)
filter<-read.csv( file="~/Documents/FROG/all_lists_filtered.csv")
counts<-txi.kallisto.tsv$counts

N1<-c("N1O","N1W")
M2<-c("M2O","M2W")
F2<-c("F2O","F2W")
FO<-c("F2O","M2O")
W2<-c("F2W","M2W")
key<-as.data.frame(rbind(N1,
                         M2,
                         F2,
                         FO,
                         W2))
key$Comparison<-row.names(key)

# treatment
comparisons<-c("N1",
               "M2",
               "F2")
list<-list()
degs_genotype<-list()
for(set in unique(comparisons)){
  comp<-filter[(filter$Comparison==set),]
  keys<-subset(key, key$Comparison %in% comp$Comparison)
  metadata<-colData[(colData$Factor==keys$V1 | colData$Factor == keys$V2),,drop=T]
  metadata<-droplevels(metadata)
  samples<-as.character(metadata$Sample)
  genes<-as.character(comp$GeneID)
  df<-counts[genes,samples]
  case<-colnames(df)==metadata$Sample
  list[[length(list)+1]]<-case
  dds <- DESeqDataSetFromMatrix(round(df), metadata, ~ Genotype)
  dds<-DESeq(dds)
  r<-as.data.frame(results(dds, contrast=c("Genotype", "O", "W") ,format='DataFrame', tidy=TRUE))
  r$comparison<-paste(set)
  degs_genotype[[length(degs_genotype)+1]]<-r
}

comparisons<-c("FO",
               "W2")
list<-list()
degs_treatment<-list()
for(set in unique(comparisons)){
  comp<-filter[(filter$Comparison==set),]
  keys<-subset(key, key$Comparison %in% comp$Comparison)
  metadata<-colData[(colData$Factor==keys$V1 | colData$Factor == keys$V2),,drop=T]
  metadata<-droplevels(metadata)
  samples<-as.character(metadata$Sample)
  genes<-as.character(comp$GeneID)
  df<-counts[genes,samples]
  case<-colnames(df)==metadata$Sample
  list[[length(list)+1]]<-case
  dds <- DESeqDataSetFromMatrix(round(df), metadata, ~ Treatment)
  dds<-DESeq(dds)
  r<-as.data.frame(results(dds, contrast=c("Treatment", "F", "M") ,format='DataFrame', tidy=TRUE))
  r$comparison<-paste(set)
  degs_treatment[[length(degs_treatment)+1]]<-r
}
genotype_filtered<-as.data.frame(do.call(rbind.data.frame, degs_genotype))
treatment_filtered<-as.data.frame(do.call(rbind.data.frame, degs_treatment))

all<-rbind(genotype_filtered, treatment_filtered)
all<-all %>%
  mutate(Comparison=case_when(
    comparison == "N1" ~ "N1O vs N1W",
    comparison == "M2" ~ "M2O vs M2W",
    comparison == "F2" ~ "F2O vs F2W",
    comparison == "FO" ~ "F2O vs M2O",
    comparison == "W2" ~ "F2W vs M2W",
  ))
all$comparison<-NULL
all_total<-all
all<-all[(all$padj<0.05),]
all<-na.omit(all)

table(all$Comparison)
write.csv(all, file="all_separate.csv")
amelie_genes<-subset(all_total, all_total$row %in% tafrog_rnaseq_all_transcripts_id$V1)

