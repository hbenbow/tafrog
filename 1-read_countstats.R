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
colData<-read.csv("~/Documents/FROG/")
setwd("~/Documents/FROG/")
colData$Timepoint<-as.factor(colData$Timepoint)
# ==================================================================================


# if already have a txi object, load it with the metadata (colData)
load("~/Documents/bmc/Data/txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

# ==================================================================================
# read count stats chunk starts here

expressed_genes<-txi.kallisto.tsv$counts
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),]
expressed_genes<-expressed_genes[,c(49, 1:48)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:49)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Sample")
sub<-all_wheat_genes[,c(9, 2, 3, 4)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>0.5, 1,0)

# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")

table<-merge(as.data.frame(table(expressed$Factor)),colData, by.x="Var1", by.y="Factor") 

ggplot(table, aes(x=as.factor(Timepoint), y=Freq)) +
  geom_col(aes(fill=Treatment), position="dodge") +
  theme_classic() +
  facet_grid(Genotype~.) +
  theme(legend.position = "right") + 
  scale_fill_manual(values=c("grey80","grey40"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Number of expressed genes") + 
  xlab("Timepoint (dpi)")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black")) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(aes(yintercept=0))

test<-aggregate(all_wheat_genes$TPM, by=list(all_wheat_genes$Factor, all_wheat_genes$GeneID, all_wheat_genes$Rep), FUN="sum")
all_wheat_genes<-all_wheat_genes[,c(9, 1,2,3,4,5,6,7,8)]


expressed_genes_filtered<-expressed_genes_long[(expressed_genes_long$TPM>0.5),]
expressed_genes_filtered$Gene<-expressed_genes_filtered$GeneID
expressed_genes_filtered<-separate(expressed_genes_filtered, Gene, c("Gene", "Transcript"), "\\.")
length(unique(expressed_genes_filtered$GeneID))
length(unique(expressed_genes_filtered$Gene))


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype + Treatment:Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Timepoint", "Genotype", "Rep"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Genotype, shape=Treatment), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  scale_colour_manual(values=c("black","grey50"))

vst_counts<-assay(vsd)

