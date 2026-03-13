library(DESeq2)
library(tidyverse)
library(rlist)
library(readr)
library(dplyr)
####read files
##setwd
Counts1 <- as.matrix (read.csv("Gene_CM.csv", row.names = "gene_id"))

Coldata <- read.csv("ColData.csv", row.names = 1)


dds <- DESeqDataSetFromMatrix(countData = Counts1,
                              colData = Coldata,
                              design= ~ Condition)
dds <- DESeq(dds)
#remove genes/rows with less than 10 transcript count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Compare HC with prodromal PD
HCvsProdromalPD<-results(dds, contrast = c("Condition", "Healthy_Control","Prodromal_PD")) 
HCvsProdromalPD <- na.omit(HCvsProdromalPD)



#Filter for , logFold2, pvalue 
 
HCvsProdromalPD <- HCvsProdromalPD[(HCvsProdromalPD$log2FoldChange > 0.58 | HCvsProdromalPD$log2FoldChange < -0.58) & HCvsProdromalPD$padj <= 0.05, ]
############No ADAR gene found differentially expressed


###setwd
write.csv(HCvsProdromalPD, "HCvsProdromalPDGene.csv")
##########extract normalized gene counts
normalized_counts_genes <- counts(dds, normalized=TRUE)
write.csv(normalized_counts_genes, "normalised_counts_gene.csv")

#Transcript expression

#read in data

Counts2 <- as.matrix (read.csv("Transcript_CM.csv", row.names = "transcript_id"))


dds2 <- DESeqDataSetFromMatrix(countData = Counts2,
                              colData = Coldata,
                              design= ~ Condition)
dds2 <- DESeq(dds2)
#remove genes/rows with less than 10 transcript count
keep <- rowSums(counts(dds2,)) >= 10
dds2 <- dds2[keep,]

#Compare HC with Prodromal PD
HCvsProdromalPD_transcript <-results(dds2, contrast = c("Condition", "Healthy_Control","Prodromal_PD")) 
HCvsProdromalPD_transcript <- na.omit(HCvsProdromalPD_transcript)
#Filter for P_value =< 0.05, logFold2 
HCvsProdromalPD_transcript<-HCvsProdromalPD_transcript [(HCvsProdromalPD_transcript$log2FoldChange > 0.58 | HCvsProdromalPD_transcript$log2FoldChange < -0.58) & HCvsProdromalPD_transcript$padj <= 0.05, ]
print(HtvsWt_pvalue)
#save results
write.csv(HCvsProdromalPD_transcript, "HCvsProdromalPD_transcript.csv")

#retrieve normalized counts
normalized_counts_transcripts <- counts(dds2, normalized=TRUE)
write.csv(normalized_counts_transcripts, "normalized_counts_transcripts.csv")
















