
##### RNAseq analysis with DESeq2 #####


##### Install packages #####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("S4Vectors")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray")


#### Activate packages #####

 # Het zou kunnen dat je deze packages eerst zelf nog moet installeren in R

library("reshape2")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("genefilter")
library("dplyr")
library("ggrepel")
library("DESeq2")
library("ggforce")

#### Read in data ####

count_data <- read.table(file="all_readcounts.txt", header=T, sep="\t", row.names = 1)

count_data <- count_data[-c(1:4),]
count_data <- count_data[,-(31:32)] 


col_data <- read.table(file="MetaData_RNAseq_ferroptosis.txt", header=T, sep="\t")

col_data$cell_line <- as.factor(col_data$cell_line)
col_data$treatment <- as.factor(col_data$treatment)
col_data$replicate <- as.factor(col_data$replicate)


  # Inspect data

head(col_data)
str(col_data)
head(count_data)
str(count_data)


#### Create DESeq2 Data Matix ####

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~cell_line + treatment)


#### DSEQ2 Analysis #### Untreated vs RSL3

dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds, test="Wald")

res <- results(dds,name=, contrast = c("treatment", "RSL3", "untreated"))
res <- res[order(res$padj),]
mcols(res, use.names=TRUE)

# Export data as "DESeq2_results_group-RSL3vsuntreated.txt"


  # Volcano Plot

plot(res$log2FoldChange, -log10(res$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "RSL3 vs UNTR",
     xlim = c(-6,15),
     ylim = c(0,50),
     col=ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "blue", "black"))

legend(6, 45, legend = c("logFC >1 & padj <0,05", "not significant"), col = c("blue", "black"), pch=16)


#### DSEQ2 Analysis #### ferrostatine-1 vs RSL3

dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds, test="Wald")

res2 <- results(dds,name=, contrast = c("treatment", "ferrostatin_RSL3","RSL3"))
res2 <- res2[order(res2$padj),]
mcols(res2, use.names=TRUE)

# Export data as "DESeq2_results_group-RFvsRSL3.txt"


# Volcano Plot

plot(res2$log2FoldChange, -log10(res2$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "RSL3 vs UNTR",
     xlim = c(-6,15),
     ylim = c(0,50),
     col=ifelse(res2$padj < 0.05 & abs(res2$log2FoldChange) > 1, "blue", "black"))

legend(6, 45, legend = c("logFC >1 & padj <0,05", "not significant"), col = c("blue", "black"), pch=16)

