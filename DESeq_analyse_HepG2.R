###Activeren packages

library("reshape2")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("genefilter")
library("dplyr")
library("ggrepel")
library("DESeq2")


###Inlezen data

count_data <- read.table(file = "all_readcounts_opendata.tsv", header = TRUE, row.names = 1)
count_data <- count_data[ -(1:4),-10]

col_data <- read.table(file = "MetaData_RNASeq_openbare_data.txt", header = TRUE)
col_data$cell_line <- as.factor(col_data$cell_line)
col_data$treatment <- as.factor(col_data$treatment)
col_data$replicate <- as.factor(col_data$replicate)


#### Create DSeq2 Data Matix ####

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ treatment)

#### DSEQ2 Analysis ####

dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds, test="Wald")


### RESULTATEN: erastin vs untreated

res <- results(dds, contrast = c("treatment", "erastin", "untreated"))
res <- res[order(res$padj),]
mcols(res, use.names=TRUE)

# Export data: naar je eigen 'directory'
write.table(as.data.frame(res), "DESeq_opendata_erastinevsuntr.txt")



### RESULTATEN: ferrostatin-1 vs erastin

res2 <- results(dds, contrast = c("treatment", "ferrostatin", "erastin"))
res2 <- res2[order(res2$padj),]
mcols(res2, use.names=TRUE)

# Export data: naar je eigen 'directory'
write.table(as.data.frame(res2), "DESeq_opendata_ferrostvserastine.txt")

