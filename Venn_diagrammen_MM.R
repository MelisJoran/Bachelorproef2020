### Activeren packages

library(dplyr)
library(GOplot)



### SIGNIFICANTE genen

# Inlezen data

"DESeq_RSL3vsuntreatedsign" <- read.table(file = "DESeq_RSL3vsuntreatedsign.txt", header = TRUE)

"DESeq_RFvsRSL3sign" <- read.table(file = "DESeq_RFvsRSL3sign.txt", header = TRUE)

"all_epifactor" <- read.table("all_epifactor.txt",
                              header = TRUE)


#Venn lijsten

l11 <- select(DESeq_RSL3vsuntreatedsign, Gene, log2FoldChange)
names(l11) <- c("ID", "logFC")


l12 <- select(DESeq_RFvsRSL3sign, Gene, log2FoldChange)
names(l12) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")

GOVenn(data1 = l11, data2 = l12, data3 = l0, 
       title = "A. Significante genen per dataset (log2FC > 1 & padj < 0.05)",
       label=c("RSL3 vs untreated", "ferrostatine-1 vs RSL3", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c("#2980B9", "#F39C12", "#E74C3C"))


### Milde sign genen

# Inlezen data

"DESeq_RSL3vsuntreatedsignmild" <- read.table(file = "DESeq_RSL3vsuntreatedsignmild.txt", header = TRUE)

"DESeq_RFvsRSL3signmild" <- read.table(file = "DESeq_RFvsRSL3signmild.txt", header = TRUE)

"all_epifactor" <- read.table("all_epifactor.txt",
                              header = TRUE)



#Venn lijsten

l21 <- select(DESeq_RSL3vsuntreatedsignmild, Gene, log2FoldChange)
names(l21) <- c("ID", "logFC")


l22 <- select(DESeq_RFvsRSL3signmild, Gene, log2FoldChange)
names(l22) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")

GOVenn(data1 = l21, data2 = l22, data3 = l0,
       title = "B. Significante genen per dataset (log2FC > 0.5 & padj < 0.05)",
       label=c("RSL3 vs untreated", "ferrostatine-1 vs RSL3", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c( "#2980B9", "#F39C12", "#E74C3C"))


