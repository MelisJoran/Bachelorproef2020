### Activeren packages

library(dplyr)
library(GOplot)


### SIGNIFICANTE genen

# Inlezen data

DESeq_opendata_erastinevsuntrsign <- read.table("DESeq_opendata_erastinevsuntrsign.txt", header = TRUE)

DESeq_opendata_ferrostvserastinesign <- read.table("DESeq_opendata_ferrostvserastinesign.txt", header = TRUE)

"all_epifactor" <- read.table("all_epifactor.txt",
                              header = TRUE)


# Venn lijsten

l11 <- select(DESeq_opendata_erastinevsuntrsign, Gene, log2FoldChange)
names(l11) <- c("ID", "logFC")

l12 <- select(DESeq_opendata_ferrostvserastinesign, Gene, log2FoldChange)
names(l12) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")


GOVenn(data1 = l11, data2 = l12, data3 = l0, 
       title = "A. Significante genen per dataset (log2FC > 1 & padj < 0.05)",
       label=c("erastine vs untreated", "ferrostatine vs erastine", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c("#2980B9", "#F39C12", "#E74C3C"))




### MILD SIGNIFICANTE genen

# Inlezen data

DESeq_opendata_erastinevsuntrsignmild <- read.table("DESeq_opendata_erastinevsuntrsignmild.txt", header = TRUE)

DESeq_opendata_ferrostvserastinesignmild <- read.table("DESeq_opendata_ferrostvserastinesignmild.txt", header = TRUE)

"all_epifactor" <- read.table("all_epifactor.txt",
                              header = TRUE)


# Venn lijsten

l21 <- select(DESeq_opendata_erastinevsuntrsignmild, Gene, log2FoldChange)
names(l21) <- c("ID", "logFC")

l22 <- select(DESeq_opendata_ferrostvserastinesignmild, Gene, log2FoldChange)
names(l22) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")


GOVenn(data1 = l21, data2 = l22, data3 = l0, 
       title = "B. Significante genen per dataset (log2FC > 0.5 & padj < 0.05)",
       label=c("erastine vs untreated", "ferrostatine vs erastine", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c("#2980B9", "#F39C12", "#E74C3C"))







