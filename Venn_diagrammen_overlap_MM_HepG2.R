### Activeren packages

library(dplyr)
library(GOplot)

### Inlezen files

all_epifactor <- read.table("all_epifactor.txt", header = TRUE)

DESeq_genensign <- read.table("DESeq_genensign.txt", header = TRUE)

DESeq_opendata_genensign <- read.table("DESeq_opendata_genensign.txt", header = TRUE)


DESeq_genensignmild <- read.table("DESeq_genensignmild.txt", header = TRUE)

DESeq_opendata_genensignmild <- read.table("DESeq_opendata_genensignmild.txt", header = TRUE)


### Venn diagrammen

# Significante genen     (log2FC = gemiddelde van absolute waarden v logFC)

DESeq_genensign$logFC <- c((abs(DESeq_genensign$log2FoldChange.x)+abs(DESeq_genensign$log2FoldChange.y))/2)
l11 <- select(DESeq_genensign, Gene, logFC)
names(l11) <- c("ID", "logFC")

DESeq_opendata_genensign$logFC <- c((abs(DESeq_opendata_genensign$log2FoldChange.x)+abs(DESeq_opendata_genensign$log2FoldChange.y))/2)
l12 <- select(DESeq_opendata_genensign, Gene, logFC)
names(l12) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")


GOVenn(data1 = l11, data2 = l12, data3 = l0, 
       title = "Significante genen per experiment (log2FC > 1 & padj < 0.05)",
       label=c("eigen experiment", "openbare data", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c("#2980B9", "#F39C12", "#E74C3C"))







# Mild significante genen     (logFC = gemiddelde van absolute waarden v logFC)

DESeq_genensignmild$logFC <- c((abs(DESeq_genensignmild$log2FoldChange.x)+abs(DESeq_genensignmild$log2FoldChange.y))/2)
l21 <- select(DESeq_genensignmild, Gene, logFC)
names(l21) <- c("ID", "logFC")

DESeq_opendata_genensignmild$logFC <- c((abs(DESeq_opendata_genensignmild$log2FoldChange.x)+abs(DESeq_opendata_genensignmild$log2FoldChange.y))/2)
l22 <- select(DESeq_opendata_genensignmild, Gene, logFC)
names(l22) <- c("ID", "logFC")

l0 <- select(all_epifactor, HGNC.approved.symbol)
l0$log2FoldChange <- c("0")
names(l0) <- c("ID", "logFC")


GOVenn(data1 = l21, data2 = l22, data3 = l0, 
       title = "Significante genen per experiment (log2FC > 0.5 & padj < 0.05)",
       label=c("eigen experiment", "openbare data", "Epifactors database"),
       circle.col = c("#7D3C98", "#2E86C1", "#28B463"),
       lfc.col = c("#2980B9", "#F39C12", "#E74C3C"))



