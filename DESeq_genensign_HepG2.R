### Activeren packages

library(dplyr)
library(ggforce)

### Inlezen data

DESeq_opendata_erastinevsuntr <- read.table("DESeq_opendata_erastinevsuntr.txt",
                                            header = TRUE)

DESeq_opendata_ferrostvserastine <- read.table("DESeq_opendata_ferrostvserastine.txt",
                                               header = TRUE)

all_epifactor <- read.table("all_epifactor.txt", header = TRUE)



# ANALYSE: erastine vs untreated

### Volcano plot 

plot(DESeq_opendata_erastinevsuntr$log2FoldChange, -log10(DESeq_opendata_erastinevsuntr$padj),
     main= "erastine vs untreated", xlab= "LogFoldChange", ylab= "-log10(padj)",
     pch="+",
     xlim = c(-10,10), ylim = c(0,300),
     col= ifelse(abs(DESeq_opendata_erastinevsuntr$log2FoldChange) > 1 & DESeq_opendata_erastinevsuntr$padj < 0.05, "darkred", "grey"))

legend("topright", legend = c("logFC >1 & padj <0.05", "not significant"), col= c("darkred", "grey"),pch= "+")


### Lijst met significanten (L2FC>1 & padj<0.05)

"DESeq_opendata_erastinevsuntrsign" <- filter(DESeq_opendata_erastinevsuntr, abs(log2FoldChange) > 1 & padj < 0.05)

write.table(DESeq_opendata_erastinevsuntrsign, file= "DESeq_opendata_erastinevsuntrsign.txt", row.names= FALSE)



# ANALYSE: ferrostatine vs erastine

### Volcano plot  

plot(DESeq_opendata_ferrostvserastine$log2FoldChange, -log10(DESeq_opendata_ferrostvserastine$padj),
     main= "ferrostatine vs erastine", xlab= "LogFoldChange", ylab= "-log10(padj)",
     pch="+",
     xlim = c(-6,7), ylim = c(0,250),
     col= ifelse(abs(DESeq_opendata_ferrostvserastine$log2FoldChange) > 1 & DESeq_opendata_ferrostvserastine$padj < 0.05, "darkred", "grey"))

legend("topright", legend = c("logFC >1 & padj <0.05", "not significant"), col= c("darkred", "grey"),pch= "+")


### Lijst met significanten (L2FC>1 & padj<0.05)

"DESeq_opendata_ferrostvserastinesign" <- filter(DESeq_opendata_ferrostvserastine, abs(log2FoldChange) > 1 & padj < 0.05)

write.table(DESeq_opendata_ferrostvserastinesign, file= "DESeq_opendata_ferrostvserastinesign.txt", row.names= FALSE)



### Tabel met overeenkomstige sign genen

DESeq_opendata_genensign <- merge(DESeq_opendata_erastinevsuntrsign, DESeq_opendata_ferrostvserastinesign, by="Gene")

write.table(DESeq_opendata_genensign, "DESeq_opendata_genensign.txt", row.names = FALSE)


### Tabel sign genen en epifactor database

DESeq_opendata_histongenensign <- merge(DESeq_opendata_genensign, all_epifactor, by.x="Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_opendata_histongenensign, "DESeq_opendata_histongenensign.txt", row.names = FALSE)


### Tabel erastine vs untr <-> all_epifactor

DESeq_opendata_erastinevsuntrhistongenensign <-
        merge(DESeq_opendata_erastinevsuntrsign, all_epifactor,
              by.x= "Gene", by.y= "HGNC.approved.symbol")

write.table(DESeq_opendata_erastinevsuntrhistongenensign, "DESeq_opendata_erastinevsuntrhistongenensign.txt", row.names = FALSE)


### Tabel ferrostatine vs erastine  <-> all_epifactor

DESeq_opendata_ferrostvserastinehistongenensign <-
        merge(DESeq_opendata_ferrostvserastinesign, all_epifactor,
              by.x= "Gene", by.y= "HGNC.approved.symbol")

write.table(DESeq_opendata_ferrostvserastinehistongenensign, "DESeq_opendata_ferrostvserastinehistongenensign.txt", row.names = FALSE)




### MILDE ANALYSE:


### Lijst met milde significanten (L2FC>0.5 & padj<0.05)

"DESeq_opendata_erastinevsuntrsignmild" <- filter(DESeq_opendata_erastinevsuntr, abs(log2FoldChange) > 0.5 & padj < 0.05)

write.table(DESeq_opendata_erastinevsuntrsignmild, file= "DESeq_opendata_erastinevsuntrsignmild.txt", row.names= FALSE)



"DESeq_opendata_ferrostvserastinesignmild" <- filter(DESeq_opendata_ferrostvserastine, abs(log2FoldChange) > 0.5 & padj < 0.05)

write.table(DESeq_opendata_ferrostvserastinesignmild, file= "DESeq_opendata_ferrostvserastinesignmild.txt", row.names= FALSE)




### Tabel met overeenkomstige milde sign genen

DESeq_opendata_genensignmild <- merge(DESeq_opendata_erastinevsuntrsignmild, DESeq_opendata_ferrostvserastinesignmild, by="Gene")

write.table(DESeq_opendata_genensignmild, "DESeq_opendata_genensignmild.txt", row.names = FALSE)


### Tabel milde sign genen en epifactor database

DESeq_opendata_histongenensignmild <- merge(DESeq_opendata_genensignmild, all_epifactor, by.x="Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_opendata_histongenensignmild, "DESeq_opendata_histongenensignmild.txt", row.names = FALSE)


### Tabel erastine vs untr <-> all_epifactor

DESeq_opendata_erastinevsuntrhistongenensignmild <-
        merge(DESeq_opendata_erastinevsuntrsignmild, all_epifactor,
              by.x= "Gene", by.y= "HGNC.approved.symbol")

write.table(DESeq_opendata_erastinevsuntrhistongenensignmild, "DESeq_opendata_erastinevsuntrhistongenensignmild.txt", row.names = FALSE)


### Tabel ferrostatine vs erastine  <-> all_epifactor

DESeq_opendata_ferrostvserastinehistongenensignmild <-
        merge(DESeq_opendata_ferrostvserastinesignmild, all_epifactor,
              by.x= "Gene", by.y= "HGNC.approved.symbol")

write.table(DESeq_opendata_ferrostvserastinehistongenensignmild, "DESeq_opendata_ferrostvserastinehistongenensignmild.txt", row.names = FALSE)





