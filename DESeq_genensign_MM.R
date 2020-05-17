### Activeren packages

library(dplyr)
library(ggforce)

### Inlezen data


"DESeq_RSL3vsuntreated" <- read.table("DESeq2_results_group-RSL3vsuntreated.txt",
                                      header = TRUE, sep = "\t")

"DESeq_RFvsRSL3" <- read.table("DESeq2_results_group-RFvsRSL3.txt",
                               header = TRUE, sep = ",")

"all_epifactor" <- read.table("all_epifactor.txt",
                              header = TRUE)



# ANALYSE: RSL3 vs untreated

### Volcano plot  

plot(DESeq_RSL3vsuntreated$log2FoldChange, -log10(DESeq_RSL3vsuntreated$padj),
     main= "RSL3 vs untreated", xlab= "Log2FoldChange", ylab= "-log10(padj)",
     pch="+",
     xlim = c(-7,13), ylim = c(0,45),
     col= ifelse( abs(DESeq_RSL3vsuntreated$log2FoldChange) > 1 & DESeq_RSL3vsuntreated$padj < 0.05, "darkred", "grey"))

legend("topleft", legend = c("logFC >1 & padj <0.05", "not significant"), col= c("darkred", "grey"),pch= "+")



### Lijst met significanten (L2FC>1 & padj<0.05)

"DESeq_RSL3vsuntreatedsign" <- filter(DESeq_RSL3vsuntreated, abs(log2FoldChange) > 1 & padj < 0.05)

write.table(DESeq_RSL3vsuntreatedsign, file= "DESeq_RSL3vsuntreatedsign.txt", row.names= FALSE)




# ANALYSE: RF vs RSL3

### Volcano plot  

plot(DESeq_RFvsRSL3$log2FoldChange, -log10(DESeq_RFvsRSL3$padj),
     main= "RF vs RSL3", xlab= "Log2FoldChange", ylab= "-log10(padj)",
     pch="+",
     xlim = c(-10,7), ylim = c(0,25),
     col= ifelse( abs(DESeq_RFvsRSL3$log2FoldChange) > 1 & DESeq_RFvsRSL3$padj < 0.05, "darkred", "grey"))

legend("topright", legend = c("logFC >1 & padj <0.05", "not significant"), col= c("darkred", "grey"),pch= "+")


### Lijst met significanten (L2FC>1 & padj<0.05)

"DESeq_RFvsRSL3sign" <- filter(DESeq_RFvsRSL3, abs(log2FoldChange) > 1 & padj < 0.05)

write.table(DESeq_RFvsRSL3sign, file= "DESeq_RFvsRSL3sign.txt", row.names= FALSE)




# VERGELIJKEN van datasets

### Tabel met overeenkomstige sign genen

"DESeq_genensign" <- merge(DESeq_RSL3vsuntreatedsign, DESeq_RFvsRSL3sign, by="Gene")

write.table(DESeq_genensign, "DESeq_genensign.txt", row.names = FALSE)
                           
### Tabel sign genen en epifactor database

"DESeq_histongenensign" <- merge(DESeq_genensign, all_epifactor, by.x="Gene", by.y = "HGNC.approved.symbol", all = TRUE)

### Tabel RSL3vsuntreated <-> all_epifactor

"DESeq_RSL3vsuntreatedhistongenensign" <- 
  merge(DESeq_RSL3vsuntreatedsign, all_epifactor, 
        by.x = "Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_RSL3vsuntreatedhistongenensign, "DESeq_RSL3vsuntreatedhistongenensign.txt", row.names = FALSE)

### Tabel RFvsRSL3 <-> all_epifactor

"DESeq_RFvsRSL3histongenensign" <-
  merge(DESeq_RFvsRSL3sign, all_epifactor,
        by.x = "Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_RFvsRSL3histongenensign,"DESeq_RFvsRSL3histongenensign.txt", row.names = FALSE)



### MILDE FILTER: LFC > 0,5 & padj < 0,05

"DESeq_RSL3vsuntreatedsignmild" <- filter(DESeq_RSL3vsuntreated, abs(log2FoldChange) > 0.5 & padj < 0.05)

write.table(DESeq_RSL3vsuntreatedsignmild, file= "DESeq_RSL3vsuntreatedsignmild.txt", row.names= FALSE)


"DESeq_RFvsRSL3signmild" <- filter(DESeq_RFvsRSL3, abs(log2FoldChange) > 0.5 & padj < 0.05)

write.table(DESeq_RFvsRSL3signmild, file= "DESeq_RFvsRSL3signmild.txt", row.names= FALSE)



# VERGELIJKEN van datasets

### Tabel met overeenkomstige sign genen

"DESeq_genensignmild" <- merge(DESeq_RSL3vsuntreatedsignmild, DESeq_RFvsRSL3signmild, by="Gene")

write.table(DESeq_genensignmild, "DESeq_genensignmild.txt", row.names = FALSE)

### Tabel sign genen en epifactor database

"DESeq_histongenensignmild" <- merge(DESeq_genensignmild, all_epifactor, by.x="Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_histongenensignmild, "DESeq_histongenensignmild.txt", row.names = FALSE)

### Tabel RSL3vsuntreated <-> all_epifactor

"DESeq_RSL3vsuntreatedhistongenensignmild" <- 
  merge(DESeq_RSL3vsuntreatedsignmild, all_epifactor, 
        by.x = "Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_RSL3vsuntreatedhistongenensignmild, "DESeq_RSL3vsuntreatedhistongenensignmild.txt", row.names = FALSE)

### Tabel RFvsRSL3 <-> all_epifactor

"DESeq_RFvsRSL3histongenensignmild" <-
  merge(DESeq_RFvsRSL3signmild, all_epifactor,
        by.x = "Gene", by.y = "HGNC.approved.symbol")

write.table(DESeq_RFvsRSL3histongenensignmild,"DESeq_RFvsRSL3histongenensignmild.txt", row.names = FALSE)



### Volcano plot  

plot(DESeq_RSL3vsuntreated$log2FoldChange, -log10(DESeq_RSL3vsuntreated$padj),
     main= "RSL3 vs untreated", xlab= "Log2FoldChange", ylab= "-log10(padj)",
     pch="+",
     xlim = c(-7,13), ylim = c(0,45),
     col= ifelse( abs(DESeq_RSL3vsuntreated$log2FoldChange) > 0.5 & DESeq_RSL3vsuntreated$padj < 0.05, "darkred", "grey"))

legend("topleft", legend = c("logFC >0,5 & padj <0.05", "not significant"), col= c("darkred", "grey"),pch= "+")





