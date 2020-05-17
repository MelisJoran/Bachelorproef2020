### Activeren packages

library(dplyr)
library(GOplot)

#Inlezen files

GOChord_data_genes_FC <- read.table("Gochord_data_genes_FC.txt", sep = "\t", header = TRUE)

GOChord_data_GO_genes <- read.table("Gochord_data_GO_genes.txt", sep = "\t", header = TRUE)

GOChord_data_processes <- read.table("GOChord_data_processes.txt", sep = "\t")
colnames(GOChord_data_processes) <- c("processes")

### Creëren plot
chord <- chord_dat(GOChord_data_GO_genes, GOChord_data_genes_FC, GOChord_data_processes$processes)


GOChord(chord, title = "Top 40 differentieel geëxpreseerde epigenetische genen na RSL3 behandeling" , gene.order = "logFC", lfc.max = 4, lfc.min = -1, gene.size = 4,
        process.label = 10, ribbon.col = c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD'))



