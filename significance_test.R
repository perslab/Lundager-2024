library(org.Hs.eg.db)
library(tidyverse)
library(writexl)
library(gridExtra)
library(Seurat)

# Load data
data <- read.csv("/projects/perslab/people/jmg776/ad_hoc/anne/baron2016_human.avg_expr.ensembl.csv")
data$gene <- mapIds(org.Hs.eg.db, keys = data$gene, keytype = "ENSEMBL", column = "SYMBOL")   # Convert ENSEMBL ID to gene symbol

goi <- c("ACSL1", "CENPU", "TENT5C", "G6PC2", "GRB10", "SSTR1", "SSTR2", "UBE2E2")   # Genes of interest
goi %in% data$gene   # Sanity check: Are the genes in the dataframe?

# Statistical Testing
data$labels <- ifelse(data$gene %in% goi, "goi", "background")
goi_results <- setNames(data.frame(matrix(0, ncol = 3, nrow = 14), row.names = colnames(data[2:15])), c("p_value", "W_statistic", "p_value_BH_adj")) # column 2-15 designates the 14 cell types

for (cell_type in colnames(data[2:15])) {    
   mann_whitney <- wilcox.test(data[data$labels == "goi",][[cell_type]], 
                   data[data$labels == "background",][[cell_type]], 
                   alternative = "greater")
   goi_results[cell_type, "p_value"] <- mann_whitney$p.value
   goi_results[cell_type, "W_statistic"] <- mann_whitney$statistic              
   goi_results[cell_type, "p_value_BH_adj"] <- p.adjust(mann_whitney$p.value, method = "BH")

}

write_xlsx(goi_results, path = "/projects/perslab/people/jmg776/ad_hoc/anne/enrichment_results.xlsx")
