setwd("C:/Users/Admin/Downloads")
library(readr)
library(DESeq2)

file1 <- read.table("GSE247649_Rv.txt", header = T)
write.csv(file1,"GSE247049_Rv.csv")
file2<- read.table("GSE247649_Ra.txt", header = T)
write.csv(file2,"GSE247649_Ra.csv")

library(dplyr)

csv1 <- read.csv("GSE247049_Rv.csv", header = T)
csv2 <- read.csv("GSE247649_Ra.csv", header = T)

merged_data <- merge(csv1, csv2, by = "ID", all = TRUE)

write.csv(merged_data, "merged_file_GSE247649.csv", row.names = FALSE)


data_rv_ra_count <- read.csv("merged_file_GSE247649.csv", row.names = "ID", header = T)

condition <- c("Rv", "Rv", "Ra", "Ra")

col_data_rv_ra <- data.frame(condition)

rownames(col_data_rv_ra) <- colnames(data_rv_ra_count)


dds_rv_ra <- DESeqDataSetFromMatrix(countData = data_rv_ra_count,
                              colData = col_data_rv_ra,  
                              design = ~ condition)  


dds_final_rv_ra<- DESeq(dds_rv_ra, test="Wald")

res_rv_ra <- results(dds_final_rv_ra, alpha  = 0.05)
summary(res_rv_ra)

significant_genes_rv_ra <- res_rv_ra[which(res_rv_ra$pvalue < 0.05), ]


write.csv(significant_genes_rv_ra, "sig_rv_ra_GSE247649.csv")

###################################################################################
library(ggplot2)
library(ComplexHeatmap)
library(Cairo)
library(org.Hs.eg.db)
sig_rv_ra <- read.csv("sig_rv_ra_GSE247649.csv", header = TRUE, row.names = 1)
new_rv_ra <- as.data.frame(sig_rv_ra)

#new_rv_ra$ID <- mapIds(org.Hs.eg.db, keys = rownames(new_rv_ra), keytype = "ENSEMBL", column = "SYMBOL")

rv_ra_final <- counts(dds_final_rv_ra, normalized = TRUE)[rownames(new_rv_ra), ]

mat.z_rv_ra <- t(apply(rv_ra_final, 1, scale))

# Set column names of scaled matrix

colnames(mat.z_rv_ra) <- colnames(rv_ra_final)

rv_ra_heatmap <- Heatmap(mat.z_rv_ra, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Z-score", 
        row_names_gp = gpar(fontsize = 8),
        row_names_side = "left")






