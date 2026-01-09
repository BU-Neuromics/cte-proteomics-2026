#fgsea groupings
#helen pennington
#september 18, 2025

#packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(qusage)
library(readxl)
library(proxy)
library(pheatmap)
library(dynamicTreeCut)
library(ComplexHeatmap)

#combined main models file - fgsea
combined <- fread("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/combined_fgsea_mainmodels.csv")
combined <- combined[which(combined$padj < 0.05),]
combined$pathway <- paste(combined$model, combined$pathway, sep = "_")
#hierarchical clustering 
all_genes <- unique(unlist(strsplit(combined$leadingEdge, ",")))

# Function to assign signed membership
gene_membership_signed <- function(genes_str, nes, all_genes) {
  genes <- strsplit(genes_str, ",")[[1]]
  vec <- numeric(length(all_genes))
  vec[all_genes %in% genes] <- ifelse(nes > 0, 1, -1)
  vec
}

# Build matrix
pathway_gene_mat <- t(mapply(
  gene_membership_signed,
  combined$leadingEdge,
  combined$NES,
  MoreArgs = list(all_genes = all_genes)
))

rownames(pathway_gene_mat) <- combined$pathway
colnames(pathway_gene_mat) <- all_genes

row_dist <- dist(pathway_gene_mat, method = "cosine")
row_hc <- hclust(row_dist, method = "average")
row_order <- row_hc$order
row_dist <- as.matrix(row_dist)
rownames(row_dist) <- row_hc$labels
colnames(row_dist) <- row_hc$labels


col_dist <- dist(t(pathway_gene_mat), method = "cosine")
col_hc <- hclust(col_dist, method = "average")
col_order <- col_hc$order

ordered_mat <- pathway_gene_mat[row_order, col_order]

matrix_df <- as.data.frame(ordered_mat)

matrix_df <- cbind(
  NES = combined$NES[match(rownames(matrix_df), combined$pathway)],
  matrix_df
)

clusters <- cutreeDynamic(
  dendro = row_hc,
  distM = row_dist,
  deepSplit = 2,              # 0–4, controls sensitivity (2 is usually good)
  pamRespectsDendro = FALSE,  # allow small, irregular clusters
  minClusterSize = 5          # can adjust based on your number of pathways
)
names(clusters) <- row_hc$labels

clustered_pathways <- data.frame(
  pathway = names(clusters),
  cluster = clusters
)

matrix_df$pathway <- rownames(matrix_df)
matrix_df$row_number <- 1:nrow(matrix_df)
clusters_withgroups <- merge(clustered_pathways, matrix_df, by.x = "pathway", by.y = "pathway")
write.csv(clusters_withgroups,"/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/clustered_fgsea_using_function_nofaq_withproteins.csv")

#create visualization of groupings
clusters_withgroups_sorted <- clusters_withgroups[order(clusters_withgroups$row_number), ]
clusters_matrix <- as.matrix(clusters_withgroups_sorted[,4:1819])
rownames(clusters_matrix) <- clusters_withgroups_sorted$pathway

heatmap <- Heatmap(
  pathway_gene_mat,
  cluster_rows = row_hc,
  cluster_columns = col_hc,
  use_raster = TRUE,          # <- key for large matrices
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = c("blue", "white", "red"),
  row_title = "Pathways",
  column_title = "Proteins",
  heatmap_legend_param = list(
    title = "NES",
    at = c(1, 0, -1),
    labels = c("Positive", "NA", "Negative"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_direction = "vertical"
  )
)

heatmap
png(file = "my_heatmap_300dpi.png",
    width = 7,         # Width in inches
    height = 7,        # Height in inches
    units = "in",      # Specify units for width and height
    res = 300          # Set resolution to 300 dpi
)
png(file = "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/heatmap_groups_300dpi.png", width = 12, height = 6,units = "in", res = 300)
print(heatmap)
dev.off()
