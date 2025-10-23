#Heirarchical Clustering of IGFBP5
#Proteomics
#Helen Pennington
#September 23, 2025

#packages
library(SummarizedExperiment)
library(pheatmap)
#open summarized experiment adat file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
which(rowData(adat)$EntrezGeneSymbol == "IGFBP5") #c("seq.19581.15","seq.2685.21")
mat <- assay(adat)

#IGFBP5
protein_vec <- mat["seq.19581.15", ]

# correlation of all proteins with your protein, decreasing = TRUE
cors <- cor(t(mat), protein_vec, method = "spearman")
n <- length(protein_vec)
tvals <- cors * sqrt((n-2) / (1 - cors^2))
pvals <- 2 * pt(-abs(tvals), df = n-2)
res_df <- data.frame(
  gene = rownames(mat),
  cor = cors,
  pval = pvals
)
sorted_cors <- res_df[order(res_df$pval), ]
sorted_cors$padj <- p.adjust(sorted_cors$pval, method = "bonferroni")
cors_filt <-  subset(sorted_cors, abs(cor) > 0.6 & padj < 0.05)
top_protein_index <- which(rownames(mat) %in% c(cors_filt$gene))
top_proteins <- data.frame(
  seqID = rownames(rowData(adat))[top_protein_index],
  geneName = rowData(adat)$EntrezGeneSymbol[top_protein_index],
  stringsAsFactors = FALSE
)
top_proteins <- merge(top_proteins, cors_filt, by.x = "seqID", by.y = "gene")
top_proteins <- top_proteins[order(top_proteins$pval), ]
#colnames(top_proteins)[3] <- 'correlation'
head(top_proteins)
# subset matrix
sub_mat <- mat[top_proteins$seqID, ]
rownames(sub_mat) <- top_proteins$geneName[ match(rownames(sub_mat), top_proteins$seqID) ]
annotation_col <- data.frame(
  Condition = colData(adat)$CTEStage, 
  SampleName = rownames(colData(adat))
)
rownames(annotation_col) <- annotation_col$SampleName
annotation_col$SampleName <- NULL
# cluster and plot
my_breaks <- seq(-5, 5, length.out = 100)
pheatmap(sub_mat, breaks = my_breaks,scale = "row", clustering_method = "complete", annotation_col = annotation_col)

install.packages("gprofiler2")
library(gprofiler2)

res_gost <- gost(query = top_proteins$geneName,
                 organism = "hsapiens",
                 sources = c("GO:BP","GO:MF","REAC","KEGG","HP","WP")) 
# view table
res_gost$result %>% head(20)
# simple plotting
gostplot(res_gost, capped = TRUE, interactive = FALSE)
########################################################
res_df <- t(apply(mat, 1, function(x) {
  #print(length(x))
  #print(length(protein_vec))
  ct <- cor.test(x, protein_vec, method = "spearman")
  c(cor = ct$estimate, pval = ct$p.value)
}))
cors <- as.data.frame(res_df)
cors$proteins <- rownames(cors)
sorted_cors <- cors[order(abs(cors$cor.rho), decreasing = TRUE), ]
sorted_cors$padj <- p.adjust(sorted_cors$pval, method = "BH")
cors_filt <- sorted_cors[which(abs(sorted_cors$cor.rho) >= 0.70),]