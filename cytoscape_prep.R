#Cytoscape Prep
#Helen Pennington
#October 28, 2025

#packages
library(data.table)
library(dplyr)
library(tidyr)
library(lsa)
library(Matrix)
library(igraph)
library(leiden)
reticulate::py_install("numpy")
reticulate::py_install("leidenalg")

combined_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/combined_fgsea_mainmodels_sigonly.csv")
cluster_df <- fread("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/clustered_fgsea_using_function_nofaq_withproteins_grouped.csv")

##################################################################################
#Cytoscape Prep Based on hierarchical clustering (cosine)
#nodes: models and pathways
combined_fgsea$pathway <- paste(combined_fgsea$model, combined_fgsea$pathway, sep = "_")
all_genes <- unique(unlist(strsplit(combined_fgsea$leadingEdge, ",")))

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
  combined_fgsea$leadingEdge,
  combined_fgsea$NES,
  MoreArgs = list(all_genes = all_genes)
))

rownames(pathway_gene_mat) <- combined_fgsea$pathway
colnames(pathway_gene_mat) <- all_genes
mat <- as.matrix(pathway_gene_mat)
mode(mat) <- "numeric"

# 1 Compute cosine similarity between pathways
cos_sim <- cosine(t(mat))   # square matrix: n_pathways x n_pathways
rownames(cos_sim) <- rownames(mat)
colnames(cos_sim) <- rownames(mat)

# 2 apply threshold
threshold <- 0.2
cos_sim_thresh <- cos_sim
cos_sim_thresh[cos_sim_thresh <= threshold] <- 0
diag(cos_sim_thresh) <- 0

# 4. Add to node metadata
cluster_df_h <- cluster_df %>% 
  select(pathway, group)
colnames(cluster_df_h) <- c("id","cluster_h")
m <- as.matrix(cos_sim_thresh)
# 3 Get upper triangle only to avoid duplicate source-target and reverse
tri_idx <- which(upper.tri(m) & m > 0, arr.ind = TRUE)
pp_edges <- data.frame(
  source = rownames(m)[tri_idx[,1]],
  target = rownames(m)[tri_idx[,2]],
  weight = m[upper.tri(m) & m > 0],
  edge_type = "pathway_similarity",
  stringsAsFactors = FALSE
)

# 4 model pathway edges
combined_fgsea <- combined_fgsea %>% 
  select(model, pathway, NES)
mp_edges <- combined_fgsea %>%
  rename(source = model, target = pathway) %>%
  mutate(weight = 1, edge_type = "pathway") %>%
  select(source, target, weight, edge_type)

# 5 Combine into one edge lis
edges <- bind_rows(pp_edges, mp_edges) %>% distinct(source, target, .keep_all = TRUE)

# 6 node table
nodes <- unique(c(edges$source, edges$target))
nodes_df <- data.frame(id = nodes, stringsAsFactors = FALSE)

nodes_df <- nodes_df %>%
  mutate(node_type = case_when(
    id %in% unique(mp_edges$source) ~ "model",
    id %in% unique(mp_edges$target) ~ "pathway",
    TRUE ~ "other"
  ))
nodes_df <- nodes_df %>%
  left_join(cluster_df_h, by = "id")

nes_info <- combined_fgsea %>% select(pathway, NES)
nes_info$direction <- ifelse(nes_info$NES < 0, "down",
                               ifelse(nes_info$NES > 0, "up", "neutral"))
nes_info <- nes_info %>%
  rename(id = pathway)
nes_info <- nes_info %>% select(id,direction)

nodes_df <- nodes_df %>%
  left_join(nes_info, by = "id")

write.csv(edges, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/hierarchical_clustering_network_edges_thresh02_825.csv", row.names = FALSE, quote = FALSE)
write.csv(nodes_df, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/hierarchical_clustering_network_nodes_nothresh02_825.csv", row.names = FALSE, quote = FALSE)
