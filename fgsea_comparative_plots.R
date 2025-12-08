#fgsea comparative scatter plots
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

#fgsea scatterplots
#models
fgsea_rl <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvslow.csv")
rl_sig <- fgsea_rl[which(fgsea_rl$padj < 0.05),]
rl_sig[, model := "RHIvsLowCTE"]
fgsea_rh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvsHigh.csv")
rh_sig <- fgsea_rh[which(fgsea_rh$padj < 0.05),]
rh_sig[, model := "RHIvsHighCTE"]
fgsea_lh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_lowvshigh.csv")
lh_sig <- fgsea_lh[which(fgsea_lh$padj < 0.05),]
lh_sig[, model := "LowvsHighCTE"]
fgsea_totyrs <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_totyrs.csv")
totyrs_sig <- fgsea_totyrs[which(fgsea_totyrs$padj < 0.05),]
totyrs_sig[, model := "totyrs"]
fgsea_AT8 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_AT8_total.csv")
AT8_sig <- fgsea_AT8[which(fgsea_AT8$padj < 0.05),]
AT8_sig[, model := "AT8_total"]
fgsea_cds <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CDStot.csv")
cds_sig <- fgsea_cds[which(fgsea_cds$padj < 0.05),]
cds_sig[, model := "CDStot"]
fgsea_dem <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_DementiaHx.csv")
dem_sig <- fgsea_dem[which(fgsea_dem$padj < 0.05),]
dem_sig[, model := "Dementia"]
fgsea_faq <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_faqtot.csv")
faq_sig <- fgsea_faq[which(fgsea_faq$padj < 0.05),]
faq_sig[, model := "faq"]
fgsea_rhi1 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHI1.csv")
rhi1_sig <- fgsea_rhi1[which(fgsea_rhi1$padj < 0.05),]
rhi1_sig[, model := "RHIvsStage1"]
fgsea_rhi2 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHI2.csv")
rhi2_sig <- fgsea_rhi2[which(fgsea_rhi2$padj < 0.05),]
rhi2_sig[, model := "RHIvsStage2"]
fgsea_rhi3 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHI3.csv")
rhi3_sig <- fgsea_rhi3[which(fgsea_rhi3$padj < 0.05),]
rhi3_sig[, model := "RHIvsStage3"]
fgsea_rhi4 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHI4.csv")
rhi4_sig <- fgsea_rhi4[which(fgsea_rhi4$padj < 0.05),]
rhi4_sig[, model := "RHIvsStage4"]

#main models
combined <- rbindlist(list(rl_sig,rh_sig,AT8_sig,totyrs_sig,cds_sig,dem_sig), use.names = TRUE, fill = TRUE)
#write.csv(combined, "/restricted/projectnb/cteseq/projects/somascan/final_plots/combined_fgsea_mainmodels.csv")
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
#write.csv(matrix_df,"/restricted/projectnb/cteseq/projects/somascan/final_plots/clustered_fgsea_res_withmodel_nofaq.csv")

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
#write.csv(clustered_pathways,"/restricted/projectnb/cteseq/projects/somascan/final_plots/clustered_fgsea_using_function_nofaq.csv")

matrix_df$pathway <- rownames(matrix_df)
matrix_df$row_number <- 1:nrow(matrix_df)
clusters_withgroups <- merge(clustered_pathways, matrix_df, by.x = "pathway", by.y = "pathway")
#write.csv(clusters_withgroups,"/restricted/projectnb/cteseq/projects/somascan/final_plots/clustered_fgsea_using_function_nofaq_withproteins.csv")

Heatmap(
  pathway_gene_mat,
  cluster_rows = row_hc,
  cluster_columns = col_hc,
  use_raster = TRUE,          # <- key for large matrices
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = c("blue", "white", "red")
)

#final grouped results sheet from clustering
full_grouped_paths <- fread("/restricted/projectnb/cteseq/projects/somascan/final_plots/clustered_fgsea_using_function_nofaq_withproteins_groupnames.csv")
groups <- full_grouped_paths %>% 
  select(pathway, Group, NES)

head(groups)
#group_comp <- merge(groups,clustered_pathways, by.x = "Model_Pathway", by.y = "pathway")
#group_comp <- group_comp[order(group_comp$cluster), ]
#write.csv(group_comp,"/restricted/projectnb/cteseq/projects/somascan/final_plots/clustered_fgsea_functionvmanual.csv")

groups <- groups %>%
  mutate(pathway = sub("^AT8_", "AT8 ", pathway)) %>%
  separate(pathway, into = c("Model", "Pathway"), sep = "_", extra = "merge")

group_colors <- c(
  "Rho GTPase/ Mitotic Regulation" = "orange",
  "MAPK and PI3K/AKT Signaling Network" = "yellow",
  "Complement and Coagulation Cascade" = "green",
  "Ubiquitin-Proteasome System" = "lightblue",
  "Lysosomal Glycan Metabolism" = "blue",
  "Growth Factor Signaling" = "purple",
  "Translation Stress and Ribosome Dysfunction" = "pink",
  "mRNA Processing" = "cyan",
  "Innate Immune and Lysosomal Activation" = "#808080",   # Gray
  "Other" = "#000000",  # Black
  "Extracellular Matrix" = "#8B4513",  # Brown
  "Immune/ Secretory Trafficking" = "#56B4E9",  # Sky Blue
  "Lipid Metabolsim and Signaling" = "#009E73",  # Bluish Green
  "DNA Damage Response" = "#D55E00",  # Vermillion
  "Purine Metabolsim and Oxidative Stress" = "#CC79A7",  # Reddish Purple
  "Transcriptional and Epigenetic Regulation" = "#DAA520"   # Goldenrod
)

#CTE models
rhivslow_groups <- groups[which(groups$Model == "RHIvsLowCTE"),]
rhivslow_groups <- full_join(rhivslow_groups, fgsea_rl, by = c("Pathway" = "pathway"))
rhivshigh_groups <- groups[which(groups$Model == "RHIvsHighCTE"),]
rhivshigh_groups <- full_join(rhivshigh_groups, fgsea_rh, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  rhivslow_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_rl = NES.y, padj_rl = padj, leadingEdge_rl = leadingEdge, Group_rl = Group),
  rhivshigh_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_rh = NES.y, padj_rh = padj, leadingEdge_rh = leadingEdge, Group_rh = Group),
  by = "Pathway"
)

merged_df <- merged_df[which(merged_df$padj_rl < 0.05 | merged_df$padj_rh < 0.05),]
merged_df$Group_rh == merged_df$Group_rl
merged_df$Group <- ifelse(is.na(merged_df$Group_rl), merged_df$Group_rh, merged_df$Group_rl)
merged_df$Group <- ifelse(is.na(merged_df$Group_rh), merged_df$Group_rl, merged_df$Group_rh)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_rl < 0.05 & padj_rh < 0.05 ~ "Both",
    padj_rl < 0.05 & padj_rh >= 0.05 ~ "RHI vs Low CTE only",
    padj_rl >= 0.05 & padj_rh < 0.05 ~ "RHI vs High CTE only",
    TRUE ~ "Not significant"
  ))
shape_values <- c(
  "Both" = 16,                     # Circle
  "RHI vs Low CTE only" = 17,     # Triangle
  "RHI vs High CTE only" = 15,    # Square
  "Not significant" = 20          # Small dot
)
size_values <- c(
  "Both" = 10,
  "RHI vs Low CTE only" = 5,
  "RHI vs High CTE only" = 5,
  "Not significant" = 1
)
p <- ggplot(merged_df, aes(
  x = NES_rl,
  y = NES_rh,
  color = Group,
  shape = Significance,
  size  = Significance
)) +
  geom_point(alpha = 1.5) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_color_manual(values = group_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(
    axis.text   = element_text(size = 20),
    axis.title  = element_text(size = 20),
    legend.title= element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  labs(
    x = "NES (RHI to Low CTE)",
    y = "NES (RHI to High CTE)"
  )
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CTE_main_models_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CTE_main_models_fgsea_comp.rds")

#AT8 and totyrs
#AT8_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "AT8")
AT8_groups <- groups[which(groups$Model == "AT8 total"),]
#AT8_groups$Pathway <- sub("^total_", "", AT8_groups$Pathway)
AT8_groups <- full_join(AT8_groups, fgsea_AT8, by = c("Pathway" = "pathway"))
#totyrs_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "totyrs")
totyrs_groups <- groups[which(groups$Model == "totyrs"),]
totyrs_groups <- full_join(totyrs_groups, fgsea_totyrs, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  AT8_groups %>% dplyr::select(Pathway,NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_AT8 = NES.y, padj_AT8 = padj, leadingEdge_AT8 = leadingEdge, Group_AT8 = Group),
  totyrs_groups %>% dplyr::select(Pathway,NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_totyrs = NES.y, padj_totyrs = padj, leadingEdge_totyrs = leadingEdge, Group_totyrs = Group),
  by = "Pathway"
)

merged_df <- merged_df[which(merged_df$padj_AT8 < 0.05 | merged_df$padj_totyrs < 0.05),]
merged_df$Group <- ifelse(is.na(merged_df$Group_AT8), merged_df$Group_totyrs, merged_df$Group_AT8)
merged_df$Group <- ifelse(is.na(merged_df$Group_totyrs), merged_df$Group_AT8, merged_df$Group_totyrs)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_AT8 < 0.05 & padj_totyrs < 0.05 ~ "Both",
    padj_AT8 < 0.05 & padj_totyrs >= 0.05 ~ "AT8 Total only",
    padj_AT8 >= 0.05 & padj_totyrs < 0.05 ~ "Total Years of Play only",
    TRUE ~ "Not significant"
  ))
shape_values <- c(
  "Both" = 16,                     # Circle
  "AT8 Total only" = 17,     # Triangle
  "Total Years of Play only" = 15,    # Square
  "Not significant" = 20          # Small dot
)
size_values <- c(
  "Both" = 10,
  "AT8 Total only" = 5,
  "Total Years of Play only" = 5,
  "Not significant" = 1
)
p <- ggplot(merged_df, aes(x = NES_AT8, y = NES_totyrs, color = Group, shape = Significance, size = Significance)) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_color_manual(values = group_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  labs(x = "NES (AT8 Total)",
       y = "NES (Total Years of Play)")
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/AT8_totyrs_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/AT8_totyrs_fgsea_comp.rds")

#Cognitive Tests
#CDS and dementia
#cds_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "cds")
cds_groups <- groups[which(groups$Model == "CDStot"),]
cds_groups <- full_join(cds_groups, fgsea_cds, by = c("Pathway" = "pathway"))
#dem_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "dementia")
dem_groups <- groups[which(groups$Model == "Dementia"),]
dem_groups <- full_join(dem_groups, fgsea_dem, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  cds_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_cds = NES.y, padj_cds = padj, leadingEdge_cds = leadingEdge, Group_cds = Group),
  dem_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_dem = NES.y, padj_dem = padj, leadingEdge_dem = leadingEdge, Group_dem = Group),
  by = "Pathway"
)
merged_df <- merged_df[which(merged_df$padj_cds < 0.05 | merged_df$padj_dem < 0.05),]
merged_df$Group <- ifelse(is.na(merged_df$Group_cds), merged_df$Group_dem, merged_df$Group_cds)
merged_df$Group <- ifelse(is.na(merged_df$Group_dem), merged_df$Group_cds, merged_df$Group_dem)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_cds < 0.05 & padj_dem < 0.05 ~ "Both",
    padj_cds < 0.05 & padj_dem >= 0.05 ~ "CDS total only",
    padj_cds >= 0.05 & padj_dem < 0.05 ~ "Dementia only",
    TRUE ~ "Not significant"
  ))
merged_df <- merged_df[which(!merged_df$Significance == "Not significant"),]
shape_values <- c(
  "Both" = 16,                     # Circle
  "Dementia only" = 17,     # Triangle
  "CDS total only" = 15,    # Square
  "Not significant" = 20          # Small dot
)
size_values <- c(
  "Both" = 10,
  "Dementia only" = 5,
  "CDS total only" = 5,
  "Not significant" = 1
)
p <- ggplot(merged_df, aes(x = NES_cds, y = NES_dem, color = Group, shape = Significance, size = Significance)) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_color_manual(values = group_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 5))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  labs(x = "NES (CDS total)",
       y = "NES (Dementia)")
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CDS_dementia_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CDS_dementia_fgsea_comp.rds")

#CDS and FAQ
cds_groups <- groups[which(groups$Model == "CDStot"),]
cds_groups <- full_join(cds_groups, fgsea_cds, by = c("Pathway" = "pathway"))
faq_groups <- groups[which(groups$Model == "faq"),]
faq_groups <- full_join(faq_groups, fgsea_faq, by = c("Pathway" = "pathway"))
#faq_test <- faq_groups[which(faq_groups$padj < 0.05),]
#create plot
merged_df <- full_join(
  cds_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_cds = NES.y, padj_cds = padj, leadingEdge_cds = leadingEdge, Group_cds= Group),
  faq_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_faq = NES.y, padj_faq = padj, leadingEdge_faq = leadingEdge, Group_faq = Group),
  by = "Pathway"
)
merged_df <- merged_df[which(merged_df$padj_cds < 0.05 | merged_df$padj_faq < 0.05),]
merged_df$Group <- ifelse(is.na(merged_df$Group_cds), merged_df$Group_faq, merged_df$Group_cds)
merged_df$Group <- ifelse(is.na(merged_df$Group_faq), merged_df$Group_cds, merged_df$Group_faq)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_cds < 0.05 & padj_faq < 0.05 ~ "Both",
    padj_cds < 0.05 & padj_faq >= 0.05 ~ "CDS total only",
    padj_cds >= 0.05 & padj_faq < 0.05 ~ "FAQ total only",
    TRUE ~ "Not significant"
  ))
shape_values <- c(
  "Both" = 16,                     # Circle
  "FAQ total only" = 17,     # Triangle
  "CDS total only" = 15,    # Square
  "Not significant" = 20          # Small dot
)
size_values <- c(
  "Both" = 10,
  "FAQ total only" = 5,
  "CDS total only" = 5,
  "Not significant" = 1
)
ggplot(merged_df, aes(x = NES_cds, y = NES_faq, color = Group, shape = Significance, size = Significance)) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_color_manual(values = group_colors) +
  guides(
    color = guide_legend(override.aes = list(size = 5))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  labs(x = "NES (CDS total)",
       y = "NES (FAQ total)")
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CDS_FAQ_fgsea_comp_300dpi.png", plot = last_plot(), width = 12, height = 6, dpi = 300)

#Dementia and FAQ
dem_groups <- groups[which(groups$Model == "Dementia"),]
dem_groups <- full_join(dem_groups, fgsea_dem, by = c("Pathway" = "pathway"))
faq_groups <- groups[which(groups$Model == "faq"),]
faq_groups <- full_join(faq_groups, fgsea_faq, by = c("Pathway" = "pathway"))

groups_dem <- dem_groups[which(dem_groups$padj < 0.05),]
groups_faq <- faq_groups[which(faq_groups$padj < 0.05),]
#create plot
merged_df <- full_join(
  dem_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_dem = NES.y, padj_dem = padj, leadingEdge_dem = leadingEdge, Group_dem = Group),
  faq_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, Group) %>% dplyr::rename(NES_faq = NES.y, padj_faq = padj, leadingEdge_faq = leadingEdge, Group_faq = Group),
  by = "Pathway"
)
merged_df <- merged_df[which(merged_df$padj_dem < 0.05 | merged_df$padj_faq < 0.05),]
merged_df$Group <- ifelse(is.na(merged_df$Group_dem), merged_df$Group_faq, merged_df$Group_dem)
merged_df$Group <- ifelse(is.na(merged_df$Group_faq), merged_df$Group_dem, merged_df$Group_faq)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_dem < 0.05 & padj_faq < 0.05 ~ "Both",
    padj_dem < 0.05 & padj_faq >= 0.05 ~ "Dementia only",
    padj_dem >= 0.05 & padj_faq < 0.05 ~ "FAQ only",
    TRUE ~ "Not significant"
  ))
merged_df <- merged_df[which(!merged_df$Significance == "Not significant"),]
ggplot(merged_df, aes(x = NES_dem, y = NES_faq, color = Group, shape = Significance, size = Significance)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 10, "Dementia only" = 5, "FAQ only" = 5, "Not significant" = 1)) +
  scale_color_manual(values = c("orange", "yellow", "green","lightblue","blue","purple", "pink", "cyan","grey")) +
  guides(
    color = guide_legend(override.aes = list(size = 5))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  labs(x = "NES (Dementia)",
       y = "NES (FAQ)")
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/FAQ_Dementia_fgsea_comp_300dpi.png", plot = last_plot(), width = 12, height = 6, dpi = 300)

