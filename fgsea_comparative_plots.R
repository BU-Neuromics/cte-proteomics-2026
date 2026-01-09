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
fgsea_rh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvsHigh.csv")
fgsea_lh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_lowvshigh.csv")
fgsea_totyrs <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_totyrs.csv")
fgsea_AT8 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_AT8_total.csv")
fgsea_cds <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CDStot.csv")
fgsea_dem <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_DementiaHx.csv")

#final grouped results sheet from clustering
full_grouped_paths <- fread("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/clustered_fgsea_using_function_nofaq_withproteins_grouped.csv")
groups <- full_grouped_paths %>% 
  select(pathway, group, NES)

groups <- groups %>%
  mutate(pathway = sub("^AT8_", "AT8 ", pathway)) %>%
  separate(pathway, into = c("Model", "Pathway"), sep = "_", extra = "merge")

group_colors <- c(
  "Rho GTPase" = "orange",
  "MAPK and PI3K/AKT" = "yellow",
  "Complement System" = "green",
  "Proteasome" = "lightblue",
  "Lysosomal Glycan Metabolism" = "blue",
  "Growth Factor Signaling" = "purple",
  "Ribosome" = "deeppink", # Reddish Purple
  "mRNA Processing" = "cyan",
  "Other" = "grey",
  "Extracellular Matrix" = "#8B4513",  # Brown
  "DNA Damage Response" = "pink",
  "Immune/ Secretory Trafficking" = "#009E73",  # Bluish Green
  "Purine Metabolism" = "salmon"  # Reddish Purple
)

##creating plots##

#1 - CTE models
rhivslow_groups <- groups[which(groups$Model == "RHIvsLowCTE"),]
rhivslow_groups <- full_join(rhivslow_groups, fgsea_rl, by = c("Pathway" = "pathway"))
rhivshigh_groups <- groups[which(groups$Model == "RHIvsHighCTE"),]
rhivshigh_groups <- full_join(rhivshigh_groups, fgsea_rh, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  rhivslow_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_rl = NES.y, padj_rl = padj, leadingEdge_rl = leadingEdge, Group_rl = group),
  rhivshigh_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_rh = NES.y, padj_rh = padj, leadingEdge_rh = leadingEdge, Group_rh = group),
  by = "Pathway"
)

merged_df <- merged_df[which(merged_df$padj_rl < 0.05 | merged_df$padj_rh < 0.05),]
merged_df$Group_rh == merged_df$Group_rl
merged_df$group <- ifelse(is.na(merged_df$Group_rl), merged_df$Group_rh, merged_df$Group_rl)
merged_df$group <- ifelse(is.na(merged_df$Group_rh), merged_df$Group_rl, merged_df$Group_rh)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_rl < 0.05 & padj_rh < 0.05 ~ "Both",
    padj_rl < 0.05 & padj_rh >= 0.05 ~ "RHI vs Low CTE only",
    padj_rl >= 0.05 & padj_rh < 0.05 ~ "RHI vs High CTE only",
    TRUE ~ "Not significant"
  ))
shape_values <- c(
  "Both" = 21,                     # Circle
  "RHI vs Low CTE only" = 24,     # Triangle
  "RHI vs High CTE only" = 22,    # Square
  "Not significant" = 21          # Small dot
)
size_values <- c(
  "Both" = 10,
  "RHI vs Low CTE only" = 5,
  "RHI vs High CTE only" = 5,
  "Not significant" = 2
)
p <- ggplot(merged_df, aes(x = NES_rl, y = NES_rh)) +
  geom_point(aes(fill = group, shape =Significance, size = Significance),
             alpha = 0.8, color = "black") +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_fill_manual(values = group_colors) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "black",size = 5)),
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
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/CTE_main_models_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/CTE_main_models_fgsea_comp.rds")

#2 - AT8 and totyrs
AT8_groups <- groups[which(groups$Model == "AT8 total"),]
AT8_groups <- full_join(AT8_groups, fgsea_AT8, by = c("Pathway" = "pathway"))
totyrs_groups <- groups[which(groups$Model == "totyrs"),]
totyrs_groups <- full_join(totyrs_groups, fgsea_totyrs, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  AT8_groups %>% dplyr::select(Pathway,NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_AT8 = NES.y, padj_AT8 = padj, leadingEdge_AT8 = leadingEdge, Group_AT8 = group),
  totyrs_groups %>% dplyr::select(Pathway,NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_totyrs = NES.y, padj_totyrs = padj, leadingEdge_totyrs = leadingEdge, Group_totyrs = group),
  by = "Pathway"
)

merged_df <- merged_df[which(merged_df$padj_AT8 < 0.05 | merged_df$padj_totyrs < 0.05),]
merged_df$group <- ifelse(is.na(merged_df$Group_AT8), merged_df$Group_totyrs, merged_df$Group_AT8)
merged_df$group <- ifelse(is.na(merged_df$Group_totyrs), merged_df$Group_AT8, merged_df$Group_totyrs)
# Define significance status
merged_df <- merged_df %>%
  mutate(Significance = case_when(
    padj_AT8 < 0.05 & padj_totyrs < 0.05 ~ "Both",
    padj_AT8 < 0.05 & padj_totyrs >= 0.05 ~ "AT8 Total only",
    padj_AT8 >= 0.05 & padj_totyrs < 0.05 ~ "Total Years of Play only",
    TRUE ~ "Not significant"
  ))
shape_values <- c(
  "Both" = 21,                  # was 16
  "AT8 Total only" = 24,        # was 17
  "Total Years of Play only" = 22,  # was 15
  "Not significant" = 21        # was 20
)
size_values <- c(
  "Both" = 10,
  "AT8 Total only" = 5,
  "Total Years of Play only" = 5,
  "Not significant" = 2        # was 1, increased for border visibility
)
p <- ggplot(merged_df, aes(x = NES_AT8, y = NES_totyrs)) +
  geom_point(
    aes(fill = group, shape = Significance, size = Significance),  # map fill and shape
    colour = "black",       # added: black border
    alpha = 0.8
  ) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_fill_manual(values = group_colors) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,       # filled shape
        colour = "black", # border
        size = 5
      )
    ),
    shape = guide_legend(
      override.aes = list(size = 5)
    )
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
    x = "NES (AT8 Total)",
    y = "NES (Total Years of Play)"
  )
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/AT8_totyrs_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/AT8_totyrs_fgsea_comp.rds")

#3 - Cognitive Tests: CDS and dementia
cds_groups <- groups[which(groups$Model == "CDStot"),]
cds_groups <- full_join(cds_groups, fgsea_cds, by = c("Pathway" = "pathway"))
dem_groups <- groups[which(groups$Model == "Dementia"),]
dem_groups <- full_join(dem_groups, fgsea_dem, by = c("Pathway" = "pathway"))
#create plot
merged_df <- full_join(
  cds_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_cds = NES.y, padj_cds = padj, leadingEdge_cds = leadingEdge, Group_cds = group),
  dem_groups %>% dplyr::select(Pathway, NES.y, padj, leadingEdge, group) %>% dplyr::rename(NES_dem = NES.y, padj_dem = padj, leadingEdge_dem = leadingEdge, Group_dem = group),
  by = "Pathway"
)
merged_df <- merged_df[which(merged_df$padj_cds < 0.05 | merged_df$padj_dem < 0.05),]
merged_df$group <- ifelse(is.na(merged_df$Group_cds), merged_df$Group_dem, merged_df$Group_cds)
merged_df$group <- ifelse(is.na(merged_df$Group_dem), merged_df$Group_cds, merged_df$Group_dem)
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
  "Both" = 21,                     # Circle
  "Dementia only" = 24,     # Triangle
  "CDS total only" = 22,    # Square
  "Not significant" = 21          # Small dot
)
size_values <- c(
  "Both" = 10,
  "Dementia only" = 5,
  "CDS total only" = 5,
  "Not significant" = 1
)
p <- ggplot(merged_df, aes(x = NES_cds, y = NES_dem)) +
  geom_point(aes(fill = group, shape = Significance, size = Significance),
             color = "black", alpha = 0.8) +
  scale_shape_manual(values = shape_values) +
  scale_size_manual(values  = size_values) +
  scale_fill_manual(values = group_colors) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "black",size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
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
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/CDS_dementia_fgsea_comp_300dpi.png", plot = p, width = 12, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots/CDS_dementia_fgsea_comp.rds")
