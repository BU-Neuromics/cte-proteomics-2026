#leading edge genes plots
#Helen Pennington

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

full_grouped_paths <- fread("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files/clustered_fgsea_using_function_nofaq_withproteins_grouped.csv")

groups <- full_grouped_paths %>%
  select(-V1, -cluster, -NES)

# Convert wide format to long format
long_df <- groups %>%
  pivot_longer(
    cols = -c(group, pathway),
    names_to = "Protein",
    values_to = "value"
  )

# Convert to binary leading edge presence (1 if in leading edge, 0 if not)
long_df <- long_df %>%
  mutate(in_leading_edge = ifelse(value != 0, 1, 0))

# Count number of pathways in each group per protein
protein_counts <- long_df %>%
  group_by(Protein, group) %>%
  summarise(count = sum(in_leading_edge), .groups = "drop")

# Pivot back to wide format with proteins as rows and groups as columns
protein_counts_wide <- protein_counts %>%
  pivot_wider(
    names_from = group,
    values_from = count,
    values_fill = 0
  )

protein_counts_wide

# Get the group with the maximum count for each protein
protein_max_group <- protein_counts_wide %>%
  rowwise() %>%
  mutate(
    max_group = names(select(cur_data(), -Protein))[which.max(c_across(-Protein))]
  ) %>%
  ungroup() %>%
  select(Protein, max_group)

protein_max_group

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

# prep models
CTE_rl_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvslow.csv")
CTE_rl_paths <- CTE_rl_fgsea[which(CTE_rl_fgsea$padj < 0.05),]
leading_proteins_rl <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", CTE_rl_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_rl <- data.frame(protein = rownames(leading_proteins_rl),`RHI vs Low CTE Proteins` = leading_proteins_rl[,1])

CTE_rh_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvsHigh.csv")
CTE_rh_paths <- CTE_rh_fgsea[which(CTE_rh_fgsea$padj < 0.05),]
#CTE_rh_paths$leadingEdge <- as.list(CTE_rh_paths$leadingEdge)
leading_proteins_rh <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", CTE_rh_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_rh <- data.frame(protein = rownames(leading_proteins_rh),`RHI vs High CTE Proteins` = leading_proteins_rh[,1])

AT8_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_AT8_total.csv")
AT8_paths <- AT8_fgsea[which(AT8_fgsea$padj < 0.05),]
leading_proteins_AT8 <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", AT8_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_AT8 <- data.frame(protein = rownames(leading_proteins_AT8),`AT8 Total Proteins` = leading_proteins_AT8[,1])

totyrs_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_totyrs.csv")
totyrs_paths <- totyrs_fgsea[which(totyrs_fgsea$padj < 0.05),]
leading_proteins_totyrs <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", totyrs_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_totyrs <- data.frame(protein = rownames(leading_proteins_totyrs),`Total Years of Play Proteins` = leading_proteins_totyrs[,1])

cds_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CDStot.csv")
cds_paths <- cds_fgsea[which(cds_fgsea$padj < 0.05),]
leading_proteins_cds <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", cds_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_cds <- data.frame(protein = rownames(leading_proteins_cds),`CDS Total Proteins` = leading_proteins_cds[,1])

dem_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_DementiaHx.csv")
dem_paths <- dem_fgsea[which(dem_fgsea$padj < 0.05),]
leading_proteins_dem <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", dem_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_dem <- data.frame(protein = rownames(leading_proteins_dem),`Dementia Proteins` = leading_proteins_dem[,1])

#function
plot_leading_edge_single_model <- function(
    leading_df,
    protein_col = "protein",
    model_col,
    protein_max_group,
    group_colors,
    top_n = 100
) {
  
  plot_data <- leading_df %>%
    select(
      protein = all_of(protein_col),
      Appearances = all_of(model_col)
    ) %>%
    arrange(desc(Appearances)) %>%
    slice_head(n = top_n) %>%
    left_join(protein_max_group, by = c("protein" = "Protein")) %>%
    mutate(
      protein = factor(protein, levels = protein)
    )
  
  ggplot(plot_data,
         aes(x = Appearances,
             y = protein,
             fill = max_group)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = group_colors) +
    labs(
      x = "Number of Appearances",
      y = "Protein",
      fill = "Group"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 10)
    )
}

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_rl,
  model_col = "RHI.vs.Low.CTE.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/lowCTE_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/lowCTE_leading_edge.rds")

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_rh,
  model_col = "RHI.vs.High.CTE.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/highCTE_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/highCTE_leading_edge.rds")

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_AT8,
  model_col = "AT8.Total.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/AT8_total_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/AT8_total_leading_edge.rds")

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_totyrs,
  model_col = "Total.Years.of.Play.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/totyrs_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/totyrs_leading_edge.rds")

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_dem,
  model_col = "Dementia.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/Dementia_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/Dementia_leading_edge.rds")

p <- plot_leading_edge_single_model(
  leading_df = leading_proteins_cds,
  model_col = "CDS.Total.Proteins",
  protein_max_group = protein_max_group,
  group_colors = group_colors,
  top_n = 100
)
p
ggsave("/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/CDS_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots/CDS_leading_edge.rds")

