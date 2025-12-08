#leading edge genes plots
#Helen Pennington

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

#CTE models
#proteomics
CTE_rl_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvslow.csv")
CTE_rl_paths <- CTE_rl_fgsea[which(CTE_rl_fgsea$padj < 0.05),]
leading_proteins_rl <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", CTE_rl_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_rl <- data.frame(protein = rownames(leading_proteins_rl),`RHI vs Low CTE Proteins` = leading_proteins_rl[,1])

CTE_rh_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvsHigh.csv")
CTE_rh_paths <- CTE_rh_fgsea[which(CTE_rh_fgsea$padj < 0.05),]
CTE_rh_paths$leadingEdge <- as.list(CTE_rh_paths$leadingEdge)
leading_proteins_rh <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", CTE_rh_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_rh <- data.frame(protein = rownames(leading_proteins_rh),`RHI vs High CTE Proteins` = leading_proteins_rh[,1])
# Merge together (full outer joins)
leading_edge <- full_join(leading_proteins_rl, leading_proteins_rh, by = "protein")

leading_edge[is.na(leading_edge)] <- 0
rownames(leading_edge) <- leading_edge$protein
leading_edge <- leading_edge[,-1]
leading_edge$sum <- rowSums(leading_edge)
leading_edge <- arrange(leading_edge, desc(sum))
leading_edge <- leading_edge[1:100,]

leading_edge <- leading_edge %>%
  tibble::rownames_to_column("protein")
leading_edge$protein <- factor(leading_edge$protein, levels = leading_edge$protein)

# Pivot longer to tidy format
leading_edge_long <- leading_edge %>%
  pivot_longer(cols = starts_with(c("RHI.vs.Low", "RHI.vs.High")),
               names_to = "model",
               values_to = "Appearances")
head(leading_edge_long)
leading_edge_long <- leading_edge_long %>%
  mutate(model = recode(model, "RHI.vs.Low.CTE.Proteins" = "RHI vs Low CTE Proteins", "RHI.vs.High.CTE.Proteins" = "RHI vs High CTE Proteins"))
# Make grouped horizontal bar plot
p <- ggplot(leading_edge_long, aes(x = Appearances,
                              y = protein,
                              fill = model)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  scale_fill_manual(values = c("blue","orange")) +
  labs(x = "Number of Appearances", y = "Protein", fill = "Model") +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/CTE_main_models_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/CTE_main_models_leading_edge.rds")
#################################################################################

#AT8 and totyrs
#proteomics
AT8_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_AT8_total.csv")
AT8_paths <- AT8_fgsea[which(AT8_fgsea$padj < 0.05),]
leading_proteins_AT8 <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", AT8_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_AT8 <- data.frame(protein = rownames(leading_proteins_AT8),`AT8 Total Proteins` = leading_proteins_AT8[,1])

totyrs_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_totyrs.csv")
totyrs_paths <- totyrs_fgsea[which(totyrs_fgsea$padj < 0.05),]
leading_proteins_totyrs <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", totyrs_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_totyrs <- data.frame(protein = rownames(leading_proteins_totyrs),`Total Years of Play Proteins` = leading_proteins_totyrs[,1])

# Merge all four together (full outer joins)
leading_edge <- full_join(leading_proteins_AT8, leading_proteins_totyrs, by = "protein")

leading_edge[is.na(leading_edge)] <- 0
rownames(leading_edge) <- leading_edge$protein
leading_edge <- leading_edge[,-1]
leading_edge$sum <- rowSums(leading_edge)
leading_edge_AT8 <- arrange(leading_edge, desc(AT8.Total.Proteins))
leading_edge_AT8 <- leading_edge_AT8[1:100,]
leading_edge_totyrs <- arrange(leading_edge, desc(Total.Years.of.Play.Proteins))
leading_edge_totyrs <- leading_edge_totyrs[1:100,]

leading_plot <- function(leading_edge){
  leading_edge <- leading_edge %>%
    tibble::rownames_to_column("protein")
  leading_edge$protein <- factor(leading_edge$protein, levels = leading_edge$protein)
  
  # Pivot longer to tidy format
  leading_edge_long <- leading_edge %>%
    pivot_longer(cols = starts_with(c("AT8", "Total")),
                 names_to = "model",
                 values_to = "Appearances")
  print(unique(leading_edge_long$model))
  leading_edge_long <- leading_edge_long %>%
    mutate(model = recode(model, "AT8.Total.Proteins" = "AT8 Total Proteins", "Total.Years.of.Play.Proteins" = "Total Years of Play Proteins"))
  # Make grouped horizontal bar plot
  p <- ggplot(leading_edge_long, aes(x = Appearances,
                                y = protein,
                                fill = model)) +
    geom_col(position = position_dodge(width = 0.8)) +
    coord_flip() +
    scale_fill_manual(values = c("blue","orange")) +
    labs(x = "Number of Appearances", y = "Protein", fill = "Model") +
    theme_minimal(base_size = 20) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  return(p)
}
p <- leading_plot(leading_edge_AT8)
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/AT8_main_totyrs_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/AT8_main_totyrs_leading_edge.rds")
p <- leading_plot(leading_edge_totyrs)
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/totyrs_main_AT8_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/totyrs_main_AT8_leading_edge.rds")
#################################################################################

#cognitive tests
#proteomics
cds_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CDStot.csv")
cds_paths <- cds_fgsea[which(cds_fgsea$padj < 0.05),]
leading_proteins_cds <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", cds_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_cds <- data.frame(protein = rownames(leading_proteins_cds),`CDS Total Proteins` = leading_proteins_cds[,1])

dem_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_DementiaHx.csv")
dem_paths <- dem_fgsea[which(dem_fgsea$padj < 0.05),]
leading_proteins_dem <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", dem_paths$leadingEdge), ","))), decreasing = TRUE))
leading_proteins_dem <- data.frame(protein = rownames(leading_proteins_dem),`Dementia Proteins` = leading_proteins_dem[,1])

#faq_fgsea <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_faqtot.csv")
#faq_paths <- faq_fgsea[which(faq_fgsea$padj < 0.05),]
#leading_proteins_faq <- as.matrix(sort(table(unlist(strsplit(gsub("\\s+", "", faq_paths$leadingEdge), ","))), decreasing = TRUE))
#leading_proteins_faq <- data.frame(protein = rownames(leading_proteins_faq),`FAQ Total Proteins` = leading_proteins_faq[,1])

# Merge all four together (full outer joins)
leading_edge <- full_join(leading_proteins_cds, leading_proteins_dem, by = "protein")

leading_edge[is.na(leading_edge)] <- 0
rownames(leading_edge) <- leading_edge$protein
leading_edge <- leading_edge[,-1]
leading_edge$sum <- rowSums(leading_edge)
leading_edge_dementia <- arrange(leading_edge, desc(Dementia.Proteins))
leading_edge_dementia <- leading_edge_dementia[1:100,]
leading_edge_cds <- arrange(leading_edge, desc(CDS.Total.Proteins))
leading_edge_cds <- leading_edge_cds[1:100,]

leading_plot <- function(leading_edge){
  leading_edge <- leading_edge %>%
    tibble::rownames_to_column("protein")
  leading_edge$protein <- factor(leading_edge$protein, levels = leading_edge$protein)
  
  # Pivot longer to tidy format
  leading_edge_long <- leading_edge %>%
    pivot_longer(cols = starts_with(c("CDS.", "Dementia.")),
                 names_to = "model",
                 values_to = "Appearances")
  print(unique(leading_edge_long$model))
  leading_edge_long <- leading_edge_long %>%
    mutate(model = recode(model, "CDS.Total.Proteins" = "CDS Total Proteins", "Dementia.Proteins" = "Dementia Proteins"))
  p <- ggplot(leading_edge_long, aes(x = Appearances,
                                     y = protein,
                                     fill = model)) +
    geom_col(position = position_dodge(width = 0.8)) +
    coord_flip() +
    scale_fill_manual(values = c("blue","orange")) +
    labs(x = "Number of Appearances", y = "Protein", fill = "Model") +
    theme_minimal(base_size = 20) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  return(p)
}
p <- leading_plot(leading_edge_dementia)
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/dementia_main_cds_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/dementia_main_cds_leading_edge.rds")
p <- leading_plot(leading_edge_cds)
ggsave("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/cds_main_dementia_leading_edge_300dpi.png", plot = p, width = 20, height = 6, dpi = 300)
saveRDS(p, "/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/cds_main_dementia_leading_edge.rds")
