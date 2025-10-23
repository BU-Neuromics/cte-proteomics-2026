#lfc proteomics vs rna
#Helen Pennington
#Started August 28, 2025

#packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(stringr)

#RHI vs Low CTE
protein_rh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_CTE_RHIvsHigh_noF_results_coef.csv")

#RHI vs Low CTE
protein_rl <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_CTE_RHIvslow_noF_results_coef.csv")
rna_rl <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ rhivslow .csv")
rna_rl_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ rhivslow _no_shrink.csv")
protein_rl_filt <- protein_rl[which(abs(protein_rl$logFC) < 1),]
rna_rl_ns_filt <- rna_rl_ns[which(abs(rna_rl_ns$log2FoldChange) < 1),]

#Low CTE vs High CTE
protein_lh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_CTE_lowvshigh_noF_results_coef.csv")
rna_lh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ lowvshigh .csv")
rna_lh_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ lowvshigh _no_shrink.csv")
protein_lh_filt <- protein_lh[which(abs(protein_lh$logFC) < 1),]
rna_lh_ns_filt <- rna_lh_ns[which(abs(rna_lh_ns$log2FoldChange) < 1),]

#scatterplot
merged_rl <- merge(protein_rl, rna_rl_ns, by.x="genes", by.y="gene_name")
head(merged_rl)
plot(merged_rl$log2FoldChange, merged_rl$logFC, main = "RHI vs Low CTE")

merged_lh <- merge(protein_lh, rna_lh_ns, by.x="genes", by.y="gene_name")
head(merged_lh)
plot(merged_lh$log2FoldChange, merged_lh$logFC, main = "Low CTE vs High CTE")

#scatterplots filtered
merged_rl_filt <- merge(protein_rl_filt, rna_rl_ns_filt, by.x="genes", by.y="gene_name")
head(merged_rl_filt)
plot(merged_rl_filt$log2FoldChange, merged_rl_filt$logFC, main = "RHI vs Low CTE Filtered")

merged_lh_filt <- merge(protein_lh_filt, rna_lh_ns_filt, by.x="genes", by.y="gene_name")
head(merged_lh_filt)
plot(merged_lh_filt$log2FoldChange, merged_lh_filt$logFC, main = "Low CTE vs High CTE Filtered")

#AT8 total
protein_AT8 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_AT8_total_noF_results_coef.csv")
rna_AT8 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ AT8_total .csv")
rna_AT8_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ AT8_total _no_shrink.csv")
protein_AT8_filt <- protein_AT8[which(abs(protein_AT8$logFC) < 1),]
rna_AT8_ns_filt <- rna_AT8_ns[which(abs(rna_AT8_ns$log2FoldChange) < 1),]

#totyrs
protein_totyrs <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_totyrs_noF_results_coef.csv")
rna_totyrs <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ totyrs .csv")
rna_totyrs_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ totyrs _no_shrink.csv")
protein_totyrs_filt <- protein_totyrs[which(abs(protein_totyrs$logFC) < 1),]
rna_totyrs_ns_filt <- rna_totyrs_ns[which(abs(rna_totyrs_ns$log2FoldChange) < 1),]

#scatterplots filtered
merged_AT8 <- merge(protein_AT8, rna_AT8_ns, by.x="genes", by.y="gene_name")
head(merged_AT8)
plot(merged_AT8$log2FoldChange, merged_AT8$logFC, main = "AT8 Total")

merged_totyrs <- merge(protein_totyrs, rna_totyrs_ns, by.x="genes", by.y="gene_name")
head(merged_totyrs)
plot(merged_totyrs$log2FoldChange, merged_totyrs$logFC, main = "Total Years of Play")

#scatterplots filtered
merged_AT8_filt <- merge(protein_AT8_filt, rna_AT8_ns_filt, by.x="genes", by.y="gene_name")
head(merged_AT8_filt)
plot(merged_AT8_filt$log2FoldChange, merged_AT8_filt$logFC, main = "AT8 Total Filtered")

merged_totyrs_filt <- merge(protein_totyrs_filt, rna_totyrs_ns_filt, by.x="genes", by.y="gene_name")
head(merged_totyrs_filt)
plot(merged_totyrs_filt$log2FoldChange, merged_totyrs_filt$logFC, main = "Total Years of Play Filtered")

#CDStot
protein_cds <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_CDStot_noF_results_coef.csv")
rna_cds <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ CDStot .csv")
rna_cds_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ CDStot _no_shrink.csv")
protein_cds_filt <- protein_cds[which(abs(protein_cds$logFC) < 1),]
rna_cds_ns_filt <- rna_cds_ns[which(abs(rna_cds_ns$log2FoldChange) < 1),]

#scatterplot
merged_cds <- merge(protein_cds, rna_cds_ns, by.x="genes", by.y="gene_name")
head(merged_cds)
plot(merged_cds$log2FoldChange, merged_cds$logFC, main = "CDS Total")

#scatterplot filtered
merged_cds_filt <- merge(protein_cds_filt, rna_cds_ns_filt, by.x="genes", by.y="gene_name")
head(merged_cds_filt)
plot(merged_cds_filt$log2FoldChange, merged_cds_filt$logFC, main = "CDS Total Filtered")

#Dementia
protein_d <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_DementiaHx_noF_results_coef.csv")
rna_d <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ DementiaHx .csv")
rna_d_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ DementiaHx _no_shrink.csv")
protein_d_filt <- protein_d[which(abs(protein_d$logFC) < 1),]
rna_d_ns_filt <- rna_d_ns[which(abs(rna_d_ns$log2FoldChange) < 1),]

#scatterplot
merged_d <- merge(protein_d, rna_d_ns, by.x="genes", by.y="gene_name")
head(merged_d)
plot(merged_d$log2FoldChange, merged_d$logFC, main = "Dementia")

#scatterplot
merged_d_filt <- merge(protein_d_filt, rna_d_ns_filt, by.x="genes", by.y="gene_name")
head(merged_d_filt)
plot(merged_d_filt$log2FoldChange, merged_d_filt$logFC, main = "Dementia Filtered")

#faqtot
protein_faq <- fread("/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_faqtot_noF_results_coef.csv")
rna_faq <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ faqtot .csv")
rna_faq_ns <- fread("/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_ faqtot _no_shrink.csv")
protein_faq_filt <- protein_faq[which(abs(protein_faq$logFC) < 1),]
rna_faq_ns_filt <- rna_faq_ns[which(abs(rna_faq_ns$log2FoldChange) < 1),]

#scatterplot
merged_faq <- merge(protein_faq, rna_faq_ns, by.x="genes", by.y="gene_name")
head(merged_faq)
plot(merged_faq$log2FoldChange, merged_faq$logFC, main = "FAQ Total")

#scatterplot filtered
merged_faq_filt <- merge(protein_faq_filt, rna_faq_ns_filt, by.x="genes", by.y="gene_name")
head(merged_faq_filt)
plot(merged_faq_filt$log2FoldChange, merged_faq_filt$logFC, main = "FAQ Total Filtered")
#################################################################################

#paried dot plot
merged <- merge(protein_rl, protein_lh, by = "unique_names", suffixes = c("_rl", "_lh"))
#proteasome/ubiquitin proteins
merged <- merged %>% filter(genes_rl %in% c("PSMA5","PSMA7","PSMB1","PSMB2","PSMB3","PSMC5","RPN1","SMURF1","KEAP1",
                                  "ATG7","GAN","KCTD6","FBXL4","KLHL13","USP10","UBE2M","UBE2C","ANAPC10",
                                  "IFNG","UBL4A","UBTD2","UCHL5","UBB"))
#synapse proteins
#merged <- merged %>% filter(genes_rl %in% c("LRRTM2","HOMER1","HOMER2","HOMER3","SYT7","SHANK1","SLITRK1","SLITRK5","SLITRK6"))
long_df <- merged %>%
  pivot_longer(cols = c(logFC_rl, logFC_lh),
               names_to = "model", values_to = "logFC") %>%
  mutate(model = ifelse(model == "logFC_rl", "RHI vs Low CTE", "Low vs High CTE"),
  model = factor(model, levels = c("RHI vs Low CTE", "Low vs High CTE")))

# Plot
ggplot(long_df, aes(x = model, y = logFC, group = unique_names)) +
  geom_line(alpha = 0.4) +              # line connecting the two points per protein
  geom_point(size = 2) +# points for LFC values
  geom_text_repel(aes(label = genes_rl), size = 3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Model", y = "log Fold Change (LFC)",
       title = "Paired LFC comparison across models") +
  theme_minimal()

#fgsea scatterplots
#CTE models
fgsea_rl <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_CTE_RHIvslow_noF_gseaResults.csv")
fgsea_lh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_CTE_lowvshigh_noF_gseaResults.csv")

#only keep unique pathways (they are in long format)
fgsea_rl_unique <- fgsea_rl[!duplicated(fgsea_rl$pathway), ]
fgsea_lh_unique <- fgsea_lh[!duplicated(fgsea_lh$pathway), ]
#fgsea_rl_u_sig <- fgsea_rl_unique[which(fgsea_rl_unique$padj < 0.05),]
#fgsea_lh_u_sig <- fgsea_lh_unique[which(fgsea_lh_unique$padj < 0.05),]

fgsea_rl <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_RHIvslow.rds")
fgsea_lh <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_lowvshigh.rds")
fgsea_rl <- fgsea_rl[order(fgsea_rl$pval), ]
fgsea_lh <- fgsea_lh[order(fgsea_lh$pval), ]
fgsea_rl$synapse <- sapply(fgsea_rl$leadingEdge, function(x) "HOMER1" %in% x)
fgsea_rl$proteasome <- sapply(fgsea_rl$leadingEdge, function(x) "PSMB1" %in% x)
fgsea_rl$ribosome <- sapply(fgsea_rl$leadingEdge, function(x) "RPS12" %in% x)
fgsea_rl$group <- ifelse(fgsea_rl$synapse == TRUE, "Synapse", 
                         ifelse(fgsea_rl$proteasome == TRUE,"Proteasome",
                                ifelse(fgsea_rl$ribosome == TRUE,"Ribosome","Other")))
fgsea_lh$synapse <- sapply(fgsea_lh$leadingEdge, function(x) "HOMER1" %in% x)
fgsea_lh$proteasome <- sapply(fgsea_lh$leadingEdge, function(x) "PSMB1" %in% x)
fgsea_lh$ribosome <- sapply(fgsea_lh$leadingEdge, function(x) "RPS12" %in% x)
fgsea_lh$group <- ifelse(fgsea_lh$synapse == TRUE, "Synapse", 
                         ifelse(fgsea_lh$proteasome == TRUE,"Proteasome",
                                ifelse(fgsea_lh$ribosome == TRUE,"Ribosome","Other")))
#create plot
merged_df <- inner_join(
  fgsea_rl %>% dplyr::select(pathway, NES, padj, group, leadingEdge) %>% dplyr::rename(NES_rl = NES, padj_rl = padj, group_rl = group, leadingEdge_rl = leadingEdge),
  fgsea_lh %>% dplyr::select(pathway, NES, padj, group, leadingEdge) %>% dplyr::rename(NES_lh = NES, padj_lh = padj, group_lh = group, leadingEdge_lh = leadingEdge),
  by = "pathway"
)

merged_df <- merged_df[which(merged_df$padj_rl < 0.05 | merged_df$padj_lh < 0.05),]
merged_df$group_lh == merged_df$group_rl
merged_df$group <- ifelse(merged_df$group_lh == merged_df$group_rl, merged_df$group_rl, ifelse(merged_df$group_lh != "Other",merged_df$group_lh, merged_df$group_rl))
merged_df$group[which(merged_df$pathway == "KEGG_ALZHEIMERS_DISEASE" | merged_df$pathway == "KEGG_PARKINSONS_DISEASE" | merged_df$pathway == "KEGG_HUNTINGTONS_DISEASE" | merged_df$pathway == "WP_METABOLIC_EPILEPTIC_DISORDERS")] <- "Neurological Disorders"
merged_df$group[which(merged_df$pathway == "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS" | merged_df$pathway == "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES")] <- "Synapse"
# Define significance status
merged_df <- merged_df %>%
  mutate(sig_status = case_when(
    padj_rl < 0.05 & padj_lh < 0.05 ~ "Both",
    padj_rl < 0.05 & padj_lh >= 0.05 ~ "Model1 only",
    padj_rl >= 0.05 & padj_lh < 0.05 ~ "Model2 only",
    TRUE ~ "Not significant"
  ))

ggplot(merged_df, aes(x = NES_rl, y = NES_lh, color = group, shape = sig_status, size = sig_status)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 5, "Model1 only" = 2, "Model2 only" = 2, "Not significant" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(x = "NES (RHI to Low CTE)",
       y = "NES (Low CTE to High CTE)",
       title = "Comparison of NES values between RHI to Low CTE and Low CTE to High CTE")

#totyrs and AT8
fgsea_totyrs <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_totyrs_noF_gseaResults.csv")
fgsea_AT8 <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_AT8_total_CTE_noF_gseaResults.csv")

#only keep unique pathways (they are in long format)
fgsea_totyrs_unique <- fgsea_totyrs[!duplicated(fgsea_totyrs$pathway), ]
fgsea_AT8_unique <- fgsea_AT8[!duplicated(fgsea_AT8$pathway), ]
#fgsea_totyrs_u_sig <- fgsea_totyrs_unique[which(fgsea_totyrs_unique$padj < 0.05),]
#fgsea_AT8_u_sig <- fgsea_AT8_unique[which(fgsea_AT8_unique$padj < 0.05),]

fgsea_totyrs <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_totyrs.rds")
fgsea_AT8 <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_AT8.rds")

fgsea_rl <- fgsea_rl[order(fgsea_rl$pval), ]
fgsea_lh <- fgsea_lh[order(fgsea_lh$pval), ]
fgsea_rl$synapse <- sapply(fgsea_rl$leadingEdge, function(x) "HOMER1" %in% x)
fgsea_rl$proteasome <- sapply(fgsea_rl$leadingEdge, function(x) "PSMB1" %in% x)
fgsea_rl$ribosome <- sapply(fgsea_rl$leadingEdge, function(x) "RPS12" %in% x)
fgsea_rl$group <- ifelse(fgsea_rl$synapse == TRUE, "Synapse", 
                         ifelse(fgsea_rl$proteasome == TRUE,"Proteasome",
                                ifelse(fgsea_rl$ribosome == TRUE,"Ribosome","Other")))
fgsea_lh$synapse <- sapply(fgsea_lh$leadingEdge, function(x) "HOMER1" %in% x)
fgsea_lh$proteasome <- sapply(fgsea_lh$leadingEdge, function(x) "PSMB1" %in% x)
fgsea_lh$ribosome <- sapply(fgsea_lh$leadingEdge, function(x) "RPS12" %in% x)
fgsea_lh$group <- ifelse(fgsea_lh$synapse == TRUE, "Synapse", 
                         ifelse(fgsea_lh$proteasome == TRUE,"Proteasome",
                                ifelse(fgsea_lh$ribosome == TRUE,"Ribosome","Other")))

merged_df <- inner_join(
  fgsea_totyrs %>% dplyr::select(pathway, NES, padj, leadingEdge) %>% dplyr::rename(NES_totyrs = NES, padj_totyrs = padj, leadingEdge_totyrs = leadingEdge),
  fgsea_AT8 %>% dplyr::select(pathway, NES, padj, leadingEdge) %>% dplyr::rename(NES_AT8 = NES, padj_AT8 = padj, leadingEdge_AT8 = leadingEdge),
  by = "pathway"
)

merged_df <- merged_df[which(merged_df$padj_totyrs < 0.05 | merged_df$padj_AT8 < 0.05),]

ggplot(merged_df, aes(x = NES_totyrs, y = NES_AT8, label = pathway)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_text(vjust = -0.5, size = 3, check_overlap = TRUE) +  # labels each pathway
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(x = "NES (Total Years of Play)",
       y = "NES (AT8 total)",
       title = "Comparison of NES values between Total Years of Play and AT8 Total")
#################################################################################
#leading edge genes 
#CTE models
#proteomics
CTE_rl_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_RHIvslow.rds")
#CTE_rl_paths_p <- CTE_rl_fgsea_p[which(CTE_rl_fgsea_p$NES > 0 & CTE_rl_fgsea_p$pval < 0.05),]
CTE_rl_paths_p <- CTE_rl_fgsea_p[which(CTE_rl_fgsea_p$pval < 0.05),]
leading_proteins_rl <- as.matrix(sort(table(unlist(CTE_rl_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_rl <- data.frame(protein = rownames(leading_proteins_rl),`RHI vs Low CTE Proteins` = leading_proteins_rl[,1])

CTE_lh_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_lowvshigh.rds")
CTE_lh_paths_p <- CTE_lh_fgsea_p[which(CTE_lh_fgsea_p$pval < 0.05),]
leading_proteins_lh <- as.matrix(sort(table(unlist(CTE_lh_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_lh <- data.frame(protein = rownames(leading_proteins_lh),`Low vs High CTE Proteins` = leading_proteins_lh[,1])

#rna-seq
CTE_rl_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_RHIvslow_rnaseq.rds")
CTE_rl_paths_r <- CTE_rl_fgsea_r[which(CTE_rl_fgsea_r$pval < 0.05),]
leading_genes_rl <- as.matrix(sort(table(unlist(CTE_rl_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_rl <- data.frame(protein = rownames(leading_genes_rl),`RHI vs Low CTE Genes` = leading_genes_rl[,1])

CTE_lh_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_lowvshigh_rnaseq.rds")
CTE_lh_paths_r <- CTE_lh_fgsea_r[which(CTE_lh_fgsea_r$pval < 0.05),]
leading_genes_lh <- as.matrix(sort(table(unlist(CTE_lh_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_lh <- data.frame(protein = rownames(leading_genes_lh),`Low vs High CTE Genes` = leading_genes_lh[,1])

# Merge all four together (full outer joins)
leading_edge <- full_join(leading_proteins_rl, leading_proteins_lh, by = "protein") %>%
  full_join(leading_genes_rl, by = "protein") %>%
  full_join(leading_genes_lh, by = "protein")

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
  pivot_longer(cols = starts_with(c("Low.", "RHI.")),
               names_to = "model",
               values_to = "Appearances")

# Make grouped horizontal bar plot
ggplot(leading_edge_long, aes(x = Appearances,
                              y = protein,
                              fill = model)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(x = "Number of Appearances", y = "Protein/Gene", fill = "Model") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

#################################################################################

#AT8 and totyrs
#proteomics
AT8_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_AT8.rds")
AT8_paths_p <- AT8_fgsea_p[which(AT8_fgsea_p$pval < 0.05),]
leading_proteins_AT8 <- as.matrix(sort(table(unlist(AT8_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_AT8 <- data.frame(protein = rownames(leading_proteins_AT8),`AT8 Total Proteins` = leading_proteins_AT8[,1])

totyrs_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_totyrs.rds")
totyrs_paths_p <- totyrs_fgsea_p[which(totyrs_fgsea_p$pval < 0.05),]
leading_proteins_totyrs <- as.matrix(sort(table(unlist(totyrs_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_totyrs <- data.frame(protein = rownames(leading_proteins_totyrs),`Total Years of Play Proteins` = leading_proteins_totyrs[,1])

#rna-seq
AT8_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_AT8_rnaseq.rds")
AT8_paths_r <- AT8_fgsea_r[which(AT8_fgsea_r$pval < 0.05),]
leading_genes_AT8 <- as.matrix(sort(table(unlist(AT8_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_AT8 <- data.frame(protein = rownames(leading_genes_AT8),`AT8 Total Genes` = leading_genes_AT8[,1])

totyrs_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_totyrs_rnaseq.rds")
totyrs_paths_r <- totyrs_fgsea_r[which(totyrs_fgsea_r$pval < 0.05),]
leading_genes_totyrs <- as.matrix(sort(table(unlist(totyrs_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_totyrs <- data.frame(protein = rownames(leading_genes_totyrs),`Total Years of Play Genes` = leading_genes_totyrs[,1])

# Merge all four together (full outer joins)
leading_edge <- full_join(leading_proteins_AT8, leading_proteins_totyrs, by = "protein") %>%
  full_join(leading_genes_AT8, by = "protein") %>%
  full_join(leading_genes_totyrs, by = "protein")

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
  pivot_longer(cols = starts_with(c("AT8.", "Total.")),
               names_to = "model",
               values_to = "Appearances")

# Make grouped horizontal bar plot
ggplot(leading_edge_long, aes(x = Appearances,
                              y = protein,
                              fill = model)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(x = "Number of Appearances", y = "Protein/Gene", fill = "Model") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
#################################################################################

#cognitive tests
#proteomics
cds_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_cds.rds")
cds_paths_p <- cds_fgsea_p[which(cds_fgsea_p$pval < 0.05),]
leading_proteins_cds <- as.matrix(sort(table(unlist(cds_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_cds <- data.frame(protein = rownames(leading_proteins_cds),`CDS Total Proteins` = leading_proteins_cds[,1])

dem_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_dementia.rds")
dem_paths_p <- dem_fgsea_p[which(dem_fgsea_p$pval < 0.05),]
leading_proteins_dem <- as.matrix(sort(table(unlist(dem_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_dem <- data.frame(protein = rownames(leading_proteins_dem),`Dementia Proteins` = leading_proteins_dem[,1])

faq_fgsea_p <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_faq.rds")
faq_paths_p <- faq_fgsea_p[which(faq_fgsea_p$pval < 0.05),]
leading_proteins_faq <- as.matrix(sort(table(unlist(faq_paths_p$leadingEdge)), decreasing = TRUE))
leading_proteins_faq <- data.frame(protein = rownames(leading_proteins_faq),`FAQ Total Proteins` = leading_proteins_faq[,1])

#rna-seq
cds_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_cds_rnaseq.rds")
cds_paths_r <- cds_fgsea_r[which(cds_fgsea_r$pval < 0.05),]
leading_genes_cds <- as.matrix(sort(table(unlist(cds_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_cds <- data.frame(protein = rownames(leading_genes_cds),`CDS Total Genes` = leading_genes_cds[,1])

dem_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_dementia_rnaseq.rds")
dem_paths_r <- dem_fgsea_r[which(dem_fgsea_r$pval < 0.05),]
leading_genes_dem <- as.matrix(sort(table(unlist(dem_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_dem <- data.frame(protein = rownames(leading_genes_dem),`Dementia Genes` = leading_genes_dem[,1])

faq_fgsea_r <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_faq_rnaseq.rds")
faq_paths_r <- faq_fgsea_r[which(faq_fgsea_r$pval < 0.05),]
leading_genes_faq <- as.matrix(sort(table(unlist(faq_paths_r$leadingEdge)), decreasing = TRUE))
leading_genes_faq <- data.frame(protein = rownames(leading_genes_faq),`FAQ Total Genes`= leading_genes_faq[,1])

# Merge all four together (full outer joins)
leading_edge <- full_join(leading_proteins_cds, leading_proteins_dem, by = "protein") %>%
  full_join(leading_proteins_faq, by = "protein") %>%
  full_join(leading_genes_cds, by = "protein") %>%
  full_join(leading_genes_dem, by = "protein") %>%
  full_join(leading_genes_faq, by = "protein")

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
  pivot_longer(cols = starts_with(c("CDS.", "Dementia.","FAQ.")),
               names_to = "model",
               values_to = "Appearances")

leading_edge_long <- leading_edge_long %>%
  mutate(appearances_signed = ifelse(str_detect(model, "\\.Proteins$"),
                                     Appearances,
                                     -Appearances))

# Make grouped horizontal bar plot
ggplot(leading_edge_long, aes(x = protein,
                              y = appearances_signed,
                              fill = model)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(x = "Protein/Gene",
       y = "Number of Appearances",
       fill = "Model") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(labels = abs)

ggplot(leading_edge_long, aes(x = Appearances,
                    y = protein,
                    fill = model)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(x = "Number of Appearances", y = "Protein/Gene", fill = "Model") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

