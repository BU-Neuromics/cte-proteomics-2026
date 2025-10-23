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


#fgsea scatterplots
#models
fgsea_rl <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_CTE_RHIvslow.csv")
rl_sig <- fgsea_rhi3[which(fgsea_rhi3$padj < 0.05),]
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

#just cte stages
combined <- rbindlist(list(rhi1_sig,rhi2_sig,rhi3_sig,rhi4_sig), use.names = TRUE, fill = TRUE)
write.csv(combined, "/restricted/projectnb/cteseq/projects/somascan/final_plots/combined_fgsea_CTEStages.csv")

#all models
combined <- rbindlist(list(rl_sig,rh_sig,lh_sig,AT8_sig,totyrs_sig,cds_sig,dem_sig,faq_sig,rhi1_sig,rhi2_sig,rhi3_sig,rhi4_sig), use.names = TRUE, fill = TRUE)
write.csv(combined, "/restricted/projectnb/cteseq/projects/somascan/final_plots/combined_fgsea_allmodels.csv")
#CTE models
rhivslow_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "rhivslow")
rhivslow_groups <- full_join(rhivslow_groups, fgsea_rl_csv, by = "pathway")
lowvshigh_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "lowvshigh")
lowvshigh_groups <- full_join(lowvshigh_groups, fgsea_lh_csv, by = "pathway")
#create plot
merged_df <- full_join(
  rhivslow_groups %>% dplyr::select(pathway, NES.y, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_rl = NES.y, padj_rl = padj.y, leadingEdge_rl = leadingEdge.y, assigned_group_rl = assigned_group),
  lowvshigh_groups %>% dplyr::select(pathway, NES.y, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_lh = NES.y, padj_lh = padj.y, leadingEdge_lh = leadingEdge.y, assigned_group_lh = assigned_group),
  by = "pathway"
)

merged_df <- merged_df[which(merged_df$padj_rl < 0.05 | merged_df$padj_lh < 0.05),]
merged_df$assigned_group_lh == merged_df$assigned_group_rl
merged_df$group <- ifelse(is.na(merged_df$assigned_group_rl), merged_df$assigned_group_lh, merged_df$assigned_group_rl)
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
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "NES (RHI to Low CTE)",
       y = "NES (Low CTE to High CTE)",
       title = "Comparison of NES values between RHI to Low CTE and Low CTE to High CTE")

#AT8 and totyrs
AT8_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "AT8")
AT8_groups <- full_join(AT8_groups, fgsea_AT8_csv, by = "pathway")
totyrs_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "totyrs")
totyrs_groups <- full_join(totyrs_groups, fgsea_totyrs_csv, by = "pathway")
#create plot
merged_df <- full_join(
  AT8_groups %>% dplyr::select(pathway,NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_AT8 = NES, padj_AT8 = padj.y, leadingEdge_AT8 = leadingEdge.y, assigned_group_AT8 = assigned_group),
  totyrs_groups %>% dplyr::select(pathway, NES.y, padj.y, leadingEdge.y, group_assigned) %>% dplyr::rename(NES_totyrs = NES.y, padj_totyrs = padj.y, leadingEdge_totyrs = leadingEdge.y, assigned_group_totyrs = group_assigned),
  by = "pathway"
)

merged_df <- merged_df[which(merged_df$padj_AT8 < 0.05 | merged_df$padj_totyrs < 0.05),]
merged_df$assigned_group_totyrs == merged_df$assigned_group_AT8
merged_df$group <- ifelse(is.na(merged_df$assigned_group_AT8), merged_df$assigned_group_totyrs, merged_df$assigned_group_AT8)
# Define significance status
merged_df <- merged_df %>%
  mutate(sig_status = case_when(
    padj_AT8 < 0.05 & padj_totyrs < 0.05 ~ "Both",
    padj_AT8 < 0.05 & padj_totyrs >= 0.05 ~ "AT8 Total only",
    padj_AT8 >= 0.05 & padj_totyrs < 0.05 ~ "Total Years of Play only",
    TRUE ~ "Not significant"
  ))
merged_df$group[86] = "GTPase"
ggplot(merged_df, aes(x = NES_AT8, y = NES_totyrs, color = group, shape = sig_status, size = sig_status)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 5, "AT8 Total only" = 2, "Total Years of Play only" = 2, "Not significant" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "NES (AT8 Total)",
       y = "NES (Total Years of Play)",
       title = "Comparison of NES values between AT8 Total and Total Years of Play")

#Cognitive Tests
#CDS and dementia
cds_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "cds")
cds_groups <- full_join(cds_groups, fgsea_cds_csv, by = "pathway")
dem_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "dementia")
dem_groups <- full_join(dem_groups, fgsea_dem_csv, by = "pathway")
#create plot
merged_df <- full_join(
  cds_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_cds = NES, padj_cds = padj.y, leadingEdge_cds = leadingEdge.y, assigned_group_cds= assigned_group),
  dem_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_dem = NES, padj_dem = padj.y, leadingEdge_dem = leadingEdge.y, assigned_group_dem = assigned_group),
  by = "pathway"
)
merged_df <- merged_df[which(merged_df$padj_cds < 0.05 | merged_df$padj_dem < 0.05),]
merged_df$assigned_group_dem == merged_df$assigned_group_cds
merged_df$group <- ifelse(is.na(merged_df$assigned_group_cds), merged_df$assigned_group_dem, merged_df$assigned_group_cds)
merged_df$group[c(25:27,55,76,143,147,222,300,306)] <- "GTPase"
# Define significance status
merged_df <- merged_df %>%
  mutate(sig_status = case_when(
    padj_cds < 0.05 & padj_dem < 0.05 ~ "Both",
    padj_cds < 0.05 & padj_dem >= 0.05 ~ "CDS only",
    padj_cds >= 0.05 & padj_dem < 0.05 ~ "Dementia only",
    TRUE ~ "Not significant"
  ))

ggplot(merged_df, aes(x = NES_cds, y = NES_dem, color = group, shape = sig_status, size = sig_status)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 5, "CDS only" = 2, "Dementia only" = 2, "Not significant" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "NES (CDS)",
       y = "NES (Dementia)",
       title = "Comparison of NES values between Cognitive Difficulty Score (CDS) and Dementia")

#CDS and FAQ
cds_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "cds")
cds_groups <- full_join(cds_groups, fgsea_cds_csv, by = "pathway")
faq_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "faq")
faq_groups <- full_join(faq_groups, fgsea_faq_csv, by = "pathway")
#create plot
merged_df <- full_join(
  cds_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_cds = NES, padj_cds = padj.y, leadingEdge_cds = leadingEdge.y, assigned_group_cds= assigned_group),
  faq_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_faq = NES, padj_faq = padj.y, leadingEdge_faq = leadingEdge.y, assigned_group_faq = assigned_group),
  by = "pathway"
)
merged_df <- merged_df[which(merged_df$padj_cds < 0.05 | merged_df$padj_faq < 0.05),]
merged_df$assigned_group_faq == merged_df$assigned_group_cds
merged_df$group <- ifelse(is.na(merged_df$assigned_group_cds), merged_df$assigned_group_faq, merged_df$assigned_group_cds)
merged_df$group[c(25:27,51,75)] <- "GTPase"
# Define significance status
merged_df <- merged_df %>%
  mutate(sig_status = case_when(
    padj_cds < 0.05 & padj_faq < 0.05 ~ "Both",
    padj_cds < 0.05 & padj_faq >= 0.05 ~ "CDS only",
    padj_cds >= 0.05 & padj_faq < 0.05 ~ "FAQ only",
    TRUE ~ "Not significant"
  ))

ggplot(merged_df, aes(x = NES_cds, y = NES_faq, color = group, shape = sig_status, size = sig_status)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 5, "CDS only" = 2, "FAQ only" = 2, "Not significant" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "NES (CDS)",
       y = "NES (FAQ)",
       title = "Comparison of NES values between Cognitive Difficulty Score (CDS) and Function Activities Questionnarie (FAQ)")

#Dementia and FAQ
dem_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "dementia")
dem_groups <- full_join(dem_groups, fgsea_dem_csv, by = "pathway")
faq_groups <- read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/groups/fgsea_grouped_all.xlsx", sheet = "faq")
faq_groups <- full_join(faq_groups, fgsea_faq_csv, by = "pathway")
#create plot
merged_df <- full_join(
  dem_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_dem = NES, padj_dem = padj.y, leadingEdge_dem = leadingEdge.y, assigned_group_dem = assigned_group),
  faq_groups %>% dplyr::select(pathway, NES, padj.y, leadingEdge.y, assigned_group) %>% dplyr::rename(NES_faq = NES, padj_faq = padj.y, leadingEdge_faq = leadingEdge.y, assigned_group_faq = assigned_group),
  by = "pathway"
)
merged_df <- merged_df[which(merged_df$padj_dem < 0.05 | merged_df$padj_faq < 0.05),]
merged_df$assigned_group_faq == merged_df$assigned_group_dem
merged_df$group <- ifelse(is.na(merged_df$assigned_group_dem), merged_df$assigned_group_faq, merged_df$assigned_group_dem)
merged_df$group[c(307:316,319)] <- "GTPase"

# Define significance status
merged_df <- merged_df %>%
  mutate(sig_status = case_when(
    padj_dem < 0.05 & padj_faq < 0.05 ~ "Both",
    padj_dem < 0.05 & padj_faq >= 0.05 ~ "CDS only",
    padj_dem >= 0.05 & padj_faq < 0.05 ~ "FAQ only",
    TRUE ~ "Not significant"
  ))

ggplot(merged_df, aes(x = NES_dem, y = NES_faq, color = group, shape = sig_status, size = sig_status)) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c("Both" = 5, "Dementia only" = 2, "FAQ only" = 2, "Not significant" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "NES (Dementia)",
       y = "NES (FAQ)",
       title = "Comparison of NES values between Dementia and Function Activities Questionnarie (FAQ)")


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