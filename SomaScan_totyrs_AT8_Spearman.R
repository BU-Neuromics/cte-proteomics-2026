#Test for ADAT workflow using SomaScan.db
#Helen Pennington
#Labadorf Lab
#Started February 2, 2024

#packages
library(dplyr)
library(GO.db)
library(SomaDataIO)
library(SomaScan.db)
library(data.table)
library(ggplot2)

#now using CTE data
# Read in the file
adat_bbid <- fread("/usr4/spclpgm/hpenn13/LabadorfRotation/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
adat <- SomaDataIO::read_adat("/restricted/projectnb/cteseq/data/somascan_2024/CTE Somascan/HMS-24-036_2024-08-09/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")
dim(adat)
adat <- dplyr::left_join(adat,adat_bbid,by=c("SampleGroup"="SampleGroup"))

#only keep samples (get rid of controls) - removes 18 rows
adat <- adat[adat$SampleType == "Sample", ]
dim(adat)

#get rid of non-human samples - removes 298 columns
num_meta_cols <- length(adat) - length(attributes(adat)$Col.Meta$Organism) 
human_samples <- attributes(adat)$Col.Meta$Organism == "Human"
adat <- adat[, c(1:num_meta_cols, (num_meta_cols + which(human_samples))), drop = FALSE]
dim(adat)
#read in metadata
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
head(metadata)

#merge adat with metadata
adat <- dplyr::left_join(adat,metadata,by=c("BBID"="SampleName"))
adat <- adat[which(!is.na(adat$BBID)),]
dim(adat)
#########################################################################
#########################################################################

#Total Years of Play

#Data Prep and Exploration
summary(adat$totyrs)
png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/totyrs_histogram_dist.png'))
hist(adat$totyrs,
     xlab = "Total Years of Play",
     main = "Years Distribution"
)
dev.off()
totyrs <- data.frame(sample = adat$SampleGroup, totyrs = adat$totyrs)
write.csv(totyrs, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/totyrs_spearman.csv'))
adat <- dplyr::filter(adat, !is.na(BBID))
adat_totyrs <- dplyr::filter(adat, !is.na(totyrs))
dim(adat_totyrs)
analytes <- SomaDataIO::getAnalytes(adat_totyrs)


#header <- attributes(adat)$Col.Meta
#proteins <- header$Target
#rownames(adat) <- proteins

##########################################################################
#Identifying Associations

log10_adat <- log10(adat_totyrs)

# Calculate correlations for each analyte
cors_df <- lapply(analytes, function(x) {
  results <- stats::cor.test(log10_adat$totyrs, log10_adat[[x]],
                             method = "spearman", exact = FALSE)
  results <- results[c("estimate", "p.value")]
  unlist(results)
}) %>% setNames(analytes)

# Bind results together into a single dataframe
cors_df <- dplyr::bind_rows(cors_df, .id = "analytes")
cors_df$padj <- p.adjust(cors_df$p.value, method = "BH")

#order by p-value
cors_df <- cors_df %>% arrange(p.value)
cors_df <- cors_df[order(abs(cors_df$estimate.rho), decreasing = TRUE),]

#create volcano plots
#padj
cors_df <- cors_df %>%
  mutate(protein_type = case_when(estimate.rho >= 0.3 & padj <= 0.05 ~ "up",
                                  estimate.rho <= -0.3 & padj <= 0.05 ~ "down",
                                  TRUE ~ "ns"))
cors_df %>%
  count(protein_type)
#create volcano plot
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
vol_plot <- cors_df %>%
  ggplot(aes(x = estimate.rho,
             y = -log10(padj), color = protein_type)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(-0.3, 0.3),linetype = "dashed") +
  scale_color_manual(values=c("#E69F00", "#999999", "#56B4E9")) +
  labs(title = "Spearman Correlation for Total Years of Play (adjusted P-value)",
       x = "Estimated Rho",
       y = "-log10(P-value)",
       colour = "Protein Level Change")

vol_plot
ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/totyrs_adjpvalue_spearman.png'), plot = vol_plot, width = 8, height = 6)
#p.value
cors_df <- cors_df %>%
  mutate(protein_type = case_when(estimate.rho >= 0.3 & p.value <= 0.05 ~ "up",
                                  estimate.rho <= -0.3 & p.value <= 0.05 ~ "down",
                                  TRUE ~ "ns"))
cors_df %>%
  count(protein_type)
#create volcano plot
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
vol_plot <- cors_df %>%
  ggplot(aes(x = estimate.rho,
             y = -log10(p.value),color = protein_type)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(-0.3, 0.3),linetype = "dashed") +
  scale_color_manual(values=c("#E69F00", "#999999", "#56B4E9")) +
  labs(title = "Spearman Correlation for Total Years of Play (nominal P-value)",
       x = "Estimated Rho",
       y = "-log10(P-value)",
       colour = "Protein Level Change")

vol_plot
ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/totyrs_pvalue_spearman.png'), plot = vol_plot, width = 8, height = 6)
top_posCor <- cors_df %>%
  dplyr::filter(p.value < 0.05) %>% # Retain significant cors only
  dplyr::filter(estimate.rho >= abs(0.19)) %>% # Retain strong correlations
  arrange(desc(estimate.rho))

top_posCor
write.csv(cors_df, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/results_file_totyrs_spearman.csv'))

###########################################################################
#annotating results
#all available annotations
columns(SomaScan.db)

#grab top 10 analytes associated with AT8
top10_analytes <- head(top_posCor$analytes, 10L)

#retrive gene names and symbols for top 10 analytes
anno <- SomaScan.db::select(SomaScan.db, 
                            keys = top10_analytes,
                            columns = c("SYMBOL", "GENENAME", "GENETYPE"))
#perform gene ontology
go_anno <- SomaScan.db::select(SomaScan.db, 
                               keys = anno$PROBEID[1:3L],
                               columns = c("SYMBOL", "GENENAME", "GENETYPE", 
                                           "GO", "ONTOLOGY")) %>%
  dplyr::filter(ONTOLOGY == "BP")
go_terms <- AnnotationDbi::select(GO.db, 
                                  keys = go_anno$GO, 
                                  columns = c("GOID", "TERM", "DEFINITION"))
final_totyrs_df <- left_join(go_anno, go_terms, by = c("GO" = "GOID"))

final_totyrs_df
write.csv(final_totyrs_df, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/results_file_totyrs_GO.csv'))
#########################################################################
#########################################################################

# Total AT8 (tau)
summary(adat$AT8_total)
png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/AT8_total_histogram_dist.png'))
hist(adat$AT8_total,
     xlab = "Total Tau",
     main = "Distribution"
)
dev.off()
adat <- dplyr::filter(adat, !is.na(BBID))
adat_AT8_total <- dplyr::filter(adat, !is.na(AT8_total))
dim(adat_AT8_total)
AT8_total <- data.frame(sample = adat_AT8_total$SampleGroup, AT8_total = adat_AT8_total$AT8_total)
write.csv(AT8_total, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/AT8_total_spearman.csv'))
analytes <- SomaDataIO::getAnalytes(adat_AT8_total)


#header <- attributes(adat)$Col.Meta
#proteins <- header$Target
#rownames(adat) <- proteins

##########################################################################
#Identifying Associations

log10_adat <- log10(adat_AT8_total)

# Calculate correlations for each analyte
cors_df <- lapply(analytes, function(x) {
  results <- stats::cor.test(log10_adat$AT8_total, log10_adat[[x]],
                             method = "spearman", exact = FALSE)
  results <- results[c("estimate", "p.value")]
  unlist(results)
}) %>% setNames(analytes)

# Bind results together into a single dataframe
cors_df <- dplyr::bind_rows(cors_df, .id = "analytes")
cors_df$padj <- p.adjust(cors_df$p.value, method = "BH")

#order by p-value
#cors_df <- cors_df %>% arrange(p.value)
#cors_df <- cors_df[order(abs(cors_df$estimate.rho), decreasing = TRUE),]

#look at abs rho > 0.3
cors_df <- cors_df %>%
  mutate(protein_type = case_when(estimate.rho >= 0.3 & padj <= 0.05 ~ "up",
                               estimate.rho <= -0.3 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns"))
cors_df %>%
  count(protein_type)
#create volcano plot
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

vol_plot <- cors_df %>%
  ggplot(aes(x = estimate.rho,
             y = -log10(padj), color = protein_type)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.3, 0.3),
             linetype = "dashed") +
  scale_color_manual(values=c("#E69F00", "#999999", "#56B4E9")) +
  labs(title = "Spearman Correlation for Total AT8 (adjusted P-value, abs(rho) > 0.3)",
       x = "Estimated Rho",
       y = "-log10(P-value)",
       colour = "Protein Level Change")

vol_plot
ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/AT8_total_adjpvalue_spearman.png'), plot = vol_plot, width = 8, height = 6)
#look at rho abs > 0.5
cors_df <- cors_df %>%
  mutate(protein_type = case_when(estimate.rho >= 0.5 & padj <= 0.05 ~ "up",
                                  estimate.rho <= -0.5 & padj <= 0.05 ~ "down",
                                  TRUE ~ "ns"))
cors_df %>%
  count(protein_type)
#create volcano plot
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

vol_plot <- cors_df %>%
  ggplot(aes(x = estimate.rho,
             y = -log10(padj), color = protein_type)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed") +
  scale_color_manual(values=c("#E69F00", "#999999", "#56B4E9")) +
  labs(title = "Spearman Correlation for Total AT8 (adjusted P-value, abs(rho)>0.5",
       x = "Estimated Rho",
       y = "-log10(P-value)",
       colour = "Protein Level Change")

vol_plot

# Isolate top positive correlations
top_posCor <- cors_df %>%
  dplyr::filter(padj < 0.05) %>% # Retain significant cors only
  dplyr::filter(estimate.rho >= abs(0.3)) %>% # Retain strong correlations
  arrange(desc(estimate.rho))

top_posCor
write.csv(cors_df, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/results_file_AT8_total_spearman.csv'))
########################################################################
#annotating results
#all available annotations
columns(SomaScan.db)

#grab top 10 analytes associated with AT8
top10_analytes <- head(top_posCor$analytes, 10L)

#retrive gene names and symbols for top 10 analytes
anno <- SomaScan.db::select(SomaScan.db, 
                            keys = top10_analytes,
                            columns = c("SYMBOL", "GENENAME", "GENETYPE"))
#perform gene ontology
go_anno <- SomaScan.db::select(SomaScan.db, 
                               keys = anno$PROBEID[1:3L],
                               columns = c("SYMBOL", "GENENAME", "GENETYPE", 
                                           "GO", "ONTOLOGY")) %>%
  dplyr::filter(ONTOLOGY == "BP")
go_terms <- AnnotationDbi::select(GO.db, 
                                  keys = go_anno$GO, 
                                  columns = c("GOID", "TERM", "DEFINITION"))
final_AT8_df <- left_join(go_anno, go_terms, by = c("GO" = "GOID"))

final_AT8_df
write.csv(final_AT8_df, paste0('/restricted/projectnb/cteseq/projects/somascan/results/spearman/results_file_AT8_total_GO.csv'))
