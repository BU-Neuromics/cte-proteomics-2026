#Helen Pennington
#Labadorf Rotation 
#Started February 7, 2025

#Spearman Correlation to Compare adat protein data to expression data 
#CTE Samples

#packages 
library(dplyr)
library(GO.db)
library(SomaDataIO)
library(SomaScan.db)
library(data.table)
library(ggplot2)
library(EnsDb.Hsapiens.v79)
library(AnnotationDbi)
library(SummarizedExperiment)
library(limma)
library(DESeq2)
library(tidyr)
library(purrr)
library(tidyverse)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan_analysis/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
dim(adat)
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan_analysis/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
#adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7313 (-283)
adat <- adat[,which(colData(adat)$SampleGroup != "K-39")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-122")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-84")]
#remove non-human samples
adat <- adat[which(rowData(adat)$Organism == "Human"),] #7301
#removes same rows
adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7285
adat <- adat[which(rowData(adat)$EntrezGeneID != ""),] #7285
dim(adat)

##################################################################################

#set unique column names that have significance
#colnames(adat_df) <- attributes(adat)$Col.Meta$EntrezGeneSymbol
#adat_df <- make.unqiue(colnames(adat_df))
#adat_df$SampleGroup <- adat$SampleGroup
rownames <- rna$V1
rna <- rna[,-1]
rna <- t(rna)
index <- as.integer(which(rownames(rna) %in% metadata$Core_ID))
rna_meta <- rna[index,]
rna_meta <- rna_meta[order(rownames(rna_meta)),]
rna_meta <- t(rna_meta)
rownames(rna_meta) <- rownames
rownames(metadata) <- metadata$Core_ID
metadata <- metadata[which(rownames(metadata) %in% colnames(rna_meta)),]
rownames(metadata) <- metadata$Core_ID
metadata <- metadata[order(rownames(metadata)),]
rownames(metadata) <- metadata$Core_ID

all(rownames(metadata) == colnames(rna_meta))
rna_meta <- round(rna_meta)

#read in expression data
ddsA <- DESeqDataSetFromMatrix(countData = rna_meta,
                               colData = metadata,
                               design= ~ 1)
vst_rna <- vst(ddsA, blind=FALSE)

head(vst_rna)

vst_rna <- assay(vst_rna)

#merge with adat
vst_rna <- as.data.frame(vst_rna)
vst_rna$gene_name <- rownames(vst_rna)
vst_rna$gene_name <- sub("\\..*", "", vst_rna$gene_name)  # Keep only the main ID
vst_rna <- vst_rna[which(!is.na(vst_rna$gene_name)),]
length(vst_rna$gene_name)
edb <- EnsDb.Hsapiens.v79
# Convert Gene IDs to Gene Names (5071 geneids did not have gene names - num=56462)
vst_rna$gene_name <- mapIds(edb, keys = vst_rna$gene_name, keytype = "GENEID", column = "GENENAME", multiVals = "first")
vst_rna <- vst_rna[which(!is.na(vst_rna$gene_name)),]
vst_rna$gene_id <- rownames(vst_rna)
head(vst_rna)

#convert adat data into data frame 
adat_df <- as.data.frame(assays(adat)$counts)
colnames(adat_df) <- colData(adat)$Core_ID
rownames(adat_df) <- rowData(adat)$SeqId
adat_df$SeqID <- rowData(adat)$SeqId
adat_df$geneName <- rowData(adat)$EntrezGeneSymbol

#inner join adat and rna files
adat_rna <- dplyr::inner_join(adat_df, vst_rna,by=c("geneName"="gene_name"), relationship = "many-to-many",suffix = c(".p", ".r"))

#create rowData file to incorporate into lm model
rowData <- data.frame(Core_ID = colData(adat)$Core_ID, Batch = colData(adat)$Batch, RIN = colData(adat)$RIN)

#combine all data needed for lm model using tidyverse
# Pivot to long format
adat_rna_long <- adat_rna %>%
  pivot_longer(-c(gene_id,geneName,SeqID), names_to = "Sample_Type", values_to = "Expression") %>%
  separate(Sample_Type, into = c("Sample", "Type"), sep = "\\.") %>%
  pivot_wider(names_from = Type, values_from = Expression)
head(adat_rna_long)
adat_rna_long <- adat_rna_long[which(is.na(adat_rna_long$`NA`)),]
adat_rna_long <- adat_rna_long[,-6]

# Merge with metadata
adat_rna_long <- adat_rna_long %>%
  left_join(rowData, by = c("Sample" = "Core_ID"))

# Run lm() separately for each Gene-Protein combination
final_results <- adat_rna_long %>%
  group_by(SeqID, gene_id) %>%  # Group by Gene and Protein
  nest() %>%
  mutate(
    model = map(data, ~ lm(p ~ r + Batch + RIN, data = .x)),  # Fit model
    summary = map(model, summary),
    coef_table = map(summary, ~ as.data.frame(coef(.x))),
    coef_table = map2(SeqID, coef_table, ~ mutate(.y, SeqID = .x, Variable = rownames(.y)))
  ) %>%
  ungroup() %>%
  unnest(coef_table)  # Flatten nested structure

final_results <- as_tibble(final_results)
final_results <- final_results[,-1]
final_results <- final_results %>% unnest(coef_table)

#prep file for histogram (beta 1 - slope)
hist_results <- as.data.frame(final_results[,c(1,5,8,9,10)])
hist_results <- hist_results[which(hist_results$Variable == "r"),]
head(hist_results)

hist(hist_results$Estimate,
     main="Range in Beta values (Slope) for Protein vs Expression CTE data",
     ylab="Beta 1 (Slope)",
     xlab="Expression and Protein Combinations",
     col="darkmagenta",
     freq=FALSE
)

#prep file for histogram (beta 0 - intercept)
hist_results <- as.data.frame(final_results[,c(1,5,8,9,10)])
hist_results <- hist_results[which(hist_results$Variable == "(Intercept)"),]
head(hist_results)

hist(hist_results$Estimate,
     main="Range in Beta values (Intercept) for Protein vs Expression CTE data",
     ylab="Beta 0 (Intercept)",
     xlab="Expression and Protein Combinations",
     col="darkmagenta",
     freq=FALSE
)