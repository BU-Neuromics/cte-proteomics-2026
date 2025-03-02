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
library(EnhancedVolcano)
library(DESeq2)

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
#remove samples that do not have a AT8_total value
adat <- adat[,which(!is.na(colData(adat)$AT8_total))] #186
#adat <- adat[,-173]
dim(adat)
adat <- adat[,which(colData(adat)$PathAD != 1)]

#counts <- assays(adat)$counts
#rownames(counts) <- rowData(adat)$EntrezGeneSymbol
#colnames(counts) <- colData(adat)$SampleGroup

#Look at duplicates
protein_table <- table(rowData(adat)$EntrezGeneSymbol)
protein_counts <- as.data.frame(protein_table)
colnames(protein_counts) <- c("protein", "frequency")
filtered_counts <- subset(protein_counts, frequency > 1)
write.csv(filtered_counts,'/restricted/projectnb/cteseq/projects/somascan_analysis/duplicate_proteins.csv')

# Plot distribution
ggplot(filtered_counts, aes(x = frequency)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of Duplicate Counts (Excluding Single Occurrences)",
    x = "Number of Duplicates",
    y = "Frequency"
  ) +
  theme_minimal()

frequency_table <- table(protein_counts$frequency)
frequency_table