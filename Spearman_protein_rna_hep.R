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

# Read in the adat file and BBID file
adat_bbid <- fread("/usr4/spclpgm/hpenn13/LabadorfRotation/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
adat <- SomaDataIO::read_adat("/restricted/projectnb/cteseq/data/somascan_2024/CTE Somascan/HMS-24-036_2024-08-09/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")
dim(adat)

#get rid of non-human samples - removes 298 columns
num_meta_cols <- length(adat) - length(attributes(adat)$Col.Meta$Organism) 
human_samples <- attributes(adat)$Col.Meta$Organism == "Human"
adat <- adat[, c(1:num_meta_cols, (num_meta_cols + which(human_samples))), drop = FALSE]
adat <- log10(adat)
dim(adat)

#convert to data frame 
adat_df <- as.data.frame(adat[1:212,31:7331])
#adat_df[nrow(adat_df) + 1,] <- attributes(adat)$Col.Meta$EntrezGeneSymbol
colnames(adat_df) <- attributes(adat)$Col.Meta$EntrezGeneSymbol
adat_df <- make.unqiue(colnames(adat_df))
adat_df$SampleGroup <- adat$SampleGroup

#merge files together to get BBID for each sample
adat_df <- dplyr::left_join(adat_df,adat_bbid,by=c("SampleGroup"="SampleGroup"))
adat_df <- adat_df[which(!is.na(adat_df$BBID)),]
dim(adat_df)

#read in metadata
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
head(metadata)
dim(metadata)

#merge adat with metadata
adat <- dplyr::left_join(adat,metadata,by=c("BBID"="SampleName"))
rownames(adat_df) <- adat_df$Core_ID
dim(adat_df)
print(adat_df[1:5,7280:7301])

#read in expression data
rna <- as.data.frame(fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv"))
rna$V1 <- sub("\\..*", "", rna$V1)  # Keep only the main ID
dim(rna)
rna <- rna[which(!is.na(rna$V1)),]
dim(rna)
edb <- EnsDb.Hsapiens.v79
# Convert Gene IDs to Gene Names
rna$V1 <- mapIds(edb, keys = rna$V1, keytype = "GENEID", column = "GENENAME", multiVals = "first")
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rna$V1, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
dim(geneIDs1)
dim(rna)
rna <- merge(rna,geneIDs1,by.x="V1",by.y="GENEID", all=TRUE)
head(rna)

#get transpose of rna data for easier usage
rownames(rna) <- rna$V1
rna <- rna[,-1]
trna <- transpose(rna)
#trna <- as.data.frame(rna)
colnames(trna)<-rownames(rna)
trna$sample_id <- colnames(rna)
print(trna[1:5,61530:61534])

#merge with protein data
#adat_rna <- dplyr::left_join(adat,trna,by=c("Core_ID"="sample_id"))

#Spearman correlation code
#Calculate correlations for each analyte
common_genes <- intersect(rownames(adat_df), rownames(trna))

# Subset data to common genes
adat_df <- adat_df[common_genes, ]
trna <- trna[common_genes, ]

# Function to calculate Spearman correlation for each gene
compute_spearman <- function(gene) {
  cor.test(adat_df[gene, ], trna[gene, ], method = "spearman")$estimate
}

# Apply function to each gene
spearman_results <- sapply(common_genes, compute_spearman)

# Convert to data frame for easier visualization
spearman_df <- data.frame(Gene = common_genes, Spearman_Correlation = spearman_results)

# View results
head(spearman_df)
