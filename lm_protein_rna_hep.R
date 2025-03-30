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
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
dim(adat)
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan/data/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
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
dim(vst_rna) #56462
length(unique(vst_rna$gene_name)) #54580
#convert adat data into data frame 
dim(adat) #7285
adat_df <- as.data.frame(assays(adat)$counts)
colnames(adat_df) <- colData(adat)$Core_ID
rownames(adat_df) <- rowData(adat)$SeqId
adat_df$SeqID <- rowData(adat)$SeqId
adat_df$geneName <- rowData(adat)$EntrezGeneSymbol
adat_df <- adat_df %>%
  separate_rows(geneName, sep = "\\|")
dim(adat_df) #7393
length(unique(adat_df$geneName)) #6417
length(unique(adat_df$SeqID)) #7285
#inner join adat and rna files
adat_rna <- dplyr::inner_join(adat_df, vst_rna,by=c("geneName"="gene_name"), relationship = "many-to-many",suffix = c(".p", ".r"))
dim(adat_rna) #6926 -> 7107
length(unique(adat_rna$SeqID)) #6988
length(unique(adat_rna$gene_id)) #6140
length(unique(adat_rna$geneName)) #6119
adat_df$geneName[c(which(!adat_df$geneName %in% adat_rna$geneName))]
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
  ungroup()

final_results <- as_tibble(final_results)
final_results <- final_results[,-1]
final_results <- final_results %>% unnest(coef_table)
final_results_df <- as.data.frame(final_results[,-c(2:4)])
write.csv(final_results_df, paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/RNA_protein_comparisons.csv'))

#prep file for histogram (beta 1 - slope)
hist_results1 <- as.data.frame(final_results[,c(1,5,8,9,10)])
hist_results1 <- hist_results1[which(hist_results1$Variable == "r"),]
head(hist_results1)
hist_results1 <- arrange(hist_results1, `Pr(>|t|)`)
hist_results1 <- dplyr::mutate(hist_results1,
                              padj=p.adjust(`Pr(>|t|)`, method = "fdr"))
write.csv(hist_results1, paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/RNA_protein_Beta1RNA_comparisons.csv'))
pivot_final_results <- pivot_wider(final_results,
            id_cols = c("gene_id","SeqID"),
            names_from = c("Variable"),
            values_from = c("Estimate"))
plot(pivot_final_results$`(Intercept)`,pivot_final_results$r)

png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/lm_beta1RNAhistogram_RNAprotein.png'))
hist(hist_results1$Estimate,
     main="Range in Beta values (Slope) for Protein vs Expression CTE data",
     ylab="Beta 1 (Slope)",
     xlab="Expression and Protein Combinations",
     col="darkmagenta",
     freq=FALSE,
     breaks=100
)
dev.off()

#prep file for histogram (beta 0 - intercept)
hist_results0 <- as.data.frame(final_results[,c(1,5,8,9,10)])
hist_results0 <- hist_results0[which(hist_results0$Variable == "(Intercept)"),]
head(hist_results0)
write.csv(hist_results0, paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/RNA_protein_Beta0_comparisons.csv'))

png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/lm_beta0histogram_RNAprotein.png'))
hist(hist_results0$Estimate,
     main="Range in Beta values (Intercept) for Protein vs Expression CTE data",
     ylab="Beta 0 (Intercept)",
     xlab="Expression and Protein Combinations",
     col="darkmagenta",
     freq=FALSE
)
dev.off()
################################################################################
#Follow-up
vst_rna$median <- apply(vst_rna, 1, median, na.rm=T)
rna_median <- data.frame(gene_id = rownames(vst_rna), median = vst_rna$median)
results_rnaMedian_merge <- merge(rna_median,hist_results1, by = "gene_id")
png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/MedianScatterplot_followUp_RNAbeta1.png'))
plot(results_rnaMedian_merge$Estimate,results_rnaMedian_merge$median)
dev.off()
#which beta r are significant (FDR < 0.05)
sig_hist_results <- hist_results1[which(hist_results1$`Pr(>|t|)` < 0.05),]
nrow(sig_hist_results)

#scatterplots
Scatterplots <- function(index_r, index_p, vst_rna, adat_df){
  rna <- as.data.frame(t(vst_rna[index_r,]))
  protein <- as.data.frame(t(adat_df[index_p,]))
  name_r <- colnames(rna)
  name_p <- colnames(protein)
  rna$samples <- rownames(rna)
  protein$samples <- rownames(protein)
  scatterplot_rp <- merge(rna, protein, by = "samples")
  colnames(scatterplot_rp) <- c("samples","rna","protein")
  scatter <- ggplot(scatterplot_rp, aes(x = as.numeric(rna), y = as.numeric(protein))) +
    geom_point() +
    labs(title = paste(name_r,' and ',name_p,' RNA/Protein Comparison'), x = "RNA", y = "Protein") +
    geom_smooth(method = "lm")
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/lm/proteinRNA_mostsig_',name_r,'_',name_p,'_scatterplot.png'), plot = scatter, width = 8, height = 6)
}
#smallest p value
index_r <- which(rownames(vst_rna) == "ENSG00000134184.13")
index_p <- which(rownames(adat_df) == "15395-15")
Scatterplots(index_r, index_p, vst_rna, adat_df)

which(rownames(assays(adat)$counts) == "seq.15395.15")
separation <- as.data.frame(colData(adat))
separation <- t(separation)
counts <- assays(adat)$counts[1663,]
counts <- t(as.data.frame(counts))
plot(separation$V1,separation$Optic.nerve)

separation <- merge(counts,separation, by = "row.names")
dim(separation)

index_r <- which(rownames(vst_rna) == "ENSG00000164308.17")
index_p <- which(rownames(adat_df) == "8960-3")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000134184.13")
index_p <- which(rownames(adat_df) == "7239-9")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000197705.10")
index_p <- which(rownames(adat_df) == "24903-7")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000133433.11")
index_p <- which(rownames(adat_df) == "11273-176")
Scatterplots(index_r, index_p, vst_rna, adat_df)

#tau
index_r <- which(vst_rna$gene_name == "MAPT")
index_p <- which(adat_df$geneName == "MAPT")
Scatterplots(index_r, index_p, vst_rna, adat_df)

#APOE
index_r <- which(rownames(vst_rna) == "ENSG00000130203.10")
index_p <- which(adat_df$SeqID == "2418-55")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000130203.10")
index_p <- which(rownames(adat_df) == "2937-10")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000130203.10")
index_p <- which(rownames(adat_df) == "2938-55")
Scatterplots(index_r, index_p, vst_rna, adat_df)

index_r <- which(rownames(vst_rna) == "ENSG00000130203.10")
index_p <- which(rownames(adat_df) == "5312-49")
Scatterplots(index_r, index_p, vst_rna, adat_df)

dim(hist_results)
length(unique(hist_results$gene_id))
length(unique(hist_results$SeqID))
head(vst_rna)
which(vst_rna$gene_name == "MAPT")
which(hist_results$gene_id == "ENSG00000130203.10")

duplicates_geneid <- subset(as.data.frame(table(hist_results$gene_id)),Freq > 1)
duplicates_seqid <- subset(as.data.frame(table(hist_results$SeqID)),Freq >1)
nrow(duplicates_geneid)
nrow(duplicates_seqid)

adatRNA_dups <- data.frame(SeqID = adat_rna$SeqID,gene_id = adat_rna$gene_id,geneName = adat_rna$geneName)
adatRNA_dups <- adatRNA_dups[order(adatRNA_dups$geneName), ] 
adatRNA_dups <- adatRNA_dups[which(duplicated(adatRNA_dups$geneName) | duplicated(adatRNA_dups$geneName, fromLast = TRUE)),]
adatRNA_dups3 <- adatRNA_dups %>%
  dplyr::group_by(geneName) %>%
  dplyr::filter(n() > 2) %>%
  dplyr::ungroup()
adatRNA_dups3 <- adatRNA_dups3[which(duplicated(adatRNA_dups3$gene_id) | duplicated(adatRNA_dups3$gene_id, fromLast = TRUE)),]
adatRNA_dups3 <- adatRNA_dups3[which(duplicated(adatRNA_dups3$SeqID) | duplicated(adatRNA_dups3$SeqID, fromLast = TRUE)),]
print(adatRNA_dups3, n=75)
write.csv(adatRNA_dups3, '/restricted/projectnb/cteseq/projects/somascan_analysis/RNAvsProteinAnalysis/RNA_Protein_DuplicateGroups.csv')
