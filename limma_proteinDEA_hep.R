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

#differential expression analysis (Limma)
#create design matrix and run limma
genes <- data.frame(unique_names = rownames(rowData(adat)), genes = rowData(adat)$EntrezGeneSymbol)
head(genes)
AT8_model <- model.matrix(~ agedeath + AT8_total, data=colData(adat))
fit_AT8 <- limma::lmFit(assays(adat)$counts, AT8_model)
elimma_AT8 <- eBayes(fit_AT8)
#create results files
AT8_limma_res <- topTable(elimma_AT8,number = 100, adjust="BH", coef = "AT8_total",sort.by="P")
AT8_limma_res$unique_names <- rownames(AT8_limma_res)
AT8_limma_res <- merge(AT8_limma_res, genes, by.x = "unique_names", by.y="unique_names")

head(AT8_limma_res)

#volcano plot
fit_AT8$unique_names <- rownames(fit_AT8)
fit_AT8 <- merge(fit_AT8, genes, by.x = "unique_names", by.y="unique_names")
volcanoplot(elimma_AT8,coef="AT8_total",highlight=20,names=fit_AT8$genes, col="red",main="AT8_total") 

#spot check significant protein C1QTNF5 = seq.7810.20
protein_index <- which(rownames(adat) == "seq.7810.20")
protein_index
#check to verify that the row number corresponds to the correct gene
rowData(adat)$EntrezGeneSymbol[protein_index]
#grab counts for C1QTNF5
protein_counts <- assays(adat)$counts[protein_index,]

#grab AT8 total for C1QTNF5
protein_AT8 <- colData(adat)$AT8_total

#check lengths 
length(protein_counts) == length(protein_AT8)

C1QTNF5_AT8 <- data.frame(counts = protein_counts, AT8 = protein_AT8)

#create scatter plot
ggplot(C1QTNF5_AT8, aes(x=AT8, y=counts)) + 
  geom_point()
#######################

#spot check significant protein PDCD6IP = seq.18174.79
protein_index <- which(rownames(adat) == "seq.18174.79")
protein_index
#check to verify that the row number corresponds to the correct gene
rowData(adat)$EntrezGeneSymbol[protein_index]
#grab counts for C1QTNF5
protein_counts <- assays(adat)$counts[protein_index,]

#grab AT8 total for C1QTNF5
protein_AT8 <- colData(adat)$AT8_total

#check lengths 
length(protein_counts) == length(protein_AT8)

PDCD6IP_AT8 <- data.frame(counts = protein_counts, AT8 = protein_AT8)

#create scatter plot
ggplot(PDCD6IP_AT8, aes(x=AT8, y=counts)) + 
  geom_point()

#####################
#Dkk-4
protein_index <- which(rowData(adat)$EntrezGeneSymbol == "DKK4")
protein_index
#grab counts for C1QTNF5
protein_counts <- assays(adat)$counts[c(protein_index),]

#grab AT8 total for C1QTNF5
protein_AT8 <- colData(adat)$AT8_total

#check lengths 
length(protein_counts) == 2*length(protein_AT8)

DKK4_AT8 <- data.frame(counts1 = protein_counts[1,],counts2 = protein_counts[2,], AT8 = protein_AT8)

#create scatter plot
ggplot(DKK4_AT8, aes(x=AT8, y=counts1)) + 
  geom_point()
ggplot(DKK4_AT8, aes(x=AT8, y=counts2)) + 
  geom_point()