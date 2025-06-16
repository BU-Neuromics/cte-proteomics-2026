#boxplot spot checks to look at progression of CTE
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(SummarizedExperiment)
library(fgsea)
library(qusage)
library(gt)
library(tidyr)
library(viridis)
library(hrbrthemes)
library(clusterProfiler)
library(DOSE)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan/data/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
C2_genesets <- read.gmt("/restricted/projectnb/cteseq/projects/somascan/data/data/c2.all.v2024.1.Hs.symbols.gmt")
adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7313 (-283)
adat <- adat[,which(colData(adat)$SampleGroup != "K-39")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-122")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-84")]
#remove non-human samples
adat <- adat[which(rowData(adat)$Organism == "Human"),] #7301
#removes same rows
#adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7285
adat <- adat[which(rowData(adat)$EntrezGeneID != ""),] #7285
AT8 <- colData(adat)$AT8_total
AT8 <- AT8[which(AT8 > 0)]
AT8min <- min(AT8, na.rm=T)
colData(adat)$AT8_total <- log(colData(adat)$AT8_total + 0.002628708)
colData(adat)$hitsperyear <- colData(adat)$chii_g / colData(adat)$totyrs

#FGG spot check
protein_index <- which(rownames(adat) == "seq.4989.7")
#ACHE spot check
protein_index <- which(rownames(adat) == "seq.10980.11")
#protein_index <- which(rownames(adat) == "seq.15553.22") - not good

#UBE2C 
protein_index <- which(rownames(adat) == "seq.12556.7")
protein_index <- which(rownames(adat) == "seq.20142.42")

#MITD1
protein_index <- which(rownames(adat) == "seq.23266.62")

#CAPNS2
protein_index <- which(rownames(adat) == "seq.22063.20")

protein_counts <- assays(adat)$counts[protein_index,]
#grab CTE data
protein <- colData(adat)$CTEStage
#check lengths 
length(protein_counts) == length(protein)
counts_CTE <- data.frame(counts = protein_counts, CTE = protein)
#write.csv(counts_CTE, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',protein_index[i,1],protein_index[i,2],'_CTE_followup_unlogged_file.csv'))
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE and CAPNS2')) +
  xlab("CTE")
CTE_box
for(i in 1:nrow(protein_index)){
  index <- protein_index[i,3]
  rowData(adat)$EntrezGeneSymbol[index] == protein_index[i,1]
  protein_counts <- assays(adat_CTE)$counts[index,]
  #grab CTE data
  protein_CTE <- colData(adat_CTE)$CTE
  #check lengths 
  length(protein_counts) == length(protein_CTE)
  counts_CTE <- data.frame(counts = protein_counts, CTE = protein_CTE)
  write.csv(counts_CTE, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',protein_index[i,1],protein_index[i,2],'_CTE_followup_unlogged_file.csv'))
  #create boxplots
  CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste0('Boxplot of CTE and ',protein_index[i,1])) +
    xlab("CTE")
  CTE_box
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',protein_index[i,1],protein_index[i,2],'_CTE_followup_boxplot.png'), plot = CTE_box, width = 8, height = 6)
}