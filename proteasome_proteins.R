#protein clusters
#Helen Pennington
#Started June 5, 2025

#packages
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(SummarizedExperiment)
library(pheatmap)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan/data/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
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

lowvshigh <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_Group_de_gseaResults.csv")
RHIvslow <- fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_Group_de_RHIl_gseaResults.csv")


#proteins <- lowvshigh$leadingEdge[1:40]
#RHIvslow$leadingEdge[1:11]
#which(RHIvslow$leadingEdge[1:11] %in% lowvshigh$leadingEdge[1:40])
proteins <- c("PSMB1","PSMB2","PSMA7","PSMB5","PSMA2","PSMB4","PSMA4","PSMB3","PSMB6","PSMA1","UBB","UBC")
proteins <- c("PSMB1","PSMB2","PSMA7","PSMB5","PSMA2","PSMB4")
proteins <- c("CD247","ACHE","MITD1","COPS5","SFRP1","MARS1","MED11","MEGF10")
proteins <- c("HDHD2","EMC2","CAPNS2","RAD51C","GNAQ","SBDS","OAF","HSPB1","MAP3K3","WNT3A","EMC4","PGM2","HCE001796","PCYT1A")
index <- rep(NA, length(proteins))
i=1
for (i in 1:length(proteins)) {
  print(i)
  index[i] <- which(rowData(adat)$EntrezGeneSymbol == proteins[i])
  print(index[i])
}
print(index)
adat_filt <- adat[index,]
#adat_filt <- colData(adat_filt)$CTEStage
adat_df <- as.matrix(assays(adat_filt)$counts)
rownames(adat_df) <- rowData(adat_filt)$EntrezGeneSymbol
anno <- data.frame(Stage=colData(adat_filt)$CTEStage, names = colnames(adat_filt))
rownames(anno) <- anno$names
pheatmap(adat_df, annotation_col = anno)
adat_df$proteins <- rowData(adat_filt)$EntrezGeneSymbol

adat_dfm <- melt(adat_df)
head(adat_dfm)
colnames(adat_dfm) <- c("y", "x", "value")
ggplot(adat_dfm, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") + # Customize colors
  labs(x = "Column", y = "Row", fill = "Value") + # Add labels
  coord_fixed() # Ensure square tiles


