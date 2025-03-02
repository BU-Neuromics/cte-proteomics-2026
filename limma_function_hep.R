#Helen Pennington
#Labadorf Rotation 
#Started March 1, 2025

#Spearman Correlation to Compare adat protein data to expression data 
#CTE Samples

#packages 
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(SummarizedExperiment)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan_analysis/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan_analysis/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7313 (-283)
adat <- adat[,which(colData(adat)$SampleGroup != "K-39")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-122")]
adat <- adat[,which(colData(adat)$SampleGroup != "K-84")]
#remove non-human samples
adat <- adat[which(rowData(adat)$Organism == "Human"),] #7301
#removes same rows
#adat <- adat[-grep("Internal Use Only", rowData(adat)$TargetFullName), ] #7285
adat <- adat[which(rowData(adat)$EntrezGeneID != ""),] #7285

#differential expression analysis (Limma)
#create design matrix and run limma
limma_model <- function(model,adat, name, coefficient) {
  fit <- limma::lmFit(assays(adat)$counts, model)
  elimma <- eBayes(fit)
  #create results files
  limma_res <- topTable(elimma,number = 100, adjust="BH", coef = coefficient,sort.by="P")
  limma_res$unique_names <- rownames(limma_res)
  limma_res <- merge(limma_res, genes, by.x = "unique_names", by.y="unique_names")
  
  #volcano plot
  fit$unique_names <- rownames(fit)
  fit <- merge(fit, genes, by.x = "unique_names", by.y="unique_names")
  png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/volcano_plots/limma_DEA_',name,'_volcanoplot.png'))
  volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  dev.off()
  write.csv(limma_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_',name,'_results.csv'))
  return(limma_res)
}

#~ agedeath + AT8_total
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8_PathAD <- adat_AT8[,which(colData(adat_AT8)$PathAD != 1)]
genes <- data.frame(unique_names = rownames(rowData(adat_AT8_PathAD)), genes = rowData(adat_AT8_PathAD)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + AT8_total, data=colData(adat_AT8_PathAD))
AT8_res <- limma_model(AT8_model,adat_AT8_PathAD, "AT8_total","AT8_total")

#~ age_at_death + totyrs
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))]
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + totyrs, data=colData(adat_totyrs))
totyrs_res <- limma_model(totyrs_model,adat_totyrs, "totyrs","totyrs")

#`~ age_at_death + PathAD + totyrs`
adat_totyrs_PathAD <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))]
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs_PathAD)), genes = rowData(adat_totyrs_PathAD)$EntrezGeneSymbol)
totyrs_PathAD_model <- model.matrix(~ agedeath + PathAD + totyrs, data=colData(adat_totyrs_PathAD))
totyrs_PathAD_res <- limma_model(totyrs_PathAD_model,adat_totyrs_PathAD, "totyrs + PathAD","totyrs")

#~ age_at_death + CTE
adat_CTE <- adat[,which(!is.na(colData(adat)$CTE))]
genes <- data.frame(unique_names = rownames(rowData(adat_CTE)), genes = rowData(adat_CTE)$EntrezGeneSymbol)
CTE_model <- model.matrix(~ agedeath + CTE, data=colData(adat_CTE))
CTE_res <- limma_model(CTE_model,adat_CTE, "CTE","CTE")

#~ age_at_death + CTEStage
adat_CTEStage <- adat[,which(!is.na(colData(adat)$CTEStage))]
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage)), genes = rowData(adat_CTEStage)$EntrezGeneSymbol)
CTEStage_model <- model.matrix(~ agedeath + CTEStage, data=colData(adat_CTEStage))
CTEStage_res <- limma_model(CTEStage_model,adat_CTEStage, "CTEStage","CTEStage")

#~ age_at_death + CTELowvsHigh
adat_CTEonly <- adat[,which(colData(adat)$CTE == 1)]
adat_CTEonly_lowvshigh <- adat_CTEonly[,which(!is.na(colData(adat_CTEonly)$Group_de))]
colData(adat_CTEonly_lowvshigh)$Group_de <- as.integer(ifelse(colData(adat_CTEonly_lowvshigh)$Group_de == "low", 0, ifelse(colData(adat_CTEonly_lowvshigh)$Group_de == "high", 1, colData(adat_CTEonly_lowvshigh)$Group_de))) 
genes <- data.frame(unique_names = rownames(rowData(adat_CTEonly_lowvshigh)), genes = rowData(adat_CTEonly_lowvshigh)$EntrezGeneSymbol)
CTEonly_lowvshigh_model <- model.matrix(~ agedeath + Group_de, data=colData(adat_CTEonly_lowvshigh))
CTEonly_lowvshigh_res <- limma_model(CTEonly_lowvshigh_model,adat_CTEonly_lowvshigh, "Group_de","Group_de")