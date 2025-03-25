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
library(fgsea)
library(qusage)
library(gt)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan_analysis/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan_analysis/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
C2_genesets <- read.gmt("/restricted/projectnb/cteseq/projects/somascan_analysis/data/c2.all.v2024.1.Hs.symbols.gmt")
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
limma_model <- function(model,adat, name, coefficient,genes) {
  fit <- limma::lmFit(assays(adat)$counts, model)
  elimma <- eBayes(fit)
  #create results files
  limma_res <- topTable(elimma,number = 100, adjust="BH", coef = coefficient,sort.by="P")
  limma_res$unique_names <- rownames(limma_res)
  limma_res <- merge(limma_res, genes, by.x = "unique_names", by.y="unique_names")
  limma_res_all <- topTable(elimma,number = 1000, adjust="BH", coef = coefficient,sort.by="P")
  limma_res_all$unique_names <- rownames(limma_res_all)
  limma_res_all<- merge(limma_res_all, genes, by.x = "unique_names", by.y="unique_names")
  
  #volcano plot
  fit$unique_names <- rownames(fit)
  fit <- merge(fit, genes, by.x = "unique_names", by.y="unique_names")
  png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/volcano_plots/limma_DEA_',name,'_volcanoplot.png'))
  volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  dev.off()
  volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  write.csv(limma_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_',name,'_results.csv'))
  return(limma_res_all)
}
fgsea_function <- function(limma_res_all,C2_genesets,name){
  #sort results file to 
  # Function: Adjacency matrix to list -------------------------
  matrix_to_list <- function(pws){
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  # Get all the genes in your dataset and assign them to my_genes 
  my_genes <- unique(limma_res_all$genes)
  hidden <- unique(unlist(C2_genesets))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(C2_genesets)),
                nrow = length(hidden), ncol = length(C2_genesets))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% C2_genesets[[i]])
  }
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(my_genes, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again using the function we previously defined
  final_list <- matrix_to_list(mat)
  genes_sorted <- limma_res_all[order(abs(limma_res_all$logFC), decreasing = TRUE), ]
  genes_sorted_unique <- genes_sorted[!duplicated(genes_sorted$genes), ]
  plot(genes_sorted_unique$logFC)
  stats <- genes_sorted_unique$logFC
  names(stats) <- genes_sorted_unique$genes
  fgseaRes <- fgsea::fgsea(pathways = final_list, 
                    stats    = stats,
                    minSize  = 15,
                    maxSize  = 500)
  head(fgseaRes[order(pval), ])
  #top 10 pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  #png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/gseaPlots/limma_DEA_',name,'_gseaPlot.png'))
  gseaTable <- plotGseaTable(final_list[topPathways], stats, fgseaRes, 
                gseaParam=0.5)
  #dev.off()
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/gseaPlots/limma_DEA_',name,'_gseaPlot.png'), plot = gseaTable, width = 8, height = 6)
  #final_results <- as_tibble(final_results)
  #final_results <- final_results[,-1]
  #final_results <- final_results %>% unnest(coef_table)
  #write.csv(fgseaRes, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_',name,'_gseaResults.csv'))
  return(fgseaRes)
}

significance_table <- matrix(data=NA,nrow = 7,ncol = 5)
significance_table <- as.data.frame(significance_table)
colnames(significance_table) <- c("Model","Nominal Significant Proteins","Adjusted Significant Proteins","Nominal Significant Pathways","Adjusted Significant Pathways")
significance_table$Model <- c("age_at_death + AT8_total","age_at_death + PathAD + AT8_total","age_at_death + totyrs","age_at_death + PathAD + totyrs","age_at_death + CTE","age_at_death + CTEStage","age_at_death + CTELowvsHigh")
significance_table

#~ agedeath + AT8_total
#model <- AT8_model
#adat <- adat_AT8_PathAD
#name <- "AT8_total"
#coefficient <- "AT8_total"
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8_PathAD <- adat_AT8[,which(colData(adat_AT8)$PathAD != 1)]
genes <- data.frame(unique_names = rownames(rowData(adat_AT8_PathAD)), genes = rowData(adat_AT8_PathAD)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + AT8_total, data=colData(adat_AT8_PathAD))
AT8_res_all <- limma_model(AT8_model,adat_AT8_PathAD, "AT8_total","AT8_total",genes)
significance_table[1,2] <- length(which(AT8_res_all$P.Value < 0.05))
significance_table[1,3] <- length(which(AT8_res_all$adj.P.Val < 0.05))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total")
write.csv(AT8_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_AT8_model_gseaResults.csv'))
significance_table[1,4] <- length(which(AT8_fgsea_res$pval < 0.05))
significance_table[1,5] <- length(which(AT8_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + AT8_total`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8_PathAD <- adat_AT8[,which(!is.na(colData(adat_AT8)$PathAD))]
genes <- data.frame(unique_names = rownames(rowData(adat_AT8_PathAD)), genes = rowData(adat_AT8_PathAD)$EntrezGeneSymbol)
AT8_PathAD_model <- model.matrix(~ agedeath + PathAD + AT8_total, data=colData(adat_AT8_PathAD))
AT8_PathAD_res_all <- limma_model(AT8_PathAD_model,adat_AT8_PathAD, "AT8_PathAD_total","AT8_total",genes)
significance_table[2,2] <- length(which(AT8_PathAD_res_all$P.Value < 0.05))
significance_table[2,3] <- length(which(AT8_PathAD_res_all$adj.P.Val < 0.05))
AT8_PathAD_fgsea_res <- fgsea_function(AT8_PathAD_res_all,C2_genesets,"AT8_PathAD_total")
write.csv(AT8_PathAD_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_AT8_PathAD_model_gseaResults.csv'))
significance_table[2,4] <- length(which(AT8_PathAD_fgsea_res$pval < 0.05))
significance_table[2,5] <- length(which(AT8_PathAD_fgsea_res$padj < 0.05))

#~ age_at_death + totyrs
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))]
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + totyrs, data=colData(adat_totyrs))
totyrs_res_all <- limma_model(totyrs_model,adat_totyrs, "totyrs","totyrs",genes)
significance_table[3,2] <- length(which(totyrs_res_all$P.Value < 0.05))
significance_table[3,3] <- length(which(totyrs_res_all$adj.P.Val < 0.05))
totyrs_fgsea_res <- fgsea_function(totyrs_res_all,C2_genesets,"totyrs")
write.csv(totyrs_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_totyrs_model_gseaResults.csv'))
significance_table[3,4] <- length(which(totyrs_fgsea_res$pval < 0.05))
significance_table[3,5] <- length(which(totyrs_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + totyrs`
adat_totyrs_PathAD <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))]
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs_PathAD)), genes = rowData(adat_totyrs_PathAD)$EntrezGeneSymbol)
totyrs_PathAD_model <- model.matrix(~ agedeath + PathAD + totyrs, data=colData(adat_totyrs_PathAD))
totyrs_PathAD_res_all <- limma_model(totyrs_PathAD_model,adat_totyrs_PathAD, "totyrs + PathAD","totyrs",genes)
significance_table[4,2] <- length(which(totyrs_PathAD_res_all$P.Value < 0.05))
significance_table[4,3] <- length(which(totyrs_PathAD_res_all$adj.P.Val < 0.05))
totyrs_PathAD_fgsea_res <- fgsea_function(totyrs_PathAD_res_all,C2_genesets,"totyrs + PathAD")
write.csv(totyrs_PathAD_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_totyrs_PathAD_model_gseaResults.csv'))
significance_table[4,4] <- length(which(totyrs_PathAD_fgsea_res$pval < 0.05))
significance_table[4,5] <- length(which(totyrs_PathAD_fgsea_res$padj < 0.05))

#~ age_at_death + CTE
adat_CTE <- adat[,which(!is.na(colData(adat)$CTE))]
genes <- data.frame(unique_names = rownames(rowData(adat_CTE)), genes = rowData(adat_CTE)$EntrezGeneSymbol)
CTE_model <- model.matrix(~ agedeath + CTE, data=colData(adat_CTE))
head(CTE_model)
CTE_res_all <- limma_model(CTE_model,adat_CTE, "CTE","CTE",genes)
significance_table[5,2] <- length(which(CTE_res_all$P.Value < 0.05))
significance_table[5,3] <- length(which(CTE_res_all$adj.P.Val < 0.05))
CTE_fgsea_res <- fgsea_function(CTE_res_all,C2_genesets,"CTE")
CTE_fgsea_res <- CTE_fgsea_res[order(CTE_fgsea_res$padj), ] 
CTE_fgsea_res <- CTE_fgsea_res[1:83,]
pathway_lists <- frequency(CTE_fgsea_res$leadingEdge)
pathway_dups <- subset(as.data.frame(table(CTE_fgsea_res$leadingEdge)),Freq >1)
head(CTE_fgsea_res$leadingEdge)
write.csv(CTE_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_CTE_model_gseaResults.csv'))
significance_table[5,4] <- length(which(CTE_fgsea_res$pval < 0.05))
significance_table[5,5] <- length(which(CTE_fgsea_res$padj < 0.05))

#~ age_at_death + CTEStage
adat_CTEStage <- adat[,which(!is.na(colData(adat)$CTEStage))]
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage)), genes = rowData(adat_CTEStage)$EntrezGeneSymbol)
CTEStage_model <- model.matrix(~ agedeath + CTEStage, data=colData(adat_CTEStage))
CTEStage_res_all <- limma_model(CTEStage_model,adat_CTEStage, "CTEStage","CTEStage",genes)
significance_table[6,2] <- length(which(CTEStage_res_all$P.Value < 0.05))
significance_table[6,3] <- length(which(CTEStage_res_all$adj.P.Val < 0.05))
CTEStage_fgsea_res <- fgsea_function(CTEStage_res_all,C2_genesets,"CTEStage")
CTEStage_fgsea_res <- CTEStage_fgsea_res[order(CTEStage_fgsea_res$padj), ] 
head(CTEStage_fgsea_res$leadingEdge)
write.csv(CTEStage_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_CTEStage_model_gseaResults.csv'))
significance_table[6,4] <- length(which(CTEStage_fgsea_res$pval < 0.05))
significance_table[6,5] <- length(which(CTEStage_fgsea_res$padj < 0.05))

#~ age_at_death + CTELowvsHigh
adat_CTEonly <- adat[,which(colData(adat)$CTE == 1)]
adat_CTEonly_lowvshigh <- adat_CTEonly[,which(!is.na(colData(adat_CTEonly)$Group_de))]
colData(adat_CTEonly_lowvshigh)$Group_de <- as.integer(ifelse(colData(adat_CTEonly_lowvshigh)$Group_de == "low", 0, ifelse(colData(adat_CTEonly_lowvshigh)$Group_de == "high", 1, colData(adat_CTEonly_lowvshigh)$Group_de))) 
genes <- data.frame(unique_names = rownames(rowData(adat_CTEonly_lowvshigh)), genes = rowData(adat_CTEonly_lowvshigh)$EntrezGeneSymbol)
CTEonly_lowvshigh_model <- model.matrix(~ agedeath + Group_de, data=colData(adat_CTEonly_lowvshigh))
CTEonly_lowvshigh_res_all <- limma_model(CTEonly_lowvshigh_model,adat_CTEonly_lowvshigh, "Group_de","Group_de",genes)
significance_table[7,2] <- length(which(CTEonly_lowvshigh_res_all$P.Value < 0.05))
significance_table[7,3] <- length(which(CTEonly_lowvshigh_res_all$adj.P.Val < 0.05))
CTEonly_lowvshigh_fgsea_res <- fgsea_function(CTEonly_lowvshigh_res_all,C2_genesets,"Group_de")
write.csv(CTEonly_lowvshigh_fgsea_res, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/limma_DEA_CTEonly_lowvshigh_model_gseaResults.csv'))
significance_table[7,4] <- length(which(CTEonly_lowvshigh_fgsea_res$pval < 0.05))
significance_table[7,5] <- length(which(CTEonly_lowvshigh_fgsea_res$padj < 0.05))

significance_table %>%
  gt()
write.csv(significance_table, paste0('/restricted/projectnb/cteseq/projects/somascan_analysis/limma_results/significant_limma_Proteins_Pathways_table.csv'))
