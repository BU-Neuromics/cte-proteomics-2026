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

#differential expression analysis (Limma)
#create design matrix and run limma
limma_model <- function(model,adat, name, coefficient,genes) {
  fit <- limma::lmFit(assays(adat)$counts, model)
  elimma <- eBayes(fit)
  
  #create results files
  #all results for all covariates
  limma_res_all <- topTable(elimma,number = Inf,adjust="BH",sort.by="F")
  limma_res_all$unique_names <- rownames(limma_res_all)
  limma_res_all <- merge(limma_res_all, genes, by.x = "unique_names", by.y="unique_names")
  limma_res_all <- arrange(limma_res_all, adj.P.Val)
  
  #all results for coef as coefficient
  coef_names <- colnames(elimma$coefficients)
  # Extract logFC and p-values for each coefficient
  limma_res_all_coef <- lapply(coef_names, function(cn) {
    tt <- topTable(elimma, coef = cn, number = Inf, sort.by = "none")
    tt <- tt[, c("logFC", "P.Value", "adj.P.Val")]
    colnames(tt) <- paste(cn, c("logFC", "P.Value", "adj.P.Val"), sep = "_")
    return(tt)
  })
  
  # Combine all into one data frame
  limma_res_all_coef_combined <- do.call(cbind, limma_res_all_coef)
  limma_res_all_coef_combined$unique_names <- rownames(limma_res_all_coef_combined)
  limma_res_all_coef_combined <- merge(limma_res_all_coef_combined, genes, by.x = "unique_names", by.y="unique_names")
  full_results <- merge(limma_res_all,limma_res_all_coef_combined, by = "unique_names")
  full_results$genes <- full_results$genes.x
  full_results <- full_results[,which(!colnames(full_results) == "genes.x" & !colnames(full_results) == "genes.y")]
  full_results <- arrange(full_results, adj.P.Val)
  #file to return for fgsea
  limma_res_all_coef <- topTable(elimma,number = Inf, adjust="BH",coef = coefficient,sort.by="P")
  limma_res_all_coef$unique_names <- rownames(limma_res_all_coef)
  limma_res_all_coef <- merge(limma_res_all_coef, genes, by.x = "unique_names", by.y="unique_names")
  limma_res_all_coef <- arrange(limma_res_all_coef, adj.P.Val)
  
  #volcano plot
  fit$unique_names <- rownames(fit)
  fit <- merge(fit, genes, by.x = "unique_names", by.y="unique_names")
  png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/covariates/limma_DEA_',name,'_volcanoplot.png'))
  volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  dev.off()
  volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  
  write.csv(limma_res_all_coef, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_',name,'_results_coef.csv'))
  write.csv(full_results, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/limma_DEA_',name,'_results_all.csv'))
  return(limma_res_all_coef)
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
  C2_genesets <- split(C2_genesets$gene, C2_genesets$term)
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
  fgseaRes <- arrange(fgseaRes, padj)
  #top 10 pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  gseaTable <- plotGseaTable(final_list[topPathways], stats, fgseaRes, 
                             gseaParam=0.5)
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseaPlots/covariates/limma_DEA_',name,'_gseaPlot.png'), plot = gseaTable, width = 8, height = 6)
  final_results <- as_tibble(fgseaRes)
  final_results <- final_results %>% unnest(leadingEdge)
  write.csv(final_results, paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_',name,'_gseaResults.csv'))
  return(fgseaRes)
}

significance_table <- matrix(data=NA,nrow = 11,ncol = 5)
significance_table <- as.data.frame(significance_table)
colnames(significance_table) <- c("Model","Nominal Significant Proteins",
                                  "Adjusted Significant Proteins","Nominal Significant Pathways",
                                  "Adjusted Significant Pathways")
significance_table$Model <- c("age_at_death + PathAD + PathLBD + PathFTD + AT8_total",
                              "age_at_death + PathAD + PathLBD + PathFTD + AT8_total*CTE",
                              "age_at_death + PathAD + PathLBD + PathFT + totyrs",
                              "age_at_death + PathAD + PathLBD + PathFT + totyrs*CTE",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTEStage",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTELowvsHigh",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTERHIvslow",
                              "age_at_death + PathAD + PathLBD + PathFTD + hitsperyear",
                              "age_at_death + PathAD + PathLBD + PathFTD + AT8_total + hitsperyear",
                              "age_at_death + PathAD + PathLBD + PathFTD + chii_g",
                              "age_at_death + PathAD + PathLBD + PathFTD + AT8_total + chii_g")
significance_table

#`~ age_at_death + PathAD + PathLBD + PathFTD + AT8_total`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8 <- adat_AT8[,which(!is.na(colData(adat_AT8)$PathAD))] #181
colData(adat_AT8)$PathLBD <- as.integer(ifelse(colData(adat_AT8)$PathLBD == 2, 1, 0)) 
genes <- data.frame(unique_names = rownames(rowData(adat_AT8)), genes = rowData(adat_AT8)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total, data=colData(adat_AT8))
AT8_res_all <- limma_model(AT8_model,adat_AT8, "AT8_total","AT8_total",genes)
significance_table[1,2] <- length(which(AT8_res_all$P.Value < 0.05))
significance_table[1,3] <- length(which(AT8_res_all$adj.P.Val < 0.05))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total")
significance_table[1,4] <- length(which(AT8_fgsea_res$pval < 0.05))
significance_table[1,5] <- length(which(AT8_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + PathLBD + PathFTD + AT8_total*CTE`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8 <- adat_AT8[,which(!is.na(colData(adat_AT8)$PathAD))] #181
colData(adat_AT8)$PathLBD <- as.integer(ifelse(colData(adat_AT8)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_AT8)), genes = rowData(adat_AT8)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total*CTE, data=colData(adat_AT8))
AT8_res_all <- limma_model(AT8_model,adat_AT8, "AT8_total_CTE","AT8_total:CTE",genes)
significance_table[2,2] <- length(which(AT8_res_all$P.Value < 0.05))
significance_table[2,3] <- length(which(AT8_res_all$adj.P.Val < 0.05))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total_CTE")
significance_table[2,4] <- length(which(AT8_fgsea_res$pval < 0.05))
significance_table[2,5] <- length(which(AT8_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + PathLBD + PathFT + totyrs`
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))]
adat_totyrs <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))]
colData(adat_totyrs)$PathLBD <- as.integer(ifelse(colData(adat_totyrs)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + totyrs, data=colData(adat_totyrs))
totyrs_res_all <- limma_model(totyrs_model,adat_totyrs, "totyrs","totyrs",genes)
significance_table[3,2] <- length(which(totyrs_res_all$P.Value < 0.05))
significance_table[3,3] <- length(which(totyrs_res_all$adj.P.Val < 0.05))
totyrs_fgsea_res <- fgsea_function(totyrs_res_all,C2_genesets,"totyrs")
significance_table[3,4] <- length(which(totyrs_fgsea_res$pval < 0.05))
significance_table[3,5] <- length(which(totyrs_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + PathLBD + PathFT + totyrs*CTE`
#interaction model
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))]
adat_totyrs <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))]
colData(adat_totyrs)$PathLBD <- as.integer(ifelse(colData(adat_totyrs)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + totyrs*CTE, data=colData(adat_totyrs))
totyrs_res_all <- limma_model(totyrs_model,adat_totyrs, "totyrs_CTE","totyrs:CTE",genes)
significance_table[4,2] <- length(which(totyrs_res_all$P.Value < 0.05))
significance_table[4,3] <- length(which(totyrs_res_all$adj.P.Val < 0.05))
totyrs_fgsea_res <- fgsea_function(totyrs_res_all,C2_genesets,"totyrs_CTE")
significance_table[4,4] <- length(which(totyrs_fgsea_res$pval < 0.05))
significance_table[4,5] <- length(which(totyrs_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTEStage
adat_CTEStage <- adat[,which(!is.na(colData(adat)$CTEStage))]
adat_CTEStage <- adat_CTEStage[,which(!is.na(colData(adat_CTEStage)$PathAD))]
colData(adat_CTEStage)$PathLBD <- as.integer(ifelse(colData(adat_CTEStage)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage)), genes = rowData(adat_CTEStage)$EntrezGeneSymbol)
CTEStage_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + CTEStage, data=colData(adat_CTEStage))
CTEStage_res_all <- limma_model(CTEStage_model,adat_CTEStage, "CTEStage","CTEStage",genes)
significance_table[5,2] <- length(which(CTEStage_res_all$P.Value < 0.05))
significance_table[5,3] <- length(which(CTEStage_res_all$adj.P.Val < 0.05))
CTEStage_fgsea_res <- fgsea_function(CTEStage_res_all,C2_genesets,"CTEStage")
significance_table[5,4] <- length(which(CTEStage_fgsea_res$pval < 0.05))
significance_table[5,5] <- length(which(CTEStage_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTELowvsHigh
#no CTE = removed
#low CTE = 0
#high CTE = 1
#remember to add back in other steps if we remove CTE stage model
adat_highlow <- adat_CTEStage[,which(!(colData(adat_CTEStage)$CTEStage) == 0)]
colData(adat_highlow)$Group_de <- as.integer(ifelse(colData(adat_highlow)$CTEStage == 1, 0, ifelse(colData(adat_highlow)$CTEStage == 2, 0, ifelse(colData(adat_highlow)$CTEStage == 3, 1, ifelse(colData(adat_highlow)$CTEStage == 4, 1, colData(adat_highlow)$Group_de)))))
genes <- data.frame(unique_names = rownames(rowData(adat_highlow)), genes = rowData(adat_highlow)$EntrezGeneSymbol)
CTE_lowvshigh_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + Group_de, data=colData(adat_highlow))
CTE_lowvshigh_res_all <- limma_model(CTE_lowvshigh_model,adat_highlow, "CTE_lowvshigh","Group_de",genes)
significance_table[6,2] <- length(which(CTE_lowvshigh_res_all$P.Value < 0.05))
significance_table[6,3] <- length(which(CTE_lowvshigh_res_all$adj.P.Val < 0.05))
CTE_lowvshigh_fgsea_res <- fgsea_function(CTE_lowvshigh_res_all,C2_genesets,"CTE_lowvshigh")
CTE_lowvshigh_fgsea_res <- arrange(CTE_lowvshigh_fgsea_res, padj)
significance_table[6,4] <- length(which(CTE_lowvshigh_fgsea_res$pval < 0.05))
significance_table[6,5] <- length(which(CTE_lowvshigh_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTERHIvslow
#no CTE = 0
#low CTE = 1
#high CTE = removed
#remember to add back in other steps if we remove CTE stage model
adat_RHIlow <- adat_CTEStage[,which(!(colData(adat_CTEStage)$Group_de) == "high")]
colData(adat_RHIlow)$Group_de <- as.integer(ifelse(colData(adat_RHIlow)$CTEStage == 1, 1, ifelse(colData(adat_RHIlow)$CTEStage == 2, 1, ifelse(colData(adat_RHIlow)$CTEStage == 0, 0, colData(adat_RHIlow)$Group_de))))
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage)), genes = rowData(adat_CTEStage)$EntrezGeneSymbol)
CTE_RHIvslow_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + Group_de, data=colData(adat_RHIlow))
CTE_RHIvslow_res_all <- limma_model(CTE_RHIvslow_model,adat_RHIlow, "CTE_RHIvslow","Group_de",genes)
significance_table[7,2] <- length(which(CTE_RHIvslow_res_all$P.Value < 0.05))
significance_table[7,3] <- length(which(CTE_RHIvslow_res_all$adj.P.Val < 0.05))
CTE_RHIvslow_fgsea_res <- fgsea_function(CTE_RHIvslow_res_all,C2_genesets,"CTE_RHIvslow")
CTE_RHIvslow_fgsea_res <- arrange(CTE_RHIvslow_fgsea_res, padj)
significance_table[7,4] <- length(which(CTE_RHIvslow_fgsea_res$pval < 0.05))
significance_table[7,5] <- length(which(CTE_RHIvslow_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + hitsperyear
adat_hits <- adat[,which(!is.na(colData(adat)$hitsperyear))]
adat_hits <- adat_hits[,which(!is.na(colData(adat_hits)$PathAD))]
colData(adat_hits)$PathLBD <- as.integer(ifelse(colData(adat_hits)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_hits)), genes = rowData(adat_hits)$EntrezGeneSymbol)
hits_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + hitsperyear, data=colData(adat_hits))
hits_res_all <- limma_model(hits_model,adat_hits, "hitsperyear","hitsperyear",genes)
significance_table[8,2] <- length(which(hits_res_all$P.Value < 0.05))
significance_table[8,3] <- length(which(hits_res_all$adj.P.Val < 0.05))
hits_fgsea_res <- fgsea_function(hits_res_all,C2_genesets,"hitsperyear")
significance_table[8,4] <- length(which(hits_fgsea_res$pval < 0.05))
significance_table[8,5] <- length(which(hits_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + hitsperyear + AT8_total
adat_hits <- adat[,which(!is.na(colData(adat)$hitsperyear))]
adat_hits <- adat_hits[,which(!is.na(colData(adat_hits)$PathAD))]
adat_hits <- adat_hits[,which(!is.na(colData(adat_hits)$AT8_total))]
colData(adat_hits)$PathLBD <- as.integer(ifelse(colData(adat_hits)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_hits)), genes = rowData(adat_hits)$EntrezGeneSymbol)
hits_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total + hitsperyear, data=colData(adat_hits))
hits_res_all <- limma_model(hits_model,adat_hits, "hitsperyear_with_tau","hitsperyear",genes)
significance_table[9,2] <- length(which(hits_res_all$P.Value < 0.05))
significance_table[9,3] <- length(which(hits_res_all$adj.P.Val < 0.05))
hits_fgsea_res <- fgsea_function(hits_res_all,C2_genesets,"hitsperyear_with_tau")
significance_table[9,4] <- length(which(hits_fgsea_res$pval < 0.05))
significance_table[9,5] <- length(which(hits_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + chii_g
adat_ghits <- adat[,which(!is.na(colData(adat)$chii_g))]
adat_ghits <- adat_ghits[,which(!is.na(colData(adat_ghits)$PathAD))]
colData(adat_ghits)$PathLBD <- as.integer(ifelse(colData(adat_ghits)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_ghits)), genes = rowData(adat_ghits)$EntrezGeneSymbol)
ghits_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + chii_g, data=colData(adat_ghits))
ghits_res_all <- limma_model(ghits_model,adat_ghits, "chii_g","chii_g",genes)
significance_table[10,2] <- length(which(ghits_res_all$P.Value < 0.05))
significance_table[10,3] <- length(which(ghits_res_all$adj.P.Val < 0.05))
ghits_fgsea_res <- fgsea_function(ghits_res_all,C2_genesets,"chii_g")
significance_table[10,4] <- length(which(ghits_fgsea_res$pval < 0.05))
significance_table[10,5] <- length(which(ghits_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + chii_g + AT8_total
adat_ghits <- adat[,which(!is.na(colData(adat)$chii_g))]
adat_ghits <- adat_ghits[,which(!is.na(colData(adat_ghits)$PathAD))]
adat_ghits <- adat_ghits[,which(!is.na(colData(adat_ghits)$AT8_total))]
colData(adat_ghits)$PathLBD <- as.integer(ifelse(colData(adat_ghits)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_ghits)), genes = rowData(adat_ghits)$EntrezGeneSymbol)
ghits_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total + chii_g, data=colData(adat_ghits))
ghits_res_all <- limma_model(ghits_model,adat_ghits, "chii_g_with_tau","chii_g",genes)
significance_table[11,2] <- length(which(ghits_res_all$P.Value < 0.05))
significance_table[11,3] <- length(which(ghits_res_all$adj.P.Val < 0.05))
ghits_fgsea_res <- fgsea_function(ghits_res_all,C2_genesets,"chii_g_with_tau")
significance_table[11,4] <- length(which(ghits_fgsea_res$pval < 0.05))
significance_table[11,5] <- length(which(ghits_fgsea_res$padj < 0.05))

significance_table %>%
  gt()
write.csv(significance_table, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/significant_limma_Proteins_Pathways_table.csv'))

#impact of AD samples - no apparent impact on distribution of age of death
#with AD samples
hist(colData(adat)$agedeath)
mean(colData(adat)$agedeath, na.rm = TRUE)
adat_noAD <- adat[,which(colData(adat)$PathAD != 1)]
#check to make sure filtering worked
colData(adat_noAD)$PathAD
hist(colData(adat_noAD)$agedeath)
mean(colData(adat_noAD)$agedeath)