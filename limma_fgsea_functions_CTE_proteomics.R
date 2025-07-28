#Helen Pennington
#Labadorf Rotation 
#Started March 1, 2025
#Limma Function + fgsea Analysis

#packages 
library(dplyr)
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

# Read in the adat file and pathways file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
C2_genesets <- read.gmt("/restricted/projectnb/cteseq/projects/somascan/data/c2.all.v2024.1.Hs.symbols.gmt")

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
    tt <- topTable(elimma, coef = cn, number = Inf, sort.by = "none",adjust="BH")
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

significance_table <- matrix(data=NA,nrow = 13,ncol = 5)
significance_table <- as.data.frame(significance_table)
colnames(significance_table) <- c("Model","Nominal Significant Proteins",
                                  "Adjusted Significant Proteins","Nominal Significant Pathways",
                                  "Adjusted Significant Pathways")
significance_table$Model <- c("age_at_death + PathAD + PathLBD + PathFTD + AT8_total",
                              "age_at_death + PathAD + PathLBD + PathFTD + AT8_total*CTE",
                              "age_at_death + PathAD + PathLBD + PathFT + totyrs",
                              "age_at_death + PathAD + PathLBD + PathFT + totyrs*CTE",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTE",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTEStage",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTELowvsHigh",
                              "age_at_death + PathAD + PathLBD + PathFTD + CTERHIvslow",
                              "age_at_death + PathAD + PathLBD + PathFTD + maxaggsum",
                              "age_at_death + PathAD + PathLBD + PathFTD + CDStot",
                              "age_at_death + PathAD + PathLBD + PathFTD + DementiaHx",
                              "age_at_death + PathAD + PathLBD + PathFTD + faqtot",
                              "age_at_death + PathAD + PathLBD + PathFTD + chii_g")
significance_table

#`~ age_at_death + PathAD + PathLBD + PathFTD + AT8_total`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8 <- adat_AT8[,which(!is.na(colData(adat_AT8)$PathAD))] #181
colData(adat_AT8)$PathLBD <- as.integer(ifelse(colData(adat_AT8)$PathLBD == 2, 1, 0)) 
genes <- data.frame(unique_names = rownames(rowData(adat_AT8)), genes = rowData(adat_AT8)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total, data=colData(adat_AT8))
AT8_res_all <- limma_model(AT8_model,adat_AT8, "AT8_total","AT8_total",genes)
significance_table[1,2] <- length(which(AT8_res_all$P.Value < 0.1))
significance_table[1,3] <- length(which(AT8_res_all$adj.P.Val < 0.1))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total")
significance_table[1,4] <- length(which(AT8_fgsea_res$pval < 0.1))
significance_table[1,5] <- length(which(AT8_fgsea_res$padj < 0.1))

#`~ age_at_death + PathAD + PathLBD + PathFTD + AT8_total*CTE`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #183
adat_AT8 <- adat_AT8[,which(!is.na(colData(adat_AT8)$PathAD))] #181
colData(adat_AT8)$PathLBD <- as.integer(ifelse(colData(adat_AT8)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_AT8)), genes = rowData(adat_AT8)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + AT8_total*CTE, data=colData(adat_AT8))
AT8_res_all <- limma_model(AT8_model,adat_AT8, "AT8_total_CTE","AT8_total:CTE1",genes)
significance_table[2,2] <- length(which(AT8_res_all$P.Value < 0.05))
significance_table[2,3] <- length(which(AT8_res_all$adj.P.Val < 0.05))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total_CTE")
significance_table[2,4] <- length(which(AT8_fgsea_res$pval < 0.05))
significance_table[2,5] <- length(which(AT8_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + PathLBD + PathFT + totyrs`
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))] #196
adat_totyrs <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))] #195
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
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))] #196
adat_totyrs <- adat_totyrs[,which(!is.na(colData(adat_totyrs)$PathAD))] #195
colData(adat_totyrs)$PathLBD <- as.integer(ifelse(colData(adat_totyrs)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + totyrs*CTE, data=colData(adat_totyrs))
totyrs_res_all <- limma_model(totyrs_model,adat_totyrs, "totyrs_CTE","totyrs:CTE1",genes)
significance_table[4,2] <- length(which(totyrs_res_all$P.Value < 0.05))
significance_table[4,3] <- length(which(totyrs_res_all$adj.P.Val < 0.05))
totyrs_fgsea_res <- fgsea_function(totyrs_res_all,C2_genesets,"totyrs_CTE")
significance_table[4,4] <- length(which(totyrs_fgsea_res$pval < 0.05))
significance_table[4,5] <- length(which(totyrs_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTE
adat_CTE <- adat[,which(!is.na(colData(adat)$CTE))] #206
adat_CTE <- adat_CTE[,which(!is.na(colData(adat_CTE)$PathAD))] #205
colData(adat_CTE)$PathLBD <- as.integer(ifelse(colData(adat_CTE)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_CTE)), genes = rowData(adat_CTE)$EntrezGeneSymbol)
CTE_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + CTE, data=colData(adat_CTE))
CTE_res_all <- limma_model(CTE_model,adat_CTE, "CTE","CTE1",genes)
significance_table[5,2] <- length(which(CTE_res_all$P.Value < 0.05))
significance_table[5,3] <- length(which(CTE_res_all$adj.P.Val < 0.05))
CTE_fgsea_res <- fgsea_function(CTE_res_all,C2_genesets,"CTE")
significance_table[5,4] <- length(which(CTE_fgsea_res$pval < 0.05))
significance_table[5,5] <- length(which(CTE_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTEStage
adat_CTEStage <- adat[,which(!is.na(colData(adat)$CTEStage))] #205
adat_CTEStage <- adat_CTEStage[,which(!is.na(colData(adat_CTEStage)$PathAD))] #204
colData(adat_CTEStage)$PathLBD <- as.integer(ifelse(colData(adat_CTEStage)$PathLBD == 2, 1, 0))
adat_CTEStage_num <- adat_CTEStage
colData(adat_CTEStage_num)$CTEStage <- as.numeric(colData(adat_CTEStage)$CTEStage)
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage_num)), genes = rowData(adat_CTEStage_num)$EntrezGeneSymbol)
CTEStage_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + CTEStage, data=colData(adat_CTEStage_num))
CTEStage_res_all <- limma_model(CTEStage_model,adat_CTEStage_num, "CTEStage","CTEStage",genes)
significance_table[6,2] <- length(which(CTEStage_res_all$P.Value < 0.05))
significance_table[6,3] <- length(which(CTEStage_res_all$adj.P.Val < 0.05))
CTEStage_fgsea_res <- fgsea_function(CTEStage_res_all,C2_genesets,"CTEStage")
significance_table[6,4] <- length(which(CTEStage_fgsea_res$pval < 0.05))
significance_table[6,5] <- length(which(CTEStage_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTELowvsHigh
#no CTE = removed
#low CTE = 0
#high CTE = 1
#DO NOT RUN THIS WITHOUT CREATING adat_CTEStage (code block above)
adat_highlow <- adat_CTEStage[,which(!(colData(adat_CTEStage)$CTEStage) == 0)] #146
adat_highlow <- adat_highlow[,which(!is.na(colData(adat_highlow)$Group_de))]
colData(adat_highlow)$Group_de <- as.integer(ifelse(colData(adat_highlow)$CTEStage == 2, 0, ifelse(colData(adat_highlow)$CTEStage == 3, 0, ifelse(colData(adat_highlow)$CTEStage == 4, 1, ifelse(colData(adat_highlow)$CTEStage == 5, 1, colData(adat_highlow)$Group_de)))))
genes <- data.frame(unique_names = rownames(rowData(adat_highlow)), genes = rowData(adat_highlow)$EntrezGeneSymbol)
CTE_lowvshigh_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + Group_de, data=colData(adat_highlow))
CTE_lowvshigh_res_all <- limma_model(CTE_lowvshigh_model,adat_highlow, "CTE_lowvshigh","Group_de",genes)
significance_table[7,2] <- length(which(CTE_lowvshigh_res_all$P.Value < 0.05))
significance_table[7,3] <- length(which(CTE_lowvshigh_res_all$adj.P.Val < 0.05))
CTE_lowvshigh_fgsea_res <- fgsea_function(CTE_lowvshigh_res_all,C2_genesets,"CTE_lowvshigh")
CTE_lowvshigh_fgsea_res <- arrange(CTE_lowvshigh_fgsea_res, padj)
significance_table[7,4] <- length(which(CTE_lowvshigh_fgsea_res$pval < 0.05))
significance_table[7,5] <- length(which(CTE_lowvshigh_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTERHIvslow
#no CTE = 0
#low CTE = 1
#high CTE = removed
##DO NOT RUN THIS WITHOUT CREATING adat_CTEStage (code block 2 above)
adat_RHIlow <- adat_CTEStage[,which(!(colData(adat_CTEStage)$CTEStage) == 3 & !(colData(adat_CTEStage)$CTEStage) == 4)] #90
adat_RHIlow <- adat_RHIlow[,which(!is.na(colData(adat_RHIlow)$Group_de))]
colData(adat_RHIlow)$Group_de <- as.integer(ifelse(colData(adat_RHIlow)$CTEStage == 1, 1, ifelse(colData(adat_RHIlow)$CTEStage == 2, 1, ifelse(colData(adat_RHIlow)$CTEStage == 0, 0, colData(adat_RHIlow)$Group_de))))
genes <- data.frame(unique_names = rownames(rowData(adat_CTEStage)), genes = rowData(adat_CTEStage)$EntrezGeneSymbol)
CTE_RHIvslow_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + Group_de, data=colData(adat_RHIlow))
CTE_RHIvslow_res_all <- limma_model(CTE_RHIvslow_model,adat_RHIlow, "CTE_RHIvslow","Group_de",genes)
significance_table[8,2] <- length(which(CTE_RHIvslow_res_all$P.Value < 0.05))
significance_table[8,3] <- length(which(CTE_RHIvslow_res_all$adj.P.Val < 0.05))
CTE_RHIvslow_fgsea_res <- fgsea_function(CTE_RHIvslow_res_all,C2_genesets,"CTE_RHIvslow")
CTE_RHIvslow_fgsea_res <- arrange(CTE_RHIvslow_fgsea_res, padj)
significance_table[8,4] <- length(which(CTE_RHIvslow_fgsea_res$pval < 0.05))
significance_table[8,5] <- length(which(CTE_RHIvslow_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + maxaggsum
adat_maxaggsum <- adat[,which(!is.na(colData(adat)$maxaggsum))] #170
adat_maxaggsum <- adat_maxaggsum[,which(!is.na(colData(adat_maxaggsum)$agedeath))] #169
genes <- data.frame(unique_names = rownames(rowData(adat_maxaggsum)), genes = rowData(adat_maxaggsum)$EntrezGeneSymbol)
maxaggsum_model <- model.matrix(~ agedeath + maxaggsum, data=colData(adat_maxaggsum))
maxaggsum_res_all <- limma_model(maxaggsum_model,adat_maxaggsum, "maxaggsum","maxaggsum",genes)
significance_table[9,2] <- length(which(maxaggsum_res_all$P.Value < 0.05))
significance_table[9,3] <- length(which(maxaggsum_res_all$adj.P.Val < 0.05))
maxaggsum_fgsea_res <- fgsea_function(maxaggsum_res_all,C2_genesets,"maxaggsum")
significance_table[9,4] <- length(which(maxaggsum_fgsea_res$pval < 0.05))
significance_table[9,5] <- length(which(maxaggsum_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CDStot
adat_cds <- adat[,which(!is.na(colData(adat)$CDStot))] #134
adat_cds <- adat_cds[,which(!is.na(colData(adat_cds)$PathAD))] #133
colData(adat_cds)$PathLBD <- as.integer(ifelse(colData(adat_cds)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_cds)), genes = rowData(adat_cds)$EntrezGeneSymbol)
cds_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + CDStot, data=colData(adat_cds))
cds_res_all <- limma_model(cds_model,adat_cds, "CDStot","CDStot",genes)
significance_table[10,2] <- length(which(cds_res_all$P.Value < 0.05))
significance_table[10,3] <- length(which(cds_res_all$adj.P.Val < 0.05))
cds_fgsea_res <- fgsea_function(cds_res_all,C2_genesets,"CDStot")
significance_table[10,4] <- length(which(cds_fgsea_res$pval < 0.05))
significance_table[10,5] <- length(which(cds_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + Dementia
adat_dementia <- adat[,which(!is.na(colData(adat)$DementiaHx))] #207
adat_dementia <- adat_dementia[,which(!is.na(colData(adat_dementia)$agedeath))] #206
genes <- data.frame(unique_names = rownames(rowData(adat_dementia)), genes = rowData(adat_dementia)$EntrezGeneSymbol)
dementia_model <- model.matrix(~ agedeath + DementiaHx, data=colData(adat_dementia))
dementia_res_all <- limma_model(dementia_model,adat_dementia, "DementiaHx","DementiaHx1",genes)
significance_table[11,2] <- length(which(dementia_res_all$P.Value < 0.05))
significance_table[11,3] <- length(which(dementia_res_all$adj.P.Val < 0.05))
dementia_fgsea_res <- fgsea_function(dementia_res_all,C2_genesets,"DementiaHx")
significance_table[11,4] <- length(which(dementia_fgsea_res$pval < 0.05))
significance_table[11,5] <- length(which(dementia_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + faqtot
adat_faq <- adat[,which(!is.na(colData(adat)$faqtot))] #180
adat_faq <- adat_faq[,which(!is.na(colData(adat_faq)$PathAD))] #178
colData(adat_faq)$PathLBD <- as.integer(ifelse(colData(adat_faq)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_faq)), genes = rowData(adat_faq)$EntrezGeneSymbol)
faq_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + faqtot, data=colData(adat_faq))
faq_res_all <- limma_model(faq_model,adat_faq, "faqtot","faqtot",genes)
significance_table[12,2] <- length(which(faq_res_all$P.Value < 0.05))
significance_table[12,3] <- length(which(faq_res_all$adj.P.Val < 0.05))
faq_fgsea_res <- fgsea_function(faq_res_all,C2_genesets,"faqtot")
significance_table[12,4] <- length(which(faq_fgsea_res$pval < 0.05))
significance_table[12,5] <- length(which(faq_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + chii_g
adat_ghits <- adat[,which(!is.na(colData(adat)$chii_g))] #172
colData(adat_ghits)$PathLBD <- as.integer(ifelse(colData(adat_ghits)$PathLBD == 2, 1, 0))
genes <- data.frame(unique_names = rownames(rowData(adat_ghits)), genes = rowData(adat_ghits)$EntrezGeneSymbol)
ghits_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + chii_g, data=colData(adat_ghits))
ghits_res_all <- limma_model(ghits_model,adat_ghits, "chii_g","chii_g",genes)
significance_table[13,2] <- length(which(ghits_res_all$P.Value < 0.05))
significance_table[13,3] <- length(which(ghits_res_all$adj.P.Val < 0.05))
ghits_fgsea_res <- fgsea_function(ghits_res_all,C2_genesets,"chii_g")
significance_table[13,4] <- length(which(ghits_fgsea_res$pval < 0.05))
significance_table[13,5] <- length(which(ghits_fgsea_res$padj < 0.05))

significance_table %>%
  gt()
write.csv(significance_table, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/covariates/significant_limma_Proteins_Pathways_table.csv'))
