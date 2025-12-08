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
library(EnhancedVolcano)
library(webshot2)

# Read in the adat file and pathways file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
C2_genesets <- read.gmt("/restricted/projectnb/cteseq/projects/somascan/data/c2.cp.v2025.1.Hs.symbols.gmt")
amelia.out1 <- readRDS("/restricted/projectnb/cteseq/projects/somascan/imputation/Stable_Releases/2025-12-02/full_imputations.rds")
colData(adat)$PMI <- amelia.out1$imputations[[1]]$PMI
colData(adat)$PathAD <- amelia.out1$imputations[[1]]$PathAD
colData(adat)$PathLBD <- factor(colData(adat)$PathLBD)
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
  limma_res_all_coef <- arrange(limma_res_all_coef, P.Value)
  
  #volcano plot
  fit$unique_names <- rownames(fit)
  fit <- merge(fit, genes, by.x = "unique_names", by.y="unique_names")
  #png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_',name,'_noF_volcanoplot.png'))
  #volcanoplot(elimma,coef=coefficient,highlight=20,names=fit$genes, col="red",main=name)
  #dev.off()
  p_cutoff <- 0.05
  fc_cutoff <- 0.05  # visual threshold you mentioned
  
  # compute -log10 pvalue for y-axis
  limma_res_all_coef <- limma_res_all_coef %>%
    mutate(
      negLogP = -log10(P.Value),
      Significance = case_when(
        P.Value <= p_cutoff & abs(logFC) >= fc_cutoff ~ "p<0.05 & |logFC|>=1",
        P.Value <= p_cutoff ~ "p<0.05",
        TRUE ~ "not-significant"
      )
    )
  
  # pick top 10 most significant genes (by P.Value)
  top10_genes <- limma_res_all_coef %>%
    arrange(P.Value) %>%
    slice_head(n = 10) %>%
    pull(genes)
  
  # base plot
  p <- ggplot(limma_res_all_coef, aes(x = logFC, y = negLogP)) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "not-significant" = "grey30",
      "p<0.05" = "royalblue",
      "p<0.05 & |logFC|>=1" = "red2"
    )) +
    theme_minimal(base_size = 20) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(P.Value)",
      caption = paste0("Total genes = ", nrow(limma_res_all_coef),
                       ", p-cutoff = ", p_cutoff)
    ) +
    # vertical lines showing your lfc cutoffs for visual reference
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "darkgrey")
  
  # add labels for top10, with segments (connectors)
  vplot <- p + geom_text_repel(
    data = limma_res_all_coef %>% filter(genes %in% top10_genes),
    aes(label = genes),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.5,
    segment.color = "grey40"
  )
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_',name,'_ggplot_volcanoplot.png'), plot = vplot, width = 10, height = 6)
  saveRDS(vplot, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_',name,'_ggplot_volcanoplot.rds'))
  write.csv(limma_res_all_coef, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/limma_DEA_',name,'_noF_results_coef.csv'))
  write.csv(full_results, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/limma_DEA_',name,'_noF_results_all.csv'))
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
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseaPlots/PMI_imp1/limma_DEA_',name,'_noF_gseaPlot.png'), plot = gseaTable, width = 8, height = 6)
  #final_results <- as_tibble(fgseaRes)
  #final_results <- final_results %>% unnest(leadingEdge)
  #write.csv(final_results, paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/covariates/limma_DEA_',name,'_noF_gseaResults.csv'))
  fgseaRes_csv <- fgseaRes %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))
  write.csv(fgseaRes_csv, paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/fgsea_compact_',name,'.csv'))
  return(fgseaRes)
}

significance_table <- matrix(data=NA,nrow = 7,ncol = 5)
significance_table <- as.data.frame(significance_table)
colnames(significance_table) <- c("Model","Nominal Significant Proteins",
                                  "Adjusted Significant Proteins","Nominal Significant Pathways",
                                  "Adjusted Significant Pathways")
significance_table$Model <- c("AT8 Total", "Total Years of Play", "RHI to Low CTE", 
                              "RHI to High CTE", "Cognitive Difficulty Score Total", 
                              "Dementia", "Functional Activities Questionnaire Total")
significance_table

#`~ age_at_death + PathAD + PathLBD + PathFTD + AT8_total`
adat_AT8 <- adat[,which(!is.na(colData(adat)$AT8_total))] #181
dim(adat_AT8)
genes <- data.frame(unique_names = rownames(rowData(adat_AT8)), genes = rowData(adat_AT8)$EntrezGeneSymbol)
AT8_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + AT8_total, data=colData(adat_AT8))
AT8_res_all <- limma_model(AT8_model,adat_AT8, "AT8_total","AT8_total",genes)
significance_table[1,2] <- length(which(AT8_res_all$P.Value < 0.05))
significance_table[1,3] <- length(which(AT8_res_all$adj.P.Val < 0.05))
AT8_fgsea_res <- fgsea_function(AT8_res_all,C2_genesets,"AT8_total")
significance_table[1,4] <- length(which(AT8_fgsea_res$pval < 0.05))
significance_table[1,5] <- length(which(AT8_fgsea_res$padj < 0.05))

#`~ age_at_death + PathAD + PathLBD + PathFT + totyrs`
adat_totyrs <- adat[,which(!is.na(colData(adat)$totyrs))] #195
dim(adat_totyrs)
genes <- data.frame(unique_names = rownames(rowData(adat_totyrs)), genes = rowData(adat_totyrs)$EntrezGeneSymbol)
totyrs_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + totyrs, data=colData(adat_totyrs))
totyrs_res_all <- limma_model(totyrs_model,adat_totyrs, "totyrs","totyrs",genes)
significance_table[2,2] <- length(which(totyrs_res_all$P.Value < 0.05))
significance_table[2,3] <- length(which(totyrs_res_all$adj.P.Val < 0.05))
totyrs_fgsea_res <- fgsea_function(totyrs_res_all,C2_genesets,"totyrs")
significance_table[2,4] <- length(which(totyrs_fgsea_res$pval < 0.05))
significance_table[2,5] <- length(which(totyrs_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTELowvsHigh
adat_highlow <- adat[,which(!(colData(adat)$Group_de) == 0)] #147
dim(adat_highlow)
colData(adat_highlow)$Group_de <- factor(colData(adat_highlow)$Group_de, levels = c(1,2))
genes <- data.frame(unique_names = rownames(rowData(adat_highlow)), genes = rowData(adat_highlow)$EntrezGeneSymbol)
CTE_lowvshigh_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + Group_de, data=colData(adat_highlow))
CTE_lowvshigh_res_all <- limma_model(CTE_lowvshigh_model,adat_highlow, "CTE_lowvshigh","Group_de2",genes)
CTE_lowvshigh_fgsea_res <- fgsea_function(CTE_lowvshigh_res_all,C2_genesets,"CTE_lowvshigh")

#RHI vs Stages
#stage 1
adat_rhi1 <- adat[,which(!(colData(adat)$CTEStage) == 4 & !(colData(adat)$CTEStage) == 2 & !(colData(adat)$CTEStage) == 3)] #71
colData(adat_rhi1)$CTEStage <- factor(colData(adat_rhi1)$CTEStage, levels = c(0,1))
genes <- data.frame(unique_names = rownames(rowData(adat_rhi1)), genes = rowData(adat_rhi1)$EntrezGeneSymbol)
CTE_rhi1_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + CTEStage, data=colData(adat_rhi1))
CTE_rhi1_res_all <- limma_model(CTE_rhi1_model,adat_rhi1, "CTE_RHI1","CTEStage1",genes)
CTE_rhi1_fgsea_res <- fgsea_function(CTE_rhi1_res_all,C2_genesets,"CTE_RHI1")
#stage 2
adat_rhi2 <- adat[,which(!(colData(adat)$CTEStage) == 1 & !(colData(adat)$CTEStage) == 4 & !(colData(adat)$CTEStage) == 3)] #77
colData(adat_rhi2)$CTEStage <- factor(colData(adat_rhi2)$CTEStage, levels = c(0,2))
genes <- data.frame(unique_names = rownames(rowData(adat_rhi2)), genes = rowData(adat_rhi2)$EntrezGeneSymbol)
CTE_rhi2_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + CTEStage, data=colData(adat_rhi2))
CTE_rhi2_res_all <- limma_model(CTE_rhi2_model,adat_rhi2, "CTE_RHI2","CTEStage2",genes)
CTE_rhi2_fgsea_res <- fgsea_function(CTE_rhi2_res_all,C2_genesets,"CTE_RHI2")
#stage 3
adat_rhi3 <- adat[,which(!(colData(adat)$CTEStage) == 1 & !(colData(adat)$CTEStage) == 2 & !(colData(adat)$CTEStage) == 4)] #105
colData(adat_rhi3)$CTEStage <- factor(colData(adat_rhi3)$CTEStage, levels = c(0,3))
genes <- data.frame(unique_names = rownames(rowData(adat_rhi3)), genes = rowData(adat_rhi3)$EntrezGeneSymbol)
CTE_rhi3_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + CTEStage, data=colData(adat_rhi3))
CTE_rhi3_res_all <- limma_model(CTE_rhi3_model,adat_rhi3, "CTE_RHI3","CTEStage3",genes)
CTE_rhi3_fgsea_res <- fgsea_function(CTE_rhi3_res_all,C2_genesets,"CTE_RHI3")
#stage 4
adat_rhi4 <- adat[,which(!(colData(adat)$CTEStage) == 1 & !(colData(adat)$CTEStage) == 2 & !(colData(adat)$CTEStage) == 3)] #125
colData(adat_rhi4)$CTEStage <- factor(colData(adat_rhi4)$CTEStage, levels = c(0,4))
genes <- data.frame(unique_names = rownames(rowData(adat_rhi4)), genes = rowData(adat_rhi4)$EntrezGeneSymbol)
CTE_rhi4_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + CTEStage, data=colData(adat_rhi4))
CTE_rhi4_res_all <- limma_model(CTE_rhi4_model,adat_rhi4, "CTE_RHI4","CTEStage4",genes)
CTE_rhi4_fgsea_res <- fgsea_function(CTE_rhi4_res_all,C2_genesets,"CTE_RHI4")

#RHI vs High CTE
adat_highrhi <- adat[,which(!(colData(adat)$Group_de) == 1)] #172
colData(adat_highrhi)$Group_de <- factor(colData(adat_highrhi)$Group_de)
genes <- data.frame(unique_names = rownames(rowData(adat_highrhi)), genes = rowData(adat_highrhi)$EntrezGeneSymbol)
CTE_rhivshigh_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + Group_de, data=colData(adat_highrhi))
CTE_rhivshigh_res_all <- limma_model(CTE_rhivshigh_model,adat_highrhi, "CTE_RHIvsHigh","Group_de2",genes)
significance_table[4,2] <- length(which(CTE_rhivshigh_res_all$P.Value < 0.05))
significance_table[4,3] <- length(which(CTE_rhivshigh_res_all$adj.P.Val < 0.05))
CTE_rhivshigh_fgsea_res <- fgsea_function(CTE_rhivshigh_res_all,C2_genesets,"CTE_RHIvsHigh")
significance_table[4,4] <- length(which(CTE_rhivshigh_fgsea_res$pval < 0.05))
significance_table[4,5] <- length(which(CTE_rhivshigh_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CTERHIvslow
adat_RHIlow <- adat[,which(!(colData(adat)$Group_de) == 2)] #91
colData(adat_RHIlow)$Group_de <- factor(colData(adat_RHIlow)$Group_de)
genes <- data.frame(unique_names = rownames(rowData(adat_RHIlow)), genes = rowData(adat_RHIlow)$EntrezGeneSymbol)
CTE_RHIvslow_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + Group_de, data=colData(adat_RHIlow))
CTE_RHIvslow_res_all <- limma_model(CTE_RHIvslow_model,adat_RHIlow, "CTE_RHIvslow","Group_de1",genes)
significance_table[3,2] <- length(which(CTE_RHIvslow_res_all$P.Value < 0.05))
significance_table[3,3] <- length(which(CTE_RHIvslow_res_all$adj.P.Val < 0.05))
CTE_RHIvslow_fgsea_res <- fgsea_function(CTE_RHIvslow_res_all,C2_genesets,"CTE_RHIvslow")
significance_table[3,4] <- length(which(CTE_RHIvslow_fgsea_res$pval < 0.05))
significance_table[3,5] <- length(which(CTE_RHIvslow_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + CDStot
adat_cds <- adat[,which(!is.na(colData(adat)$CDStot))] #133
genes <- data.frame(unique_names = rownames(rowData(adat_cds)), genes = rowData(adat_cds)$EntrezGeneSymbol)
cds_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + CDStot, data=colData(adat_cds))
cds_res_all <- limma_model(cds_model,adat_cds, "CDStot","CDStot",genes)
significance_table[5,2] <- length(which(cds_res_all$P.Value < 0.05))
significance_table[5,3] <- length(which(cds_res_all$adj.P.Val < 0.05))
cds_fgsea_res <- fgsea_function(cds_res_all,C2_genesets,"CDStot")
significance_table[5,4] <- length(which(cds_fgsea_res$pval < 0.05))
significance_table[5,5] <- length(which(cds_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + Dementia
adat_dementia <- adat[,which(!is.na(colData(adat)$DementiaHx))] #205
genes <- data.frame(unique_names = rownames(rowData(adat_dementia)), genes = rowData(adat_dementia)$EntrezGeneSymbol)
dementia_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + DementiaHx, data=colData(adat_dementia))
dementia_res_all <- limma_model(dementia_model,adat_dementia, "DementiaHx","DementiaHx1",genes)
significance_table[6,2] <- length(which(dementia_res_all$P.Value < 0.05))
significance_table[6,3] <- length(which(dementia_res_all$adj.P.Val < 0.05))
dementia_fgsea_res <- fgsea_function(dementia_res_all,C2_genesets,"DementiaHx")
significance_table[6,4] <- length(which(dementia_fgsea_res$pval < 0.05))
significance_table[6,5] <- length(which(dementia_fgsea_res$padj < 0.05))

#~ age_at_death + PathAD + PathLBD + PathFTD + faqtot
adat_faq <- adat[,which(!is.na(colData(adat)$faqtot))] #179
genes <- data.frame(unique_names = rownames(rowData(adat_faq)), genes = rowData(adat_faq)$EntrezGeneSymbol)
faq_model <- model.matrix(~ agedeath + PathAD + PathLBD + PathFTD + PMI + faqtot, data=colData(adat_faq))
faq_res_all <- limma_model(faq_model,adat_faq, "faqtot","faqtot",genes)
significance_table[7,2] <- length(which(faq_res_all$P.Value < 0.05))
significance_table[7,3] <- length(which(faq_res_all$adj.P.Val < 0.05))
faq_fgsea_res <- fgsea_function(faq_res_all,C2_genesets,"faqtot")
significance_table[7,4] <- length(which(faq_fgsea_res$pval < 0.05))
significance_table[7,5] <- length(which(faq_fgsea_res$padj < 0.05))

#number of samples included in each model
num_table <- matrix(data=NA,nrow = 7,ncol = 2)
num_table <- as.data.frame(num_table)
colnames(num_table) <- c("Model","Number of Samples Included Out of 205")
num_table$Model <- c("AT8 Total", "Total Years of Play", "RHI to Low CTE", 
                     "RHI to High CTE", "Cognitive Difficulty Score Total", 
                     "Dementia", "Functional Activities Questionnaire Total")
model_list <- list(adat_AT8, adat_totyrs, adat_RHIlow, adat_highrhi, adat_cds, 
                   adat_dementia, adat_faq)
for (i in 1:nrow(num_table)){
  num_table$'Number of Samples Included Out of 205'[i] <- ncol(model_list[[i]])
}

final_table <- cbind(
  Model = significance_table$Model,
  'Number of Samples Included Out of 204' = num_table$'Number of Samples Included Out of 205',
  significance_table[, c(2,4:ncol(significance_table))] # Columns from the target index onwards
)
final_table <- final_table[-nrow(final_table), ]
final_table %>%
  gt()
gtsave(sig_table, "/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/significant_limma_Proteins_Pathways_shortenedtable_noF.png")
gtsave(sig_table, "/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/significant_limma_Proteins_Pathways_shortenedtable_noF.html")
write.csv(significance_table, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/significant_limma_Proteins_Pathways_shortenedtable_noF.csv'))

#Stage 1 - 9
#Stage 2 - 18
#Stage 3 - 42
#Stage 4 - 56
#RHI - 54

#full table
num_table <- matrix(data=NA,nrow = 12,ncol = 2)
num_table <- as.data.frame(num_table)
colnames(num_table) <- c("Model","Number of Samples Included Out of 205")
num_table$Model <- c("AT8 Total", "Total Years of Play", "RHI to Low CTE", 
                     "RHI to High CTE", "Cognitive Difficulty Score Total", 
                     "Dementia", "Functional Activities Questionnaire Total", 
                     "RHI vs Stage 1", "RHI vs Stage 2", "RHI vs Stage 3", 
                     "RHI vs Stage 4", "CTE low vs CTE high")
model_list <- list(adat_AT8, adat_totyrs, adat_RHIlow, adat_highrhi, adat_cds, 
                   adat_dementia, adat_faq, adat_rhi1, adat_rhi2, adat_rhi3, 
                   adat_rhi4, adat_highlow)
for (i in 1:nrow(num_table)){
  num_table$'Number of Samples Included Out of 205'[i] <- ncol(model_list[[i]])
}

final_table <- rbind(num_table,significance_table)
num_table <- num_table %>%
  gt()
write.csv(num_table, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/limma_model_num_samps_fulltable_noF.csv'))
gt_two_column_layout(list(sum_table, sig_table))
#main model table
num_table <- matrix(data=NA,nrow = 7,ncol = 2)
num_table <- as.data.frame(num_table)
colnames(num_table) <- c("Model","Number of Samples Included Out of 205")
num_table$Model <- c("AT8 Total", "Total Years of Play", "RHI to Low CTE", 
                     "RHI to High CTE", "Cognitive Difficulty Score Total", 
                     "Dementia", "Functional Activities Questionnaire Total")
model_list <- list(adat_AT8, adat_totyrs, adat_RHIlow, adat_highrhi, adat_cds, 
                   adat_dementia, adat_faq)
for (i in 1:nrow(num_table)){
  num_table$'Number of Samples Included Out of 205'[i] <- ncol(model_list[[i]])
}
num_table %>%
  gt()
write.csv(num_table, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results/PMI_imp1/limma_model_num_samps_maintable_noF.csv'))


