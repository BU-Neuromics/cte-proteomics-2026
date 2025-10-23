#Helen Pennington
#Labadorf Lab 
#Started August 19, 2025
#RNA Sequencing Validation of SomaScan Proteomics Results

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
library(qusage)
library(fgsea)

################################################################################
#Data Inputs
################################################################################

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
C2_genesets <- read.gmt("/restricted/projectnb/cteseq/projects/somascan/data/c2.cp.v2025.1.Hs.symbols.gmt")
#61533 transcripts vs 7285 proteins

################################################################################
#Data Filtering and Preparation
################################################################################

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
rna_meta <- as.data.frame(rna_meta)
rna_meta$gene_name <- rownames
rna_meta$gene_name <- sub("\\..*", "", rna_meta$gene_name)  # Keep only the main ID
length(rna_meta$gene_name)
edb <- EnsDb.Hsapiens.v79
# Convert Gene IDs to Gene Names (5071 geneids did not have gene names - num=56462)
rna_meta$gene_name <- mapIds(edb, keys = rna_meta$gene_name, keytype = "GENEID", column = "GENENAME", multiVals = "first")
rna_meta <- rna_meta[which(!is.na(rna_meta$gene_name)),]
head(rna_meta)

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

#merge adat_df and rna_meta to only look at genes who's proteins are evaluated in adat
overlap_genes <- intersect(rna_meta$gene_name, adat_df$geneName)
rna_counts_filtered <- rna_meta[rna_meta$gene_name %in% overlap_genes, ]
genes <- data.frame(gene_name = rna_counts_filtered$gene_name, gene_id = rownames(rna_counts_filtered))
rna_counts_filtered$gene_name <- NULL

#prep metadata
rownames(metadata) <- metadata$Core_ID
metadata <- metadata[which(rownames(metadata) %in% colnames(rna_meta)),]
rownames(metadata) <- metadata$Core_ID
metadata <- metadata[order(rownames(metadata)),]

#removing samples that were not used in somascan analyses
#adat <- adat[,which(colData(adat)$BBID != "K-0666")]
#metadata <- metadata[which(metadata$SampleName %in% colData(adat)$BBID),]
#rna_meta <- rna_meta[,which(colnames(rna_meta) %in% metadata$Core_ID)]
rownames(metadata) <- metadata$Core_ID
all(rownames(metadata) == colnames(rna_counts_filtered))
rna_meta <- round(rna_counts_filtered)

#prep
metadata_pp <- metadata[which(!is.na(metadata$PathAD)),] 
rna_meta_pp <- rna_meta[,which(!is.na(metadata$PathAD))]
metadata_p <- metadata_pp[which(!is.na(metadata_pp$PathFTD)),] 
rna_meta_p <- rna_meta_pp[,which(!is.na(metadata_pp$PathFTD))]
metadata_p$PathLBD <- as.integer(ifelse(metadata_p$PathLBD == 2, 1, 0))
metadata_p$agedeath <- scale(as.numeric(metadata_p$agedeath))
metadata_p$PathAD   <- factor(metadata_p$PathAD)
metadata_p$PathLBD  <- factor(metadata_p$PathLBD)
metadata_p$PathFTD  <- factor(metadata_p$PathFTD)
rownames(metadata_p) <- metadata_p$Core_ID
all(rownames(metadata_p) == colnames(rna_meta_p))
################################################################################
#functions
################################################################################

#volcano plot
volcano_plot <- function(res, name_title, name_path, genes){
  #volcano plot
  resLFC_df <- as.data.frame(res)
  resLFC_df <- merge(resLFC_df, genes, by.x = 0, by.y = "gene_id")
  resLFC_df$threshold <- resLFC_df$padj < 0.05
  
  top_genes <- resLFC_df[order(resLFC_df$padj), ][1:10, ]
  
  deseq2_volcano <- ggplot(resLFC_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = threshold), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text(data = top_genes, aes(label = gene_name), vjust = 1.5, size = 3) +
    theme_minimal() +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste("Differential Expression Analysis of RNA-Sequencing Data for ",name_title, " Model"))
  ggsave(paste('/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_',name_path,'_no_shrink_volcano_plot.png'), plot = deseq2_volcano, width = 8, height = 6)
  deseq2_volcano
  resLFC_df <- resLFC_df[which(!is.na(resLFC_df$log2FoldChange)),]
  write.csv(resLFC_df, paste('/restricted/projectnb/cteseq/projects/somascan/results/deseq2/deseq2_results_',name_path,'_no_shrink.csv'))
  return(resLFC_df)
}

#fgsea
fgsea_function <- function(resLFC_df,C2_genesets,name){
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
  my_genes <- unique(resLFC_df$gene_name)
  #C2_genesets <- split(C2_genesets$gene, C2_genesets$term)
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
  genes_sorted <- resLFC_df[order(abs(resLFC_df$log2FoldChange), decreasing = TRUE), ]
  genes_sorted_unique <- genes_sorted[!duplicated(genes_sorted$gene_name), ]
  plot(genes_sorted_unique$log2FoldChange)
  stats <- genes_sorted_unique$log2FoldChange
  names(stats) <- genes_sorted_unique$gene_name
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
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseaPlots/rna_seq/deseq2_',name,'_no_shrink_gseaPlot.png'), plot = gseaTable, width = 8, height = 6)
  final_results <- as_tibble(fgseaRes)
  final_results <- final_results %>% unnest(leadingEdge)
  write.csv(final_results, paste0('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/rna_seq/deseq2_',name,'_no_shrink_gseaResults.csv'))
  return(fgseaRes)
}

################################################################################
#models
################################################################################

#AT8 model
metadata_AT8 <- metadata_p[which(!is.na(metadata_p$AT8_total)),] 
rna_meta_AT8 <- rna_meta_p[,which(!is.na(metadata_p$AT8_total))]
metadata_AT8$AT8_total <- scale(as.numeric(metadata_AT8$AT8_total))
rownames(metadata_AT8) <- metadata_AT8$Core_ID
all(rownames(metadata_AT8) == colnames(rna_meta_AT8))
metadata_AT8 <- dplyr::select(as.data.frame(metadata_AT8), PathAD, PathLBD, PathFTD, AT8_total, agedeath, Core_ID)
apply(metadata_AT8, 2, function(x) length(unique(x)))
ddsA <- DESeqDataSetFromMatrix(countData = rna_meta_AT8,
                               colData = metadata_AT8,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + AT8_total)
#vst normalization and pca plot
vst_rna <- vst(ddsA, blind = TRUE)
plotPCA(vst_rna, intgroup = "AT8_total")

dds <- DESeq(ddsA)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="AT8_total")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="AT8_total", type="apeglm")
resLFC_df <- volcano_plot(res, "AT8 Total", "AT8_total",genes)
fgsea_res_AT8 <- fgsea_function(resLFC_df,C2_genesets,"AT8_total")
saveRDS(fgsea_res_AT8, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_AT8_rnaseq.rds")

#totyrs model
metadata_yrs <- metadata_p[which(!is.na(metadata_p$totyrs)),] 
rna_meta_yrs <- rna_meta_p[,which(!is.na(metadata_p$totyrs))]
metadata_yrs$totyrs <- scale(as.numeric(metadata_yrs$totyrs))
rownames(metadata_yrs) <- metadata_yrs$Core_ID
all(rownames(metadata_yrs) == colnames(rna_meta_yrs))
metadata_yrs <- dplyr::select(as.data.frame(metadata_yrs), PathAD, PathLBD, PathFTD, totyrs, agedeath, Core_ID)
apply(metadata_yrs, 2, function(x) length(unique(x)))
ddsy <- DESeqDataSetFromMatrix(countData = rna_meta_yrs,
                               colData = metadata_yrs,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + totyrs)
vst_rna <- vst(ddsy, blind = TRUE)
plotPCA(vst_rna, intgroup = "totyrs")
dds <- DESeq(ddsy)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="totyrs")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="totyrs", type="apeglm")
resLFC_df <- volcano_plot(res, "Total Years of Play", "totyrs",genes)
fgsea_res_totyrs <- fgsea_function(resLFC_df,C2_genesets,"totyrs")
saveRDS(fgsea_res_totyrs, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_totyrs_rnaseq.rds")

#dementia model
metadata_dem <- metadata_p[which(!is.na(metadata_p$DementiaHx)),] 
rna_meta_dem <- rna_meta_p[,which(!is.na(metadata_p$DementiaHx))]
metadata_dem$DementiaHx  <- factor(metadata_dem$DementiaHx)
rownames(metadata_dem) <- metadata_dem$Core_ID
all(rownames(metadata_dem) == colnames(rna_meta_dem))
metadata_dem <- dplyr::select(as.data.frame(metadata_dem), PathAD, PathLBD, PathFTD, DementiaHx, agedeath, Core_ID)
apply(metadata_dem, 2, function(x) length(unique(x)))
ddsd <- DESeqDataSetFromMatrix(countData = rna_meta_dem,
                               colData = metadata_dem,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + DementiaHx)
vst_rna <- vst(ddsd, blind = TRUE)
plotPCA(vst_rna, intgroup = "DementiaHx")
dds <- DESeq(ddsd)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="DementiaHx_1_vs_0")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="DementiaHx_1_vs_0", type="apeglm")
resLFC_df <- volcano_plot(res, "Dementia", "DementiaHx",genes)
fgsea_res_dem <- fgsea_function(resLFC_df,C2_genesets,"DementiaHx")
saveRDS(fgsea_res_dem, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_dementia_rnaseq.rds")

#CDS model
metadata_cds <- metadata_p[which(!is.na(metadata_p$CDStot)),] 
rna_meta_cds <- rna_meta_p[,which(!is.na(metadata_p$CDStot))]
metadata_cds$CDStot <- scale(as.numeric(metadata_cds$CDStot))
rownames(metadata_cds) <- metadata_cds$Core_ID
all(rownames(metadata_cds) == colnames(rna_meta_cds))
metadata_cds <- dplyr::select(as.data.frame(metadata_cds), PathAD, PathLBD, PathFTD, CDStot, agedeath, Core_ID)
apply(metadata_cds, 2, function(x) length(unique(x)))
ddsc <- DESeqDataSetFromMatrix(countData = rna_meta_cds,
                               colData = metadata_cds,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + CDStot)
dds <- DESeq(ddsc)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="CDStot")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="CDStot", type="apeglm")
resLFC_df <- volcano_plot(res, "Cognitive Difficulty Score", "CDStot",genes)
fgsea_res_cds <- fgsea_function(resLFC_df,C2_genesets,"CDStot")
saveRDS(fgsea_res_cds, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_cds_rnaseq.rds")

#faq model
metadata_faq <- metadata_p[which(!is.na(metadata_p$faqtot)),] 
rna_meta_faq <- rna_meta_p[,which(!is.na(metadata_p$faqtot))]
metadata_faq$totyrs <- scale(as.numeric(metadata_faq$faqtot))
rownames(metadata_faq) <- metadata_faq$Core_ID
all(rownames(metadata_faq) == colnames(rna_meta_faq))
metadata_faq <- dplyr::select(as.data.frame(metadata_faq), PathAD, PathLBD, PathFTD, faqtot, agedeath, Core_ID)
apply(metadata_faq, 2, function(x) length(unique(x)))
ddsf <- DESeqDataSetFromMatrix(countData = rna_meta_faq,
                               colData = metadata_faq,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + faqtot)
dds <- DESeq(ddsf)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="faqtot")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="faqtot", type="apeglm")
resLFC_df <- volcano_plot(res, "faqtot", "faqtot",genes)
fgsea_res_faq <- fgsea_function(resLFC_df,C2_genesets,"faqtot")
saveRDS(fgsea_res_faq, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_faq_rnaseq.rds")

#rhi vs low cte model
metadata_rl <- metadata_p[which(!(metadata_p$CTEStage) == 3 & !(metadata_p$CTEStage) == 4),] #90
rna_meta_rl <- rna_meta_p[,which(!(metadata_p$CTEStage) == 3 & !(metadata_p$CTEStage) == 4)]
metadata_rl$Group_de <- as.integer(ifelse(metadata_rl$CTEStage == 1, 1, 
                                          ifelse(metadata_rl$CTEStage == 2, 1, 
                                                 ifelse(metadata_rl$CTEStage == 0, 0,
                                                               metadata_rl$Group_de))))
metadata_rl <- metadata_rl[which(!is.na(metadata_rl$Group_de)),] 
rna_meta_rl <- rna_meta_rl[,which(!is.na(metadata_rl$Group_de))]
metadata_rl$Group_de  <- factor(metadata_rl$Group_de)
rownames(metadata_rl) <- metadata_rl$Core_ID
all(rownames(metadata_rl) == colnames(rna_meta_rl))
metadata_rl <- dplyr::select(as.data.frame(metadata_rl), PathAD, PathLBD, PathFTD, Group_de, agedeath, Core_ID)
apply(metadata_rl, 2, function(x) length(unique(x)))
ddsr <- DESeqDataSetFromMatrix(countData = rna_meta_rl,
                               colData = metadata_rl,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + Group_de)
dds_rl <- DESeq(ddsr)
resultsNames(dds_rl) # lists the coefficients
res_rl <- results(dds_rl, name="Group_de_1_vs_0")
# or to shrink log fold changes association with condition:
#res_rl <- lfcShrink(dds_rl, coef="Group_de_1_vs_0", type="apeglm")
resLFC_rl <- volcano_plot(res_rl, "RHI vs CTE Low", "rhivslow",genes)
fgsea_res_rl <- fgsea_function(resLFC_rl,C2_genesets,"rhivslow")
saveRDS(fgsea_res_rl, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_RHIvslow_rnaseq.rds")

CTE_rl_paths <- fgsea_res_rl[which(fgsea_res_rl$NES > 0 & fgsea_res_rl$pval < 0.05),]
leading_proteins_rl <- sort(table(unlist(CTE_rl_paths$leadingEdge)), decreasing = TRUE)
length(leading_proteins_rl)

#low vs high cte model
metadata_hl <- metadata_p[which(!(metadata_p$CTEStage) == 0),] #146
rna_meta_hl <- rna_meta_p[,which(!(metadata_p$CTEStage) == 0)]
metadata_hl$Group_de <- as.integer(ifelse(metadata_hl$CTEStage == 1, 0, 
                                          ifelse(metadata_hl$CTEStage == 2, 0, 
                                                 ifelse(metadata_hl$CTEStage == 3, 1, 
                                                        ifelse(metadata_hl$CTEStage == 4, 1, 
                                                               metadata_hl$Group_de)))))
metadata_hl <- metadata_hl[which(!is.na(metadata_hl$Group_de)),] 
rna_meta_hl <- rna_meta_hl[,which(!is.na(metadata_hl$Group_de))]
metadata_hl$Group_de <- factor(metadata_hl$Group_de)
rownames(metadata_hl) <- metadata_hl$Core_ID
all(rownames(metadata_hl) == colnames(rna_meta_hl))
metadata_hl <- dplyr::select(as.data.frame(metadata_hl), PathAD, PathLBD, PathFTD, Group_de, agedeath, Core_ID)
apply(metadata_hl, 2, function(x) length(unique(x)))
ddsh <- DESeqDataSetFromMatrix(countData = rna_meta_hl,
                               colData = metadata_hl,
                               design= ~ agedeath + PathAD + PathLBD + PathFTD + Group_de)
dds_lh <- DESeq(ddsh)
resultsNames(dds_lh) # lists the coefficients
res_lh <- results(dds_lh, name="Group_de_1_vs_0")
# or to shrink log fold changes association with condition:
#res_lh <- lfcShrink(dds_lh, coef="Group_de_1_vs_0", type="apeglm")
resLFC_lh <- volcano_plot(res_lh, "Low CTE vs High CTE", "lowvshigh",genes)
fgsea_res_lh <- fgsea_function(resLFC_lh,C2_genesets,"lowvshigh")
saveRDS(fgsea_res_lh, file = "/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_lowvshigh_rnaseq.rds")

CTE_lh_paths <- fgsea_res_lh[which(fgsea_res_lh$NES < 0 & fgsea_res_lh$pval < 0.05),]
leading_proteins_lh <- sort(table(unlist(CTE_lh_paths$leadingEdge)), decreasing = TRUE)
length(leading_proteins_lh)


