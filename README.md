# hp_rotation
# SomaScan_totyrs_AT8_Spearman.R 
preliminary analysis using SomaScan workflow to determine spearman correlations based on total years of play and AT8 total. Analyses include all 212 samples, not removing outlier samples or samples not included in metadata. Proteins are also not fully filtered. Outputs here: /restricted/projectnb/cteseq/projects/somascan/results/spearman

# Creating_SummarizedExperiment_hep.R 
Adat data fully filtered and converted to summarized experiment format. All following analyses done using this data file found here: /restricted/projectnb/cteseq/projects/somascan/data and pca plot to identify outliers done in this file and found here: /restricted/projectnb/cteseq/projects/somascan/results

# identifying_adat_duplicates_hep.R 
Idenitfies duplicates in SomaScan protein data and produces barplot and table. Outputs here: /restricted/projectnb/cteseq/projects/somascan/results/lm

# limma_proteinDEA_hep.R
Runs limma DEA for AT8_total + age_at_death and spot checks top significant proteins and a top significant gene from previous analysis to verify accuracy of our method before expanding to other models. Outputs here: /restricted/projectnb/cteseq/projects/somascan/results/limma/AT8_spotchecks

# limma_function_hep.R
Runs limma DEA and fgsea for all seven models and produces final table with number of significant proteins and pathways for each model. Output files/table: /restricted/projectnb/cteseq/projects/somascan/results/limma/limma_results Output plots: /restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots Output fgsea plots and files: /restricted/projectnb/cteseq/projects/somascan/results/fgsea

# lm_protein_rna_hep.R
Runs lm model for each gene_id vs each SeqID for each gene in SomaScan and RNA-sequencing datasets. Final file with all beta value outputs, histograms of beta value distributions, follow-up scatterplots of significant/relavent genes, scatterplot of median RNA level vs beta1 (RNA comparison), and duplicate geneID and SeqID analysis are all here. Outputs here: /restricted/projectnb/cteseq/projects/somascan/results/lm
