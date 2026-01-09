# CTE_proteomics
# Creating_SummarizedExperiment_hep.R 
Adat data fully filtered and converted to summarized experiment format. All following analyses done using this data file found here: /restricted/projectnb/cteseq/projects/somascan/data and pca plot to identify outliers done in this file and found here: /restricted/projectnb/cteseq/projects/somascan/results. Summarized experiment file found here: /restricted/projectnb/cteseq/projects/somascan/data

# SomaScan_Amelia_Imputation.R
Imputation on metadata to fill in missing values for Post Mortem Inverval (PMI) and Alzheimer's Disease (PathAD) values. Metadata file found here in the column data of summarized experiment file found here: /restricted/projectnb/cteseq/projects/somascan/data and imputation file found here: /restricted/projectnb/cteseq/projects/somascan/imputation/Stable_Releases/2025-12-02

# limma_fgsea_functions_CTE_proteomics.R
Imputed values filled in, limma analyses for all models, and fgsea analyses for all models. File with all limma analysis results and file with all fgsea analysis results found here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files. Volcano plots detailing limma results for all models found here: /restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1 

# fgsea_grouping.R
fgsea significant pathways for all models hierarchically clustered based on cosine similarity scores of leading edge genes and assigned group numbers using cutreeDynamic. File saved here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files. Visualization of groupings shown in heatmap saved here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper

# fgsea_comparative_plots.R
Manually grouped pathways based on hierarchical clustering and cutreeDynamic groupings was inputted with file saved here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files along with fgsea results file saved in the same folder. Scatterplots created to compare significant pathway group differences between models. Final plots found here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots

# leading_edge_genes.R
Leading edge genes that appeared in significnat pathways with adjusted pvalues less than 0.05 where included in barplots for each model to display how many times the top 100 genes appeared in significant pathways for that model. The bar plots were also color coded to sort them into their hierarchically clustered higher order group. These bar plots were saved here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots

# cytoscape_prep.R
The fgsea results file and hierarchical clustering with manual groups file both found here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/final_files were used to create nodes and edges files to input into Cytoscape also found in the same final_files folder. These files were used to generate a network of all models, their significant pathways, and connect pathways with similar leading edge genes (those having greater than 0.2 cosine similarity scores). 

# Final_figures_patchwork
Three of the five final figures were created by combining several plots using the Patchwork package. The volcano plots generated in limma_fgsea_functions_CTE_proteomics.R and stored here: /restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1, the fgsea comparitive plots generated in fgsea_comparative_plots.R and stored here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_fgsea_comp_plots, and the leading edge plots generated in leading_edge_genes.R and stored here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/publication_ready_leading_edge_plots were used in these figures. Final Patchwork generated figures were saved here: /restricted/projectnb/cteseq/projects/somascan/proteomics_paper/Patchworks
