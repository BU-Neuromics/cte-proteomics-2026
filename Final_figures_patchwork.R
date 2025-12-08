#Final Plots with Patchwork
#Helen Pennington
#December 8, 2025

#packages
library(patchwork)

#volcano plots
vol_AT8 <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_AT8_total_ggplot_volcanoplot.rds")
vol_totyrs <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_totyrs_ggplot_volcanoplot.rds")
vol_low <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_CTE_RHIvslow_ggplot_volcanoplot.rds")
vol_high <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_CTE_RHIvsHigh_ggplot_volcanoplot.rds")
vol_dem <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_DementiaHx_ggplot_volcanoplot.rds")
vol_cds <- readRDS("/restricted/projectnb/cteseq/projects/somascan/results/limma/volcano_plots/PMI_imp1/limma_DEA_CDStot_ggplot_volcanoplot.rds")

#leading edge
le_cte <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/CTE_main_models_leading_edge.rds")
le_AT8_totyrs <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/AT8_main_totyrs_leading_edge.rds")
le_totyrs_AT8 <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/totyrs_main_AT8_leading_edge.rds")
le_dem_cds <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/dementia_main_cds_leading_edge.rds")
le_cds_dem <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_leading_edge_plots/cds_main_dementia_leading_edge.rds")

#fgsea_comp
fgsea_cte <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CTE_main_models_fgsea_comp.rds")
fgsea_AT8_totyrs <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/AT8_totyrs_fgsea_comp.rds")
fgsea_dem_cds <- readRDS("/restricted/projectnb/cteseq/projects/somascan/final_plots/publication_ready_fgsea_comp_plots/CDS_dementia_fgsea_comp.rds")

#Figure 1
tag_theme <- theme(
  plot.tag = element_text(size = 30, face = "bold"),
  plot.tag.position = c(0, 1),   # top-left inside
  plot.tag.hjust = 0,
  plot.tag.vjust = 1
)

fig1 <- (
  (vol_AT8 | vol_totyrs) /
    le_AT8_totyrs /
    le_totyrs_AT8 /
    fgsea_AT8_totyrs
) +
  plot_annotation(tag_levels = "A") &
  tag_theme
pdf("/restricted/projectnb/cteseq/projects/somascan/final_plots/Patchworks/patchwork_output_fig1.pdf", width = 20, height = 20)
print(fig1)
dev.off()

#Figure 2
tag_theme <- theme(
  plot.tag = element_text(size = 30, face = "bold"),
  plot.tag.position = c(0, 1),   # top-left inside
  plot.tag.hjust = 0,
  plot.tag.vjust = 1
)

fig2 <- (
  (vol_low | vol_high) /
    le_cte /
    fgsea_cte
) +
  plot_annotation(tag_levels = "A") &
  tag_theme
pdf("/restricted/projectnb/cteseq/projects/somascan/final_plots/Patchworks/patchwork_output_fig2.pdf", width = 20, height = 20)
print(fig2)
dev.off()

#Figure 3
tag_theme <- theme(
  plot.tag = element_text(size = 30, face = "bold"),
  plot.tag.position = c(0, 1),   # top-left inside
  plot.tag.hjust = 0,
  plot.tag.vjust = 1
)

fig3 <- (
  (vol_dem | vol_cds) /
    le_dem_cds /
    le_cds_dem/
    fgsea_dem_cds
) +
  plot_annotation(tag_levels = "A") &
  tag_theme
pdf("/restricted/projectnb/cteseq/projects/somascan/final_plots/Patchworks/patchwork_output_fig3.pdf", width = 20, height = 20)
print(fig3)
dev.off()

