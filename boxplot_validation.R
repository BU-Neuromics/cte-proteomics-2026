#boxplot spot checks to look at progression of CTE
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(gt)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")

###############################################################################
#Single Protein Checks - can use function to produce a boxplot of interest
#FGG spot check
protein_index <- which(rownames(adat) == "seq.4989.7")

#ACHE spot check
protein_index <- which(rownames(adat) == "seq.10980.11")
#protein_index <- which(rownames(adat) == "seq.15553.22") - not good

#UBE2C 
protein_index <- which(rownames(adat) == "seq.12556.7")
protein_index <- which(rownames(adat) == "seq.20142.42")

#MITD1
protein_index <- which(rownames(adat) == "seq.23266.62")

#CAPNS2
protein_index <- which(rownames(adat) == "seq.22063.20")

################################################################################
#List of Proteins Check for Dillon - prep for a list of proteins to produce many boxplots in loop
#NRF2 does not exist in proteomics data set
box_plot_table <- matrix(data=NA,nrow = 23,ncol = 3)
box_plot_table <- as.data.frame(box_plot_table)
colnames(box_plot_table) <- c("protein_names","seqID_names","protein_index")
box_plot_table$protein_names <- c("PSMA5","PSMA7","PSMB1","PSMB2","PSMB3","RPN1","SMURF1","KEAP1",
                   "GAN","KCTD6","FBXL4","KLHL13","USP10","UBQLN2","MITD1",
                   "UBB","UBC","IGFBP5","NPTX2","C3","C5","C9","CFB")
for(i in 1:length(box_plot_table$protein_names)){
  print(i)
  print(box_plot_table$protein_names[i])
  box_plot_table$protein_index[i] <- which(rowData(adat)$EntrezGeneSymbol == box_plot_table$protein_names[i])
  print(box_plot_table$protein_index[i])
  box_plot_table$seqID_names[i] <- rownames(rowData(adat))[box_plot_table$protein_index[i]]
}

box_plot_table
################################################################################
#loop to create many boxplots at once
for(i in 1:length(box_plot_table$protein_index)){
  index <- box_plot_table$protein_index[i]
  protein_counts <- assays(adat)$counts[index,]
  #grab CTE data
  protein_CTE <- colData(adat)$CTEStage
  #check lengths 
  length(protein_counts) == length(protein_CTE)
  counts_CTE <- data.frame(counts = protein_counts, CTE = protein_CTE)
  write.csv(counts_CTE, paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',box_plot_table$protein_names[i],"_",box_plot_table$seqID_names[i],'_CTE_followup_file.csv'))
  #create boxplots
  CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste0('Boxplot of CTE and ',box_plot_table$protein_names[i])) +
    xlab("CTE")
  CTE_box
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',box_plot_table$protein_names[i],"_",box_plot_table$seqID_names[i],'_CTE_followup_boxplot.png'), plot = CTE_box, width = 8, height = 6)
}
#finalized table in a nice form
box_plot_table %>%
  gt()
################################################################################
#function to look at any boxplot
protein_counts <- assays(adat)$counts[protein_index,]
#grab CTE data
protein <- colData(adat)$CTEStage
#check lengths 
length(protein_counts) == length(protein)
counts_CTE <- data.frame(counts = protein_counts, CTE = protein)
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE and CAPNS2')) +
  xlab("CTE")
CTE_box