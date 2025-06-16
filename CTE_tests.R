#Investigation of tests
#Helen Pennington
#Started June 5, 2025

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
library(Metrics)

# Read in the adat file and metadata file
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
rna <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/all_counts.csv")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan/data/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
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
colData(adat)$Group_de <- ifelse(colData(adat)$CTEStage == 1, "low", ifelse(colData(adat)$CTEStage == 2, "low", ifelse(colData(adat)$CTEStage == 3, "high", ifelse(colData(adat)$CTEStage == 4, "high", 0))))
colData(adat)$hitsperyear <- colData(adat)$chii_g / colData(adat)$totyrs

#test
CDS <- head(colData(adat)$CDStot) #46
adat_tests <- adat[,which(!is.na(colData(adat)$CDStot))]
FAQ <- colData(adat)$faqtot #54
adat_tests <- adat_tests[,which(!is.na(colData(adat_tests)$faqtot))]
dem <- colData(adat)$DementiaHx #52
length(which(!is.na(colData(adat_tests)$DementiaHx)))
agg <- colData(adat)$maxaggsum #45
length(which(!is.na(colData(adat_tests)$maxaggsum)))
adat_tests <- adat_tests[,which(!is.na(colData(adat_tests)$maxaggsum))]
adat_tests <- adat_tests[,which(!is.na(colData(adat_tests)$AT8_total))]
adat_tests <- adat_tests[,which(!is.na(colData(adat_tests)$totyrs))]
adat_tests_df <- as.data.frame(colData(adat)[c(40,45,46,52,54,100,150,156)])
adat_tests2_df <- as.data.frame(colData(adat)[c(49,157,40,156, 150,100)])
adat_tests2_df$Group_de <- as.integer(ifelse(adat_tests2_df$CTEStage == 1, 2, ifelse(adat_tests2_df$CTEStage == 2, 2, ifelse(adat_tests2_df$CTEStage == 3, 3, ifelse(adat_tests2_df$CTEStage == 4, 3, ifelse(adat_tests2_df$CTEStage == 0, 1, adat_tests2_df$Group_de))))))

#removing NAs brough samples from 207 to 98
pairs(adat_tests_df[c(1:6,8)], main = "CTE Proteasome Metadata Exploration",
      pch = 21, bg = c("grey","lightblue","blue")[adat_tests_df$Group_de])
adat_tests2_df <- adat_tests2_df[which(!is.na(adat_tests2_df$Group_de)),]
pairs(adat_tests2_df[c(1:4)], main = "CTE head impacts and years of play exploration",
      pch = 21, bg = c("grey","lightblue","blue")[adat_tests2_df$Group_de])
#comparisons
comp_1 <- "AT8_total"
comp_2 <- "totyrs"
name_1 <- "Tau Level"
name_2 <- "Total Years of Play"
comparisons_lm <- function(adat_res,comp_1, comp_2, name_1, name_2){
  col_data <- colData(adat_res)
  model <- lm(col_data[[comp_1]] ~ col_data[[comp_2]])
  plot <- ggplot(data.frame(col_data[[comp_1]], col_data[[comp_2]]), aes(x = col_data[[comp_1]], y = col_data[[comp_2]])) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue") +
    labs(title = paste0("Comparison of ", name_1, " and ", name_2), x = paste0(name_1), y = paste0(name_2))
  ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/modeling/regression_plot_',comp_1,'_',comp_2,'.png'),
         plot = last_plot(),width = 8, height = 5)
  print(paste("slope: ",model$coefficients[2]))
  model.res = resid(model)
  plot(col_data[[comp_2]], model.res,
       ylab="Residuals", xlab=paste0(name_2),
       main=paste0(name_1, "and ", name_2)) 
  abline(0, 0)
  print(paste("r-squared: ",summary(model)$r.squared))
  print(paste("rmse: ",rmse(col_data[[comp_1]], col_data[[comp_2]])))
  return(plot)
}

#tau and totyrs
adat_res <- adat[,which(!is.na(colData(adat)$AT8_total))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$totyrs))]
lm_plot <- comparisons_lm(adat_res,"AT8_total","totyrs","Tau Level","Total Years of Play")
lm_plot
#tau and CTE
adat_res <- adat[,which(!is.na(colData(adat)$AT8_total))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$CTEStage))]
lm_plot <- comparisons_lm(adat_res,"AT8_total","CTEStage","Tau Level","CTE Stage")
lm_plot
mean(colData(adat[,which(colData(adat)$CTEStage == 0)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$AT8_total, na.rm = TRUE)

#AT8 and CTE
#grab CTE data
protein_CTE <- colData(adat_res)$CTEStage
AT8_counts <- colData(adat_res)$AT8_total
#check lengths 
length(AT8_counts) == length(protein_CTE)
counts_CTE <- data.frame(counts = AT8_counts, CTE = protein_CTE)
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE Stage and AT8_total')) +
  xlab("CTE Stage")
CTE_box
ggsave(paste0('/restricted/projectnb/cteseq/projects/somascan/results/limma/CTE_followup/',
              protein_index[i,1],protein_index[i,2],'_CTEstage_followup_boxplot.png'), 
       plot = CTE_box, width = 8, height = 6)

#tau and hits/totyrs
adat_res <- adat[,which(!is.na(colData(adat)$AT8_total))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$hitsperyear))]
lm_plot <- comparisons_lm(adat_res,"AT8_total","hitsperyear","Tau Level","Average Hits per Year")
lm_plot

adat_res <- adat[,which(!is.na(colData(adat)$AT8_total))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$hitsperyear))]
lm_plot <- comparisons_lm(adat_res,"AT8_total","hitsperyear","Tau Level","Average Hits per Year")
lm_plot
mean(colData(adat[,which(colData(adat)$CTEStage == 0)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$AT8_total, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$AT8_total, na.rm = TRUE)


#totyrs and CTE
adat_res <- adat[,which(!is.na(colData(adat)$totyrs))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$CTEStage))]
lm_plot <- comparisons_lm(adat_res,"totyrs","CTEStage","Total Years of Play","CTE Stage")
lm_plot
mean(colData(adat[,which(colData(adat)$CTEStage == 0)])$totyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$totyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$totyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$totyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$totyrs, na.rm = TRUE)
#grab CTE data
protein_CTE <- colData(adat_res)$CTEStage
totyrs_counts <- colData(adat_res)$totyrs * colData(adat_res)$agedeath
#check lengths 
length(totyrs_counts) == length(protein_CTE)
counts_CTE <- data.frame(counts = totyrs_counts, CTE = protein_CTE)
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE Stage and Total Years of Play')) +
  xlab("CTE Stage")
CTE_box

lm_plot <- comparisons_lm(adat_res,"totyrs","agedeath","Total Years of Play","Age of Death")
lm_plot
plot(colData(adat)$totyrs)
adat_totyrs_10 <- adat[,which(colData(adat)$totyrs == 10)]
plot(colData(adat_totyrs_10)$agedeath)
colData(adat_totyrs_10)$BBID
colData(adat_totyrs_10)$footyrs


#footyrs and CTE
adat_res <- adat[,which(!is.na(colData(adat)$footyrs))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$CTEStage))]
lm_plot <- comparisons_lm(adat_res,"footyrs","CTEStage","Total Years of Football","CTE Stage")
lm_plot
mean(colData(adat[,which(colData(adat)$CTEStage == 0)])$footyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$footyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$footyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$footyrs, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$footyrs, na.rm = TRUE)
#grab CTE data
protein_CTE <- colData(adat_res)$CTEStage
footyrs_counts <- colData(adat_res)$footyrs * colData(adat_res)$agedeath
#check lengths 
length(footyrs_counts) == length(protein_CTE)
counts_CTE <- data.frame(counts = footyrs_counts, CTE = protein_CTE)
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE Stage and Total Years of Football')) +
  xlab("CTE Stage")
CTE_box

mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$agedeath - colData(adat[,which(colData(adat)$CTEStage == 1)])$AFE, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$Group_de == "low")])$agedeath - colData(adat[,which(colData(adat)$Group_de == "low")])$AFE, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$agedeath - colData(adat[,which(colData(adat)$CTEStage == 2)])$AFE, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$agedeath - colData(adat[,which(colData(adat)$CTEStage == 3)])$AFE, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$agedeath - colData(adat[,which(colData(adat)$CTEStage == 4)])$AFE, na.rm = TRUE)


#cognition and CTE
adat_res <- adat[,which(!is.na(colData(adat)$CDStot))]
adat_res <-adat_res[,which(!is.na(colData(adat_res)$CTEStage))]
lm_plot <- comparisons_lm(adat_res,"CDStot","CTEStage","Cognitive Difficulty Scale","CTE Stage")
lm_plot
mean(colData(adat[,which(colData(adat)$CTEStage == 0)])$CDStot, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 1)])$CDStot, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 2)])$CDStot, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 3)])$CDStot, na.rm = TRUE)
mean(colData(adat[,which(colData(adat)$CTEStage == 4)])$CDStot, na.rm = TRUE)
#grab CTE data
protein_CTE <- colData(adat_res)$CTEStage
CDStot_counts <- colData(adat_res)$CDStot 
#check lengths 
length(CDStot_counts) == length(protein_CTE)
counts_CTE <- data.frame(counts = CDStot_counts, CTE = protein_CTE)
#create boxplots
CTE_box <- ggplot(counts_CTE, aes(x=CTE, y=counts, group = CTE)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(paste0('Boxplot of CTE Stage and Cognitive Difficulty Scale')) +
  xlab("CTE Stage")
CTE_box