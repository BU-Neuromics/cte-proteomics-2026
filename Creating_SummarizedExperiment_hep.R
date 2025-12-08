#Coverting adat file to summarized experiment
#Labadorf rotation
#Helen Pennington
#February 13, 2025

#packages 
library(dplyr)
library(GO.db)
library(SomaDataIO)
library(SomaScan.db)
library(data.table)
library(SummarizedExperiment)
library(ggfortify)
library(ggrepel)
library(readxl)

# Read in the adat file and BBID file
adat <- SomaDataIO::read_adat("/restricted/projectnb/cteseq/data/somascan_2024/CTE Somascan/HMS-24-036_2024-08-09/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan/data/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
metadata <- readRDS("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.rds")
metadata_2 <- read_excel("/restricted/projectnb/cteseq/projects/somascan/data/PMI_somascan_decimals.xlsx")
metadata_2$ID_new <- sub("([A-Za-z]+)([0-9]+)", "\\1-\\2", metadata_2$`SOMASCAN ID`)
adat <- log2(adat)
analytes <- SomaDataIO::getAnalytes(adat)
num_metadata <- length(adat) - length(analytes)
analyte_start <- num_metadata+1

#separate into parts
#count data
counts <- as.data.frame(adat[,analyte_start:length(adat)])
dim(counts)
counts <- t(counts)
dim(counts)

#row data
rdata <- as.data.frame(attributes(adat)$Col.Meta)
rownames(rdata) <- attributes(adat)$names[analyte_start:length(adat)]

#column data
cdata <- as.data.frame(adat[,1:num_metadata])
cdata.rownames <- rownames(cdata)
#merge files together to get BBID for each sample
cdata <- dplyr::left_join(cdata,adat_bbid,by=c("SampleGroup"="SampleGroup"))
rownames(cdata) <- cdata.rownames
#only keep samples with BBID
cdata <- cdata[which(!is.na(cdata$BBID)),] #248->212
#merge metadata into cdata
cdata$Identifyer <- rownames(cdata)
which(!cdata$BBID %in% metadata$SampleName)
#only keep samples with metadata
cdata <- dplyr::inner_join(cdata, metadata,by=c("BBID"="SampleName")) #212->210
rownames(cdata) <- cdata$Identifyer

#count data filtering
counts <- counts[,which(colnames(counts)%in%rownames(cdata))]

#check that everything got separated correctly
dim(cdata)
dim(rdata)
dim(counts)

which(rownames(counts) != rownames(rdata))
which(colnames(counts) != rownames(cdata))

#build summarized experiment data frame

counts <- as.matrix(counts)
colData <- DataFrame(cdata)
rowData <- DataFrame(rdata)

adat_summarized_experiment <- SummarizedExperiment(assays=list(counts=counts),
                                                   rowData=rowData, colData=colData)
colData(adat_summarized_experiment)
rowData(adat_summarized_experiment)

#pca plot (at n=212) = not merged with metadata yet
pca_res <- prcomp(t(assays(adat_summarized_experiment)$counts), scale. = TRUE)
scores <- as.data.frame(pca_res$x)
z_pc1 <- scale(scores$PC1)
z_pc2 <- scale(scores$PC2)

outliers <- rownames(scores)[abs(z_pc1) > 3]
outliers


saveRDS(pca_res, file = "/restricted/projectnb/cteseq/projects/somascan/results/pcaplot_212samp_file.rds")
png(filename=paste0('/restricted/projectnb/cteseq/projects/somascan/results/pcaplot_212samp.png'))
autoplot(
  pca_res,
  data = colData(adat_summarized_experiment),
  shape = FALSE,   # turn off autoplot’s shapes
  label = FALSE    # turn off autoplot’s default rowname labels
) +
  geom_point(
    aes(x = PC1, y = PC2),
    shape = 21, size = 4, color = "black", fill = "black"
  ) +
  geom_text_repel(
    aes(x = PC1, y = PC2, label = 1:nrow(colData(adat_summarized_experiment))),
    size = 3
  )
autoplot(pca_res, data = colData(adat_summarized_experiment), label = TRUE, shape = FALSE,label.label = 1:nrow(colData(adat_summarized_experiment)), label.size = 3)
dev.off()

#added filtering
adat_summarized_experiment <- adat_summarized_experiment[-grep("Internal Use Only", rowData(adat_summarized_experiment)$TargetFullName), ] #7313 (-283)
#remove outliers from pca (210 samples to 207 samples)
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$SampleGroup != "K-39")]
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$SampleGroup != "K-122")]
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$SampleGroup != "K-84")]
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$SampleGroup != "K-580")]
#remove non-human samples
adat_summarized_experiment <- adat_summarized_experiment[which(rowData(adat_summarized_experiment)$Organism == "Human"),] #7301
adat_summarized_experiment <- adat_summarized_experiment[which(rowData(adat_summarized_experiment)$EntrezGeneID != ""),] #7285
colData(adat_summarized_experiment)$PathLBD <- as.integer(ifelse(colData(adat_summarized_experiment)$PathLBD == 2, 1, 0)) 
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$BBID != "K-0666")] #206
adat_summarized_experiment <- adat_summarized_experiment[,which(colData(adat_summarized_experiment)$BBID != "SLI-171")] #205
#adat_summarized_experiment <- adat_summarized_experiment[,which(!is.na(colData(adat_summarized_experiment)$PathAD))] #204
#use the minimum AT8 total value to log normalize the AT8 column in metadata
AT8 <- colData(adat_summarized_experiment)$AT8_total
AT8 <- AT8[which(AT8 > 0)]
AT8min <- min(AT8, na.rm=T)
colData(adat_summarized_experiment)$AT8_total <- log(colData(adat_summarized_experiment)$AT8_total + 0.002628708)

#investigating metadata
metadata_adat <- colData(adat_summarized_experiment)
metadata_adat$BBID[c(11, 19, 23, 30, 39, 49, 59, 69)] <- c("SLI-016", "SLI-020", "SLI-096", "SLI-023", "SLI-073", "SLI-079", "SLI-084", "SLI-085")
head(metadata_2)
#metadata_2 <- metadata_2[which(metadata_2$ID_new %in% metadata_adat$BBID),]
metadata_adat_m <- merge(metadata_adat,metadata_2, by.x = "BBID", by.y = "ID_new")
colData(adat_summarized_experiment) <- metadata_adat_m
colData(adat_summarized_experiment)$PMI <- as.numeric(colData(adat_summarized_experiment)$PMI)
colData(adat_summarized_experiment)$Group_de <- ifelse(colData(adat_summarized_experiment)$CTEStage == 0, 0, ifelse(colData(adat_summarized_experiment)$CTEStage == 1, 1,ifelse(colData(adat_summarized_experiment)$CTEStage == 2, 1,ifelse(colData(adat_summarized_experiment)$CTEStage == 3, 2,ifelse(colData(adat_summarized_experiment)$CTEStage == 4, 2,colData(adat_summarized_experiment)$Group_de)))))
colData(adat_summarized_experiment)$Group_de[168] <- 1
colData(adat_summarized_experiment)$Group_de <- as.numeric(colData(adat_summarized_experiment)$Group_de)
#adat_summarized_experiment <- adat_summarized_experiment[,which(!is.na(colData(adat_summarized_experiment)$PMI))]
saveRDS(adat_summarized_experiment, file = "/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
