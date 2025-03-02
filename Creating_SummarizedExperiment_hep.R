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

# Read in the adat file and BBID file
adat <- SomaDataIO::read_adat("/restricted/projectnb/cteseq/data/somascan_2024/CTE Somascan/HMS-24-036_2024-08-09/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")
adat_bbid <- fread("/restricted/projectnb/cteseq/projects/somascan_analysis/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
dim(adat)
adat <- log2(adat)
analytes <- SomaDataIO::getAnalytes(adat)
num_metadata <- length(adat) - length(analytes)
analyte_start <- num_metadata+1
dim(adat)

#separate into parts
#count data
counts <- as.data.frame(adat[,analyte_start:length(adat)])
dim(counts)
counts <- t(counts)
dim(counts)

#row data
rdata <- as.data.frame(attributes(adat)$Col.Meta)
rownames(rdata) <- attributes(adat)$names[analyte_start:length(adat)]
dim(rdata)
#rdata <- t(rdata)
#dim(rdata)

#column data
cdata <- as.data.frame(adat[,1:num_metadata])
cdata.rownames <- rownames(cdata)
#merge files together to get BBID for each sample
cdata <- dplyr::left_join(cdata,adat_bbid,by=c("SampleGroup"="SampleGroup"))
rownames(cdata) <- cdata.rownames
cdata <- cdata[which(!is.na(cdata$BBID)),]
#merge metadata into cdata
cdata$Identifyer <- rownames(cdata)
#cdata.rownames <- rownames(cdata)
#index_rm <- which(!(cdata$BBID %in% metadata$SampleName) == TRUE)
cdata <- dplyr::inner_join(cdata, metadata,by=c("BBID"="SampleName"))
#cdata.rownames <- cdata.rownames[-c(index_rm)]
rownames(cdata) <- cdata$Identifyer
dim(cdata)

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
library(ggfortify)
pca_res <- prcomp(t(assays(adat_summarized_experiment)$counts), scale. = TRUE)

autoplot(pca_res, data = colData(adat_summarized_experiment), label = TRUE, shape = FALSE,label.label = "SampleGroup", label.size = 3)

saveRDS(adat_summarized_experiment, file = "/restricted/projectnb/cteseq/projects/somascan_analysis/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
