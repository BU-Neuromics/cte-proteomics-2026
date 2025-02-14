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
adat <- log10(adat)
analytes <- SomaDataIO::getAnalytes(adat)
num_metadata <- length(adat) - length(analytes)
analyte_start <- num_metadata+1

#separate into parts
#column data
cdata <- as.data.frame(attributes(adat)$Col.Meta)
rownames(cdata) <- attributes(adat)$names[analyte_start:length(adat)]
head(cdata)

#row data
rdata <- as.data.frame(adat[,1:num_metadata])
rdata.rownames <- rownames(rdata)
#merge files together to get BBID for each sample
rdata <- dplyr::left_join(rdata,adat_bbid,by=c("SampleGroup"="SampleGroup"))
rownames(rdata) <- rdata.rownames
rdata <- rdata[which(!is.na(rdata$BBID)),]
#merge metadata into rdata
rdata.rownames <- rownames(rdata)
rdata <- merge(rdata, metadata,by.x="BBID", by.y="SampleName")
rdata.rownames <- rdata.rownames[-c(117, 175)]
rownames(rdata) <- rdata.rownames
dim(rdata)
head(rdata)

#count data
counts <- as.data.frame(adat[,analyte_start:length(adat)])
counts <- counts[which(rownames(counts)%in%rownames(rdata)),]
counts[1,1:20]

#check that everything got separated correctly
dim(cdata)
dim(rdata)
dim(counts)

which(rownames(rdata) != rownames(counts))
which(rownames(cdata) != colnames(counts))

#build summarized experiment data frame

counts <- as.matrix(counts)
colData <- DataFrame(cdata)
rowData <- DataFrame(rdata)

adat_summarized_experiment <- SummarizedExperiment(assays=list(counts=counts),
                                                   rowData=rowData, colData=colData)
colData(adat_summarized_experiment)
rowData(adat_summarized_experiment)

saveRDS(adat_summarized_experiment, file = "/restricted/projectnb/cteseq/projects/somascan_analysis/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
