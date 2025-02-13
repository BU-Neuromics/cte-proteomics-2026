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
adat_bbid <- fread("/usr4/spclpgm/hpenn13/LabadorfRotation/BBIDs_adatfile_LabadorfRotaiton_hp.csv")
adat <- SomaDataIO::read_adat("/restricted/projectnb/cteseq/data/somascan_2024/CTE Somascan/HMS-24-036_2024-08-09/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")
dim(adat)
analytes <- SomaDataIO::getAnalytes(adat)
num_metadata <- length(adat) - length(analytes)
analyte_start <- num_metadata+1

#separate into parts
#row data
cdata <- as.data.frame(attributes(adat)$Col.Meta)
rownames(cdata) <- attributes(adat)$names[analyte_start:length(adat)]
head(cdata)

#column data
rdata <- as.data.frame(adat[,1:num_metadata])
head(rdata)

#count data
counts <- as.data.frame(adat[,analyte_start:length(adat)])
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
