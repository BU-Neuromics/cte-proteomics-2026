#Amelia just for covariates

library("Amelia")
library("dplyr")
library("sjlabelled") ## Remove labels from variables
library(tidyr)
library(ggplot2)
library(broom)
library(data.table)
library(haven)
library(SummarizedExperiment)
library(S4Vectors)
library(limma)
library(EnhancedVolcano)
library(clusterProfiler)
library(fgsea)
library(gt)
library(Amelia)

### Run Amelia ###
filepath <- "/restricted/projectnb/cteseq/projects/somascan/imputation"
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
col_data <- as.data.frame(colData(adat))
col_data_filt <- subset(col_data, select = -c(PlateId,ScannerID,PlatePosition,SlideId,Subarray,SampleId,PlateRunDate,Group,
                                              SampleType,PercentDilution,SampleMatrix,Barcode,Barcode2d,SampleName,SampleNotes,AliquotingNotes,
                                              SampleDescription,AssayNotes,TimePoint,ExtIdentifier,SsfExtId,SampleGroup,SiteId,SubjectID,CLI,
                                              RMA,SampleNumber,StudyId,HybControlNormScale,RowCheck,NormScale_10,Identifyer,Core_ID,Batch,Year,
                                              nphemo,npavas,npwmr,nphipscl,npbraak,npneur,npadnc,npdiff,npamy,npold,npold1,npold2,npold3,npold4,
                                              npoldd,npoldd1,npoldd2,npoldd3,npoldd4,nparter,nppath,nppath6,nplbod,nptdpb,nptdpc,nptdpe,npftdtau,
                                              nppick,npftdt2,npcort,npprog,npftdt5,npftdt8,npftdt9,npftdt10,npftdtdp,npoftd,PathPrion, micdorfront, 
                                              micinfpar, micalc, PathMND,suicide,TMEM106B_dom,TMEM106B_invrec,disdur))
colSums(is.na(col_data_filt))
cds <- col_data_filt$CDStot
max_nas <- 30
col_data_filt <- col_data_filt[, colSums(is.na(col_data_filt)) <= max_nas]
col_data_filt$CDStot <- cds

noms=c("DementiaHx","ParknismHx","PathAD","PathLBD","PathFTD","CTE",
       "sport","apoe_de","cod","csparCTE","rs1990622","rs3173615","race")

ords=c("CTEStage","apoe","Group_de")

idvars=c("BBID", "SOMASCAN.ID")
cov_data_subset_clean <- remove_all_labels(col_data_filt)
# Run imputation
amelia.out1 <- amelia(
  cov_data_subset_clean, 
  m=1, 
  p2s=2, 
  noms=noms, 
  ords=ords, 
  idvars=idvars, 
  incheck=T,
  empri=0.005*nrow(cov_data_subset_clean),
  parallel = "multicore"
  # ncpus = 30
)
compare.density(amelia.out1, "PMI")
compare.density(amelia.out1, "totyrs")
compare.density(amelia.out1, "AT8_total")
compare.density(amelia.out1, "CDStot")
compare.density(amelia.out1, "DementiaHx")
compare.density(amelia.out1, "faqtot")

st.filepath <- paste(filepath,"/Stable_Releases/",Sys.Date(),sep="")
dir.create(st.filepath)
saveRDS(amelia.out1, file = paste(st.filepath,"/full_imputations.rds",sep=""))
write.amelia(obj = amelia.out1, file.stem = paste(st.filepath,"/outdata",sep=""))

amelia.out1 <- readRDS("/restricted/projectnb/cteseq/projects/somascan/imputation/Stable_Releases/2025-10-14/full_imputations.rds")
compare.density(amelia.out1, "PMI")
compare.density(amelia.out1, "totyrs")
compare.density(amelia.out1, "AT8_total")
compare.density(amelia.out1, "CDStot")
compare.density(amelia.out1, "DementiaHx")
compare.density(amelia.out1, "faqtot")
