library(mice)
adat <- readRDS("/restricted/projectnb/cteseq/projects/somascan/data/HMS-24-036_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP_summarizedexperiment.rds")
colData(adat)$Group_de <- ifelse(colData(adat)$CTEStage == 0, "none", ifelse(colData(adat)$CTEStage == 1, "low",ifelse(colData(adat)$CTEStage == 2, "low",ifelse(colData(adat)$CTEStage == 3, "high",ifelse(colData(adat)$CTEStage == 4, "high",colData(adat)$Group_de)))))
colData(adat)$Group_de[169] <- "low"
colData(adat)$Group_de
col_data <- as.data.frame(colData(adat))
col_data_filt <- subset(col_data, select = -c(PlateId,ScannerID,PlatePosition,SlideId,Subarray,SampleId,PlateRunDate,Group,
                                              SampleType,PercentDilution,SampleMatrix,Barcode,Barcode2d,SampleName,SampleNotes,AliquotingNotes,
                                              SampleDescription,AssayNotes,TimePoint,ExtIdentifier,SsfExtId,SampleGroup,SiteId,SubjectID,CLI,
                                              RMA,SampleNumber,StudyId,HybControlNormScale,RowCheck,NormScale_10,Identifyer,Core_ID,Batch,Year))
colSums(is.na(col_data_filt))
cds <- col_data_filt$CDStot
max_nas <- 30
col_data_filt <- col_data_filt[, colSums(is.na(col_data_filt)) <= max_nas]
col_data_filt$CDStot <- cds
# Your data
dat <- col_data_filt

# Step 1: Create method vector (how each variable should be imputed)
meth <- make.method(dat)
meth[] <- ""   # default: no imputation

# Only impute the two variables of interest
meth["PMI"] <- "norm"      # numeric (could use "pmm" instead for better realism)
meth["PathAD"] <- "logreg" # binary categorical (0/1)

# Step 2: Create predictor matrix
pred <- make.predictorMatrix(dat)
pred[,] <- 0                # start with no predictors

# Let *all other variables* predict PMI and PathAD
pred["PMI", ] <- 1
pred["PathAD", ] <- 1

# Remove self-prediction
pred["PMI", "PMI"] <- 0
pred["PathAD", "PathAD"] <- 0

# Exclude identifiers (so they don’t predict)
pred[ , c("SOMASCAN.ID")] <- 0
pred[ , c("BBID")] <- 0
# Step 3: Run multiple imputation
imp <- mice(dat, m = 5, method = meth, predictorMatrix = pred, seed = 123)
head(imp)
stripplot(imp, PMI ~ .imp)
densityplot(imp, ~ PMI)
densityplot(imp, ~ PMI)
