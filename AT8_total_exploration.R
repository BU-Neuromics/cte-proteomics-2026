#Investigating AT8
#Helen Pennington
#May 15, 2025

#packages 
library(data.table)
library(dplyr)

#files
metadata <- fread("/restricted/projectnb/cteseq/projects/challenge-project-2024/merged_cte_meta.csv")
old_metadata <- fread("/restricted/projectnb/cteseq/projects/somascan/data/CTE_sample_metainfo_updated_02032021.csv")
head(metadata)
head(old_metadata)
length(which(metadata$Core_ID %in% old_metadata$ID))
metadata_m <- merge(metadata, old_metadata, by.x = "Core_ID", by.y = "ID")
metadata_m <- metadata_m %>% select(c("Core_ID","AT8","AT8Log","totyrs.y","AT8sulcus","AT8sulcusLog","AT8crest","AT8crestLog","AT8_total","totyrs.x"))
metadata_m$AT8Log2 <- log2(metadata_m$AT8)
metadata_m$AT8Log10 <- log10(metadata_m$AT8)
metadata_m$AT8Loge <- log(metadata_m$AT8)
head(metadata_m)
write.csv(metadata_m, "/restricted/projectnb/cteseq/projects/somascan/data/AT8_total_exploration.csv")

