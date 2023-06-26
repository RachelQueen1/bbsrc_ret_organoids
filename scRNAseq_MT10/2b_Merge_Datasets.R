library(dplyr)
library(Seurat)

### merge datasets
sObj_List1 <- readRDS(paste0("rObjects/sObj_List_", "run1", ".rds"))
sObj_List2 <- readRDS(paste0("rObjects/sObj_List_", "run2", ".rds"))
sObj_List <- c(sObj_List1, sObj_List2)
saveRDS(sObj_List, paste0("rObjects/sObj_List_", "all", ".rds"))


sObj_Filtered1 <- readRDS(paste0("rObjects/sObj_Filtered_", "run1", ".rds"))
sObj_Filtered2 <- readRDS(paste0("rObjects/sObj_Filtered_", "run2", ".rds"))


sObj_Filtered <- c(sObj_Filtered1, sObj_Filtered2)
saveRDS(sObj_Filtered, paste0("rObjects/sObj_Filtered_", "all", ".rds"))

sampleInfo1 <- readRDS(paste0("rObjects/sampleInfo_", "run1", ".rds"))
sampleInfo2 <- readRDS(paste0("rObjects/sampleInfo_", "run2", ".rds"))
sampleInfo <- rbind(sampleInfo1, sampleInfo2)
saveRDS(sampleInfo, paste0("rObjects/sampleInfo_", "all", ".rds"))
