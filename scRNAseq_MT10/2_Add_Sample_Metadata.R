library(dplyr)
library(Seurat)


print("Add Metadata ...")
metadata <- read.csv(paste0("csvFiles/metadata"))


# Run 1
runNo <- "run1"

sampleInfo <- readRDS(paste0("rObjects/sampleInfo_", runNo, ".rds"))
sampleInfo <- sampleInfo %>% left_join(metadata)

sObj_Filtered <- readRDS(paste0("rObjects/sObj_Filtered_", runNo, ".rds"))

## get file names for each seurat object
nameList <- function(sObj){
  sObjName <- sObj$orig.ident %>% levels()
  return(sObjName)}

names(sObj_Filtered) <- sapply(sObj_Filtered, nameList)

## make sure that seurat list is the same order as sample info
sObj_Filtered <- sObj_Filtered[sampleInfo$fName]


## add metadata to seurat objects
for (i in 1:length(sObj_Filtered)){
  sObj_Filtered[[i]]$batch <- sampleInfo$Batch[i]
  sObj_Filtered[[i]]$seq_batch <- sampleInfo$Seq_Batch[i]
  sObj_Filtered[[i]]$day <- sampleInfo$Day[i]
  sObj_Filtered[[i]]$cell_line <- "AD3"
}

saveRDS(sampleInfo, paste0("rObjects/sampleInfo_", runNo, ".rds"))
saveRDS(sObj_Filtered, paste0("rObjects/sObj_Filtered_", runNo, ".rds"))


# Run 2
runNo <- "run2"
sampleInfo <- readRDS(paste0("rObjects/sampleInfo_", runNo, ".rds"))
sampleInfo <- sampleInfo %>% left_join(metadata)

sObj_Filtered <- readRDS(paste0("rObjects/sObj_Filtered_", runNo, ".rds"))

## get file names for each seurat object
nameList <- function(sObj){
  sObjName <- sObj$orig.ident %>% levels()
  return(sObjName)}

names(sObj_Filtered) <- sapply(sObj_Filtered, nameList)

## make sure that seurat list is the same order as sample info
sObj_Filtered <- sObj_Filtered[sampleInfo$fName]


## add metadata to seurat objects
for (i in 1:length(sObj_Filtered)){
  sObj_Filtered[[i]]$batch <- sampleInfo$Batch[i]
  sObj_Filtered[[i]]$seq_batch <- sampleInfo$Seq_Batch[i]
  sObj_Filtered[[i]]$day <- sampleInfo$Day[i]
  sObj_Filtered[[i]]$cell_line <- "AD3"
}

saveRDS(sampleInfo, paste0("rObjects/sampleInfo_", runNo, ".rds"))
saveRDS(sObj_Filtered, paste0("rObjects/sObj_Filtered_", runNo, ".rds"))

