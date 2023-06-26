library(Spaniel)
library(scran)
library(scater)
library(batchelor)
library(dplyr)
library(Seurat)
library(SCFunctionsV3)
library(cowplot)
library(RColorBrewer)
library(SCFunctionsV3)
library(harmony)


allData <- "/data/rachel/Linda_Lako/Spatial/spaceranger_outs_2021_102"
tissue <- "AD3"

projDir <- paste0("/data/rachel/Linda_Lako/Spatial/", tissue)
outputDir <- file.path(projDir, "/rObjects/")


## find markers
seuratObjs_path <- list.files(outputDir)[grepl("_seurat", 
                                          list.files(outputDir))]


setwd(projDir)

for (sObj_path in seuratObjs_path){
  sObj <- readRDS(file.path(outputDir, sObj_path))
  sampleName <- gsub("_seurat.rds", "", sObj_path)
  markers <- FindAllMarkers(sObj, only.pos = TRUE, logfc.threshold = 0.5 )
  markers <- markers %>% arrange(cluster, desc(avg_log2FC))
  
  fName <- paste0("csvFiles/", sampleName, ".csv")
  write.csv(markers, fName)
  fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
  saveRDS(markers, fName)
}