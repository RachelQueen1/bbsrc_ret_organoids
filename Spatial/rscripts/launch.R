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

## import all files, cluster replicates and and save
# days <- c("10", "150", "20", "210", "35", "60", "90")
# for (day in days){
#   print(tissue)
#   source(file.path(projDir, "rscripts/1_load_data_AD3.R"))
# }


## find markers
source(file.path(projDir,"rscripts/2_Find_Markers.R"))



## Day 45 sample
allData <- "/data/rachel/Linda_Lako/Spatial/2022_017"
tissue <- "AD3"
projDir <- paste0("/data/rachel/Linda_Lako/Spatial/", tissue)
outputDir <- file.path(projDir, "/rObjects/")

day <-"45"
source(file.path(projDir, "rscripts/1_load_data_AD3.R"))


sampleName <- "AD3_D45"

markers <- FindAllMarkers(seuratObj_batch_corrected, only.pos = TRUE, logfc.threshold = 0.5 )
markers <- markers %>% arrange(cluster, desc(avg_log2FC))

fName <- paste0("csvFiles/", sampleName, ".csv")
write.csv(markers, fName)
fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
saveRDS(markers, fName)


thetitle <- paste(sampleName)
day <- gsub("AD3_D", "", sampleName) %>% as.numeric()
render('individual_results.rmd', output_file = paste0("AD3_results_", sampleName))
