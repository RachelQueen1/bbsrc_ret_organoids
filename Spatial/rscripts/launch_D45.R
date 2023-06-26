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
library(rmarkdown)


tissue <- "D45"
projDir <- "/data/rachel/Linda_Lako/Spatial/AD3/"
outputDir <- file.path(projDir, "/rObjects/")

## import all files, cluster replicates and and save
sampleName <- c("AD3_D45")

## input data
dataDir <- file.path("/data/rachel/Linda_Lako/Spatial/2022_017", sampleName)

## output data
sce_output <- paste0(outputDir, sampleName, "_sce.rds")
batch_corrected_seurat_output <- paste0(outputDir, sampleName, "_seurat.rds")


## samples
filepath <- list.files(dataDir)
id <- filepath %>% gsub("results_", "", .) %>% gsub("_Results", "", .)
section <- id %>% gsub(".*_", "", .)
sampleInfo <- data.frame(filePath = filepath,
                         id = id,
                         tissue = tissue, 
                         stage = sampleName, 
                         section = section)



## create sce list
sce_list <- list()
for (i in 1:length(filepath)){
  sce_list[[i]] <- Spaniel::createVisiumSCE(file.path(dataDir, filepath[i]))
  sce_list[[i]] <- addPerCellQC(sce_list[[i]], subsets = list(Mito=1:10))
  sce_list[[i]]$id <- sampleInfo$id[i]
  sce_list[[i]]$stage <- sampleInfo$stage[i]
  sce_list[[i]]$tissue <- sampleInfo$tissue[i]
  sce_list[[i]]$section <- sampleInfo$section[i]
}


names(sce_list) <-sampleInfo$filePath



## create seurat list
seurat_list <- lapply(sce_list, as.Seurat, data = NULL)
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", verbose = FALSE)
}


## batch correct
seuratObj_combined <- merge(seurat_list[[1]], seurat_list[-1])
seuratObj_batch_corrected <- seuratObj_combined %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE, 
            vars.to.regress = c("detected", 
                                "total")) 



seuratObj_batch_corrected <- seuratObj_batch_corrected %>% 
  RunPCA(pc.genes = VariableFeatures(seuratObj_batch_corrected), 
         npcs = 20, verbose = FALSE) %>%
  RunHarmony("id", plot_convergence = TRUE, assay.use = "originalexp") %>% 
  RunUMAP(reduction = "harmony", dims = 1:10, assay = "originalexp") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = c(0.5))

saveRDS(seuratObj_batch_corrected, batch_corrected_seurat_output)





### add clusters to sce object
ids <- sampleInfo$id
for (i in 1:4){
  id <- ids[i]
  sce_list[[i]]$cluster <- seuratObj_batch_corrected[,seuratObj_batch_corrected$id == id]$originalexp_snn_res.0.5
}

## save with added clusters
saveRDS(sce_list, sce_output)


## find markers
markers <- FindAllMarkers(seuratObj_batch_corrected, only.pos = TRUE, logfc.threshold = 0.5 )
markers <- markers %>% arrange(cluster, desc(avg_log2FC))

fName <- paste0("csvFiles/", sampleName, ".csv")
write.csv(markers, fName)
fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
saveRDS(markers, fName)


# thetitle <- paste(sampleName)
# stage <- sampleName
# render('../11_PCW_individual_results.rmd', output_file = paste0(sampleName))

