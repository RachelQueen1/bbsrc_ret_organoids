library(Seurat)
library(cowplot)
library(SCFunctionsV3)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(DoubletFinder)
library(pheatmap)
print("QC ...")
# Functions
### 1) how many cells
howManyCells <- function(sObj){dim(sObj)[2]}

### 2) pre doublet filtering
preDoublets <- function(seuratObj){  
  seuratObj <- NormalizeData(seuratObj)
  
  seuratObj <- ScaleData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj, 
                                    x.low.cutoff = 0.0125, 
                                    y.cutoff = 0.25, 
                                    do.plot=FALSE)
  seuratObj <- RunPCA(seuratObj, pc.genes = seuratObj@var.genes, pcs.print = 0)
  
  # set.seed(1234)
  # seuratObj <- RunTSNE(seuratObj, dims.use = 1:10, verbose=TRUE)
  return(seuratObj)
}

### 3) doublet filtering
findDoublets <- function(seuratObj){
  
  ### calculate expected number of doublets
  nExp_poi  <- round(0.15*nrow(seuratObj@meta.data))
  
  ### predict doublets
  seuratObj <- doubletFinder_v3(seuratObj, 
                                PCs = 1:10, 
                                pN = 0.25, 
                                pK = 0.01, 
                                nExp = nExp_poi, 
                                reuse.pANN = FALSE, 
                                sct=FALSE)
  
  
  seuratObj@meta.data <- seuratObj@meta.data %>% 
    rename_at(vars(starts_with("DF.classifications")), 
              funs(str_replace(., ".*", "DF.classifications"))) %>%
    rename_at(vars(starts_with("pANN")), 
              funs(str_replace(., ".*", "pANN")))
  
  return(seuratObj) 
  
}

### 4) Read and filter Data
read_filter <- function(dataDir){
  samples <- list.files(dataDir)
  sampleInfo <- data.frame(fName=samples)
  
  
  
  sObj_List <- list()
  for (i in 1:nrow(sampleInfo)){
    df <- Read10X(file.path(dataDir, sampleInfo$fName[i], "filtered_feature_bc_matrix"))
    sObj_List[[i]] <- CreateSeuratObject(counts = df, project = sampleInfo$fName[i], min.cells = 10, min.features = 0)
    sObj_List[[i]][["percent.mt"]] <- PercentageFeatureSet(sObj_List[[i]], pattern = "^MT-")
  }
  sampleInfo$beforeFiltering <- sapply(sObj_List, howManyCells)
  saveRDS(sObj_List, paste0("rObjects/sObj_List_", dataDir, ".rds"))
  
  
  # Filter Thresholds
  minCounts <- 1000
  minFeatures <- 500
  maxMit <- 20
  
  ### Filter
  sObj_Filtered <- list()
  for (i in 1:length(samples)){
    filter <- sObj_List[[i]]$nCount_RNA > minCounts & 
      sObj_List[[i]]$nFeature_RNA > minFeatures & 
      sObj_List[[i]]$percent.mt < maxMit
    
    sObj_Filtered[[i]] <- sObj_List[[i]][ , filter]
  }
  
  
  ### Find and remove doublets
  for (i in 1:length(sObj_Filtered)){
    sObj_Filtered[[i]] <-  preDoublets(sObj_Filtered[[i]])
    sObj_Filtered[[i]] <- findDoublets(sObj_Filtered[[i]])
    
    cellFilter <- sObj_Filtered[[i]]$DF.classifications == "Singlet"
    sObj_Filtered[[i]] <- sObj_Filtered[[i]][, cellFilter]
    
    sObj_Filtered[[i]] <- NormalizeData(sObj_Filtered[[i]], verbose = FALSE)
    sObj_Filtered[[i]] <- FindVariableFeatures(sObj_Filtered[[i]], 
                                               selection.method = "vst",
                                               nfeatures = 2000, 
                                               verbose = FALSE)
  }
  
  
  
  
  
  sampleInfo$afterFiltering <- sapply(sObj_Filtered, howManyCells)
  
  
  saveRDS(sampleInfo, paste0("rObjects/sampleInfo_", dataDir, ".rds"))
  saveRDS(sObj_Filtered, paste0("rObjects/sObj_Filtered_", dataDir, ".rds"))
}

# QC steps
read_filter("run1")