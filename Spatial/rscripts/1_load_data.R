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



sce_output <- ""
batch_corrected_seurat_output <- ""


sampleInfo <- data.frame(id = samples,
                         tissue = "AD3", 
                         day = "35")


## create sce list





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
