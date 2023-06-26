library(harmony)
library(Seurat)
library(cowplot)
library(SCFunctionsV3)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(DoubletFinder)
library(pheatmap)
library(ggplot2)

## read data
sObj_Mapped <- readRDS("rObjects/AD3_individual_reference_map.rds")





## cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## create batch corrected object
seuratObj <- merge(sObj_Mapped[[1]], sObj_Mapped[-1]) %>% 
  CellCycleScoring(s.features = s.genes, 
                   g2m.features = g2m.genes, 
                   set.ident = FALSE) %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE, 
            vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")) %>% 
  RunPCA(pc.genes = seuratObj_ds@var.genes, npcs = 20, verbose = FALSE) %>%
  RunHarmony(c("batch", "seq_batch", "orig.ident"), plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = c(0.2,1,2.2))


#find_markers.R
markers <- FindAllMarkers(seuratObj, logfc.threshold = 0.7)## DE

markers <- markers %>% dplyr::arrange(cluster, desc(avg_log2FC))

fName <- paste0("csvFiles/", "integrated_markers", ".csv")
write.csv(markers, fName)
fName <- gsub("csvFiles", "rObjects", fName) %>% gsub(".csv", ".rds", .)
saveRDS(markers, fName)

saveRDS(seuratObj, "rObjects/seuratObj_integrated.rds")





