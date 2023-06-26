library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)
howManyCells <- function(sObj){dim(sObj)[2]}
set.seed(1234)

# filter data
# add after filtering to sampleInfo
# cluster each sample


#sObj_List <- readRDS("rObjects/sObj_List.rds")
sObj_List <- readRDS("rObjects/sObj_List_2.rds")
#sampleInfo <- readRDS("rObjects/sampleInfo.rds")
sampleInfo <- readRDS("rObjects/sampleInfo_2.rds")
## filter data ---------------------------

# sObj_Filtered <- list()
# 
# for (i in 1:nrow(sampleInfo)){
#   sObj_Filtered[[i]] <- subset(
#     x = sObj_List[[i]],
#     subset = peak_region_fragments > 3000 &
#       peak_region_fragments < 20000 &
#       pct_reads_in_peaks > 0.20 &
#       blacklist_ratio < 0.05 &
#       nucleosome_signal < 4 &
#       TSS.enrichment > 2
#   )
# 
# }
# 
# 
# sampleInfo$AfterFiltering <- sapply(sObj_Filtered, howManyCells)
# #saveRDS(sampleInfo, "rObjects/SampleInfo.rds")
# saveRDS(sampleInfo, "rObjects/SampleInfo_2.rds")

## normalise, dimension reduction, and cluster
#sObj_Filtered <- readRDS("rObjects/sObj_cluster_2.rds")

### analysis of shared peaks
sObj_Filtered <- readRDS("rObjects/shared_sObj_List_D45_D60.rds")


for (i in 1:nrow(sampleInfo)){
  sObj_Filtered[[i]] <- RunTFIDF(sObj_Filtered[[i]])
  sObj_Filtered[[i]] <- FindTopFeatures(sObj_Filtered[[i]], min.cutoff = 'q0')
  sObj_Filtered[[i]] <- RunSVD(sObj_Filtered[[i]])
  sObj_Filtered[[i]] <- RunUMAP(object = sObj_Filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_Filtered[[i]] <- FindNeighbors(object = sObj_Filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_Filtered[[i]] <- FindClusters(object = sObj_Filtered[[i]], verbose = FALSE, algorithm = 3)
}
saveRDS(sObj_Filtered, "rObjects/sObj_cluster_2.rds")

## add gene activities ---------------------------------
# sObj_Filtered <- readRDS("rObjects/sObj_cluster.rds")
# sObj_Filtered <- readRDS("rObjects/sObj_cluster_2.rds")
for (i in 1:nrow(sampleInfo)){
  gene.activities <- GeneActivity(sObj_Filtered[[i]])
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  sObj_Filtered[[i]][['RNA']] <- CreateAssayObject(counts = gene.activities)
  sObj_Filtered[[i]] <- NormalizeData(
    object = sObj_Filtered[[i]],
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(sObj_Filtered[[i]]$nCount_RNA)
  )
}

# saveRDS(sObj_Filtered, "rObjects/sObj_cluster.rds")
saveRDS(sObj_Filtered, "rObjects/sObj_cluster_2.rds")




### DE genes  ------------------------------------------------
#sObj_Filtered <- readRDS("rObjects/sObj_cluster.rds")
sObj_Filtered <- readRDS("rObjects/sObj_cluster_2.rds")
names(sObj_Filtered) <- sapply(sObj_Filtered, function(x){x$day %>% 
    unique()}) %>% 
  unname()



markers <- list()
for (i in 1:length(sObj_Filtered)){
  day <- names(sObj_Filtered)[i]
  DefaultAssay(sObj_Filtered[[i]]) <- "RNA"
  markers[[day]] <- FindAllMarkers(sObj_Filtered[[i]],
                            only.pos = TRUE,
                            logfc.threshold = 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  
}
#saveRDS(markers, "rObjects/individual_markers.rds")
saveRDS(markers, "rObjects/individual_markers_2.rds")
library(openxlsx)
#write.xlsx(markers, "csvFiles/individual_markers.xls")
write.xlsx(markers, "csvFiles/individual_markers_2.xls")


