## ---Libraries ------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

##--- Functions --------
howManyCells <- function(sObj){dim(sObj)[2]}

# ##--- Setup -----------
set.seed(1234)
sampleInfo <- readRDS("rObjects/sampleInfo.rds")
sampleInfo <- sampleInfo %>% filter(!Day %in% c("D0", "D2"))

## add in d45 and d60
sampleInfo2 <- readRDS("/data/rachel/Linda_Lako/Retina/AD3/ATAC_analysis/rObjects/SampleInfo_2.rds")[colnames(sampleInfo)]

sampleInfo<- rbind(sampleInfo, sampleInfo2)
saveRDS(sampleInfo, "rObjects/sampleInfo_usedIntegration.rds")



# ---- Load shared list
sObj_List <- readRDS("rObjects/shared_sObj_List.rds")


# ##--- Original Data -------
sObj_Orig <- readRDS("rObjects/sObj_annotated.rds")
#sapply(sObj_Orig, function(x){unique(x$CellType)}) %>% unlist() %>% unique()
nonRetina <- c("OV", "neuronal progenitors", "OSE Prog", "OSE", "optic stalk",
               "OSE prog", "OV epithelium", "RPE", "optic stalk/optic nerve?",
               "neural crest derived cells", "neural progenitors",
               "Optic stalk/OPTIC NERVE", "NRPC/RPCs", "NRPCs", "NRPCs ", "PSCs",
               "corneal endothelium", "fibroblasts", "Corneal stroma", "Microglia", "NRPCs")

sObj_Orig_Retina <- lapply(sObj_Orig, function(x){x[,!x$CellType %in% nonRetina]})


# ## --- Filter --------
sObj_filtered <- list()
ids <- intersect(names(sObj_List), names(sObj_Orig_Retina))

for (id in ids){
print(id)
sObj_orig <- sObj_Orig_Retina[[id]]
sObj_shared <- sObj_List[[id]]

cellsUse <- intersect(colnames(sObj_orig), colnames(sObj_shared))
sObj_orig <- sObj_orig[,cellsUse]
sObj_shared <- sObj_shared[,cellsUse]
sObj_shared$CellType <- sObj_orig$CellType
sObj_filtered[[id]] <- sObj_shared
}

# ## update cell numbers
# rownames(sampleInfo) <- sampleInfo$ID
# sampleInfo$shared_BeforeFiltering <- sapply(sObj_List, howManyCells)
# sampleInfo$shared_AfterFiltering <- NA
# after <- sapply(sObj_filtered, howManyCells)
# sampleInfo[names(after), "shared_AfterFiltering"] <- after
# 
# saveRDS(sampleInfo, "rObjects/shared_SampleInfo_retina.rds")



### day 45 and day 60
sObj_45_60_shared <- readRDS("rObjects/shared_sObj_List_D45_D60.rds")
sObj_45_60_shared_Retina <- lapply(sObj_45_60_shared, function(x){x[,!x$CellType %in% nonRetina]})





# # ## --- Normalise, Dimension reduction, cluster
for (i in 1:length(sObj_filtered)){
  sObj_filtered[[i]] <- RunTFIDF(sObj_filtered[[i]])
  sObj_filtered[[i]] <- FindTopFeatures(sObj_filtered[[i]], min.cutoff = 'q0')
  sObj_filtered[[i]] <- RunSVD(sObj_filtered[[i]])
  sObj_filtered[[i]] <- RunUMAP(object = sObj_filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_filtered[[i]] <- FindNeighbors(object = sObj_filtered[[i]], reduction = 'lsi', dims = 2:30)
  sObj_filtered[[i]] <- FindClusters(object = sObj_filtered[[i]], verbose = FALSE, algorithm = 3)
  }


for (i in 1:length(sObj_45_60_shared_Retina)){
  sObj_45_60_shared_Retina[[i]] <- RunTFIDF(sObj_45_60_shared_Retina[[i]])
  sObj_45_60_shared_Retina[[i]] <- FindTopFeatures(sObj_45_60_shared_Retina[[i]], min.cutoff = 'q0')
  sObj_45_60_shared_Retina[[i]] <- RunSVD(sObj_45_60_shared_Retina[[i]])
  sObj_45_60_shared_Retina[[i]] <- RunUMAP(object = sObj_45_60_shared_Retina[[i]], reduction = 'lsi', dims = 2:30)
  sObj_45_60_shared_Retina[[i]] <- FindNeighbors(object = sObj_45_60_shared_Retina[[i]], reduction = 'lsi', dims = 2:30)
  sObj_45_60_shared_Retina[[i]] <- FindClusters(object = sObj_45_60_shared_Retina[[i]], verbose = FALSE, algorithm = 3)
  sObj_45_60_shared_Retina[[i]][["peaks"]] <- sObj_45_60_shared_Retina[[i]][["ATAC"]]
  }




# # ## --- Add Gene Expression  --------
for (i in 1:length(sObj_filtered)){
  sObj <- sObj_filtered[[i]]
  gene.activities <- GeneActivity(sObj)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  sObj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  sObj <- NormalizeData(
    object = sObj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(sObj$nCount_RNA)
  )
  sObj_filtered[[i]] <- sObj
}


### merge data sets
sObj_filtered <- c(sObj_filtered, sObj_45_60_shared_Retina)

saveRDS(sObj_filtered , "rObjects/shared_sObj_filtered_retina_with45_60.rds")

# ## --integrate samples -----
sObj_merge <- merge(sObj_filtered[[1]], sObj_filtered[-1])

topFeat <- FindTopFeatures(object = sObj_merge[['peaks']][])
close <- ClosestFeature(sObj_merge, regions = rownames(sObj_merge))

retinaDir <- "/data/rachel/Linda_Lako/Retina/"
sObj_RNA <- sObj_Retina <- readRDS(paste0(retinaDir, "developmental_retina_and_eye_integrated/rObjects/seuratObj_retina_eye_mt10_noHB.rds"))
vf <- VariableFeatures(sObj_RNA)
featuresUse <- close[close$gene_name %in%  vf , "query_region"]


# ## Normalization and linear dimensional reduction
sObj_merge <- RunTFIDF(sObj_merge)
sObj_merge <- FindTopFeatures(sObj_merge, min.cutoff = 'q0')
sObj_merge <- RunSVD(sObj_merge, features = featuresUse)

elbowDat <- ElbowPlot(sObj_merge, reduction = "lsi", ndims = 50)
corVals <- DepthCor(sObj_merge, n = 50)
componentUse <- corVals$data$Component[abs(corVals$data$nCount_peaks) < 0.4 &
                                         elbowDat$data$stdev > 1.4]



# ## non linear dimension reduction and clustering (no batch correction)
sObj_merge <- RunUMAP(object = sObj_merge,
                      reduction = 'lsi',
                      dims = componentUse,
                      n.neighbors = 30,
                      min.dist = 0.25,
                      n.epochs = 5000
                      )

DimPlot(sObj_merge, group.by = "CellType", label = TRUE)

sObj_merge <- FindNeighbors(object = sObj_merge, reduction = 'lsi', dims = componentUse)
sObj_merge <- FindClusters(object = sObj_merge, verbose = FALSE, algorithm = 3)

saveRDS(sObj_merge, "rObjects/shared_sObj_merge_retina_with45_60.rds")

# ### DA peaks  ------------------------------------------------
#sObj_merge <- readRDS("rObjects/shared_sObj_merge.rds")
DA_peaks <- FindAllMarkers(sObj_merge,
                           only.pos = FALSE,
                           min.pct = 0.05,
                           test.use = 'LR',
                           latent.vars = 'peak_region_fragments'
)

saveRDS(DA_peaks, "rObjects/shared_DA_peaks_retina_with45_60.rds" )


# # ### Add motif information ---------------------------------
DefaultAssay(sObj_merge) <- "peaks"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

sObj_merge <- AddMotifs(
  object = sObj_merge,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

sObj_merge <- RunChromVAR(
  object = sObj_merge,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(sObj_merge,
        "rObjects/shared_sObj_merge_retina_with45_60.rds")

# ### Find enriched motifs  ---------------------------------
#sObj_merge <- readRDS("rObjects/shared_sObj_merge_retina_with45_60.rds")

DefaultAssay(sObj_merge) <- "peaks"
DA_peaks <- readRDS("rObjects/shared_DA_peaks_retina_with45_60.rds")

# get top differentially accessible peaks
DA_peak.sig <- DA_peaks[DA_peaks$p_val < 0.005, ]

noMotifs <- DA_peak.sig$cluster %>% table() %>% data.frame()
colnames(noMotifs) <- c("cluster", "number")
noMotifs <- noMotifs[noMotifs$number > 10,]
clusters <- noMotifs$cluster


DA_peak.sig <- DA_peak.sig[DA_peak.sig$cluster %in% clusters, ]
top.da.peak <- rownames(DA_peak.sig)


# find peaks open in in each cluster
open.peaks <- AccessiblePeaks(sObj_merge, idents = clusters)

# match the overall GC content in the peak set
meta.feature <- GetAssayData(sObj_merge, assay = "peaks", slot = "meta.features")

open <- meta.feature[open.peaks, ]
query <- meta.feature[top.da.peak, ]#
query <- query[!is.na(query$count), ]

## background
peaks.matched <- MatchRegionStats(
  meta.feature = open,
  query.feature = query,
  n = 10000
)


motif_list <- list()

for (cluster in clusters){
  # get top differentially accessible peaks for the cluster
  DA_peaks_cluster <- DA_peak.sig[DA_peak.sig$cluster == cluster,]
  top.da.peak <- rownames(DA_peaks_cluster)
  top.da.peak <- top.da.peak[top.da.peak %in% rownames(sObj_merge)]

  motif_list[[cluster]] <- FindMotifs(
    object = sObj_merge,
    features = top.da.peak,
    background = peaks.matched
  )

}
library(openxlsx)
write.xlsx(motif_list, "csvFiles/shared_enriched_motifs_retina_with45_60.xls")
saveRDS(motif_list, "rObjects/shared_enriched_motifs_retina_with45_60.rds")



### DE genes  ------------------------------------------------
#sObj_merge <- readRDS("rObjects/shared_sObj_merge_retina_with45_60.rds")
DefaultAssay(sObj_merge) <- "RNA"

markers <- FindAllMarkers(sObj_merge,
                          only.pos = TRUE,
                          logfc.threshold = 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
saveRDS(markers, "rObjects/shared_markers_retina_with45_60.rds")
write.csv(markers, "csvFiles/shared_markers_retina_with45_60.csv")
