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
#library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

##--- Functions --------
howManyCells <- function(sObj){dim(sObj)[2]}

sObj <- readRDS("rObjects/shared_sObj_merge_retina_with45_60.rds")
## Remove space in MG
sObj$CellType <- sObj$CellType %>% as.character()
sObj$CellType[grepl("MG", sObj$CellType)] <- "MG"
sObj$CellType[grepl("RPCs", sObj$CellType)] <- "RPCs"


cellorder <- c("pCM", "RPCs", "T1","T3", "RGCs", "HCs + ACs", "cones", "rods", "BPs",  "MG")
sObj <- SetIdent(sObj, value = "CellType")
DefaultAssay(sObj) <- "chromvar"

differential.activity <- FindAllMarkers(
  object = sObj,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

write.csv(differential.activity, "csvFiles/chromvar_enriched_motifs.csv")