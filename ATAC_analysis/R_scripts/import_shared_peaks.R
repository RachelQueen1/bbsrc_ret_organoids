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
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

##--Sample Info -----------
dataDir <- "shared"
sampleInfo <- readRDS("rObjects/sampleInfo.rds")
sampleInfo <- sampleInfo %>% filter(!Day %in% c("D0", "D2"))
sampleInfo$fName <- paste0(sampleInfo$fName, "_Shared")


# ##--Load data ------------
sObj_List <- list()

for (i in 1:nrow(sampleInfo)){
  c_data <- Read10X_h5(filename = (file.path(dataDir, sampleInfo$fName[i], "/filtered_peak_bc_matrix.h5")))
  m_data <- read.csv(file.path(dataDir, sampleInfo$fName[i], "/singlecell.csv"),  header = TRUE,  row.names = 1)
  chrom_assay <- CreateChromatinAssay(counts = c_data,
                                      sep = c(":", "-"),
                                      genome = "hg38",
                                      fragments = (file.path(dataDir,
                                                             sampleInfo$fName[i],
                                                             "/fragments.tsv.gz")) ,
                                      min.cells = 10, min.features = 200)
  sObj_List[[i]] <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = m_data)
  Annotation(sObj_List[[i]]) <- annotations

  ## add sample metadata
  sObj_List[[i]]$day <- sampleInfo[i,"Day"]
  sObj_List[[i]]$batch <- "batch1"

  ## add QC metrics
  sObj_List[[i]] <- NucleosomeSignal(object = sObj_List[[i]])
  sObj_List[[i]] <- TSSEnrichment(object = sObj_List[[i]], fast = FALSE)
  sObj_List[[i]]$pct_reads_in_peaks <- sObj_List[[i]]$peak_region_fragments / sObj_List[[i]]$passed_filters * 100
  sObj_List[[i]]$blacklist_ratio <- sObj_List[[i]]$blacklist_region_fragments / sObj_List[[i]]$peak_region_fragments
  }

names(sObj_List) <- sampleInfo$Day


saveRDS(sObj_List, "rObjects/shared_sObj_List.rds")


##--end------













