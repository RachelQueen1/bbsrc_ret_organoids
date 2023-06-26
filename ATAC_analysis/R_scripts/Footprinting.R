library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(TFBSTools)
library(future)
plan()
plan("multicore", workers = 8)
#increase the maximum memory usage
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 50 Gb RAM

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

motifs <- sapply(pfm, function(x){x@name})
motif_df <- data.frame(motif = names(motifs),
                       name = unname(motifs)
)
top5 <- readRDS("rObjects/top5_motifs.rds")
sObj <- readRDS("rObjects/sObj_with_motifs_fp.rds")

# retina_markers <- readRDS("/data/rachel/Linda_Lako/Retina/markers/Retina_Markers.rds")
# retina_markers_df <-
# retina_markers <- unlist(retina_markers) %>% unname() %>% intersect(Annotation(sObj)$gene_name)
# motifs <- intersect(retina_markers, motif_df$name)
footprints <- GetAssayData(object = sObj, assay = "peaks", 
                           slot = "positionEnrichment") %>% names()




fp_motifs <-setdiff(top5$name, footprints)

sObj <- Footprint(
  object = sObj,
  motif.name = fp_motifs,
  in.peaks = TRUE,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# sObj <- Footprint(
#   object = sObj,
#   motif.name = "OTX2",
#   in.peaks = TRUE,
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )


saveRDS(sObj, "rObjects/sObj_with_motifs_fp.rds")

#ggsave(PlotFootprint(sObj, features = motifs[[2]]), "footprint.tiff")

#PlotFootprint(sObj_with_motifswithout_correction_fp, features = motifs[[15]])
