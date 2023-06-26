library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(cowplot)


# Read in individual atac seq samples
# Create seurat object
# Add metadata
# Compute QC metrics
# create sample info with numbers of cells before filtering



set.seed(1234)

#dataDir <- "AD3_ATAC_results/"

dataDir <- "results_AD3_2022_006/"
samples <- list.files(dataDir) 


howManyCells <- function(sObj){dim(sObj)[2]}

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC' 


sampleInfo <- data.frame(fName=samples)
sampleInfo$Tissue <- "Organoids"
sampleInfo$Day <- sampleInfo$fName %>% 
  gsub("AD3_hiPSCs_", "", .) %>% 
  gsub("AD3_Organoids_", "", .) %>%
  gsub("_.*", "", .)

sObj_List <- list()

for (i in 1:nrow(sampleInfo)){
  c_data <- Read10X_h5(filename = (file.path(dataDir, sampleInfo$fName[i], "filtered_peak_bc_matrix.h5")))
  m_data <- read.csv(file.path(dataDir, sampleInfo$fName[i], "singlecell.csv"),  header = TRUE,  row.names = 1) 
  chrom_assay <- CreateChromatinAssay(counts = c_data, 
                                      sep = c(":", "-"),
                                      genome = "hg38", 
                                      fragments = (file.path(dataDir, 
                                                             sampleInfo$fName[i], 
                                                             "fragments.tsv.gz")) ,
                                      min.cells = 10, min.features = 200)
  sObj_List[[i]] <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = m_data)
  Annotation(sObj_List[[i]]) <- annotations
  
  ### add sample info
  sObj_List[[i]]$sample_id <- sampleInfo$fName[i]
  sObj_List[[i]]$tissue <- sampleInfo$Tissue[i]
  sObj_List[[i]]$day <- sampleInfo$Day[i]
  
  ### compute QC metrics
  sObj_List[[i]] <- NucleosomeSignal(object = sObj_List[[i]])
  sObj_List[[i]] <- TSSEnrichment(object = sObj_List[[i]], fast = FALSE)
  
  sObj_List[[i]]$high.tss <- ifelse(sObj_List[[i]]$TSS.enrichment > 2, 'High', 'Low')
  sObj_List[[i]]$nucleosome_group <- ifelse(sObj_List[[i]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  
  sObj_List[[i]]$pct_reads_in_peaks <- sObj_List[[i]]$peak_region_fragments / 
    sObj_List[[i]]$passed_filters * 100
  sObj_List[[i]]$blacklist_ratio <- sObj_List[[i]]$blacklist_region_fragments / 
    sObj_List[[i]]$peak_region_fragments
}


sampleInfo$BeforeFiltering <- sapply(sObj_List, howManyCells)


  
#saveRDS(sObj_List, "rObjects/sObj_List.rds")
#saveRDS(sampleInfo, "rObjects/sampleInfo.rds")  
  
saveRDS(sObj_List, "rObjects/sObj_List_2.rds")
saveRDS(sampleInfo, "rObjects/sampleInfo_2.rds")

sObj_Filtered <- list()

for (i in 1:nrow(sampleInfo)){
  sObj_Filtered[[i]] <- subset(
    x = sObj_List[[i]],
    subset = peak_region_fragments > 3000 &
      peak_region_fragments < 20000 &
      pct_reads_in_peaks > 0.20 &
      blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )

}


sampleInfo$AfterFiltering <- sapply(sObj_Filtered, howManyCells)
#saveRDS(sampleInfo, "rObjects/SampleInfo.rds")
saveRDS(sampleInfo, "rObjects/SampleInfo_2.rds")


