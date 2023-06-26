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


sObj_Filtered <- readRDS("rObjects/sObj_Filtered_all.rds")
reference <- readRDS("/data/rachel/Linda_Lako/Retina/fetal_and_embryonic/rObjects/seuratObj_retina_eye_annotated.rds")

reference_noUMAP <- reference
reference_noUMAP@reductions$umap <- NULL
reference_model <- RunUMAP(reference_noUMAP, reduction = "harmony", dims = 1:10, return.model = TRUE)

saveRDS(reference_model, 
        "/data/rachel/Linda_Lako/Retina/fetal_and_embryonic/rObjects/seuratObj_retina_eye_annotated_with_model.rds")



#p1 <- DimPlot(reference_model, group.by = "annotation") + ggtitle("for mapping")
#p2 <- DimPlot(reference, group.by = "annotation") + ggtitle("original umap")
 
#p1 + p2 + ggsave("images/reference_model_umap.png")

i <- 1





  
  
query <- merge(sObj_Filtered[[1]], sObj_Filtered[-1])
anchors <- FindTransferAnchors(reference = reference_model, 
                               query = query,
                               dims = 1:20, 
                               reference.reduction = "harmony")

query <- MapQuery(anchorset = anchors, 
                  reference = reference_model, 
                  query = query,
                  refdata = reference_model$annotation, 
                  reference.reduction = "harmony", 
                  reduction.model = "umap",
                  transferdata.args = list(k.weight = 10, dims = 1:20)
                  )


predictions <- TransferData(anchorset = anchors,
                            refdata = reference_model$annotation,
                            dims = 1:20,
                            k.weight = 10)


query$cell_type <- predictions$predicted.id
query <- AddMetaData(query, metadata = predictions)


DimPlot(query, group.by = "cell_type") + DimPlot(reference_model, group.by = "annotation")
saveRDS(query, "rObjects/AD3_all_reference_map.rds")



namesToRemove <- colnames(sObj_Filtered[[1]][[]])[13:29]

removePredictions <- function(sObj, namesToRemove){
  sObj@meta.data <- sObj@meta.data[, setdiff(colnames(sObj[[]]), namesToRemove)]
  return(sObj)}


sObj_Filtered <- lapply(sObj_Filtered, 
                        removePredictions, 
                        namesToRemove)


sObj_Mapped <- list()
for (i in 1:length(sObj_Filtered)){
  
  print(paste("sample", i,  "anchors..."))
    query <- sObj_Filtered[[i]]
    anchors <- FindTransferAnchors(reference = reference_model, 
                                   query = query,
                                   dims = 1:20, 
                                   reference.reduction = "harmony",
                                   k.filter = NA)
    
    print(paste("sample", i,  "mapping..."))
    query <- MapQuery(anchorset = anchors,
                      reference = reference_model,
                      query = query,
                      refdata = reference_model$annotation,
                      reference.reduction = "harmony",
                      reduction.model = "umap",
                      transferdata.args = list(k.weight = 10, dims = 1:20)
    )
    
    print(paste("sample", i,  "transfer data..."))
    predictions <- TransferData(anchorset = anchors,
                                refdata = reference_model$annotation,
                                dims = 1:20,
                                k.weight = 10)
    query <- AddMetaData(query, metadata = predictions)
    
    sObj_Mapped[[i]] <- query
}

#saveRDS(sObj_Mapped, "rObjects/AD3_individual_reference_map.rds")

sObj_Mapped <- readRDS("rObjects/AD3_individual_reference_map.rds")
p1 <- DimPlot(reference_model, group = "annotation", label = TRUE)

i = 1
p2 <- DimPlot(sObj_Mapped[[i]], reduction = "ref.umap", group.by = "predicted.id", label = TRUE)
umap <- p2$data
umap$sample <- i
umap$barcode <- paste0(rownames(umap), "_", i)

for (i in 2:length(sObj_Mapped)){
  p3 <- DimPlot(sObj_Mapped[[i]], reduction = "ref.umap", 
                group.by = "predicted.id", 
                label = TRUE)
  tmp <- p3$data
  tmp$sample <- i
  tmp$barcode <- paste0(rownames(tmp), "_", i)
  umap <- rbind(umap, tmp)
}


ggplot(umap, aes("refUMAP_1", "refUMAP_2")) + geom_point()

sObj_Mapped_Merged <- merge(sObj_Mapped[[1]], sObj_Mapped[-1])
saveRDS(sObj_Mapped_Merged, "rObjects/sObj_Mapped_Merged.rds")
saveRDS(sObj_Mapped, "rObjects/sObj_Mapped.rds")

# umap_embeddings <- umap[, 1:2]
# rownames(umap_embeddings) <- umap$barcode
# seurat_embeddings <- CreateDimReducObject(embeddings = as.matrix(umap_embeddings))
# 
# seurat_object@reductions$pca@cell.embeddings
# 
# DimReduc
# 
# p1 + p2

# sObj_Mapped_Merged[["ref.umap"]] <- seurat_embeddings
# 
# DimPlot(sObj_Mapped_Merged, reduction = "ref.umap")
## create a merged umap and harmony table



DimPlot
