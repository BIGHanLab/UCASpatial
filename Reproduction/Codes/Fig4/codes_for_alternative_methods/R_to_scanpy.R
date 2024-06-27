library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratDisk)

setwd("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/scanpy_data/")

sc_UCASpatial <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/")
# remove NA meta.data in sc_UCASpatial
for(i in colnames(sc_UCASpatial@meta.data))
{
  if(length(sc_UCASpatial@meta.data[[i]][is.na(sc_UCASpatial@meta.data[[i]])]) != 0)
    sc_UCASpatial@meta.data[[i]] <- NULL
}

DefaultAssay(sc_UCASpatial) <- "RNA"
sc_UCASpatial@assays$RNA@meta.features$'GeneID-2' <- rownames(sc_UCASpatial)
sc_UCASpatial@assays$RNA@meta.features$SYMBOL <- rownames(sc_UCASpatial)
sc_UCASpatial@assays$RNA@data <- sc_UCASpatial@assays$RNA@counts 

sc_UCASpatial.loom <- as.loom(x = sc_UCASpatial, filename = "sc_UCASpatial.loom", verbose = T)
sc_UCASpatial.loom
sc_UCASpatial.loom$close_all()
write.csv(sc_UCASpatial@meta.data,'sc_UCASpatial_meta.data.csv')
write.csv(sc_UCASpatial@assays$RNA@meta.features,'sc_UCASpatial_meta.features.csv')


