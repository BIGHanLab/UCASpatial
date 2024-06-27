## 01 Initialize the library and working directory.
library(Seurat)

setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/')

## 02 Load in the scdata and matched simulated data.
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")
simu_TC <- readRDS("../../Data/Fig2/simu_TC.rds")
simu_TM <- readRDS("../../Data/Fig2/simu_TM.rds")
simu_TS <- readRDS("../../Data/Fig2/simu_TS.rds")


# RCTD
library(RCTD)
library(Matrix) 
library(ggplot2) 
library(ggpubr) 
library(gridExtra) 
library(reshape2) 
library(Seurat)
library(spacexr)
# RCTD
RunRCTD <- function(scdata,simulatData_list,clust_vr = "subclass",n = 50){
  Seurat::Idents(scdata) <- scdata@meta.data[[clust_vr]]
  cell_types <- as.data.frame(scdata@meta.data[[clust_vr]]) %>%
    mutate(subclass = scdata@meta.data[[clust_vr]])
  rownames(cell_types) <- colnames(scdata)
  coords <-
    as.data.frame(
      cbind(
        c(1:n),
        c(1:n)
      )
    )
  rownames(coords) <- paste("Simu",1:n,sep = "")
  levels(Idents(scdata))
  table(scdata$subclass)
  selected_ct <- names(table(cell_types$subclass))[table(cell_types$subclass) > 25]
  # levels(cell_types$subclass)[levels(cell_types$subclass)=="ILC1/NK"] <- "ILC1_NK"
  selected_cell <- WhichCells(scdata,idents = selected_ct)
  
  cell_types <- cell_types[rownames(cell_types)%in%selected_cell,]
  
  table(cell_types$subclass)
  
  cell_types_subclass = as.factor(cell_types$subclass)
  names(cell_types_subclass) <- rownames(cell_types)
  cell_types_subclass <- factor(cell_types_subclass,levels = levels(droplevels(cell_types_subclass)))
  
  reference <-
    Reference(scdata@assays$RNA@counts[,colnames(scdata) %in% selected_cell],
              cell_types = cell_types_subclass)
  
  result_RCTD <- lapply(simulatData_list, function(x){
    test_spot_counts <- as.matrix(x[[1]])
    counts <- test_spot_counts
    colnames(counts) <- paste("Simu",1:n,sep = "")
    spatialRNA <- SpatialRNA(coords = coords, counts = counts)
    myRCTD <-
      create.RCTD(spatialRNA,
                  reference,
                  max_cores = 2,
                  test_mode = F) # here puck is the SpatialRNA object, and reference is the Reference object.
    myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
    return(myRCTD_full)
  })
  return(result_RCTD)
}


simulatData_TME <- c(simu_TC,simu_TM,simu_TS)

RCTD_TME_highres.rds <- RunRCTD(sc.silico,simulatData_TME,clust_vr = "high_res_ident")
saveRDS(RCTD_TME_highres.rds,"RCTD_TME_highres.rds")
