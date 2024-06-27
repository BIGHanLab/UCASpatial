## 01 Initialize the library and working directory.

library(Seurat)
library(SPOTlight)

setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/')

## 02 Load in the scdata and matched simulated data.
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")
simu_TC <- readRDS("../../Data/Fig2/simu_TC.rds")
simu_TM <- readRDS("../../Data/Fig2/simu_TM.rds")
simu_TS <- readRDS("../../Data/Fig2/simu_TS.rds")


# SPOTlight
RunSPOTlight <- function(scdata,simulatData_list,clust_vr = "cluster_last_version"){
  Seurat::Idents(scdata) <- scdata@meta.data[[clust_vr]]
  marker_genes <- Seurat::FindAllMarkers(object = scdata,
                                         assay = "RNA",
                                         slot = "data",
                                         min.pct = 0.2,
                                         only.pos = TRUE,
                                         logfc.threshold = 0.25)
  
  marker_genes %>% dplyr::count(cluster)
  spotlight_ls <- lapply(simulatData_list, function(x){
    result <- spotlight_deconvolution(se_sc = scdata,
                                      counts_spatial = x[[1]],
                                      clust_vr = clust_vr,
                                      cluster_markers = marker_genes,
                                      cl_n = 100,
                                      hvg = 3000,
                                      ntop = NULL,
                                      transf = "uv",
                                      method = "nsNMF",
                                      min_cont = 0,
                                      assay = "RNA",
                                      slot = "counts")
    return(result)
  }) 
  return(spotlight_ls)
}

simulatData_TME <- c(simu_TC,simu_TM,simu_TS)

SPOTlight_TME_highres.rds <- RunSPOTlight(sc.silico,simulatData_TME,clust_vr = "high_res_ident")
saveRDS(SPOTlight_TME_highres.rds,"SPOTlight_TME_highres.rds")
