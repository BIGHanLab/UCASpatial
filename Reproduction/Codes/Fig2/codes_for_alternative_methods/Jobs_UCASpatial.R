## 01 Initialize the library and working directory.

library(Seurat)
source('/data/xy/scripts/UCASpatial_final_v1.R')

setwd('/data/xy/Spatial_transcriptome/UCASpatial/Reproduction/Codes/Fig2/')

## 02 Load in the scdata and matched simulated data.
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")
simu_TC <- readRDS("../../Data/Fig2/simu_TC.rds")
simu_TM <- readRDS("../../Data/Fig2/simu_TM.rds")
simu_TS <- readRDS("../../Data/Fig2/simu_TS.rds")

# UCASpatial

RunUCASpatial <- function(scdata,simulatData_list,clust_vr,meta.assay = 'integrated',
                      meta.resolution = 100,meta.purity = 0.95,
                      ent.filter.threshold = 0.5,cos.filter.threshold = 0.05,weight.filter.threshold = 0.2){
  Seurat::Idents(scdata) <- scdata@meta.data[[clust_vr]]
  UCASpatial_deconv_result <- lapply(simulatData_list, function(x){
    st_vis = CreateSeuratObject(counts = x[[1]],assay = 'Spatial')
    result <- UCASpatial_deconv(
      sc_ref = scdata,
      st_vis = st_vis,
      clust_vr = clust_vr,
      meta.assay = meta.assay,
      meta.resolution = meta.resolution,
      meta.purity = meta.purity,
      cos.filter = cos.filter.threshold,
      ent.filter.threshold = ent.filter.threshold,
      weight.filter.threshold = weight.filter.threshold
    )
    return(result)
  }) 
  return(UCASpatial_deconv_result)
}

simulatData_TME <- c(simu_TC,simu_TM,simu_TS)

UCASpatial_TME_highres.rds <- RunUCASpatial(sc.silico,simulatData_TME,clust_vr = "high_res_ident")
saveRDS(UCASpatial_TME_highres.rds,"UCASpatial_TME_highres.rds")


