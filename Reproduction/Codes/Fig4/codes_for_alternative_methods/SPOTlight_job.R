library(Seurat)
library(SPOTlight)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/SPOTlight/')

sc.use <-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/integrated_sc_fil_v7.rds")


# Kras_Pt8_Tm <- readRDS("/data/xy/Spatial_transcriptome/Kras_data/Single_cell_data/20210912/Kras_Pt8_Tm_sub_latest.rds")
WT_Pt3_Tm <- readRDS('/data/xy/Spatial_transcriptome/Kras_data/00_Latest_Data/WT_Pt3_Tm_220812.rds')
KRAS_Pt2_TmTc <- readRDS("/data/xy/Spatial_transcriptome/Kras_data/00_Latest_Data/KRAS_Pt2_TmTc_220812.rds")
pt37_3 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt37_3_region.RDS')
pt37_7 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt37_7_region.RDS')
# load("/data/xy/Spatial_transcriptome/Kras_data/RowSpatialData/WT_Pt6_TcTm.rdata")

NC_CRC1 <- readRDS("/data/huangzr/Spatial/Public/2022_nc_crc/Analysis/NC_CRC1_v1.RDS")
NC_CRC2 <- readRDS("/data/huangzr/Spatial/Public/2022_nc_crc/Analysis/NC_CRC2_v1.RDS")
NC_CRC3 <- readRDS("/data/huangzr/Spatial/Public/2022_nc_crc/Analysis/NC_CRC3_v1.RDS")
NC_CRC4 <- readRDS("/data/huangzr/Spatial/Public/2022_nc_crc/Analysis/NC_CRC4_v1.RDS")
# st.add1 <- readRDS('/data/huangzr/Spatial/Kras/version4/data_add/CRC_add_1.RDS')

st.list <- list(KRAS_Pt2_TmTc,WT_Pt3_Tm,pt37_3,pt37_7,
                NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4)

Idents(sc.use) <- sc.use$UCASpatial_clus_v7
SPOTlight_marker_KRAS <- Seurat::FindAllMarkers(
  object = sc.use,
  assay = "RNA",
  slot = "data",
  min.pct = 0.2,
  only.pos = TRUE,
  logfc.threshold = 0.25
)
saveRDS(SPOTlight_marker_KRAS,'SPOTlight_marker.rds')

SPOTlight_marker_KRAS %>% dplyr::count(cluster)

KRAS_SPOTlight_result <- lapply(st.list, function(x){
  result <- spotlight_deconvolution(
    se_sc = sc.use,
    counts_spatial = x@assays$Spatial@counts,
    clust_vr = 'UCASpatial_clus_v7',
    cluster_markers = SPOTlight_marker_KRAS,
    cl_n = 100,
    hvg = 3000,
    ntop = NULL,
    transf = "uv",
    method = "nsNMF",
    min_cont = 0,
    assay = "RNA",
    slot = "counts"
  )
  return(result)
})
saveRDS(KRAS_SPOTlight_result,'KRAS_SPOTlight_result.rds')