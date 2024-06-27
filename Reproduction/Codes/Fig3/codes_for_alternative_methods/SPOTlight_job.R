library(Seurat)
library(SPOTlight)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig3/codes_for_alternative_methods/')

total_coord_norm <- readRDS("../../../Data/Fig3/MRL_coord.rds")
sc.Reg <- readRDS('../../../Data/Fig3/MRL.scRNA.rds')
load('../../../Data/Fig3/MRL.st.primary.Rdata')
st.Reg <- ear.integrated.filter
sc.use <-sc.Reg

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

KRAS_SPOTlight_result <- spotlight_deconvolution(
    se_sc = sc.use,
    counts_spatial = st.Reg@assays$Spatial@counts,
    clust_vr = 'Minor_class_reduction',
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
saveRDS(KRAS_SPOTlight_result,'KRAS_SPOTlight_result.rds')