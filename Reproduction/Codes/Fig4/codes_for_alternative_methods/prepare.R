library(Seurat)



sc.use <-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20231117_Spatial_Annotation_Data/sc.silico.RDS")

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
st.add1 <- readRDS('/data/huangzr/Spatial/Kras/version4/data_add/CRC_add_1.RDS')

st.list <- list(KRAS_Pt2_TmTc,WT_Pt3_Tm,pt37_3,pt37_7,
                NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4,st.add1)
