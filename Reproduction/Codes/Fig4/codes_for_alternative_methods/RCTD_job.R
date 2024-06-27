.libPaths(new = '/data/xy/R/R_package/R-3.6.1/R_package/')
library(RCTD)
library(Seurat)
library(ggpubr) 
library(gridExtra) 
library(reshape2) 
library(spacexr)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/')

scdata <-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/integrated_sc_fil_v7.rds")
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

st.list <- list(KRAS_Pt2_TmTc,WT_Pt3_Tm,pt37_3,pt37_7,
                NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4)

clust_vr <- 'UCASpatial_clus_v7'
scdata@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                     x = scdata@meta.data[, clust_vr],
                                     perl = TRUE)

Seurat::Idents(scdata) <- scdata$UCASpatial_clus_v7
cell_types <- as.data.frame(scdata@meta.data[[clust_vr]]) %>%
  mutate(UCASpatial_clus_v7 = scdata@meta.data[[clust_vr]])
rownames(cell_types) <- colnames(scdata)

levels(Idents(scdata))
table(scdata$UCASpatial_clus_v7)
selected_ct <- names(table(cell_types$UCASpatial_clus_v7))[table(cell_types$UCASpatial_clus_v7) > 25]
# levels(cell_types$UCASpatial_clus_v7)[levels(cell_types$UCASpatial_clus_v7)=="ILC1/NK"] <- "ILC1_NK"
selected_cell <- WhichCells(scdata,idents = selected_ct)

cell_types <- cell_types[rownames(cell_types)%in%selected_cell,]

table(cell_types$UCASpatial_clus_v7)

cell_types_UCASpatial_clus_v7 = as.factor(cell_types$UCASpatial_clus_v7)
names(cell_types_UCASpatial_clus_v7) <- rownames(cell_types)
cell_types_UCASpatial_clus_v7 <- factor(cell_types_UCASpatial_clus_v7,levels = levels(droplevels(cell_types_UCASpatial_clus_v7)))

reference <-
  Reference(scdata@assays$RNA@counts[,colnames(scdata) %in% selected_cell],
            cell_types = cell_types_UCASpatial_clus_v7)

# sliceD1
STData_list1 <- list(KRAS_Pt2_TmTc,WT_Pt3_Tm)
# sliceC1
STData_list2 <- list(pt37_3,pt37_7)
# image
# STData_list3 <- list(NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4,st.add1)
STData_list3 <- list(NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4)

result_RCTD1 <- lapply(STData_list1, function(x){
  coords <- as.data.frame(cbind(x@images$sliceD1@coordinates$row,x@images$sliceD1@coordinates$col))
  colnames(coords) <- c('row','col')
  rownames(coords) <- colnames(x)
  test_spot_counts <- x@assays$Spatial@counts
  counts <- test_spot_counts
  spatialRNA <- SpatialRNA(coords = coords, counts = counts)
  myRCTD <-
    create.RCTD(spatialRNA,
                reference,
                max_cores = 2,
                test_mode = F) # here puck is the SpatialRNA object, and reference is the Reference object.
  myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
  return(myRCTD_full)
})
saveRDS(result_RCTD1,'result_RCTD1.rds')

result_RCTD2 <- lapply(STData_list2, function(x){
  coords <- as.data.frame(cbind(x@images$sliceC1@coordinates$row,x@images$sliceC1@coordinates$col))
  colnames(coords) <- c('row','col')
  rownames(coords) <- colnames(x)
  test_spot_counts <- x@assays$Spatial@counts
  counts <- test_spot_counts
  spatialRNA <- SpatialRNA(coords = coords, counts = counts)
  myRCTD <-
    create.RCTD(spatialRNA,
                reference,
                max_cores = 2,
                test_mode = F) # here puck is the SpatialRNA object, and reference is the Reference object.
  myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
  return(myRCTD_full)
})
saveRDS(result_RCTD2,'result_RCTD2.rds')


result_RCTD3 <- lapply(STData_list3, function(x){
  coords <- as.data.frame(cbind(x@images$image@coordinates$row,x@images$image@coordinates$col))
  colnames(coords) <- c('row','col')
  rownames(coords) <- colnames(x)
  test_spot_counts <- x@assays$Spatial@counts
  counts <- test_spot_counts
  spatialRNA <- SpatialRNA(coords = coords, counts = counts)
  myRCTD <-
    create.RCTD(spatialRNA,
                reference,
                max_cores = 2,
                test_mode = F) # here puck is the SpatialRNA object, and reference is the Reference object.
  myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')
  return(myRCTD_full)
})
saveRDS(result_RCTD3,'result_RCTD3.rds')
