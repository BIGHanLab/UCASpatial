library(Seurat)
source('/data/xy/scripts/eWEIDE_final_v1.R')

#### data prepare -------
pt2 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt2_region.RDS')
pt3 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt3_region.RDS')
# pt6 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt6_region.RDS')
pt37_3 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt37_3_region.RDS')
pt37_7 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/Pt37_7_region.RDS')

NC_CRC1 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/CRC1_region.RDS')
NC_CRC2 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/CRC2_region.RDS')
NC_CRC3 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/CRC3_region.RDS')
NC_CRC4 <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/eSpanno_deconv/boundary/CRC4_region.RDS')

pt37_3@images$sliceD1 <- pt37_3@images$sliceC1
pt37_3@images$sliceC1 <- NULL
pt37_7@images$sliceD1 <- pt37_7@images$sliceC1
pt37_7@images$sliceC1 <- NULL
NC_CRC1@images$sliceD1 <- NC_CRC1@images$image
NC_CRC1@images$image <- NULL
NC_CRC2@images$sliceD1 <- NC_CRC2@images$image
NC_CRC2@images$image <- NULL
NC_CRC3@images$sliceD1 <- NC_CRC3@images$image
NC_CRC3@images$image <- NULL
NC_CRC4@images$sliceD1 <- NC_CRC4@images$image
NC_CRC4@images$image <- NULL

st_list <- list(pt2,pt3,pt37_3,pt37_7,NC_CRC1,NC_CRC2,NC_CRC3,NC_CRC4)
st_name <- c('pt2','pt3','pt37_3','pt37_7','NC_CRC1','NC_CRC2','NC_CRC3','NC_CRC4')

#### Visulization -------
setwd('/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/Figures/')
cell_types_all <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                       x=levels(integrated_sc_fil$UCASpatial_clus_v7),perl = TRUE)

# UCASpatial
{
  UCASpatial_v7_res200_pur0.9 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/UCASpatial_test/UCASpatial_v8/UCASpatial_v7_res200_pur0.9.rds")
  
  UCASpatial_results_200.9_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x[[2]]
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('UCASpatial_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },UCASpatial_v7_res200_pur0.9, st_list,st_name)
  saveRDS(UCASpatial_results_200.9_and_draw_pdf,
          'UCASpatial_results.rds')
}

# SPOTlight
{
  SPOTlight_result <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/SPOTlight/KRAS_SPOTlight_result.rds")
  
  deconv.list <- list(SPOTlight_result[[1]][[2]],SPOTlight_result[[2]][[2]],
                      SPOTlight_result[[3]][[2]],SPOTlight_result[[4]][[2]],
                      SPOTlight_result[[5]][[2]],SPOTlight_result[[6]][[2]],
                      SPOTlight_result[[7]][[2]],SPOTlight_result[[8]][[2]])

  SPOTlight_results_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('SPOTlight_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },deconv.list, st_list,st_name)
  saveRDS(SPOTlight_results_and_draw_pdf,
          'SPOTlight_results.rds')
}

# RCTD
{
  deconv_pt2 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD1.rds")[[1]]
  deconv_pt3 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD1.rds")[[2]]
  deconv_pt37_3 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD2.rds")[[1]]
  deconv_pt37_7 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD2.rds")[[2]]
  deconv_CRC1 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD3.rds")[[1]]
  deconv_CRC2 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD3.rds")[[2]]
  deconv_CRC3 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD3.rds")[[3]]
  deconv_CRC4 <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/RCTD/result_RCTD3.rds")[[4]]
  
  deconv.list <- list(as.matrix(normalize_weights(deconv_pt2@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_pt3@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_pt37_3@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_pt37_7@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_CRC1@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_CRC2@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_CRC3@results[["weights"]])),
                      as.matrix(normalize_weights(deconv_CRC4@results[["weights"]])))
  
  RCTD_results_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    # rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('RCTD_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },deconv.list, st_list,st_name)
  saveRDS(RCTD_results_and_draw_pdf,
          'RCTD_results.rds')
  
}


# C2L
{
  deconv_pt2 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/KRAS_Pt2_TmTc_C2L.csv",row.names = 1)
  deconv_pt3 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/WT_Pt3_Tm_C2L.csv",row.names = 1)
  deconv_pt37 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/Pt3_Tc_Pt7_Tm_C2L.csv",row.names = 1)
  deconv_CRC1 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/NC_CRC1_C2L.csv",row.names = 1)
  deconv_CRC2 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/NC_CRC2_C2L.csv",row.names = 1)
  deconv_CRC3 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/NC_CRC3_C2L.csv",row.names = 1)
  deconv_CRC4 <- read.csv("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/cell2location/result/NC_CRC4_C2L.csv",row.names = 1)
  
  deconv.list <- list(deconv_pt2,deconv_pt3,deconv_pt37[rownames(deconv_pt37_3),],
                      deconv_pt37[rownames(deconv_pt37_7),],
                      deconv_CRC1,deconv_CRC2,deconv_CRC3,deconv_CRC4)
  
  deconv.list <- lapply(deconv.list, function(x){
    colnames(x) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                  x=colnames(x))
    return(x)
  })
  cell_types_all <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                         x=levels(integrated_sc_fil$UCASpatial_clus_v7),perl = TRUE)
  
  C2L_results_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('C2L_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },deconv.list, st_list,st_name)
  saveRDS(C2L_results_and_draw_pdf,
          'C2L_results.rds')
}

# CARD
{
  deconv_pt2 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/KRAS_Pt2_TmTc_220812_CARD_res.txt",sep = '\t')
  deconv_pt3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/WT_Pt3_Tm_220812_CARD_res.txt",sep = '\t')
  deconv_pt37_3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/Pt37_3_region_CARD_res.txt",sep = '\t')
  deconv_pt37_7 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/Pt37_7_region_CARD_res.txt",sep = '\t')
  deconv_CRC1 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/NC_CRC1_v1_CARD_res.txt",sep = '\t')
  deconv_CRC2 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/NC_CRC2_v1_CARD_res.txt",sep = '\t')
  deconv_CRC3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/NC_CRC3_v1_CARD_res.txt",sep = '\t')
  deconv_CRC4 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/CARD/NC_CRC4_v1_CARD_res.txt",sep = '\t')
  
  deconv.list <- list(deconv_pt2,deconv_pt3,deconv_pt37_3,deconv_pt37_7,
                      deconv_CRC1,deconv_CRC2,deconv_CRC3,deconv_CRC4)
  
  CARD_results_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('CARD_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },deconv.list, st_list,st_name)
  saveRDS(CARD_results_and_draw_pdf,
          'CARD_results.rds')
}

# stereoscope
{
  deconv_pt2 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ KRAS_Pt2_TmTc_220812 _decon_df.txt",sep = '\t')
  deconv_pt3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ WT_Pt3_Tm_220812 _decon_df.txt",sep = '\t')
  deconv_pt37_3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ Pt37_3_region _decon_df.txt",sep = '\t')
  deconv_pt37_7 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ Pt37_7_region _decon_df.txt",sep = '\t')
  deconv_CRC1 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ CRC1_ST_SeuratObject _decon_df.txt",sep = '\t')
  deconv_CRC2 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ CRC2_ST_SeuratObject _decon_df.txt",sep = '\t')
  deconv_CRC3 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ CRC3_ST_SeuratObject _decon_df.txt",sep = '\t')
  deconv_CRC4 <- read.table("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/stereoscope/ CRC4_ST_SeuratObject _decon_df.txt",sep = '\t')
  
  deconv.list <- list(deconv_pt2,deconv_pt3,deconv_pt37_3,deconv_pt37_7,
                      deconv_CRC1,deconv_CRC2,deconv_CRC3,deconv_CRC4)
  
  stereoscope_results_and_draw_pdf <- mapply(function(x,y,z){
    test <- y
    decon_matr <- x
    test@meta.data <- test@meta.data[,1:3]
    decon_Pt7_mtrx <- as.matrix(decon_matr)
    decon_Pt7_mtrx <- decon_Pt7_mtrx[, colnames(decon_Pt7_mtrx) != "res_ss"]
    decon_Pt7_mtrx <- decon_Pt7_mtrx/rowSums(decon_Pt7_mtrx)
    rownames(decon_Pt7_mtrx) <- colnames(test)
    decon_Pt7_df <- decon_Pt7_mtrx %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
    test@meta.data <- test@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_Pt7_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    Annotation_assay <- CreateAssayObject(t(decon_Pt7_mtrx))
    test@assays$Annotation <- Annotation_assay
    # cell_types_all <- colnames(decon_Pt7_mtrx)
    ggsave(plot = Seurat::SpatialFeaturePlot(
      object = test,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03,ncol=9),filename = paste('stereoscope_',z,'cut0.03.pdf',sep = ''),
      height = 12,width = 20)
    return(test)
  },deconv.list, st_list,st_name)
  saveRDS(stereoscope_results_and_draw_pdf,
          'stereoscope_results.rds')
}





