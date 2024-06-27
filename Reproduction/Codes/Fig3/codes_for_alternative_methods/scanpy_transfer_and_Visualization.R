library(Seurat)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(pacman)
library(tidyr)
library(dplyr)
library(readr)
p_unload(Seurat)
p_unload(SeuratDisk)
p_unload(SeuratObject)
p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig3/codes_for_alternative_methods/')
source('../../source/UCASpatial_final_v1.R')

total_coord_norm <- readRDS("../../../Data/Fig3/MRL_coord.rds")
sc.Reg <- readRDS('../../../Data/Fig3/MRL.scRNA.rds')
load('../../../Data/Fig3/MRL.st.primary.Rdata')
st.Reg <- ear.integrated.filter


source('../../source/my_spatial_featureplot.R')




#### Step 1: transfer scanpy data ####

library(SeuratDisk)
setwd("../../../Data/Fig3/scanpy_data/")

# sc_ref porocess
{
  for(i in colnames(sc.Reg@meta.data))
  {
    if(length(sc.Reg@meta.data[[i]][is.na(sc.Reg@meta.data[[i]])]) != 0)
      sc.Reg@meta.data[[i]] <- NULL
  }
  
  DefaultAssay(sc.Reg) <- "RNA"
  sc.Reg@assays$RNA@meta.features$'GeneID-2' <- rownames(sc.Reg)
  sc.Reg@assays$RNA@meta.features$SYMBOL <- rownames(sc.Reg)
  sc.Reg@assays$RNA@data <- sc.Reg@assays$RNA@counts 
  
  sc.Reg.loom <- as.loom(x = sc.Reg, filename = "sc.Reg.loom", verbose = T)
  sc.Reg.loom
  sc.Reg.loom$close_all()
  write.csv(sc.Reg@meta.data,'sc.Reg_meta.data.csv')
  write.csv(sc.Reg@assays$RNA@meta.features,'sc.Reg_meta.features.csv')
}

# st_vis process
{
  name <- 'Reg'
  st.scanpy <- st.Reg
  st.scanpy$in_tissue <- 1
  st.scanpy$array_row <- total_coord_norm$row
  st.scanpy$array_col <- total_coord_norm$col
  st.scanpy$sample <- name
  
  DefaultAssay(st.scanpy) <- "Spatial"
  st.scanpy@assays$Spatial@meta.features$gene_ids <- rownames(st.scanpy)
  st.scanpy@assays$Spatial@meta.features$SYMBOL <- rownames(st.scanpy)
  st.scanpy@assays$Spatial@data <- st.scanpy@assays$Spatial@counts
  SeuratDisk::SaveH5Seurat(st.scanpy, filename = paste(name,"_scanpy_Spatial.h5Seurat",sep = ''),overwrite = T)
  SeuratDisk::Convert(paste(name,"_scanpy_Spatial.h5Seurat",sep = ''), dest = "h5ad",overwrite = T)
  
}

#### Visulization ####

st.vis <- ear.integrated.filter
DefaultAssay(st.vis) <- 'Spatial'
st.vis@assays$eWEIDE <- NULL
st.vis@assays$cell_type <- NULL
st.vis@assays$cell_subtype <- NULL

# Visulization UCASpatial
{
  UCASpatial_reg <- readRDS("../../../Data/Fig3/UCASpatial/UCASpatial_regene.rds")
  decon_mtrx <- UCASpatial_reg[[2]][,-21]
  rownames(decon_mtrx) <- names(st.vis$orig.ident)
  UCASpatial <- st.vis
  UCASpatial@meta.data <- UCASpatial@meta.data[1]
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  UCASpatial@meta.data <- UCASpatial@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  ct.all <- colnames(decon_mtrx)
  
  pdf('UCASpatial_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = UCASpatial,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
  
}

# Visulization SPOTlight
{
  SPOTlight_result <- readRDS("../../../Data/Fig3/SPOTlight/KRAS_SPOTlight_result.rds")
  decon_mtrx <- SPOTlight_result[[2]][,-21]
  rownames(decon_mtrx) <- names(st.vis$orig.ident)
  SPOTlight <- st.vis
  SPOTlight@meta.data <- SPOTlight@meta.data[1]
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  SPOTlight@meta.data <- SPOTlight@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  ct.all <- colnames(decon_mtrx)
  
  pdf('SPOTlight_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = SPOTlight,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
}

# Visulization RCTD
{
  result_RCTD1 <- readRDS("../../../Data/Fig3/RCTD/result_RCTD.rds")
  decon_mtrx <- as.matrix(result_RCTD1@results[["weights"]])
  
  # rownames(decon_mtrx) <- names(st.vis$orig.ident)
  RCTD <- st.vis
  RCTD@meta.data <- RCTD@meta.data[1]
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  RCTD@meta.data <- RCTD@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  RCTD@meta.data[is.na(RCTD@meta.data)] <- 0
  ct.all <- gsub(colnames(decon_mtrx),pattern = '-| |[+]',replacement = '.')
  pdf('RCTD_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = RCTD,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
}

# Visulization cell2location
{
  C2L_pt37 <- read.csv("../../../Data/Fig3/cell2location/result/Reg_C2L.csv",row.names = 1)
  C2L_pt37 <- C2L_pt37/rowSums(C2L_pt37)
  rownames(C2L_pt37) <- names(st.vis$orig.ident)
  C2L <- st.vis
  C2L@meta.data <- C2L@meta.data[1]
  C2L_decon_df <- C2L_pt37 %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  C2L@meta.data <- C2L@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(C2L_decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  ct.all <- colnames(C2L_pt37)
  
  pdf('cell2location_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = C2L,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
  
}

# Visulization CARD
{
  CARD_pt37 <- read.table("../../../Data/Fig3/CARD/CRAD_regeneration_decon_result.txt",sep = '\t')
  rownames(CARD_pt37) <- names(st.vis$orig.ident)
  CARD <- st.vis
  CARD@meta.data <- CARD@meta.data[1]
  CARD_decon_df <- CARD_pt37 %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  CARD@meta.data <- CARD@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(CARD_decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  ct.all <- colnames(C2L_pt37)
  pdf('CARD_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = CARD,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
  
}

# Visulization stereoscope
{
  stereoscope_pt37 <- read.table("../../../Data/Fig3/stereoscope/stereo_decon_df.txt",sep = '\t')
  rownames(stereoscope_pt37) <- names(st.vis$orig.ident)
  stereoscope <- st.vis
  stereoscope@meta.data <- stereoscope@meta.data[1]
  stereoscope_decon_df <- stereoscope_pt37 %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  stereoscope@meta.data <- stereoscope@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(stereoscope_decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  # saveRDS(stereoscope,'/data/xy/Spatial_transcriptome/eWEIDE/Final_version_data/deconv_results/stereoscope_pt37.rds')
  ct.all <- colnames(stereoscope_pt37)
  pdf('stereoscope_regeneration.pdf',width= 11.5, height= 7.88)
  lapply(ct.all, function(x){
    my_spatial_featureplot(spatial_data = stereoscope,data_coord = total_coord_norm,
                           features = x,min.cutoff = 0.03)
  })
  dev.off()
  
}




