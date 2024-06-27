library(Seurat)
library(RColorBrewer)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig4/')

### Data integration ####
CRC_obj<-readRDS("../../Data/Fig4_5/CRC_sc_list.RDS")
all.genes<-readRDS("../../Data/Fig4_5/all.genes.RDS")
features <- SelectIntegrationFeatures(object.list = CRC_obj,nfeatures = 2000)
CRC_ICB.anchors <- FindIntegrationAnchors(object.list = CRC_obj, anchor.features = features)
CRC_ICB.combined <- IntegrateData(anchorset = CRC_ICB.anchors,normalization.method = "SCT", verbose = TRUE, features.to.integrate = all_genes)
saveRDS(CRC_ICB.combined,'integrated_data.RDS')

integrated_sc<-readRDS("../../Data/Fig4_5/integrated_data.RDS")
DefaultAssay(integrated_sc) <- "integrated"

integrated_sc <- ScaleData(integrated_sc, verbose = FALSE)
integrated_sc <- RunPCA(integrated_sc, npcs = 50, verbose = FALSE)
ElbowPlot(integrated_sc)
#VizDimLoadings(integrated_sc, dims = 25:30, reduction = "pca")
DefaultAssay(integrated_sc) <- "integrated"
integrated_sc <- RunUMAP(integrated_sc, reduction = "pca", dims = 1:20)
integrated_sc <- FindNeighbors(integrated_sc, reduction = "pca", dims = 1:20)
integrated_sc <- FindClusters(integrated_sc, resolution = 1)

### first_around_annotation ####
integrated_sc@meta.data$first_around_anno<-"Remain"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 0] <- "CD8_gd_NK"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 1] <- "CD4"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 2] <- "B"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 3] <- "Treg"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 4] <- "Plasma"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 5] <- "Plasma"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 6] <- "Mast"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 7] <- "Epithelial c1"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 8] <- "APOE-Fibroblast"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 9] <- "Mac-like"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 10] <- "Th1-like"

integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 11] <- "Epithelial c2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 12] <- "Epithelial c3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 13] <- "Myofib"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 14] <- "Endo"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 15] <- "Fibroblast c3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 16] <- "OGN Fib"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 17] <- "Mac+ cDC2 Neu"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 18] <- "Goblet-like"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 19] <- "Prolif T"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 20] <- "Fib (pericyte-like)"

integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 21] <- "pericyte"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 22] <- "Endo c2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 23] <- "CD45-CD11b-CD14+ cell + prolif Myeloid"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 24] <- "cDC2"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 25] <- "ILC"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 26] <- "Glial"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 27] <- "Endo_3"
integrated_sc@meta.data$first_around_anno[integrated_sc@meta.data$seurat_clusters == 28] <- "Fib c4"

DimPlot(integrated_sc,group.by = 'first_around_anno',label = T,repel = T)

### detailed annotation ####

# Fib
{
  Fib <- integrated_sc[,integrated_sc$first_around_anno %in% c("APOE-Fibroblast","Fib (pericyte-like)","Fib c4","Fibroblast c3","Myofib","OGN Fib","pericyte")]
  DefaultAssay(Fib)<-"RNA"
  Idents(Fib)<-Fib$first_around_anno
  table(Fib$first_around_anno)
  Fib_marker<-FindAllMarkers(Fib,logfc.threshold = 0.3,only.pos = T)
  View(Fib_marker)
  Fib@meta.data$second_around_anno<-"Remain"
  Fib@meta.data$second_around_anno[Fib@meta.data$first_around_anno %in% c("APOE-Fibroblast","OGN Fib")]<-"Fibroblast"
  Fib@meta.data$second_around_anno[Fib@meta.data$first_around_anno %in% c("Fib (pericyte-like)","Fib c4")]<-"Smooth Muscle"
  Fib@meta.data$second_around_anno[Fib@meta.data$first_around_anno == "Fibroblast c3"]<-"POSTN Myofib"
  Fib@meta.data$second_around_anno[Fib@meta.data$first_around_anno == "Myofib"]<-"FAP Myofib"
  Fib@meta.data$second_around_anno[Fib@meta.data$first_around_anno == "pericyte"]<-"Pericyte"
  
}

# T cell
{
  Tcell<-integrated_sc[,integrated_sc$first_around_anno %in% c("CD4","Treg","Th1-like","CD8_gd_NK")]
  Kras_Tcell<-Tcell[,Tcell@meta.data$orig.ident %in% c("KRAS_Pt2_immune","KRAS_Pt7_immune","WT_Pt3_immune")]
  Ye_Tcell<-Tcell[,!(Tcell@meta.data$orig.ident %in% c("KRAS_Pt2_immune","KRAS_Pt7_immune","WT_Pt3_immune"))]
  CRC_T_list<-list(Kras_Tcell,Ye_Tcell)
  features <- SelectIntegrationFeatures(object.list = CRC_T_list,nfeatures = 2000)
  CRC_ICB.anchors <- FindIntegrationAnchors(object.list = CRC_T_list, anchor.features = features)
  CRC_ICB.combined <- IntegrateData(anchorset = CRC_ICB.anchors)
  T_integrated<-CRC_ICB.combined
  
  library(RColorBrewer)
  T_integrated<-readRDS("../../Data/Fig4_5/T_integrated_data.RDS")
  DefaultAssay(T_integrated) <- "integrated"
  
  T_integrated <- ScaleData(T_integrated, verbose = FALSE)
  T_integrated <- RunPCA(T_integrated, npcs = 50, verbose = FALSE)
  ElbowPlot(T_integrated)
  T_integrated <- RunUMAP(T_integrated, reduction = "pca", dims = 1:15)
  T_integrated <- FindNeighbors(T_integrated, reduction = "pca", dims = 1:15)
  T_integrated <- FindClusters(T_integrated, resolution = 0.4)
  
  
  T_integrated@meta.data$second_around_anno_new_2<-"Remain"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 0] <- "CD8 Tem"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 1] <- "CD4 Tm"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 2] <- "CD4 Treg"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 3] <- "CD4 Tm"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 4] <- "CD4 Tn"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 5] <- "CD8 Trm"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 6] <- "CD8 Tex"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 7] <- "CD8 Trm"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 8] <- "NK"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 9] <- "Th1/Th17"
  T_integrated@meta.data$second_around_anno_new_2[T_integrated@meta.data$seurat_clusters == 10] <- "Tfh"
  p1<-DimPlot(T_integrated,group.by = "second_around_anno_new_2",cols = c(brewer.pal(12,"Paired"),brewer.pal(12,"Set1"),brewer.pal(12,"Set2")),label = T)
  p1
  saveRDS(T_integrated,file = "T_new.RDS")
  T_new<-readRDS("../../Data/Fig4_5/T_new.RDS")
  
}

# myeloid
{
  reintegrated_sc_myeloid <- readRDS("../../Data/Fig4_5/myeloid_integrated_data.RDS")
  DefaultAssay(reintegrated_sc_myeloid) <- 'integrated'
  
  #### recluster
  sc_CRC_mereintegrated_sc_myeloidrge <- FindVariableFeatures(reintegrated_sc_myeloid, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(reintegrated_sc_myeloid)
  reintegrated_sc_myeloid <- ScaleData(reintegrated_sc_myeloid, features = all.genes)
  reintegrated_sc_myeloid <- RunPCA(reintegrated_sc_myeloid, features = VariableFeatures(object = sc_CRC_merge))
  
  reintegrated_sc_myeloid <- FindNeighbors(reintegrated_sc_myeloid, dims = 1:30)
  reintegrated_sc_myeloid <- FindClusters(reintegrated_sc_myeloid, resolution = seq(0.1,1,0.1))
  reintegrated_sc_myeloid <- RunUMAP(reintegrated_sc_myeloid, dims = 1:30)
  
  DimPlot(reintegrated_sc_myeloid,group.by = 'integrated_snn_res.0.4')
  
  FeaturePlot(reintegrated_sc_myeloid,
              features = c('FOLR2','CD163L1','FCN1','THBS1',
                           'CD1C','FCER1A', 'CLEC10A','SPP1','FN1',
                           'CLEC9A','XCR1','LAMP3','CCR7',
                           'CSF3R','S100A8','PDPN','CXCL9','CXCL10','IFNAR1'))
  
  Idents(reintegrated_sc_myeloid) <- reintegrated_sc_myeloid$integrated_snn_res.0.4
  mk_mye_v2 <- FindAllMarkers(reintegrated_sc_myeloid,min.pct = 0.2,logfc.threshold = 0.25,min.diff.pct = 0.1,
                              only.pos = T)
  
  table(reintegrated_sc_myeloid$myeloid_Minor,reintegrated_sc_myeloid$integrated_snn_res.0.4)
  
  DefaultAssay(reintegrated_sc_myeloid) <- 'RNA'
  VlnPlot(reintegrated_sc_myeloid,group.by = 'integrated_snn_res.0.4',pt.size = 0,
          features = c('PTPRC','MKI67','LAMP3','S100A8','ITGAM','CCR7','CLEC9A',
                       'ITGAX','CD1C','SPP1','CSF3R','LILRB2'))
  
  reintegrated_sc_myeloid$myeloid_Minor <- droplevels(reintegrated_sc_myeloid$integrated_snn_res.0.4)
  
  levels(reintegrated_sc_myeloid$myeloid_Minor)[1] <- 'C1QC+ Mac'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[2] <- 'Monocyte'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[3] <- 'cDC2'      
  levels(reintegrated_sc_myeloid$myeloid_Minor)[4] <- 'CD45- myeloid'      #--> filter
  levels(reintegrated_sc_myeloid$myeloid_Minor)[5] <- 'FN1+ Mac'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[6] <- 'PDPN+ Mac' 
  levels(reintegrated_sc_myeloid$myeloid_Minor)[7] <- 'proliferation myeloid' #--> filter
  levels(reintegrated_sc_myeloid$myeloid_Minor)[8] <- 'cDC1'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[9] <- 'Neutrophil'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[10] <- 'cDC3'
  levels(reintegrated_sc_myeloid$myeloid_Minor)[11] <- 'pDC'
  
  levels(reintegrated_sc_myeloid$myeloid_Minor)
  DimPlot(reintegrated_sc_myeloid,group.by = 'myeloid_Minor',label = T)
  DotPlot(reintegrated_sc_myeloid,group.by = 'myeloid_Minor',features = c('MKI67','PTPRC'))
  
  
  # filter prolif myeloid and CD45- myeloid
  reintegrated_sc_myeloid_fil <- subset(reintegrated_sc_myeloid,subset= 
                                          integrated_snn_res.0.4 %in% c(0,1,2,4,5,7:10))
  reintegrated_sc_myeloid_fil <- RunUMAP(reintegrated_sc_myeloid_fil, dims = 1:30)
  DimPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4')
  
  mk_mye_v2 <- FindAllMarkers(reintegrated_sc_myeloid_fil,min.pct = 0.2,logfc.threshold = 0.25,min.diff.pct = 0.1,
                              only.pos = T)
  
  VlnPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4',pt.size = 0,
          features = c('PTPRC','MKI67','LAMP3','S100A8','ITGAM','CCR7','CLEC9A',
                       'ITGAX','CD1C','SPP1','CSF3R','LILRB2'))
  
  Dotplot(DotPlot(reintegrated_sc_myeloid_fil,group.by = 'integrated_snn_res.0.4',
                  features = c('SELENOP','IGF1','FOLR2','CD163L1','FCN1','THBS1',
                               'CD1C','FCER1A', 'CLEC10A','SPP1','FN1',
                               'CLEC9A','XCR1','LAMP3','CCR7','GZMB','IL3RA' ,'LILRB4' ,'CLEC4C',
                               'CSF3R','S100A8','PDPN','CXCL9','CXCL10','IFNAR1'))$data)
  
  reintegrated_sc_myeloid_fil$myeloid_Minor <- droplevels(reintegrated_sc_myeloid_fil$integrated_snn_res.0.4)
  
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[1] <- 'C1QC+ Mac'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[2] <- 'Monocyte'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[3] <- 'cDC2'      
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[4] <- 'FN1+ Mac'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[5] <- 'PDPN+ Mac' 
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[6] <- 'cDC1'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[7] <- 'Neutrophil'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[8] <- 'cDC3'
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)[9] <- 'pDC'
  
  levels(reintegrated_sc_myeloid_fil$myeloid_Minor)
}

# endothelial cell
{
  integrated_sc_endo <- subset(integrated_sc,subset = first_around_anno %in% 
                                 c("Endo",'Endo c2','Endo_3'))
  integrated_sc_endo$Endo_major <- "Endothelium"
}

# epithelial cell
{
  
  integrated_sc_epi <- subset(integrated_sc,subset = first_around_anno %in% 
                                c("Epithelial c1","Epithelial c2",
                                  'Epithelial c3','Goblet-like'))
  
  DefaultAssay(integrated_sc_epi) <- 'integrated'
  integrated_sc_epi <- FindNeighbors(integrated_sc_epi, dims = 1:30)
  integrated_sc_epi <- FindClusters(integrated_sc_epi, resolution = 0.05)
  integrated_sc_epi <- RunUMAP(integrated_sc_epi, dims = 1:30)
  
  table(integrated_sc_epi$integrated_snn_res.0.05)
  DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.05')
  
  p1 <- DimPlot(integrated_sc_epi,group.by = 'integrated_snn_res.0.05')
  
  p1+DimPlot(integrated_sc_epi,group.by = 'orig.ident')
  
  DefaultAssay(integrated_sc_epi) <- 'RNA'
  Idents(integrated_sc_epi) <- integrated_sc_epi$integrated_snn_res.0.05
  mk_epi <- FindAllMarkers(integrated_sc_epi,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
  
  FeaturePlot(integrated_sc_epi,features = c('nCount_RNA','nFeature_RNA'),max.cutoff = 1000)
  p1+FeaturePlot(integrated_sc_epi,features = c('LEFTY1','TM4SF1','SPINK4','MUC2','REG1B','REG1A'))
  
  integrated_sc_epi$Epi_Minor_2 <- integrated_sc_epi$integrated_snn_res.0.05
  table(integrated_sc_epi$Epi_Minor_2)
  
  levels(integrated_sc_epi$Epi_Minor_2)[1] <- 'Epithelium c1' ## DPEP1
  levels(integrated_sc_epi$Epi_Minor_2)[2] <- 'LEFTY1+ Epithelium c2' ## LEFTY1
  levels(integrated_sc_epi$Epi_Minor_2)[3] <- 'TM4SF1+ Epithelium c3' ## NNMT
  levels(integrated_sc_epi$Epi_Minor_2)[4] <- 'MUC2+ Epithelium c4' ## LEFTY1
  levels(integrated_sc_epi$Epi_Minor_2)[5] <- 'REG1A+ Epithelium c5' ## NNMT
  
  DimPlot(integrated_sc_epi,group.by = 'Epi_Minor_2',label = T,repel = T)+
    FeaturePlot(integrated_sc_epi,features = c('LEFTY1','TM4SF1','SPINK4','MUC2','REG1B','REG1A'))
  saveRDS(integrated_sc_epi,'integrated_sc_epi.rds')
}


### Combination
{
  integrated_sc$UCASpatial_clus_v7 <- integrated_sc$first_around_anno
  
  ### T
  common_cells <- intersect(colnames(integrated_sc), colnames(T_integrated))
  integrated_sc$UCASpatial_clus_v7 <- as.character(integrated_sc$UCASpatial_clus_v7)
  integrated_sc$UCASpatial_clus_v7[common_cells] <- as.character(T_integrated$second_around_anno_new_2)
  table(integrated_sc$UCASpatial_clus_v7)
  
  ### myeloid
  common_cells <- intersect(colnames(integrated_sc), colnames(reintegrated_sc_myeloid_fil))
  integrated_sc$UCASpatial_clus_v7 <- as.character(integrated_sc$UCASpatial_clus_v7)
  integrated_sc$UCASpatial_clus_v7[common_cells] <- as.character(reintegrated_sc_myeloid_fil$myeloid_Minor)
  table(integrated_sc$UCASpatial_clus_v7)
  
  ### endo
  common_cells <- intersect(colnames(integrated_sc), colnames(integrated_sc_endo))
  integrated_sc$UCASpatial_clus_v7 <- as.character(integrated_sc$UCASpatial_clus_v7)
  integrated_sc$UCASpatial_clus_v7[common_cells] <- as.character(integrated_sc_endo$Endo_major)
  table(integrated_sc$UCASpatial_clus_v7)
  
  ### Epi
  common_cells <- intersect(colnames(integrated_sc), colnames(integrated_sc_epi))
  integrated_sc$UCASpatial_clus_v7 <- as.character(integrated_sc$UCASpatial_clus_v7)
  integrated_sc$UCASpatial_clus_v7[common_cells] <- as.character(integrated_sc_epi$Epi_Minor_2)
  table(integrated_sc$UCASpatial_clus_v7)
  
  ### Fibro/Myofibro
  integrated_sc@meta.data$UCASpatial_clus_v7[integrated_sc@meta.data$first_around_anno %in% c("APOE-Fibroblast","OGN Fib")]<-"Fibroblast"
  integrated_sc@meta.data$UCASpatial_clus_v7[integrated_sc@meta.data$first_around_anno %in% c("Fib (pericyte-like)","Fib c4")]<-"Smooth Muscle"
  integrated_sc@meta.data$UCASpatial_clus_v7[integrated_sc@meta.data$first_around_anno == "Fibroblast c3"]<-"POSTN+ Myofib"
  integrated_sc@meta.data$UCASpatial_clus_v7[integrated_sc@meta.data$first_around_anno == "Myofib"]<-"FAP+ Myofib"
  integrated_sc@meta.data$UCASpatial_clus_v7[integrated_sc@meta.data$first_around_anno == "pericyte"]<-"Pericyte"
  table(integrated_sc$UCASpatial_clus_v7)
  
  
  DimPlot(integrated_sc,group.by = 'UCASpatial_clus_v7')
  names(table(integrated_sc$UCASpatial_clus_v7))
  integrated_sc$UCASpatial_clus_v7 <- factor(integrated_sc$UCASpatial_clus_v7,
                                                    levels = c("CD8 Tem","CD8 Tex","CD8 Trm",
                                                               "CD4 Tn","CD4 Tm","Tfh","Th1/Th17",'Treg',
                                                               "NK","ILC",
                                                               "B","Plasma",
                                                               "Mast",
                                                               "Monocyte","Neutrophil",
                                                               "C1QC+ Mac",'FN1+ Mac','PDPN+ Mac',
                                                               "cDC1","cDC2","cDC3","pDC",
                                                               "Endothelial","Glial",
                                                               "Fibroblast","FAP+ Myofib","POSTN+ Myofib",
                                                               "Smooth Muscle","Pericyte",
                                                               'Epithelial c1','Epithelial c2','Epithelial c3',
                                                               'Epithelial c4','Epithelial c5'))
  DimPlot(integrated_sc,group.by = 'UCASpatial_clus_v7',label = T,repel = T)
  saveRDS(integrated_sc,'integrated_sc_fil_v7.rds')
}


### Fig 4b visualization ####
{
  integrated_sc <- readRDS('../../Data/Fig4_5/integrated_sc_fil_v7.rds')
  p <- DimPlot(integrated_sc,group.by = "UCASpatial_clus_v7")
  head(p$data)
  all_cluster.umap <- p$data
  head(all_cluster.umap)
  colnames(all_cluster.umap) <- c('UMAP_1','UMAP_2','Subset')
  
  colors=c("#F9D4C7","#F5B09B","#F58773","#F16A54","#E83F3B","#C7242D","#922629","#5A1312","#F7A06D","#E2662C","#D3B67E","#EBDAB6",
           "#983D2B","#E1DEE9","#D0E1EE","#B6CEE4","#8CBED2","#62A2C6","#448AB6","#2A70AB","#1B5191","#2A3864","#82522C","#B47C3C","#CCCFE2",
           "#A6ABC9","#7779B0","#64569C","#403066","#B7E0DA","#71C0B6","#408F8A","#1F655B","#143930")
  
  library(ggthemes)
  library(ggrastr)
  ggplot(all_cluster.umap,aes(x=UMAP_1,y=UMAP_2,fill=Subset))+
    geom_point_rast(shape=21,color='black',size=0.5,stroke = NA)+
    theme_base(base_size = 16)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          panel.background	 = element_blank())+
    guides(fill=guide_legend(ncol=2,override.aes =list(size = 3)))+
    scale_fill_manual(values = colors)+
    coord_fixed(ratio = 0.8)+
    ggtitle(' ')+
    labs(x='UAMP 1',y='UMAP 2',fill=NULL)
  ggsave('Fig2b.pdf')
}
