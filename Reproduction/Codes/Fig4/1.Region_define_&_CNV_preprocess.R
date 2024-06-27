# Region define&CNV preprocess

library(AnnotationDbi)
library(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig4/')
UCASpatial_list <- readRDS("../../Data/Fig4_5/UCASpatial_results_200.9_v7.rds")
st_list <- c("CRC1","CRC3a","CRC3b","CRC2","CRC5","CRC4","CRC6","CRC7")
UCASpatial_list
# Epithelial and non_epithelial region define
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Data/Fig4_5/clone_region/')
Epithelial_cutoff=0.35
data_name=st_list
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list=AnnotationDbi::select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
list<-list[,c(1,3)]
i=1
for (i in 1:8){
  Add_st<-UCASpatial_list[[i]]
  new_data_name=st_list[i]
  Add_st@meta.data$Epithelial_sum <- colSums(Add_st@assays$Annotation@data[c("Epithelium.c1","LEFTY1..Epithelium.c2","TM4SF1..Epithelium.c3",
                                                                             "MUC2..Epithelium.c4","REG1A..Epithelium.c5"),])
  #plot(density(Add_st@meta.data$Epithelial_sum))
  #SpatialPlot(Add_st,features = c("Epithelial_sum"))
  Add_st@meta.data$Epithelial_class <-"Low_Epithelial" 
  Add_st@meta.data$Epithelial_class[Add_st@meta.data$Epithelial_sum > Epithelial_cutoff] <-"High_Epithelial" 
  #SpatialPlot(Add_st,group.by = c("Epithelial_class"))
  
  # Epithelial Region filter
  Epi_id <- rownames(Add_st@meta.data)[Add_st@meta.data$Epithelial_class=="High_Epithelial"]
  
  Epi_cor<-Add_st@images$sliceD1@coordinates[Epi_id,c(2,3)]
  Epi_cor$id<-paste(Epi_cor$row,Epi_cor$col,sep="_")
  Epi_filter<-as.vector(rep(FALSE,nrow(Epi_cor)))
  for (i in 1:nrow(Epi_cor)){
    Epi_filter[i]<-sum(paste(Epi_cor$row[i]+1,Epi_cor$col[i]-1,sep="_") %in% Epi_cor$id,
                       paste(Epi_cor$row[i]-1,Epi_cor$col[i]-1,sep="_") %in% Epi_cor$id,
                       paste(Epi_cor$row[i],Epi_cor$col[i]-2,sep="_") %in% Epi_cor$id,
                       paste(Epi_cor$row[i],Epi_cor$col[i]+2,sep="_") %in% Epi_cor$id,
                       paste(Epi_cor$row[i]+1,Epi_cor$col[i]+1,sep="_") %in% Epi_cor$id,
                       paste(Epi_cor$row[i]-1,Epi_cor$col[i]+1,sep="_") %in% Epi_cor$id)>=2
  }
  Epi_cor<-Epi_cor[Epi_filter,]
  Add_st@meta.data$Epithelial_class[rownames(Add_st@meta.data) %in% rownames(Epi_cor)]<-"High_Epithelial_True"
  #SpatialPlot(Add_st,group.by = c("Epithelial_class"))
  Add_st@meta.data$Region<-"Remaining"
  Add_st@meta.data$Region[Add_st@meta.data$Epithelial_class == "High_Epithelial_True"] <- "Epithelial"
  SpatialPlot(Add_st,group.by = c("Region"))
  ggsave(paste0("clone_region/",new_data_name,"_region.pdf"),width = 10,height = 8)
  saveRDS(Add_st,file = paste0("clone_region/",new_data_name,"_region.RDS"))
  
  #prepare for CNV calling
  
  count_matrix<-as.data.frame(Add_st@assays$Spatial@counts)
  count_matrix$SYMBOL <- rownames(count_matrix)
  count_matrix<-merge(list,count_matrix,by="SYMBOL")
  count_matrix<-count_matrix[!duplicated(count_matrix$ENSEMBL),]
  rownames(count_matrix) <- count_matrix$ENSEMBL
  count_matrix<-count_matrix[,c(-1,-2)]
  colnames(count_matrix) <- colnames(count_matrix)
  
  #export counts matrix
  write.table(count_matrix, paste0("clone_region/",new_data_name,"_count_matrix.tsv"), sep = "\t")
  
  # get annotation matrix
  Region_annotation <-data.frame(Barcode=rownames(Add_st@meta.data),Celltype=Add_st@meta.data$Region)
  
  #export annotation matrix
  write.table(Region_annotation, paste0("clone_region/Region_",new_data_name,".tsv"), 
              sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
}

