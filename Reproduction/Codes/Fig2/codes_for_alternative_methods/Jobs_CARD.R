########CARD for CRC
library(CARD)
library(Seurat)
library(magrittr)
library(data.table)

##### 20231113
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/')
sc_data <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")


#sc data
sc_count <- sc_data@assays$RNA@counts
sc_meta <- data.frame(cellID = colnames(sc_data),cellType = sc_data$high_res_ident,
                      sample = "sample",row.names = colnames(sc_data))


#st data
st_files <- list.files("../../Data/Fig2","simu",full.names = T)
lapply(st_files,function(x){
  st <- readRDS(x)
  dir_name <- paste0("./res/",gsub(".rds","",basename(x)))
  dir.create(dir_name)
  for(i in 1:length(st)){
    spatial_count <- st[[i]]$topic_profiles
    spatial_location <- data.frame(x = rep(1:10,5),y = rep(1:5,10))
    # spatial_location <- as.data.frame(spatial_location)
    rownames(spatial_location) <- colnames(spatial_count)
    # colnames(spatial_location) <- c("x","y")
    ##object
    card_obj <- createCARDObject(sc_count = sc_count,sc_meta = sc_meta,
                                 spatial_count = spatial_count,spatial_location = spatial_location,
                                 ct.varname = "cellType",
                                 sample.varname = "sample",
                                 ct.select = unique(sc_meta$cellType),
                                 minCountGene = 10,minCountSpot = 2)
    ##Deconvolution 
    card_obj <- CARD_deconvolution(card_obj)
    
    write.table(card_obj@Proportion_CARD,
                paste0(dir_name,"/",i,"_CARD_res.txt"),
                quote = F,sep = "\t",row.names = T,col.names = T)
    
  }
})