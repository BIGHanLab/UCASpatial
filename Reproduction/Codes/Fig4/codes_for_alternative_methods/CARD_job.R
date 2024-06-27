########CARD for CRC
library(CARD)
library(Seurat)
library(magrittr)
library(data.table)

##### 240329 new sc data
sc_data <- readRDS("./integrated_sc_fil_v7.rds")

# sc data
sc_count <- sc_data@assays$RNA@counts
sc_meta <- data.frame(cellID = colnames(sc_data),cellType = sc_data$UCASpatial_clus_v7,
                      sample = "sample",row.names = colnames(sc_data))
rm(sc_data)

# st data
st_files <- list.files("./","rds",ignore.case = T,full.names = T) %>% 
  grep("sc.silico.RDS|integrated_sc_fil_v7",.,invert = T,value = T)
lapply(st_files,function(x){
  st_data <- readRDS(x)
  dir_name <- paste0("./res/240329_new_scdata/",gsub(".RDS",ignore.case = T,"",basename(x)))
  dir.create(dir_name)
  
  #st data
  spatial_count <- st_data@assays$Spatial@counts
  spatial_location <- cbind(st_data@images$sliceD1@coordinates$row,
                            st_data@images$sliceD1@coordinates$col)
  spatial_location <- as.data.frame(spatial_location)
  rownames(spatial_location) <- colnames(st_data)
  colnames(spatial_location) <- c("x","y")
  
  spatial_location$x <- -spatial_location$x
  spatial_location$y <- -spatial_location$y
  
  ##object
  card_obj <- createCARDObject(sc_count = sc_count,sc_meta = sc_meta,
                               spatial_count = spatial_count,spatial_location = spatial_location,
                               ct.varname = "cellType",
                               sample.varname = "sample",
                               ct.select = unique(sc_meta$cellType),
                               minCountGene = 10,minCountSpot = 2)
  ##Deconvolution 
  card_obj <- CARD_deconvolution(card_obj)
  
  print(c(dim(card_obj@Proportion_CARD),dim(st_data)))
  
  write.table(card_obj@Proportion_CARD,
              paste0(dir_name,"/",gsub(".RDS",ignore.case = T,"",basename(x)),"_CARD_res.txt"),
              quote = F,sep = "\t",row.names = T,col.names = T)
  
})

###
x <- "./Pt37_3_region.RDS"
x <- "./Pt37_7_region.RDS"
x <- "./WT_Pt3_Tm_220812.rds"
st_data <- readRDS(x)
dir_name <- paste0("./res/240329_new_scdata/",gsub(".RDS",ignore.case = T,"",basename(x)))
dir.create(dir_name)

#st data
spatial_count <- st_data@assays$Spatial@counts
spatial_location <- cbind(st_data@images$sliceD1@coordinates$row,
                          st_data@images$sliceD1@coordinates$col)
spatial_location <- as.data.frame(spatial_location)
rownames(spatial_location) <- colnames(st_data)
colnames(spatial_location) <- c("x","y")

spatial_location$x <- -spatial_location$x
spatial_location$y <- -spatial_location$y

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
            paste0(dir_name,"/",gsub(".RDS",ignore.case = T,"",basename(x)),"_CARD_res.txt"),
            quote = F,sep = "\t",row.names = T,col.names = T)

