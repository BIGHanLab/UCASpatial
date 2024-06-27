########CARD for CRC
library(CARD)
library(Seurat)
library(magrittr)
library(data.table)

########CARD for regeneration data
load("../../../Data/Fig3/MRL.st.primary.Rdata")
st_data <- ear.integrated.filter
rm(ear.integrated.filter)
sc_data <- readRDS("../../../Data/Fig3/MRL.scRNA.rds")

#st data  
spatial_count <- st_data@assays$Spatial@counts
total_coord_norm <- readRDS("../../../Data/Fig3/MRL_coord.rds")
spatial_location <- data.frame(x = total_coord_norm$row,y = total_coord_norm$col,
                               row.names = colnames(st_data))
# spatial_location <- as.data.frame(spatial_location)
# rownames(spatial_location) <- colnames(st_data)
# colnames(spatial_location) <- c("x","y")

# spatial_location$x <- -spatial_location$x
# spatial_location$y <- -spatial_location$y

#sc data
sc_count <- sc_data@assays$RNA@counts
sc_meta <- data.frame(cellID = colnames(sc_data),cellType = sc_data$Minor_class_reduction,
                      sample = "sample",row.names = colnames(sc_data))
##object
card_obj <- createCARDObject(sc_count = sc_count,sc_meta = sc_meta,
                             spatial_count = spatial_count,spatial_location = spatial_location,
                             ct.varname = "cellType",
                             sample.varname = "sample",
                             ct.select = unique(sc_meta$cellType),
                             minCountGene = 10,minCountSpot = 2)
##Deconvolution 
card_obj <- CARD_deconvolution(card_obj)

write.table(card_obj@Proportion_CARD,"./CRAD_regeneration_decon_result.txt",
            quote = F,sep = "\t",row.names = T,col.names = T)

