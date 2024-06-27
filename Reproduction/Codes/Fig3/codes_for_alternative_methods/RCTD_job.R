library(RCTD)
library(Seurat)
library(ggpubr) 
library(gridExtra) 
library(reshape2) 
library(spacexr)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig3/codes_for_alternative_methods/')

total_coord_norm <- readRDS("../../../Data/Fig3/MRL_coord.rds")
sc.Reg <- readRDS('../../../Data/Fig3/MRL.scRNA.rds')
load('../../../Data/Fig3/MRL.st.primary.Rdata')
st.Reg <- ear.integrated.filter

scdata <- sc.Reg
clust_vr <- 'Minor_class_reduction'
Seurat::Idents(scdata) <- scdata$Minor_class_reduction
cell_types <- as.data.frame(scdata@meta.data[[clust_vr]]) %>%
  mutate(Minor_class_reduction = scdata@meta.data[[clust_vr]])
rownames(cell_types) <- colnames(scdata)

levels(Idents(scdata))
table(scdata$Minor_class_reduction)
selected_ct <- names(table(cell_types$Minor_class_reduction))[table(cell_types$Minor_class_reduction) > 25]
# levels(cell_types$Minor_class_reduction)[levels(cell_types$Minor_class_reduction)=="ILC1/NK"] <- "ILC1_NK"
selected_cell <- WhichCells(scdata,idents = selected_ct)

cell_types <- cell_types[rownames(cell_types)%in%selected_cell,]

table(cell_types$Minor_class_reduction)

cell_types_Minor_class_reduction = as.factor(cell_types$Minor_class_reduction)
names(cell_types_Minor_class_reduction) <- rownames(cell_types)
cell_types_Minor_class_reduction <- factor(cell_types_Minor_class_reduction,levels = levels(droplevels(cell_types_Minor_class_reduction)))

reference <-
  Reference(scdata@assays$RNA@counts[,colnames(scdata) %in% selected_cell],
            cell_types = cell_types_Minor_class_reduction)


x <- ear.integrated.filter
coords <- as.data.frame(cbind(total_coord_norm$row,total_coord_norm$col))
colnames(coords) <- c('row', 'col')
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
saveRDS(myRCTD_full,'result_RCTD.rds')


