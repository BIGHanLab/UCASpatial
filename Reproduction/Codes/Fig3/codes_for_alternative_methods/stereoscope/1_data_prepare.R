##### prepare data
library(Seurat)
library(data.table)

##### regeneration data
setwd("/data/zhangff/ST/data/regeneration/")
##sp counts
load("./ear.integrated.filter.2.spot_20231109.Rdata")
st <- ear.integrated.filter
st_count_data <- st@assays$Spatial@counts %>% t() %>% as.data.frame()
st_count_data <- cbind(rownames(st_count_data),st_count_data)
colnames(st_count_data)[1] <- "cell"
write.table(st_count_data,"regeneration_sp_count.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)

##sc data
sc <- readRDS("./MRL_WT_D7_all_cluster_singlet_20230406.RDS")
sc_count_data <- sc@assays$RNA@counts %>% t() %>% as.data.frame()
sc_count_data <- cbind(rownames(sc_count_data),sc_count_data)
colnames(sc_count_data)[1] <- "cell"

write.table(sc_count_data,"./regeneration_sc_count.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)

meta_data <- data.frame(cell = rownames(sc_count_data),bio_celltype = sc$Minor_class_reduction)

write.table(meta_data,"./regeneration_meta.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)