##### prepare data
library(Seurat)
library(data.table)

##### 240329 new sc data
setwd("/data/zhangff/ST/data/20231119_trueData/")
##sp counts
sc_data <- readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/integrated_sc_fil_v7.rds")
sc_count_data <- sc_data@assays$RNA@counts %>% t() %>% as.data.frame()
sc_count_data <- cbind(rownames(sc_count_data),sc_count_data)
colnames(sc_count_data)[1] <- "cell"

write.table(sc_count_data,"./240329_sc_count.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)

meta_data <- data.frame(cell = rownames(sc_count_data),bio_celltype = sc_data$UCASpatial_clus_v7)

write.table(meta_data,"./240329_meta.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)


