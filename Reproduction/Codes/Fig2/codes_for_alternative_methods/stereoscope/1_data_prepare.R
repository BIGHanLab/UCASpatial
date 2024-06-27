##### prepare data
library(Seurat)
library(data.table)

###### simulated data
setwd("/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/")
##sp counts
files <- list.files("../../Data/Fig2/","simu.*rds$",full.names = T)
lapply(files,function(x){
  st <- readRDS(x)
  dir_name <- paste0("./sp_count/",gsub(".rds","",basename(x)))
  dir.create(dir_name)
  for(i in 1:length(st)){
    st_count_data <- st[[i]]$topic_profiles %>% t() %>% as.data.frame()
    st_count_data <- cbind(rownames(st_count_data),st_count_data)
    colnames(st_count_data)[1] <- "cell"
    write.table(st_count_data,
                paste0(dir_name,"/",i,"_sp_count.tsv"),
                sep = "\t",col.names = T,row.names = F,quote = F)
  }
})

#sc
sc_data <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")
sc_count_data <- sc_data@assays$RNA@counts %>% t() %>% as.data.frame()
sc_count_data <- cbind(rownames(sc_count_data),sc_count_data)
colnames(sc_count_data)[1] <- "cell"

write.table(sc_count_data,"./simu_sc_count.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)

meta_data <- data.frame(cell = rownames(sc_count_data),bio_celltype = sc_data$low_res_ident)
write.table(meta_data,"./meta/low_res_ident_meta.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)
meta_data <- data.frame(cell = rownames(sc_count_data),bio_celltype = sc_data$med_res_ident)
write.table(meta_data,"./meta/med_res_ident_meta.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)
meta_data <- data.frame(cell = rownames(sc_count_data),bio_celltype = sc_data$high_res_ident)
write.table(meta_data,"./meta/high_res_ident_meta.tsv",
            sep = "\t",col.names = T,row.names = F,quote = F)