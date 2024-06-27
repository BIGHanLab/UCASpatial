##### stereoscope results
library(Seurat)
library(data.table)

###### simulated data
setwd("/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/")
###res
res_list <- list.dirs("../../Data/Fig2/stereoscope_result/stereoscope_res/",recursive = T) %>% 
  grep("sp_count",.,value = T) %>% 
  list.files(.,"tsv",full.names = T,recursive = T)

lapply(res_list,function(x){
  decon_df <- fread(x,data.table = F)
  rownames(decon_df) <- decon_df$V1
  decon_df$V1 <- NULL
  file <- x %>% strsplit(.,"/") %>% lapply(.,function(x){x[4:6] %>% paste0(.,collapse = "_")}) %>% unlist()
  write.table(decon_df,paste("./stereoscope_res/res_df/",file,"_decon_df.txt"),
              row.names = T,col.names = T,sep = "\t",quote = F)
})
