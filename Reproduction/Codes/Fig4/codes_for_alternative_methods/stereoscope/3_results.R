##### stereoscope results
library(Seurat)
library(data.table)

##### 240329 new sc data
###res
res_list <- list.dirs("./stereo_res/240329_new_scdata/") %>% 
  grep("sp_count",.,value = T) %>% 
  list.files(.,"tsv",full.names = T,recursive = T)

lapply(res_list,function(x){
  decon_df <- fread(x,data.table = F)
  rownames(decon_df) <- decon_df$V1
  decon_df$V1 <- NULL
  file <- x %>% strsplit(.,"/") %>% lapply(.,function(x){x[5]}) %>% unlist()
  write.table(decon_df,paste("./stereo_res/240329_new_scdata/res_df/",file,"_decon_df.txt"),
              row.names = T,col.names = T,sep = "\t",quote = F)
})


