##### stereoscope results
library(Seurat)
library(data.table)

##### regeneration data
setwd("/data/zhangff/ST/data/regeneration/")
## res
stereo_decon_df <- fread("./stereo_res/regeneration/regeneration_sp_count/W.2024-01-09043158.852231.tsv",data.table = F)
rownames(stereo_decon_df) <- stereo_decon_df$V1
stereo_decon_df$V1 <- NULL

write.table(stereo_decon_df,"./stereo_decon_df.txt",row.names = T,col.names = T,sep = "\t",quote = F)
