library(Seurat)
setwd('/data/xy/Spatial_transcriptome/UCASpatial/Reproduction/Codes/Fig3/codes_for_alternative_methods/')
source('../../source/UCASpatial_final_v1.R')

## scRNA
sc.data <- readRDS('../../../Data/Fig3/MRL.scRNA.rds')
load('../../../Data/Fig3/MRL.st.primary.Rdata')



## spatial
UCASpatial_regene <- UCASpatial_deconv(sc_ref = sc.data,downsample_n = 1000,meta.resolution = 200,meta.filter = 0.95,
                               cos.filter.threshold = 0.1,ent.filter.threshold = 0.5,weight.filter.threshold = 0.2,
                               cos.filter = T,output_path = '/data/zhangyw/project/regeneration/spatial_20210214/figure/20230330_ewide/v5/res200_filter0.95',
                               st_vis = ear.integrated.filter,clust_vr = "Minor_class_reduction")

save(UCASpatial_regene,file = './UCASpatial_regene.Rdata')  
saveRDS(UCASpatial_regene,'UCASpatial_regene.rds')