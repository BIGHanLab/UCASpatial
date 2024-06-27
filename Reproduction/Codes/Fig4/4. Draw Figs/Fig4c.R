library(ggthemes)
theme_base()
theme_set(theme_cowplot())
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig4/')

SPOTlight<-readRDS("../../Data/Fig4/processed_data/SPOTlight_results.rds")
colnames(SPOTlight[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[1]]@meta.data))))))))))))
colnames(SPOTlight[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[2]]@meta.data))))))))))))
colnames(SPOTlight[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[3]]@meta.data))))))))))))
colnames(SPOTlight[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[4]]@meta.data))))))))))))
colnames(SPOTlight[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[5]]@meta.data))))))))))))
colnames(SPOTlight[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[6]]@meta.data))))))))))))
colnames(SPOTlight[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[7]]@meta.data))))))))))))
colnames(SPOTlight[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(SPOTlight[[8]]@meta.data))))))))))))

UCASpatial<-readRDS("../../Data/Fig4/processed_data/UCASpatial_results.rds")
colnames(UCASpatial[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[1]]@meta.data))))))))))))
colnames(UCASpatial[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[2]]@meta.data))))))))))))
colnames(UCASpatial[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[3]]@meta.data))))))))))))
colnames(UCASpatial[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[4]]@meta.data))))))))))))
colnames(UCASpatial[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[5]]@meta.data))))))))))))
colnames(UCASpatial[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[6]]@meta.data))))))))))))
colnames(UCASpatial[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[7]]@meta.data))))))))))))
colnames(UCASpatial[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(UCASpatial[[8]]@meta.data))))))))))))

C2L<-readRDS("../../Data/Fig4/processed_data/C2L_results.rds")
colnames(C2L[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[1]]@meta.data))))))))))))
colnames(C2L[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[2]]@meta.data))))))))))))
colnames(C2L[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[3]]@meta.data))))))))))))
colnames(C2L[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[4]]@meta.data))))))))))))
colnames(C2L[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[5]]@meta.data))))))))))))
colnames(C2L[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[6]]@meta.data))))))))))))
colnames(C2L[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[7]]@meta.data))))))))))))
colnames(C2L[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(C2L[[8]]@meta.data))))))))))))

CARD<-readRDS("../../Data/Fig4/processed_data/CARD_results.rds")
colnames(CARD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[1]]@meta.data))))))))))))
colnames(CARD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[2]]@meta.data))))))))))))
colnames(CARD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[3]]@meta.data))))))))))))
colnames(CARD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[4]]@meta.data))))))))))))
colnames(CARD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[5]]@meta.data))))))))))))
colnames(CARD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[6]]@meta.data))))))))))))
colnames(CARD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[7]]@meta.data))))))))))))
colnames(CARD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(CARD[[8]]@meta.data))))))))))))

RCTD<-readRDS("../../Data/Fig4/processed_data/RCTD_results.rds")
colnames(RCTD[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[1]]@meta.data))))))))))))
colnames(RCTD[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[2]]@meta.data))))))))))))
colnames(RCTD[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[3]]@meta.data))))))))))))
colnames(RCTD[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[4]]@meta.data))))))))))))
colnames(RCTD[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[5]]@meta.data))))))))))))
colnames(RCTD[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[6]]@meta.data))))))))))))
colnames(RCTD[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[7]]@meta.data))))))))))))
colnames(RCTD[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(RCTD[[8]]@meta.data))))))))))))


stereoscope<-readRDS("../../Data/Fig4/processed_data/stereoscope_results.rds")
colnames(stereoscope[[1]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[1]]@meta.data))))))))))))
colnames(stereoscope[[2]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[2]]@meta.data))))))))))))
colnames(stereoscope[[3]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[3]]@meta.data))))))))))))
colnames(stereoscope[[4]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[4]]@meta.data))))))))))))
colnames(stereoscope[[5]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[5]]@meta.data))))))))))))
colnames(stereoscope[[6]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[6]]@meta.data))))))))))))
colnames(stereoscope[[7]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[7]]@meta.data))))))))))))
colnames(stereoscope[[8]]@meta.data)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(stereoscope[[8]]@meta.data))))))))))))



p1<-SpatialPlot(UCASpatial[[3]],features = "GZMK",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FEFAF8',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#FEFAF8")+
  labs(fill ='GZMK')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p2<-SpatialPlot(UCASpatial[[3]],features = "CXCL13",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#F9FCFE',brewer.pal(n = 9, name = "Blues")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#F9FCFE")+
  labs(fill ='CXCL13')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(UCASpatial[[3]],features = "XCL1",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#F9F8FB',brewer.pal(n = 9, name = "Purples")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "#F9F8FB")+
  labs(fill ='XCL1')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(UCASpatial[[3]],features = "FAP",image.alpha = 0,stroke = 0,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FBF2DE',brewer.pal(n = 9, name = "#E2BC3C")))(n=100),
                       na.value = "#FBF2DE")+
  labs(fill ='FAP')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4
p1+p2+p3+p4


p_list <- list(p1,p2,p3,p4)
marker_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


DotPlot(sc_new,features = c("GZMK","CXCL13","KLRD1"))



p1<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
            min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(UCASpatial[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(UCASpatial[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('#FBF2DE',brewer.pal(n = 9, name = "#E2BC3C")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p_list <- list(p1,p2,p3,p4)
UCASpatial_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')




p1<-SpatialPlot(CARD[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p1
p2<-SpatialPlot(CARD[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(CARD[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(CARD[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))


p_list <- list(p1,p2,p3,p4)
CARD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')



p1<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(SPOTlight[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(SPOTlight[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))




p_list <- list(p1,p2,p3,p4)
SPOTlight_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')



p1<-SpatialPlot(C2L[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(C2L[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(C2L[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(C2L[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p1+p2+p3


p_list <- list(p1,p2,p3,p4)
C2L_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialPlot(stereoscope[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(stereoscope[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(stereoscope[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(stereoscope[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p_list <- list(p1,p2,p3,p4)
stereoscope_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


p1<-SpatialPlot(RCTD[[3]],features = "CD8+ Tem",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tem`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))

p2<-SpatialPlot(RCTD[[3]],features = "CD8+ Tex",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Tex`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p3<-SpatialPlot(RCTD[[3]],features = "CD8+ Trm",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`CD8+ Trm`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p4<-SpatialPlot(RCTD[[3]],features = "FAP+ Myofib",image.alpha = 0,stroke = 0.01,pt.size.factor = 2,
                min.cutoff = 0.03)+
  scale_fill_gradient2(low = "white", high =  "gold2",
                       limits = c(0.03, quantile(UCASpatial[[3]]@meta.data$`FAP+ Myofib`, 1)),
                       #limits = c(0.03,0.5),
                       na.value = "white",
                       labels = scales::label_percent())+
  labs(fill ='Proportion')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))



p_list <- list(p1,p2,p3,p4)
RCTD_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')


plot_grid(marker_p,UCASpatial_p,CARD_p,C2L_p,RCTD_p,SPOTlight_p,stereoscope_p,
          align = 'hv',nrow = 1)


p1<-SpatialPlot(UCASpatial[[3]],features = "CD3G",image.alpha = 0,stroke = 0.01,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "white")+
  labs(fill ='CD3G')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p2<-SpatialPlot(UCASpatial[[3]],features = "CD8A",image.alpha = 0,stroke = 0.01,pt.size.factor = 2)+
  scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                       #limits = c(0.03,0.5),
                       na.value = "white")+
  labs(fill ='CD8A')+theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))
p_list <- list(p1,p2)
plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
