setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

# CNV co-occurance
copy_number_CNV_module_sep <- readRDS("../../Data/Fig4_5/copy_number_CNV_module_sep.RDS")
CNV_sep_2<-copy_number_CNV_module_sep[,1:27]
CNV_sep_2_01_matrix<-(CNV_sep_2 > 2.549121 | CNV_sep_2 <1.569168)
CNV_sep_2_01_matrix[is.na(CNV_sep_2_01_matrix)]<-FALSE
CNV_sep_2_filtered<-CNV_sep_2[,colSums(CNV_sep_2_01_matrix)/466 > 0.05]
CNV_cor_matrix_2<-cor(CNV_sep_2)
cor_mat_2<-cor(CNV_cor_matrix_2,use = "pairwise.complete.obs")
colnames(CNV_cor_matrix_2)<-gsub("_adj","",colnames(CNV_cor_matrix_2))
rownames(CNV_cor_matrix_2)<-gsub("_adj","",rownames(CNV_cor_matrix_2))
colnames(cor_mat_2)<-gsub("_adj","",colnames(cor_mat_2))
rownames(cor_mat_2)<-gsub("_adj","",rownames(cor_mat_2))

pval_mat <- matrix(NA, ncol(CNV_sep_2), ncol(CNV_sep_2))
for(i in 1:ncol(CNV_sep_2)) {
  for(j in 1:ncol(CNV_sep_2)) {
    temp <- cor.test(CNV_sep_2[,i], CNV_sep_2[,j])
    pval_mat[i,j] <- temp$p.value
  }
}
sig_mat <- ifelse(pval_mat < 0.001, "***", ifelse(pval_mat < 0.01, "**", ifelse(pval_mat < 0.05, "*", "")))

rownames(sig_mat) = rownames(CNV_cor_matrix_2)
colnames(sig_mat) = colnames(CNV_cor_matrix_2)
sig_mat<-sig_mat[c("chr20p","chr20q","chr7p","chr7q","chr17p"),c("chr20p","chr20q","chr7p","chr7q","chr17p")]

library(pheatmap)
pheatmap(cor_mat_2[c("chr20p","chr20q","chr7p","chr7q","chr17p"),
                          c("chr20p","chr20q","chr7p","chr7q","chr17p")],
         color = colorRampPalette(c("#3384BB","#C3DAEB","white","white","white","white","#FED9B4","#FE7F02"))(100),
         border_color = "black",display_numbers = sig_mat,cluster_rows = F,cluster_cols = F)

