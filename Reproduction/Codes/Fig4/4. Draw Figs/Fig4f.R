tumor_prop_avg<-readRDS("../../Data/Fig4_5/tumor_prop_avg.RDS")

tumor_prop_avg[,c(8,15,16,25,27)]

cell_composition_clone<-tumor_prop_avg[,c(-8,-15,-16,-25,-27)]


PCC_matrix <- matrix(NA, 29, 29)
for(i in 1:29) {
  for(j in 1:29) {
    temp <- cor.test(cell_composition_clone[,i], cell_composition_clone[,j],method = "pearson")
    PCC_matrix[i,j] <- temp$estimate
  }
}

rownames(PCC_matrix) = colnames(cell_composition_clone[,1:29])
colnames(PCC_matrix) = colnames(cell_composition_clone[,1:29])

pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                   cellwidth = 10,cellheight = 10,fontsize = 8,breaks = seq(-1, 1, length.out = 100),
                   border_color = "black",cutree_cols = 5,cutree_rows = 5,
                   clustering_method = "complete")

