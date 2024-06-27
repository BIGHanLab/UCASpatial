#signature calculation
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

TCGA_COAD_origin_tumor<-readRDS("../../Data/Fig4_5/TCGA_COAD_origin_tumor.rds")
sc_marker<-readRDS("../../Data/Fig4_5/Seurat_cluster_markers.rds")
top_genes <- sc_marker %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

for (cluster in clusters) {
  TCGA_COAD_origin_tumor[[paste0(cluster, "_signature")]] <- rowMeans(
    as.data.frame(scale(TCGA_COAD_origin_tumor[, colnames(TCGA_COAD_origin_tumor) %in% top_genes$gene[top_genes$cluster == cluster]])) %>%dplyr::select(where(~ !any(is.na(.))))
  )
}


# split gain and loss CNV
adj_cols <- grep("_adj", colnames(TCGA_COAD_origin_tumor))
adj_CNV<-TCGA_COAD_origin_tumor[,adj_cols]
mat_gain <- apply(adj_CNV, 2, function(x) ifelse(x > 2, x, NA))
mat_loss <- apply(adj_CNV, 2, function(x) ifelse(x < 2, x, NA))

mat_new <- cbind(mat_gain, mat_loss)
colnames(mat_new) <- c(paste0(colnames(adj_CNV), "_gain"), paste0(colnames(adj_CNV), "_loss"))

TCGA_COAD_origin_tumor<-cbind(TCGA_COAD_origin_tumor,mat_new)


#TCGA matrix
colnames(TCGA_COAD_origin_tumor)[adj_cols]
adj_cols <- grep("_adj_", colnames(TCGA_COAD_origin_tumor))
signature_cols <- grep("signature", colnames(TCGA_COAD_origin_tumor))

cor_matrix <- matrix(0, nrow = length(adj_cols), ncol = length(signature_cols))
rownames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[adj_cols]
colnames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[signature_cols]

for (i in 1:length(adj_cols)) {
  for (j in 1:length(signature_cols)) {
    cor_matrix[i, j] <- cor(TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), adj_cols[i]], TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), signature_cols[j]], method = "pearson")
  }
}

cor_matrix<-cor_matrix[,4:ncol(cor_matrix)]
colnames(cor_matrix)<-gsub("_signature","",colnames(cor_matrix))
rownames(cor_matrix)<-gsub("_adj","",rownames(cor_matrix))
bk = unique(c(seq(-0.4,0.4, length=100)))
cor_matrix[abs(cor_matrix) < 0.2] = 0
cor_matrix[grep("loss",rownames(cor_matrix)),]<- (-cor_matrix[grep("loss",rownames(cor_matrix)),])
p1<-pheatmap::pheatmap(cor_matrix[c("chr20q_gain","chr7p_gain","chr20p_gain","chr7q_gain","chr17p_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),1:29],
                       color = colorRampPalette(c("#3384BB","#A6CEE3","#E4F0F6","white","white","white","white","white","white","#FEECD4","#FDBF6F","#FE7F02"))(100),cluster_rows = F,
                       cluster_cols = F,breaks = bk,border_color = "#262626",angle_col = 90,cellwidth = 15,cellheight = 15,
                       gaps_col = c(10,20,25,27,29),gaps_row = c(4,5,7))


library(reshape2)
CNV_Frequency<-readRDS("../../Data/Fig4_5/CNV_Frequency.RDS")
CNV_Frequency<-data.frame(loss_frequency=colSums(CNV_module_filtered_1 < (-0.35))/nrow(CNV_module_filtered_1),gain_frequency=colSums(CNV_module_filtered_1 > (0.35))/nrow(CNV_module_filtered_1))
CRC4_CNV_freq<-CNV_Frequency[c("chr20q","chr7p","chr7q","chr20p","chr17p","chr17q","chr8p","chr8q","chr1p_2","chr14q"),]
CRC4_CNV_freq$region<-rownames(CRC4_CNV_freq)
CRC4_CNV_freq<-melt(CRC4_CNV_freq)
CRC4_CNV_freq$region<-factor(CRC4_CNV_freq$region,levels=rev(c("chr20q","chr20p","chr7p","chr7q","chr17p","chr17q","chr8q","chr8p","chr14q","chr1p_2")))
theme_set(theme_cowplot())

CRC4_CNV_freq_2<-CRC4_CNV_freq[c(11,12,14,13,5,18,17,9),]
CRC4_CNV_freq$region<-factor(CRC4_CNV_freq$region,levels=rev(c("chr20q","chr7p","chr20p","chr7q","chr17p","chr17q","chr8q","chr8p","chr14q","chr1p_2")))
ggplot(CRC4_CNV_freq_2,aes(x=value,y=region,fill=variable))+geom_bar(stat = "identity",width = 0.7)+scale_fill_manual(values = c("#97C3D8","#F6B373"))
p2<-TCGA_COAD_origin_tumor$chr20p_status
