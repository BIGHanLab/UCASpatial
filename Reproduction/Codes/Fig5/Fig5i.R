#diff TE analysis
library(Seurat)
library(cowplot)
# p_unload(Seurat)
# p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

CRC4_TE <- readRDS("../../Data/Fig4_5/CRC4_st_with_TE.RDS")
DefaultAssay(CRC4_TE)<-"TE"
CRC4_TE@meta.data<-CRC4_TE@meta.data[colnames(CRC4_TE),]
table(colnames(CRC4_TE) == rownames(CRC4_TE@meta.data))

CRC4_TE_Tumor<-CRC4_TE[,CRC4_TE$Region %in% c("Clone1","Clone2","Clone3")]
Idents(CRC4_TE_Tumor)<-CRC4_TE_Tumor$Region
diff_TE<-FindAllMarkers(CRC4_TE_Tumor,assay = "TE",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,only.pos = F)
repeats_info<-read.csv("../../Data/Fig4_5/repeats_family_info.csv",header = F)
CRC_specific_TE<-read.csv("../../Data/Fig4_5/GCP_COAD_tumor_vs_normal_Diff_TE.csv")
CRC_specific_TE$change_tumor<-"Remain"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC>0] <- "Up"
CRC_specific_TE$change_tumor[CRC_specific_TE$logFC<0] <- "Down"
CRC_specific_TE<-CRC_specific_TE$X[CRC_specific_TE$change_tumor == "Up"]
diff_CRC_TE<-diff_TE[diff_TE$gene %in% CRC_specific_TE,]
Clone3_TE<-diff_TE[diff_TE$cluster == "Clone3",]
Clone3_TE$p_val_adj <- -log10(Clone3_TE$p_val_adj)
Clone3_TE$change<-"Remain"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC > 0.25] = "Up"
Clone3_TE$change[Clone3_TE$p_val_adj > 2 & Clone3_TE$avg_logFC < (-0.25)] = "Down"
Clone3_TE$group<-"Remain"
Clone3_TE$group[Clone3_TE$gene %in% c("LTR7","LTR7Y","LTR7B","HERVH-int")]<-"HERVH"

theme_set(theme_cowplot())
library(ggrepel)
for_label <-  Clone3_TE[(Clone3_TE$change != "Remain") & (Clone3_TE$gene %in% CRC_specific_TE),]
for_label_2 <-  Clone3_TE[(Clone3_TE$group == "HERVH") & (Clone3_TE$gene %in% CRC_specific_TE),]
theme_set(theme_cowplot())

ggplot(Clone3_TE[Clone3_TE$gene %in% CRC_specific_TE,],aes(x=avg_logFC,y=p_val_adj,fill=change,color=group))+geom_point(shape=21,stroke=0.5,size=1.5)+
  scale_color_manual(values = c("black",alpha("white",0),rep(alpha("white",0),9)))+ 
  scale_fill_manual(values = c("#1E72A9","lightgrey","#D1332D"))+
  labs(x="log2(Clone3 vs Other Subclones)",y="-log2(adj.P)",title = "Tumor Specific TE in CRC")+
  theme(plot.title = element_text(size=12,hjust=0.5))+geom_vline(xintercept = 0.25,lty=2)+geom_vline(xintercept = -0.25,lty=2)+
  geom_hline(yintercept = 2,lty=2)+coord_cartesian(xlim = c(-1.5,1.5))+
  geom_text_repel(data = for_label_2[for_label_2$avg_logFC<0,],aes(label = gene),alpha = 1,color="#1F78B4",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(-0.3,-1.2),size=4)+
  geom_text_repel(data = for_label[for_label$avg_logFC>0,],aes(label = gene),alpha = 1,color="#D1332D",
                  max.overlaps = 15, segment.size = 0.5,box.padding = 1,force = 0.5,nudge_x = 0.15,xlim = c(0.3,1.2),size=4)