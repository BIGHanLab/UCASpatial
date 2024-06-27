setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

CRC4_st <- readRDS("../../Data/Fig4_5/CRC4_st.rds")
library(pacman)
library(Seurat)
library(patchwork)
library(ggrepel)
library(clusterProfiler)

# p_unload(Seurat)
# p_load(Seurat)
DefaultAssay(CRC4_st) <- "Spatial"
CRC4_st_tumor<-CRC4_st[,CRC4_st$Region %in% c("Clone1","Clone2","Clone3")]

# saveRDS(CRC4_st_tumor,"../../Data/Fig4_5/CRC4_st_tumor.rds")
CRC4_st_tumor <- readRDS("../../Data/Fig4_5/CRC4_st_tumor.rds")
Idents(CRC4_st_tumor)<-CRC4_st_tumor@meta.data$Region
diff_ST_subclone<-FindAllMarkers(CRC4_st_tumor,assay = "Spatial",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,
                                 min.cells.feature = 0)
diff_ST_subclone$p_val_adj <- -log10(diff_ST_subclone$p_val_adj)
diff_ST_subclone$change<-"Remain"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC > 0.1 & diff_ST_subclone$p_val_adj >1.30103]<-"Up"
diff_ST_subclone$change[diff_ST_subclone$avg_logFC < (-0.1) & diff_ST_subclone$p_val_adj >1.30103]<-"Down"
diff_ST_subclone$change <- factor(diff_ST_subclone$change,levels = c("Up","Down","Remain"))
label<-c("IRF7","IRF1","ISG15","IFI27","IFI16","CXCL9","CXCL10","CXCL11","B2M","PSMB9","PSMB8","TAP1","TAP2",
         "HLA-A","HLA-B","STAT1","OLFM4","FCGBP","MUC2","RPL37","RPL31","RPS5")
label_data_down<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Down",]
label_data_up<-diff_ST_subclone[diff_ST_subclone$cluster == "Clone3" & diff_ST_subclone$gene %in% label & diff_ST_subclone$change == "Up",]
theme_set(theme_cowplot())

ggplot(diff_ST_subclone[diff_ST_subclone$cluster == "Clone3",],aes(x=avg_logFC,y=p_val_adj,color=change))+geom_point()+
  scale_color_manual(values = c("#FE7F00","#1F78B4","lightgrey"))+coord_cartesian(xlim=c(-2,2))+
  labs(x="log2(Clone3 vs Other Clones)",y="-log10(adj.P)",title = "Differential Expression Genes in Clone3")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.text = element_text(hjust = 0.5))+
  geom_vline(xintercept = -0.1,lty=2)+geom_vline(xintercept = 0.1,lty=2)+
  geom_hline(yintercept = 2,lty=2)+geom_text_repel(data=label_data_down,aes(color=change, label=gene), 
                                                   point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 12, 
                                                   nudge_y = 0.1,ylim = c(5,300),xlim = c(-2,-1))+
  geom_text_repel(data=label_data_up,aes(color=change, label=gene), 
                  point.padding = unit(0.2, "lines"), box.padding = unit(0.4, "lines"),max.overlaps = 10, 
                  nudge_y = 0.1,ylim = c(5,300),xlim=c(1,2))+guides(color=F)

Sub3_up<-diff_ST_subclone[(diff_ST_subclone$change == "Up") & (diff_ST_subclone$cluster == "Clone3"),]
Sub3_down<-diff_ST_subclone[(diff_ST_subclone$change == "Down") & (diff_ST_subclone$cluster == "Clone3"),]
# Sub3_up<-readRDS("../../Data/Fig4_5/Sub3_up_gene.RDS")
#Sub3_down<-readRDS("../../Data/Fig4_5/Sub3_down_gene.RDS")
up_GO<-enrichGO(Sub3_up,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
down_GO<-enrichGO(Sub3_down,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
up_GO@result$pvalue=-log10(up_GO@result$pvalue)
down_GO@result$pvalue=-log10(down_GO@result$pvalue)
selected_GO_up<-up_GO@result[c(9,8,6,1),]
selected_GO_down<-down_GO@result[c(1133,169,105,2),]
selected_GO_down$Description<-factor(selected_GO_down$Description,levels=selected_GO_down$Description)
selected_GO_up$Description<-factor(selected_GO_up$Description,levels=selected_GO_up$Description)
theme_set(theme_cowplot())
p1<-ggplot(selected_GO_down,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Down")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,20)+scale_fill_manual(values = c(rep("#A5CDE2",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p2<-ggplot(selected_GO_up,aes(x=Description,y=pvalue,fill=Description))+
  geom_bar(stat = "identity",width=0.8)+coord_flip()+guides(fill=F)+
  labs(x="",y=expression(-log[10]~italic(P)),title = "Clone3 Up")+theme(axis.text.x =element_text(size=16),aspect.ratio = 1)+
  ylim(0,70)+scale_fill_manual(values = c(rep("#FDBE6D",6)))+
  theme(plot.title = element_text(hjust = 0.5))
p1+p2+plot_layout(ncol=1)