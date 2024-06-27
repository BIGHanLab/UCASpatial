library(Seurat)
library(cowplot)
# p_unload(Seurat)
# p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')


CRC4_TE <- readRDS("../../Data/Fig4_5/CRC4_st_with_TE.RDS")
DefaultAssay(CRC4_TE)<-"TE"
CRC4_TE@meta.data<-CRC4_TE@meta.data[colnames(CRC4_TE),]
table(colnames(CRC4_TE) == rownames(CRC4_TE@meta.data))
CRC4_TE@meta.data$HERVH<-colMeans(CRC4_TE@assays$TE@data[c("LTR7Y","LTR7","LTR7B","HERVH-int"),])
DefaultAssay(CRC4_TE) <- "Spatial"
p1<-SpatialPlot(CRC4_TE,features = "HERVH",image.alpha = 0,stroke = NA,pt.size.factor = 1.5)+
  scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "HERV-H",label.position = "right",label.hjust = 3,
                                direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),
                                title.theme = element_text(size=16)))

Idents(CRC4_TE)<-"Region"
HERVH_expression<-CRC4_TE@meta.data[,c("Region","HERVH")]
p2<-ggplot(HERVH_expression[HERVH_expression$Region %in% c("Clone1","Clone2","Clone3"),],aes(x=Region,y=HERVH,fill=Region,color=Region))+
  geom_violin()+scale_fill_manual(values = c("#A2C581","#9BBFD0","#D77566"))+
  scale_color_manual(values = c("#378B39","#1C6797","#B62D28"))+
  labs(x="Clones",y="Expression",title = "HERV-H")+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_blank())+guides(color=F)+stat_compare_means(label="p.signif",ref.group = "Clone3")
p1+p2
