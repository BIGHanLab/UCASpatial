library(pacman)
library(Seurat)
p_unload(Seurat)
p_load(Seurat)

CRC4_st<-readRDS("../../Data/Fig4_5/clone_region/CRC4_st_clone.RDS")
CRC4_st_tumor<-CRC4_st[,CRC4_st$Node != "Remaining"]
p1<-SpatialPlot(CRC4_st_tumor,group.by ="Node",cols = c("#A7CC83","#9BBFD0","#ED7B6C"),stroke = F,pt.size=1.2)

tumor_prop_avg<-readRDS("../../Data/Fig4_5/tumor_prop_avg.RDS")
tumor_prop_avg$community_3<-rowSums(tumor_prop_avg[,c(4,5,11,13,32)])
tumor_prop_avg$community_3 = tumor_prop_avg$community_3*100
p2<-ggplot(tumor_prop_avg[tumor_prop_avg$slide=="CRC4",],aes(x=clone,y=community_3,fill=clone))+geom_bar(stat="identity")+
  scale_fill_manual(values = c("#96BB84","#8CB3C4","#E27768"))+labs(x="Clones",y="Proportions (%)",title = "Community3 in CRC4")+
  theme(plot.title = element_text(hjust = 0.5))+guides(fill=F)

p3<-SpatialPlot(CRC4_st,features = c("CD8.Tex"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.2,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tex",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p4<-SpatialPlot(CRC4_st,features = c("CD8.Tem"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.2,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "CD8+ Tem",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
p5<-SpatialPlot(CRC4_st,features = c("IGF1..Mac"),image.alpha = 0,stroke = NA,ncol = 1,max.cutoff = 0.25,pt.size.factor = 1.5,min.cutoff = 0.03)+
  scale_fill_gradientn(colors =c("#8DB8D1","white","#FAC08A","#F6B373","#F0966B","#E06058","#D22F45"))+
  theme(legend.position = "right",plot.title = element_text(size=24,hjust=0.5,vjust=2))+
  guides(fill = guide_colourbar(title = "C1QC+ Mac",label.position = "right",label.hjust = 3,direction = "vertical",barheight = 18,barwidth = 1,ticks = T,ticks.colour = "black",ticks.linewidth = 0.3,
                                frame.colour = "black",frame.linewidth = 0.3,label.theme = element_text(size=16,hjust = 2),title.theme = element_text(size=16)))
library(patchwork)
p1+p2+p3+p4+p5+plot_layout(ncol = 5)
