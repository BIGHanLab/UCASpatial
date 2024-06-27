setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

TCGA_COAD_total_with_20q<-readRDS("../../Data/Fig4_5/TCGA_COAD_total_with_20q.RDS")
TCGA_COAD_total_with_20q$chr20q_status.y<-factor(TCGA_COAD_total_with_20q$chr20q_status.y,levels=c("Normal","WT","Gain"))
ggplot(TCGA_COAD_total_with_20q[!is.na(TCGA_COAD_total_with_20q$chr20q_status.y),],aes(x=chr20q_status.y,y=log2(HERVH_expression),color=chr20q_status.y))+
  geom_boxplot(size=1)+scale_color_manual(values = c("#BFBDBD","#1C6797","#DD7E2D"))+labs(x="",y="log2(RPM)",title="HERV-H")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust=0.5))+geom_signif(comparisons = list(c("Normal", "WT"),c("Normal", "Gain"),c("WT", "Gain")),y_position=c(8, 9, 8))

