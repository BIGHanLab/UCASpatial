setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')
TCGA_GSVA_tumor<-readRDS("../../Data/Fig4_5/TCGA_GSVA_tumor.rds")
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggsignif)
colnames(TCGA_GSVA_tumor) <- c("ID","chr20q","pathway","GSVA_score")
theme_set(theme_cowplot())
TCGA_GSVA_tumor$chr20q[TCGA_GSVA_tumor$chr20q == FALSE]<-"WT"
TCGA_GSVA_tumor$chr20q[TCGA_GSVA_tumor$chr20q == TRUE]<-"Gain"
TCGA_GSVA_tumor$chr20q<-factor(TCGA_GSVA_tumor$chr20q,levels=c("WT","Gain"))
theme_set(theme_cowplot())
ggplot(TCGA_GSVA_tumor[!is.na(TCGA_GSVA_tumor$chr20q),],
       aes(x=pathway,y=GSVA_score,fill=chr20q,color=chr20q))+geom_boxplot(size=0.6)+scale_fill_manual(values = c("#97C3D8","#F6B373"))+
  scale_color_manual(values = c("#568DB0","#E78844"))+labs(x="Pathway",y="GSVA Score",title = "Pathway Activity (TCGA-COAD)")+
  theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(angle = 45,hjust = 1))+stat_compare_means(label="p.signif")


