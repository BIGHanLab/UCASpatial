setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

TCGA_COAD_origin_tumor<-readRDS("../../Data/Fig4_5/TCGA_COAD_origin_tumor.rds")
p1<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=log2(CD8A),fill=chr20q_20p_status,color=chr20q_20p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p2<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=log2(CD8A),fill=chr20q_7p_status,color=chr20q_7p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p3<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=log2(CD8A),fill=chr20q_7q_status,color=chr20q_7q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p4<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=log2(CD8A),fill=chr20q_17p_status,color=chr20q_17p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p1+p2+p3+p4
