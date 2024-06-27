tumor_prop_avg<-readRDS("../../Data/Fig4_5/tumor_prop_avg.RDS")

tumor_prop_avg$community_3<-rowSums(tumor_prop_avg[,c(4,5,11,13,32)])
tumor_prop_avg$community_3 = tumor_prop_avg$community_3*100
library(ggpubr)

ggplot(tumor_prop_avg,aes(x=subTME_category,y=community_3,color=subTME_category,fill=subTME_category))+geom_boxplot(size=1)+
  scale_fill_manual(values = c("#86B1C3","#E5A173"))+scale_color_manual(values = c("#2B6D98","#DE7C3A"))+
  stat_compare_means(label="p.signif",hjust = 0.5)+
  labs(x="ClonalTIME",y="Proportion (%)",title = "Community3")+guides(fill=F,color=F)+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits=c(0, 9),)
