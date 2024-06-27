setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

#CNV profile
CNV_clonal_region<-readRDS("../../Data/Fig4_5/CNV_clonal_region.RDS")
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
clone_order<-gsub(" ","_",clone_order)
pheatmap::pheatmap(CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr13q","chr5p","chr9q","chr8q","chr12p","chr22q","chr8p_2","chr12q_1","chr12q_2",
                                       "chr8p_1","chr1p_2"),clone_order],cluster_rows = F,cluster_cols = F,border_color = "black",
                   color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),gaps_col = c(3),cellwidth = 10,cellheight = 10,fontsize = 8)


#Odds ratio calculating
CNV_clonal_region<-CNV_clonal_region[rowSums(CNV_clonal_region == "0") != 13,]
CNV_clonal_region<-CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr13q","chr5p","chr9q","chr8q","chr12p","chr22q","chr8p_2","chr12q_1","chr12q_2",
                                       "chr8p_1","chr1p_2"),clone_order]

group1 <- abs(CNV_clonal_region[, 1:3])
group2 <- abs(CNV_clonal_region[, 4:ncol(CNV_clonal_region)])
ncol(group2)

OR_df<-matrix(nrow=15,ncol=4)
rownames(OR_df)<-rownames(CNV_clonal_region)
colnames(OR_df) <- c("OR","lower","upper","p")
program <- c('CNV', 'WT')
outcome <- c('Low', 'High')
library(epitools)
for (i in 0:14){
  i=i+1
  b<-sum(group1[i,])
  c<-sum(group2[i,])
  data <- matrix(c(c, b, 10-c, 3-b), nrow=2, ncol=2, byrow=TRUE)
  dimnames(data) <- list('Program'=program, 'Outcome'=outcome)
  OR_df[i,1:4]<-c(oddsratio(data,method = "small")$measure[2],oddsratio(data,method = "small")$measure[4],
                  oddsratio(data,method = "small")$measure[6],oddsratio(data,method = "small")$p.value[6])
}

OR_df<-as.data.frame(OR_df)
#OR_df=OR_df[order(OR_df[,1],decreasing = F),]
OR_df$CNV<-rownames(OR_df)
#OR_df<-OR_df[c(1:4,6,7,9,10,11,5,8,12,13,15,14),]
OR_df$CNV <- factor(OR_df$CNV,levels = rev(OR_df$CNV))
library(ggplot2)
#OR_df$OR = OR_df$OR + 0.01
theme_set(theme_cowplot())

colnames(OR_df)[5]<-"CNV"
ggplot(OR_df, aes(x = OR, y = CNV)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = `upper`, xmin = `lower`), size = 0.6, 
                 height = 0.4,color="black") +
  geom_point(size = 2,color="black") +
  ylab("") + xlab("Odds ratio")+coord_trans(x = 'log10')+
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(0.1,1,10,100))

