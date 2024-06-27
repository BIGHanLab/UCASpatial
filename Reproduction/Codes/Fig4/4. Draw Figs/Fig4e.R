# 4. Celltype infiltration and clonalTME phenotyping
library(dplyr)

CRC1_st<-readRDS("../../Data/Fig4_5/clone_region/CRC1_st_clone.RDS")
CRC2_st<-readRDS("../../Data/Fig4_5/clone_region/CRC2_st_clone.RDS")
CRC3a_st<-readRDS("../../Data/Fig4_5/clone_region/CRC3a_st_clone.RDS")
CRC3b_st<-readRDS("../../Data/Fig4_5/clone_region/CRC3b_st_clone.RDS")
CRC4_st<-readRDS("../../Data/Fig4_5/clone_region/CRC4_st_clone.RDS")
CRC5_st<-readRDS("../../Data/Fig4_5/clone_region/CRC5_st_clone.RDS")
CRC6_st<-readRDS("../../Data/Fig4_5/clone_region/CRC6_st_clone.RDS")
CRC7_st<-readRDS("../../Data/Fig4_5/clone_region/CRC7_st_clone.RDS")

Celltype_num=34

inhouse_cohort<-rbind(CRC1_st@meta.data[,c(4:37,40)],CRC2_st@meta.data[,c(4:37,40)],CRC3a_st@meta.data[,c(4:37,40)],CRC3b_st@meta.data[,c(4:37,40)])
public_cohort<-rbind(CRC4_st@meta.data[,c(5:38,41)],CRC5_st@meta.data[,c(4:37,40)],CRC6_st@meta.data[,c(4:37,40)],CRC7_st@meta.data[,c(4:37,40)])
inhouse_cohort$slice<-c(rep("CRC1",nrow(CRC1_st@meta.data)),rep("CRC2",nrow(CRC2_st@meta.data)),rep("CRC3a",nrow(CRC3a_st@meta.data)),rep("CRC3b",nrow(CRC3b_st@meta.data)))
public_cohort$slice<-c(rep("CRC4",nrow(CRC4_st@meta.data)),rep("CRC5",nrow(CRC5_st@meta.data)),rep("CRC6",nrow(CRC6_st@meta.data)),rep("CRC7",nrow(CRC7_st@meta.data)))

public_cohort$sample<-paste(public_cohort$slice,public_cohort$Region)
inhouse_cohort$sample<-paste(inhouse_cohort$slice,inhouse_cohort$Region)
ncol(public_cohort)

total_sample_porp<-rbind(public_cohort,inhouse_cohort)
total_sample_tumor_prop<-total_sample_porp[total_sample_porp$Region %in% c("Clone1","Clone2","Clone3"),]


#total_sample_tumor_prop[,1:Celltype_num]


# cutoff: proportion > 0.02
#table(total_sample_tumor_prop$slice %in% c("CRC1","CRC2","CRC3a","CRC3b")
total_sample_tumor_prop[,1:Celltype_num] <- apply(total_sample_tumor_prop[,1:Celltype_num], 2, function(x) ifelse(x < 0.03, 0, x))

# 将每行的值除以该行的总和，使每行的总和等于1
total_sample_tumor_prop[,1:Celltype_num] <- t(apply(total_sample_tumor_prop[,1:Celltype_num], 1, function(x) x / sum(x)))


tumor_prop_avg<-aggregate(total_sample_tumor_prop[,1:Celltype_num],by = list(total_sample_tumor_prop[,Celltype_num+3]),mean)
rownames(tumor_prop_avg)<- tumor_prop_avg$Group.1
tumor_prop_avg<-tumor_prop_avg[,2:ncol(tumor_prop_avg)] %>% t() %>% as.data.frame()
tumor_prop_avg<-t(tumor_prop_avg)
#saveRDS(tumor_prop_avg,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/clonalTME_prop_avg.RDS")

#clonalTME phenotyping
colnames(tumor_prop_avg)
#tumor_prop_avg<-as.data.frame(tumor_prop_avg[,2:ncol(tumor_prop_avg)])
set.seed(1024)

km.res<-kmeans(scale(tumor_prop_avg[,c(2:6,28:30)]),centers = 2,nstart = 1024)
library(factoextra)
library(cowplot)
theme_set(theme_cowplot())
fviz_cluster(km.res, data = scale(tumor_prop_avg[,c(2:6,28:30)]),
             frame.alpha = 0, frame.level = 0.7, geom = "text",main = "clonalTME Clustering")+theme_cowplot()


tumor_prop_avg<-as.data.frame(tumor_prop_avg)
tumor_prop_avg$subTME_category<-c(rep("T cell-rich",1),rep("T cell-excluded",4),rep("T cell-rich",2),rep("T cell-excluded",6))

colnames(tumor_prop_avg)
my_colour = list(subTME_category = c(`T cell-rich` = "#E78844", `T cell-excluded`="#568DB0"),slide=c(CRC1="#B2DF8A",CRC2="#EC6749",CRC3a="#FB9A99",CRC3b="#FDBF6F",CRC4="#CAB2D6", CRC5="#FFFF99", CRC6="#DDCC35",CRC7="#B15928"),
                 clone=c(Clone1=alpha("#B7A4C7",0.5),Clone2="#B7A4C7",Clone3="#553184",Monoclone="white"))

split_names <- strsplit(rownames(tumor_prop_avg), " ")

# 将分割后的名字作为新的列添加到数据框中
tumor_prop_avg$slide <- sapply(split_names, "[", 1)
tumor_prop_avg$clone <- sapply(split_names, "[", 2)


tumor_prop_avg[,c(36,37)]
tumor_prop_avg$clone[c(3:5,12:13)]<-"Monoclone"
levels(rich_Roe$features.plot)
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)]),cluster_cols = F,cluster_rows = F,border_color = NA,
                   gaps_col = c(3,8),gaps_row = c(3,8,10,20,22,25,29),
                   scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                   annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8)
ggsave("subTME_totalcell_heatmap.pdf",width = 6,height = 15)

#celltype difference between each TME phenotype
library(tidyr)
library(dplyr)

tumor_prop_avg_df <- as.data.frame(tumor_prop_avg)
tumor_prop_avg_df$Sample <- rownames(tumor_prop_avg_df)
colnames(tumor_prop_avg)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(tumor_prop_avg))))))))))))
colnames(total_sample_tumor_prop)<-gsub("Epithelial[+]","Epithelial",gsub("Th1[+] ","Th1/",gsub("Smooth[+]","Smooth",gsub("TM4SF1[+] ","",gsub("LEFTY1[+] ","",gsub("MUC2[+] ","",gsub("REG1A[+] ","",gsub("IGF1","C1QC",gsub("lium","lial",gsub("[.]","+ ",gsub("[.][.]","+ ",colnames(total_sample_tumor_prop))))))))))))


celltype_list<-colnames(tumor_prop_avg)[c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)]

#calculate Roe
total_sample_tumor_prop$subTME<-"T cell-excluded"
total_sample_tumor_prop$subTME[total_sample_tumor_prop$sample %in% c("CRC1 Clone1","CRC4 Clone1","CRC4 Clone2")] <- "T cell-rich"
library(dplyr)
total_sample_tumor_prop %>% group_by(subTME)
Roe_total<-aggregate(total_sample_tumor_prop[,1:34],by = list(total_sample_tumor_prop[,38]),mean)
rownames(Roe_total)<- Roe_total$Group.1
Roe_total<-Roe_total[,2:ncol(Roe_total)] %>% t() %>% as.data.frame()

Calculate_Roe <- function(avg_df)
{
  Ro_e_ex <- avg_df %>% dplyr::mutate(RowSum = rowSums(avg_df))
  NRow <- Ro_e_ex$RowSum
  Ro_e_ex <- Ro_e_ex%>% t() %>% as.data.frame() %>% dplyr::mutate(ColSum = colSums(Ro_e_ex))
  NCol <- Ro_e_ex$ColSum
  Roe_df <- avg_df
  N <- Ro_e_ex['RowSum','ColSum']
  for (i in rownames(avg_df)) {
    for (j in colnames(avg_df)) {
      Roe_df[i,j] = N*Roe_df[i,j] / (NRow[i] * NCol[j])
    }
  }
  return(Roe_df)
}


Roe_v_inhouse <- Calculate_Roe(Roe_total)
celltype_list[1:29]
Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
  select(celltype_list[1:29]) %>% t() %>% as.data.frame()
# Roe_v_inhouse <- Roe_v_inhouse %>% t() %>% as.data.frame() %>% 
#   select(names(sort(kmeans_Roe_v_inhouse[["cluster"]]))) %>% t() %>% as.data.frame()

#colnames(Roe_v_inhouse) <- c('Dessert','Inflamed','Suppression')
max(Roe_v_inhouse)
min(Roe_v_inhouse)
library("tidyr")
#dp <- DotPlot(NC_CRC1,features = 'CD8A',group.by = 'Region')$data
dotplot_Roe_v_inhouse <- pivot_longer(as.data.frame(t(Roe_v_inhouse)) %>% mutate(id=colnames(Roe_v_inhouse)),
                                      cols=rownames(Roe_v_inhouse),names_to='features.plot',values_to = 'avg.exp')


dotplot_Roe_v_inhouse <- dotplot_Roe_v_inhouse %>% mutate(Group = avg.exp,
                                                          'Ro/e' = avg.exp)
dotplot_Roe_v_inhouse[,3:5][dotplot_Roe_v_inhouse[,3:5]>2] <- 2
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8] <- 2-dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<0.8]
dotplot_Roe_v_inhouse$Group[dotplot_Roe_v_inhouse$Group<1.2 & dotplot_Roe_v_inhouse$Group >0.8] <- 0.5

#dotplot_Roe_v_inhouse$id<-factor(dotplot_Roe_v_inhouse$id,levels = c("TLS","Normal Epithelial","Tumor","Tumor Boundary","Remaining"))

dotplot_Roe_v_inhouse$features.plot <- as.factor(dotplot_Roe_v_inhouse$features.plot)
dotplot_Roe_v_inhouse$features.plot <- factor(dotplot_Roe_v_inhouse$features.plot,levels = arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot)
arrange(dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id == "T cell-rich",],`Ro/e`)$features.plot


dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]

rich_Roe<-dotplot_Roe_v_inhouse[dotplot_Roe_v_inhouse$id %in% c("T cell-rich"),]
rich_Roe$group<-"Remain"
rich_Roe$group[rich_Roe$`Ro/e` > 1.2] <- "Enrich"
rich_Roe$group[rich_Roe$`Ro/e` < 0.8] <- "Excluded"
rich_Roe$`Ro/e_adj`<-rich_Roe$`Ro/e`-1
theme_set(theme_cowplot())

ggplot(rich_Roe,aes(y=features.plot, x= `Ro/e_adj`, fill=group,color=group))+
  geom_bar(stat = "identity",width = 0.05)+geom_point(size=4,shape=21,stroke=1)+
  scale_fill_manual(values = c("#FDBF6F","#A6CEE3","white"))+scale_color_manual(values = c("#FF7F00","#4A91C1","darkgrey"))+
  geom_vline(xintercept =0,lty=2)+coord_cartesian(xlim = c(-0.8,0.8))+labs(x="Ro/e",y="Subset")

ggsave("Celltype_enrich_subTME.pdf",height = 10,width = 5)

pheatmap::pheatmap(t(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),rev(levels(rich_Roe$features.plot))]),cluster_cols = F,cluster_rows = F,border_color = NA,
                   gaps_col = c(3,8),gaps_row = c(13,26),
                   scale = "row",color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),annotation_col = tumor_prop_avg[,c(35:37)],
                   annotation_colors = my_colour,show_colnames = F,cellwidth = 10,cellheight = 10,fontsize = 8)
