# ICB response
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig5/')

library(survival)
library(ggplot2)
library(survminer)
MSK_cohort <- readRDS("../../Data/Fig4_5/MSK_cohort.rds")
# ICB survival CNV
MSK_cohort$ID<-gsub("-",".",MSK_cohort$Sample.ID)
ICB_CNV_all<-read.csv("../../Data/Fig4_5/all_data_by_genes.txt",sep = "\t")
info<-ICB_CNV_all[,1:3]
ICB_CNV_all<-ICB_CNV_all[,colnames(ICB_CNV_all) %in% MSK_cohort$ID]
ICB_CNV_all<-cbind(info,ICB_CNV_all)
ICB_CNV_all$arm <- sub("(p|q).*", "\\1", ICB_CNV_all$Cytoband)

ICB_CNV_all_df<-data.frame(chr20q_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "20q",4:74]),
                           chr20p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "20p",4:74]),
                           chr17p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "17p",4:74]),
                           chr7p_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "7p",4:74]),
                           chr7q_CNV=colMeans(ICB_CNV_all[ICB_CNV_all$arm == "7q",4:74]))

ICB_CNV_all_df_new <- apply(ICB_CNV_all_df, 2, function(x) {
  # 使用ifelse函数根据条件选择值
  ifelse(x > 0.35, "gain", ifelse(x < -0.35, "loss", "WT"))
})
ICB_CNV_all_df_new<-as.data.frame(ICB_CNV_all_df_new)
ICB_CNV_all_df_new$ID <- rownames(ICB_CNV_all_df_new)

MSK_cohort_with_CNV_status<-merge(MSK_cohort,ICB_CNV_all_df_new,by="ID")
MSK_cohort_with_CNV_status$CNV_status<-paste(MSK_cohort_with_CNV_status[,52],MSK_cohort_with_CNV_status[,53],MSK_cohort_with_CNV_status[,54],
                                             MSK_cohort_with_CNV_status[,55],MSK_cohort_with_CNV_status[,56])
#survival
MSK_cohort_with_CNV_status$status<-1
MSK_cohort_with_CNV_status$status[MSK_cohort_with_CNV_status$Overall.Survival.Status.x == "1:DECEASED"]<-2
colnames(MSK_cohort_with_CNV_status)[15]<-"time"
#msk_ICB_cohort
fit <- survfit(Surv(time,status) ~ chr20q_CNV,  # 创建生存对象 
               data = MSK_cohort_with_CNV_status)
p1<-ggsurvplot(fit, data = MSK_cohort_with_CNV_status,pval = TRUE,palette = c("#FF7F00","#226C9E"))+
  labs(x="Time post first dose of ICB (Months)",title = "CRC ICB cohort (MSK-IMPACT)")
p1

#TCGA survival
TCGA_COAD_origin_tumor<-readRDS("../../Data/Fig4_5/TCGA_COAD_origin_tumor.RDS")
clinical <- read.table('../../Data/Fig4_5/COAD_clinical_filter.txt',sep="\t",header=T)
clinical$vital_status <- ifelse(clinical$vital_status=="Alive", 0, 1)
clinical$days_to_death <- gsub('--',-1,clinical$days_to_death)
clinical$days_to_last_follow_up <- gsub('--',-1,clinical$days_to_last_follow_up)
clinical$time <- ifelse(clinical$vital_status==0, clinical$days_to_last_follow_up, clinical$days_to_death)
clinical_filter <- clinical[-which(clinical$time < 10),]
clinical$time <- as.numeric(clinical$time)/365
clinical_filter$time <- as.numeric(clinical_filter$time)/365
colnames(clinical_filter)[1]<-"case"
clinical_COAD_CNV<-merge(clinical_filter,TCGA_COAD_origin_tumor,by="case")

TCGA_COAD_origin_tumor_surv <- survfit(Surv(clinical_COAD_CNV$time,event=clinical_COAD_CNV$vital_status)~chr20q_status,
                                       data=clinical_COAD_CNV)
p2<-ggsurvplot(TCGA_COAD_origin_tumor_surv, data = clinical_COAD_CNV,pval = TRUE,palette = c("#FF7F00","#226C9E"))+
  labs(x="Time (Months)",title = "CRC non-ICB cohort (TCGA-COAD)")
p2

