library(Seurat)
library(pacman)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(forcats)
library(ggradar)
p_unload(Seurat)
p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/')

source('../source/UCASpatial_final_v1.R')
source('../source/Test_Accuracy.R')
col.sc <- readRDS('../../Data/Fig2/col.sc.re.rds')
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")

plot_cols <- c(brewer.pal(name = "Set3",n=12)[c(-1:-3)],brewer.pal(name = "Set1",n=9))[-5]
ggradar_col <- c(brewer.pal(name = "Set3",n=12)[c(-1:-3)],brewer.pal(name = "Set1",n=9))

simu_TC <- readRDS("../../Data/Fig2/simu_TC.rds")
simu_TM <- readRDS("../../Data/Fig2/simu_TM.rds")
simu_TS <- readRDS("../../Data/Fig2/simu_TS.rds")


#### Test the performance ####
simulatData_TME <- c(simu_TC,simu_TM,simu_TS)
normal_cutoff <- 0.05
e_cutoff <- 0.05
## UCASpatial
{
  UCASpatial_highres <- readRDS("../../Data/Fig2/UCASpatial_result/UCASpatial_TME_highres.rds")
  
  performance_UCASpatial_highres <- mapply(function(simulatData_list,UCASpatial_result)
  {
    
    UCASpatial_deconv_rs <- UCASpatial_result[[2]][, colnames(UCASpatial_result[[2]]) != "res_ss"]
    UCASpatial_deconv_rs[UCASpatial_deconv_rs < e_cutoff] <- 0
    UCASpatial_deconv_rs <- as.matrix(UCASpatial_deconv_rs)/rowSums(as.matrix(UCASpatial_deconv_rs))
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    # test_Pearson(UCASpatial_deconv_rs,synthetic_comp)
    result <- test_MultiIndex(deconv_result = UCASpatial_deconv_rs[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },
  simulatData_TME,UCASpatial_highres)
}

## RCTD
{
  RCTD_highres <- readRDS("../../Data/Fig2/RCTD_result/RCTD_TME_highres.rds")
  
  performance_RCTD_highres <- mapply(function(simulatData_list,myRCTD_full){
    RCTD_deconv <- as.matrix(myRCTD_full@results[["weights"]])
    RCTD_deconv[RCTD_deconv < normal_cutoff] <- 0
    RCTD_deconv <- as.matrix(RCTD_deconv)/rowSums(as.matrix(RCTD_deconv))
    colnames(RCTD_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))]
    RCTD_deconv <- RCTD_deconv[,match(colnames(synthetic_comp) ,colnames(RCTD_deconv))]
    result <- test_MultiIndex(deconv_result = RCTD_deconv[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },simulatData_TME,RCTD_highres)
  
}

## SPOTlight
{
  SPOTlight_highres <- readRDS("../../Data/Fig2/SPOTlight_result/SPOTlight_TME_highres.rds")
  
  performance_SPOTlight_highres <- mapply(function(simulatData_list,spotlight_ls){
    
    spotlight_deconv <- spotlight_ls[[2]][, colnames(spotlight_ls[[2]]) != "res_ss"]
    spotlight_deconv[spotlight_deconv < normal_cutoff] <- 0
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    #test_Pearson(spotlight_deconv,synthetic_comp)
    result <- test_MultiIndex(deconv_result = spotlight_deconv[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },simulatData_TME,SPOTlight_highres)
}

## cell2location
{
  suffix_list = as.list(paste(c(rep('TC',5),rep('TM',5),rep('TS',5)),'_simu_',rep(c(1:5),3),sep = ''))
  C2L_path = "../../Data/Fig2/cell2location_result/results/"
  
  C2L_highres <- mapply(function(simulatData_list,suffix_list){
    C2L_deconv <- read.csv(paste(C2L_path,suffix_list,"_highres.csv",sep = ""),row.names = 1)
    C2L_deconv <- as.matrix(C2L_deconv)/rowSums(as.matrix(C2L_deconv))
    C2L_deconv[C2L_deconv < 0.1] <- 0
    C2L_deconv <- as.matrix(C2L_deconv)/rowSums(as.matrix(C2L_deconv))
    colnames(C2L_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(C2L_deconv))
    return(list(C2L_deconv))
  },simulatData_TME,suffix_list)
  performance_C2L_highres <- mapply(function(simulatData_list,C2L_deconv){
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(C2L_deconv))]
    C2L_deconv <- C2L_deconv[,match(colnames(synthetic_comp) ,colnames(C2L_deconv))]
    result <- test_MultiIndex(deconv_result = C2L_deconv[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },simulatData_TME,C2L_highres)
  
}

## CARD
{
  suffix_list = as.list(paste(c('simu_TC/','simu_TM/','simu_TS/'),1:5,'_CARD_res.txt',sep=''))
  CARD_path = "../../Data/Fig2/"
  res_list <- list.dirs("CARD_result") %>% 
    grep("high",.,value = T) %>% 
    list.files(.,"txt",full.names = T,recursive = T)
  CARD_highres <- mapply(function(simulatData_list,suffix_list){
    CARD_deconv <- read.table(paste(CARD_path,suffix_list,sep = ''),row.names = 1)
    CARD_deconv <- as.matrix(CARD_deconv)/rowSums(as.matrix(CARD_deconv))
    CARD_deconv[CARD_deconv < normal_cutoff] <- 0
    CARD_deconv <- as.matrix(CARD_deconv)/rowSums(as.matrix(CARD_deconv))
    colnames(CARD_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(CARD_deconv))
    return(list(CARD_deconv))
  },simulatData_TME,res_list)
  
  performance_CARD_highres <- mapply(function(simulatData_list,CARD_deconv){
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(CARD_deconv))]
    CARD_deconv <- CARD_deconv[,match(colnames(synthetic_comp) ,colnames(CARD_deconv))]
    result <- test_MultiIndex(deconv_result = CARD_deconv[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },simulatData_TME,CARD_highres)
  
}

## stereoscope
{
  suffix_list = as.list(paste(c('simu_TC/','simu_TM/','simu_TS/'),1:5,'_stereoscope_res.txt',sep=''))
  stereoscope_path = "../../Data/Fig2/"
  res_list <- list.dirs('stereoscope_result') %>% 
    list.files(.,"txt",full.names = T,recursive = T) %>%
    grep("high",.,value = T)
  
  stereoscope_highres <- mapply(function(simulatData_list,suffix_list){
    stereoscope_deconv <- read.table(paste(stereoscope_path,suffix_list,sep = ''),row.names = 1)
    stereoscope_deconv <- as.matrix(stereoscope_deconv)/rowSums(as.matrix(stereoscope_deconv))
    stereoscope_deconv[stereoscope_deconv < normal_cutoff] <- 0
    stereoscope_deconv <- as.matrix(stereoscope_deconv)/rowSums(as.matrix(stereoscope_deconv))
    colnames(stereoscope_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(stereoscope_deconv))
    return(list(stereoscope_deconv))
  },simulatData_TME,res_list)
  
  performance_stereoscope_highres <- mapply(function(simulatData_list,stereoscope_deconv){
    synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                  rowSums(simulatData_list[[2]][[3]]))
    colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
    synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(stereoscope_deconv))]
    stereoscope_deconv <- stereoscope_deconv[,match(colnames(synthetic_comp) ,colnames(stereoscope_deconv))]
    result <- test_MultiIndex(deconv_result = stereoscope_deconv[, colnames(synthetic_comp)],
                              synthetic_comp = synthetic_comp)
    return(list(result))
  },simulatData_TME,stereoscope_highres)
  
}


#### Downsteam analysis ####

# Fig. 2b Construct entropy gradient
{
  # ground truth results
  {
    groud_truth <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(groud_truth) <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
    for(i in 1:15)
    {
      synthetic_comp <- as.matrix(simulatData_TME[[i]][["cell_composition"]][[3]] /
                                    rowSums(simulatData_TME[[i]][["cell_composition"]][[3]]))
      colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
      groud_truth <- rbind(groud_truth,synthetic_comp)
    }
    groud_truth$entropy <- 0
    for (i in 1:nrow(groud_truth)) {
      groud_truth$entropy[i] <- entropy::entropy(as.numeric(groud_truth[i,]),unit = 'log2')
    }
  }
  
  # Fig. 2b RMSE&F1 score
  {
    # Calculate RMSE spot level
    {
      RMSE_SPOTlight_highres <- mapply(function(simulatData_list,spotlight_ls){
        
        spotlight_deconv <- spotlight_ls[[2]][, colnames(spotlight_ls[[2]]) != "res_ss"]
        spotlight_deconv[spotlight_deconv < normal_cutoff] <- 0
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
        result <- test_RMSE_spotlevel(spotlight_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,SPOTlight_highres)
      
      RMSE_UCASpatial_highres <- mapply(function(simulatData_list,UCASpatial_ls){
        
        UCASpatial_deconv <- UCASpatial_ls[[2]][, colnames(UCASpatial_ls[[2]]) != "res_ss"]
        # UCASpatial_deconv[UCASpatial_deconv < e_cutoff] <- 0
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
        result <- test_RMSE_spotlevel(UCASpatial_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,UCASpatial_highres)
      
      RMSE_RCTD_highres <- mapply(function(simulatData_list,myRCTD_full){
        RCTD_deconv <- as.matrix(myRCTD_full@results[["weights"]])
        # RCTD_deconv[RCTD_deconv < normal_cutoff] <- 0
        colnames(RCTD_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))]
        RCTD_deconv <- RCTD_deconv[,match(colnames(synthetic_comp) ,colnames(RCTD_deconv))]
        result <- test_RMSE_spotlevel(RCTD_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,RCTD_highres)
      
      RMSE_C2L_highres <- mapply(function(simulatData_list,C2L_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(C2L_deconv))]
        C2L_deconv <- C2L_deconv[,match(colnames(synthetic_comp) ,colnames(C2L_deconv))]
        # C2L_deconv[C2L_deconv < normal_cutoff] <- 0
        result <- test_RMSE_spotlevel(C2L_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,C2L_highres)
      
      RMSE_CARD_highres <- mapply(function(simulatData_list,CARD_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(CARD_deconv))]
        CARD_deconv <- CARD_deconv[,match(colnames(synthetic_comp) ,colnames(CARD_deconv))]
        result <- test_RMSE_spotlevel(CARD_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,CARD_highres)
      
      RMSE_stereoscope_highres <- mapply(function(simulatData_list,stereoscope_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(stereoscope_deconv))]
        stereoscope_deconv <- stereoscope_deconv[,match(colnames(synthetic_comp) ,colnames(stereoscope_deconv))]
        result <- test_RMSE_spotlevel(stereoscope_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,stereoscope_highres)
      
    }
    # evaluate 
    {
      spot_level_RMSE <- as.data.frame(t(rbind(unlist(RMSE_UCASpatial_highres),unlist(RMSE_SPOTlight_highres),
                                               unlist(RMSE_RCTD_highres),unlist(RMSE_C2L_highres),
                                               unlist(RMSE_CARD_highres),unlist(RMSE_stereoscope_highres))))
      spot_level_RMSE_ent <- cbind(spot_level_RMSE,groud_truth$entropy)
      for(i in 1:ncol(spot_level_RMSE_ent))
      {
        spot_level_RMSE_ent[,i] <- as.numeric(spot_level_RMSE_ent[,i])
      }
      colnames(spot_level_RMSE_ent) <- c("UCASpatial","SPOTlight","RCTD",'cell2location','CARD','stereoscope','entropy')
      colMeans(spot_level_RMSE_ent)
      max(spot_level_RMSE_ent[,-7])
    }
    # Plot
    {
      
      spot_level_RMSE_ent_plot <- arrange(spot_level_RMSE_ent,entropy)
      spot_level_RMSE_ent_plot$group <- c(rep('low',250),rep('medium',250),rep('high',250))
      
      spot_level_RMSE_ent_plot_l <- pivot_longer(spot_level_RMSE_ent_plot,
                                                 -c(entropy,group),names_to = "algorithm",values_to = "RMSE")
      spot_level_RMSE_ent_plot_l$group <- factor(spot_level_RMSE_ent_plot_l$group,levels = c('low','medium','high'))
      spot_level_RMSE_ent_plot_l$algorithm <- 
        factor(spot_level_RMSE_ent_plot_l$algorithm,
               levels = c("UCASpatial","SPOTlight","RCTD",'cell2location','CARD','stereoscope','entropy'))
      
      dark_cols <- c('#1F78B4','#33A02C','#FE7F00')
      light_cols <- c('#A6CEE3','#B2DF8A','#FDBF6F')
      
      
      ggboxplot(spot_level_RMSE_ent_plot_l,x = "algorithm",y = "RMSE",#facet.by = 'algorithm',
                title = "RMSE along entropy gradient",fill = "group",color = "group",lwd= 0.8)+
        scale_fill_manual(values = alpha(light_cols,alpha = 0.5))+
        scale_color_manual(values = dark_cols)+
        theme_bw(base_size = 14)+ylab('RMSE')+xlab(NULL)+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+ylim(c(0,max(spot_level_RMSE_ent_plot_l$RMSE)))#+NoLegend()
      ggsave("entropy_gradient_RMSE_boxplot_splitbyalg.pdf",width = 7,height = 4)
      
      
      
      spot_level_RMSE_ent_plot_l$algorithm <- 
        factor(spot_level_RMSE_ent_plot_l$algorithm,
               levels = c("UCASpatial","RCTD","SPOTlight",'CARD','cell2location','stereoscope'))
      
      ggboxplot(spot_level_RMSE_ent_plot_l,x = "group",y = "RMSE",#facet.by = 'algorithm',
                title = "RMSE along entropy gradient",fill = "algorithm",lwd= 0.4)+
        scale_fill_manual(values = (plot_cols)[c(1,3,2,5,4,6)])+
        scale_color_manual(values = plot_cols)+
        theme_bw(base_size = 14)+ylab('RMSE')+xlab('Chaotic Level')+
        theme(axis.text.x = element_text(size = 13),plot.title = element_text(hjust = 0.5))+
        # stat_compare_means()+
        theme(legend.position = "right",)+
        ylim(c(0,max(spot_level_RMSE_ent_plot_l$RMSE)))#+NoLegend()
      ggsave("entropy_gradient_RMSE_boxplot_splitbyent.pdf",width = 7,height = 4)
      
      
      spot_level_RMSE_ent_plot_l_summ <- spot_level_RMSE_ent_plot_l %>% 
        group_by(algorithm) %>% 
        dplyr::summarise(Mean=median(RMSE))
      spot_level_RMSE_ent_plot_l_summ$Mean_div <- spot_level_RMSE_ent_plot_l_summ$Mean-spot_level_RMSE_ent_plot_l_summ$Mean[1]
      spot_level_RMSE_ent_plot_l_summ$impr <- spot_level_RMSE_ent_plot_l_summ$Mean_div/spot_level_RMSE_ent_plot_l_summ$Mean*100
      
    }
    
    ## F1 score level (PRC)
    {
      # summary Acur,F1,... (with 10 division or not)
      {
        All_auc <- lapply(highres_list, function(x){
          Predict <- as.matrix(x$proportion_filter)
          Truth <- as.matrix(groud_truth_split$proportion)
          colnames(Predict) <- 'Value'
          colnames(Truth) <- 'Value'
          auc <- Test_Acurracy(Predict,Truth)
          return(auc)
        })
        All_auc_value <- matrix(unlist(All_auc) ,nr=5) %>% as.data.frame()
        rownames(All_auc_value) <- colnames(All_auc[[1]])
        colnames(All_auc_value) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
        # saveRDS(All_auc_value,'All_auc_value_0.05cutoff.rds')
        
        
        All_F1_spot <- lapply(highres_list, function(x){
          F1 <- as.data.frame(matrix(nrow=750,ncol=1))
          for(i in 1:750){
            Predict <- as.matrix(x$proportion_filter[((i-1)*20+1):((i)*20)])
            Truth <- as.matrix(groud_truth_split$proportion[((i-1)*20+1):((i)*20)])
            colnames(Predict) <- 'Value'
            colnames(Truth) <- 'Value'
            F1[i,] <- Test_Acurracy(Predict,Truth)[5]
          }
          return(F1)
        })
        All_F1_spot_df <- matrix(unlist(All_F1_spot) ,nr=750) %>% as.data.frame()
        colMeans(All_F1_spot_df)
        colnames(All_F1_spot_df) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      }
      # draw line chart of F1 for all methods
      {
        spot_level_F1_ent <- cbind(All_F1_spot_df,spot_level_RMSE_ent$entropy)
        colnames(spot_level_F1_ent)[7] <- 'entropy'
        spot_level_F1_ent_plot <- arrange(spot_level_F1_ent,entropy)
        spot_level_F1_ent_plot$group <- c(rep('low',250),rep('medium',250),rep('high',250))
        
        spot_level_F1_ent_plot_l <- pivot_longer(spot_level_F1_ent_plot,
                                                 -c(entropy,group),names_to = "algorithm",values_to = "F1")
        spot_level_F1_ent_plot_l$group <- factor(spot_level_F1_ent_plot_l$group,levels = c('low','medium','high'))
        spot_level_F1_ent_plot_l$algorithm <- 
          factor(spot_level_F1_ent_plot_l$algorithm,
                 levels = c("UCASpatial","SPOTlight","RCTD",'cell2location','CARD','stereoscope'))
        
        dark_cols <- c('#1F78B4','#33A02C','#FE7F00')
        light_cols <- c('#A6CEE3','#B2DF8A','#FDBF6F')
        
        
        ggboxplot(spot_level_F1_ent_plot_l,x = "algorithm",y = "F1",#facet.by = 'algorithm',
                  title = "F1 along entropy gradient",fill = "group",color = "group",lwd= 0.8)+
          scale_fill_manual(values = alpha(light_cols,alpha = 0.5))+
          scale_color_manual(values = dark_cols)+
          theme_bw(base_size = 14)+ylab('F1')+xlab(NULL)+
          theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
          theme(legend.position = "right")+ylim(c(0,max(spot_level_F1_ent_plot_l$F1)))#+NoLegend()
        ggsave("entropy_gradient_F1_boxplot_splitbyalg.pdf",width = 7,height = 4)
        
        
        
        spot_level_F1_ent_plot_l$algorithm <- 
          factor(spot_level_F1_ent_plot_l$algorithm,
                 levels = c("UCASpatial",'stereoscope',"SPOTlight","RCTD",'CARD','cell2location'))
        
        ggboxplot(spot_level_F1_ent_plot_l,x = "group",y = "F1",#facet.by = 'algorithm',
                  title = "F1 score along entropy gradient",fill = "algorithm",lwd= 0.4)+
          scale_fill_manual(values = (plot_cols)[c(1,6,2,3,5,4)])+
          scale_color_manual(values = plot_cols)+
          theme_bw(base_size = 14)+ylab('F1 score')+xlab('Chaotic Level')+
          theme(axis.text.x = element_text(size = 13),plot.title = element_text(hjust = 0.5))+
          # stat_compare_means()+
          theme(legend.position = "right",)+
          ylim(c(0,max(spot_level_F1_ent_plot_l$F1)))#+NoLegend()
        ggsave("entropy_gradient_F1_boxplot_splitbyent.pdf",width = 7,height = 4)

        
        ggbarplot(spot_level_F1_ent_plot_l,x= "algorithm",y = "F1",add =c("mean_se"),facet.by = 'group',
                  fill = "algorithm",#color = "white",
                  position = position_dodge())+
          scale_fill_manual(values = (plot_cols)[c(1,4,3,5,6,2)])+
          scale_color_manual(values = plot_cols)+
          theme_bw(base_size = 14)+ylab('F1 score')+xlab('Chaotic Level')+
          theme(axis.text.x = element_text(size = 13),plot.title = element_text(hjust = 0.5))+
          # stat_compare_means()+
          theme(legend.position = "right",)
        
        spot_level_F1_ent_plot_l_summ <- spot_level_F1_ent_plot_l %>% 
          group_by(algorithm) %>% 
          dplyr::summarise(Mean=median(F1))
        spot_level_F1_ent_plot_l_summ$Mean_div <- spot_level_F1_ent_plot_l_summ$Mean-spot_level_F1_ent_plot_l_summ$Mean[1]
        spot_level_F1_ent_plot_l_summ$impr <- spot_level_F1_ent_plot_l_summ$Mean_div/spot_level_F1_ent_plot_l_summ$Mean*100
      }
    }
  }
  
  # Fig SI. PCC
  {
    # Calculate Pearson correlation spot level
    {
      Pearson_UCASpatial_highres <- mapply(function(simulatData_list,UCASpatial_ls){
        
        UCASpatial_deconv <- UCASpatial_ls[[2]][, colnames(UCASpatial_ls[[2]]) != "res_ss"]
        # UCASpatial_deconv[UCASpatial_deconv < e_cutoff] <- 0
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
        result <- test_Pearson_spotlevel(UCASpatial_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,UCASpatial_highres)
      
      Pearson_SPOTlight_highres <- mapply(function(simulatData_list,spotlight_ls){
        
        spotlight_deconv <- spotlight_ls[[2]][, colnames(spotlight_ls[[2]]) != "res_ss"]
        # spotlight_deconv[spotlight_deconv < normal_cutoff] <- 0
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
        result <- test_Pearson_spotlevel(spotlight_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,SPOTlight_highres)
      
      Pearson_RCTD_highres <- mapply(function(simulatData_list,myRCTD_full){
        RCTD_deconv <- as.matrix(myRCTD_full@results[["weights"]])
        # RCTD_deconv[RCTD_deconv < normal_cutoff] <- 0
        colnames(RCTD_deconv) <- gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_deconv))]
        RCTD_deconv <- RCTD_deconv[,match(colnames(synthetic_comp) ,colnames(RCTD_deconv))]
        result <- test_Pearson_spotlevel(RCTD_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,RCTD_highres)
      
      Pearson_C2L_highres <- mapply(function(simulatData_list,C2L_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(C2L_deconv))]
        C2L_deconv <- C2L_deconv[,match(colnames(synthetic_comp) ,colnames(C2L_deconv))]
        result <- test_Pearson_spotlevel(C2L_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,C2L_highres)
      
      Pearson_CARD_highres <- mapply(function(simulatData_list,CARD_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(CARD_deconv))]
        CARD_deconv <- CARD_deconv[,match(colnames(synthetic_comp) ,colnames(CARD_deconv))]
        result <- test_Pearson_spotlevel(CARD_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,CARD_highres)
      
      Pearson_stereoscope_highres <- mapply(function(simulatData_list,stereoscope_deconv){
        synthetic_comp <- as.matrix(simulatData_list[[2]][[3]] /
                                      rowSums(simulatData_list[[2]][[3]]))
        colnames(synthetic_comp) <- gsub("[_/ /+/-/(/)/-]",".",colnames(synthetic_comp))
        synthetic_comp <- synthetic_comp[,colnames(synthetic_comp)%in%gsub("[_/ /+/-/(/)/-]",".",colnames(stereoscope_deconv))]
        stereoscope_deconv <- stereoscope_deconv[,match(colnames(synthetic_comp) ,colnames(stereoscope_deconv))]
        result <- test_Pearson_spotlevel(stereoscope_deconv,synthetic_comp)
        return(list(result))
      },simulatData_TME,stereoscope_highres)
      
    }
    # evaluate 
    {
      spot_level_PCC <- as.data.frame(t(rbind(unlist(Pearson_UCASpatial_highres),unlist(Pearson_SPOTlight_highres),
                                              unlist(Pearson_RCTD_highres),unlist(Pearson_C2L_highres),
                                              unlist(Pearson_CARD_highres),unlist(Pearson_stereoscope_highres))))
      spot_level_PCC_ent <- cbind(spot_level_PCC,groud_truth$entropy)
      for(i in 1:ncol(spot_level_PCC_ent))
      {
        spot_level_PCC_ent[,i] <- as.numeric(spot_level_PCC_ent[,i])
      }
      colnames(spot_level_PCC_ent) <- c("UCASpatial","SPOTlight","RCTD",'cell2location','CARD','stereoscope','entropy')
      colMeans(spot_level_PCC_ent)
    }
    # Plot
    {
      
      spot_level_PCC_ent_plot <- arrange(spot_level_PCC_ent,entropy)
      spot_level_PCC_ent_plot$group <- c(rep('low',250),rep('medium',250),rep('high',250))
      
      spot_level_PCC_ent_plot_l <- pivot_longer(spot_level_PCC_ent_plot,
                                                -c(entropy,group),names_to = "algorithm",values_to = "PCC")
      spot_level_PCC_ent_plot_l$group <- factor(spot_level_PCC_ent_plot_l$group,levels = c('low','medium','high'))
      spot_level_PCC_ent_plot_l$algorithm <- 
        factor(spot_level_PCC_ent_plot_l$algorithm,
               levels = c("UCASpatial","SPOTlight","RCTD",'cell2location','CARD','stereoscope'))
      
      dark_cols <- c('#1F78B4','#33A02C','#FE7F00')
      light_cols <- c('#A6CEE3','#B2DF8A','#FDBF6F')
      plot_cols <- c(brewer.pal(name = "Set3",n=12)[c(-1:-3)],brewer.pal(name = "Set1",n=9))[-5]
      
      spot_level_PCC_ent_plot_l$algorithm <- 
        factor(spot_level_PCC_ent_plot_l$algorithm,
               levels = c("UCASpatial",'cell2location',"RCTD",'CARD','stereoscope',"SPOTlight"))
      
      ggboxplot(spot_level_PCC_ent_plot_l,x = "group",y = "PCC",#facet.by = 'algorithm',
                title = "PCC along entropy gradient",fill = "algorithm",size = 5,lwd=0.4)+
        scale_fill_manual(values = (plot_cols)[c(1,4,3,5,6,2)])+
        scale_color_manual(values = plot_cols)+
        theme_bw(base_size = 14)+ylab('PCC')+xlab('Chaotic Level')+
        theme(axis.text.x = element_text(size = 13),plot.title = element_text(hjust = 0.5))+
        # stat_compare_means()+
        theme(legend.position = "right")+
        ylim(c(min(spot_level_PCC_ent_plot_l$PCC),max(spot_level_PCC_ent_plot_l$PCC)))#+NoLegend()
      ggsave("entropy_gradient_PCC_boxplot_splitbyent.pdf",width = 7,height = 4)
      
      ggboxplot(spot_level_PCC_ent_plot_l,x = "algorithm",y = "PCC",#facet.by = 'algorithm',
                title = "PCC along entropy gradient",fill = "group",color = "group",lwd= 0.8)+
        scale_fill_manual(values = alpha(light_cols,alpha = 0.5))+
        scale_color_manual(values = dark_cols)+
        theme_bw(base_size = 14)+ylab('PCC')+xlab(NULL)+
        theme(axis.text.x = element_text(),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+ylim(c(0,max(spot_level_PCC_ent_plot_l$PCC)))#+NoLegend()
      ggsave("entropy_gradient_PCC_boxplot_splitbyalg.pdf",width = 7,height = 4)
      
      
      spot_level_PCC_ent_plot_l_summ <- spot_level_PCC_ent_plot_l %>% 
        group_by(algorithm) %>% 
        dplyr::summarise(Mean=median(PCC))
      spot_level_PCC_ent_plot_l_summ$Mean_div <- spot_level_PCC_ent_plot_l_summ$Mean-spot_level_PCC_ent_plot_l_summ$Mean[1]
      spot_level_PCC_ent_plot_l_summ$impr <- spot_level_PCC_ent_plot_l_summ$Mean_div/spot_level_PCC_ent_plot_l_summ$Mean*100
    }
  }

}

# Fig. 2c-e cell type level line chart
{
  ## Fig. 2c-e PCC line chart
  {
    df_Pearson_l_highres <- Pearson_l_highres %>% na.omit() %>% group_by(algorithm,index) %>% 
      mutate(Mean = mean(value)) %>% dplyr::select(algorithm,index,Mean) %>% unique()
    df_Pearson_l_highres <- pivot_wider(df_Pearson_l_highres,names_from = 'index',values_from = 'Mean') %>% as.data.frame()
    
    rownames(df_Pearson_l_highres) <- as.character(df_Pearson_l_highres$algorithm)
    df_Pearson_l_highres <- df_Pearson_l_highres[,-1] 
    
    df_Pearson_l_highres[is.na(df_Pearson_l_highres)] <- 0
    # remove Tangram
    # df_Pearson_l_highres_avg <- df_Pearson_l_highres[-5,] %>% t() %>% as.data.frame()
    df_Pearson_l_highres_avg <- df_Pearson_l_highres %>% t() %>% as.data.frame()
    df_Pearson_l_highres_avg <- df_Pearson_l_highres_avg %>% 
      mutate(avg = rowMeans(.),minm= Rfast::rowMaxs(as.matrix(df_Pearson_l_highres_avg[,-1]),value=T)) %>%
      arrange(minm) %>% t() %>% as.data.frame()
    
    
    df_Pearson_l_highres_normalized <- t(t(df_Pearson_l_highres_avg) - as.data.frame(t(df_Pearson_l_highres_avg))$UCASpatial) %>%
      as.data.frame()
    
    
    df_Pearson_l_highres_normalized <- df_Pearson_l_highres_normalized %>% 
      dplyr::select(c("Stromal.c1","Stromal.c2","Epithelial.c1","Epithelial.c2",
                      "CD8.Trm","CD8.Tem",'CD8.Tex',"NK","CD4.T.CCR7","CD4.Tfh","CD4.Treg",
                      "B.GC","B.memory.like","B.naive.like","Plasma",
                      "Macrophage","Monocyte","DC","pDC","Mast"))
    
    
    # Fig. 2c example overall
    {
      example_overall <- df_Pearson_l_highres_avg %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = colnames(.)[1:20],names_to = 'Celltype',values_to = 'PCC')
      # example_overall <- example_overall %>% filter(Algorithm%in%c('UCASpatial','minm','avg'))
      example_overall$Algorithm <- as.factor(example_overall$Algorithm)
      example_overall$Algorithm <- factor(example_overall$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                               'cell2location','CARD','stereoscope',
                                                                               'avg','minm'))
      example_overall$Celltype <- as.factor(example_overall$Celltype)
      example_overall$Celltype <- factor(example_overall$Celltype,levels = colnames(df_Pearson_l_highres_avg))
      
      p0 <- ggplot(example_overall,aes(Celltype,PCC,group=Algorithm,color=Algorithm,
                                       linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.5,alpha=0.5)+ ylim(0,1)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,avg=2,CARD=3,stereoscope=3,dif = 2,minm = 2))+
        scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,avg=1,CARD=0.5,stereoscope=0.5,dif=1,minm = 1))+
        scale_color_manual(values = ggradar_col[c(1:4,6:7,16,9)])+
        theme_bw()+xlab(NULL)+ylab('PCC')+
        # ggtitle(label = 'PCC for')+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")
      # ggsave('PCC_example_CD8T.pdf',height = 2.6,width = 4)
      p0
      ggsave('PCC_overall_lineplot.pdf',height = 4,width = 5)
      }
    # Fig. 2d example CD8T/NK
    {
      # PCC CD8T/NK
      {
        example_stro_mye <- df_Pearson_l_highres[-5,] %>% dplyr::select("CD8.Tem","CD8.Tex","CD8.Trm","NK") %>% 
          mutate(Algorithm = rownames(.)) %>% 
          pivot_longer(cols = c("CD8.Tem","CD8.Tex","CD8.Trm","NK"),names_to = 'Celltype',values_to = 'PCC')
        example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
        example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                   'cell2location','CARD','stereoscope'))
        
        p3 <- ggplot(example_stro_mye,aes(Celltype,PCC,group=Algorithm,color=Algorithm,
                                          linetype=Algorithm,linewidth=Algorithm))+
          geom_point(size=1.8)+ ylim(0,0.65)+
          geom_line()+
          scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                           cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
          scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                            cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
          scale_color_manual(values = ggradar_col[-5])+
          theme_bw()+xlab(NULL)+ylab('PCC')+labs(title = 'CD8T/NK')+
          # ggtitle(label = 'PCC for')+
          theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
          theme(legend.position = "right")+NoLegend()
        p3
        ggsave('PCC_example_2_CD8T.pdf',height = 4,width = 2)
      }
      # False rate of CD8T/NK
      {
        source('/data/xy/scripts/Summary_True_False_rate.R')
        highres_list_TFR <- list(UCASpatial_highres_all,SPOTlight_highres_all,RCTD_highres_all,C2L_highres_all,
                                 CARD_highres_all,stereoscope_highres_all)
        highres_list_TFR <- lapply(highres_list_TFR, function(x){
          a <- select(x,colnames(UCASpatial_highres_all)[1:20])
          return(a)
        })
        groud_truth <- select(groud_truth,colnames(UCASpatial_highres_all)[1:20])
        FR_ct <- lapply(highres_list_TFR, function(x){
          False_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
          temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
          for(i in 1:20)
          {
            deconv_result <- as.matrix(x[,i])
            deconv_result[deconv_result<normal_cutoff] <- 0
            synthetic_comp = as.matrix(groud_truth[,i])
            colnames(deconv_result) <- 'deconv_result'
            colnames(synthetic_comp) <- 'synthetic_comp'
            temp <- True_False_rate(deconv_result,synthetic_comp)
            False_rate[i] <- temp[3]+temp[4]
          }
          return(False_rate)
        })
        FR_ct_all <- matrix(unlist(FR_ct) ,nr=20) %>% as.data.frame()
        
        rownames(FR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
        colnames(FR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
        
        p_FR_ct_all <- pivot_longer(as.data.frame(t(FR_ct_all)) %>% mutate(algorithm=rownames(.)),
                                    cols = rownames(FR_ct_all),
                                    names_to = 'cell_type',values_to = 'False_Counts')
        p_FR_ct_all_example <- filter(p_FR_ct_all,cell_type %in% c("CD8.Tem","CD8.Tex","CD8.Trm","NK"))
        p_FR_ct_all_example$algorithm <- factor(p_FR_ct_all_example$algorithm,
                                                levels = c('UCASpatial','stereoscope','cell2location','RCTD','CARD','SPOTlight'))
        ggbarplot(p_FR_ct_all_example,x='algorithm',y='False_Counts',
                  fill = 'cell_type',color = 'cell_type')+
          xlab(NULL)+ylab('False Positive + False Negative')+
          theme_bw(base_size = 20)+
          theme(legend.position = 'top',legend.justification = c(1,0))+labs(fill = "Cell States",color= "Cell States")+
          theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values = RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],guide=guide_legend(ncol = 2))+
          scale_color_manual(values = alpha(RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],alpha = 0))+
          scale_y_continuous(position = "right")
        ggsave('PCC_example_FalseCounts_2_barplot.pdf',height = 8,width = 4)
        
        TR_ct <- lapply(highres_list_TFR, function(x){
          True_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
          temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
          for(i in 1:20)
          {
            deconv_result <- as.matrix(x[,i])
            deconv_result[deconv_result<normal_cutoff] <- 0
            synthetic_comp = as.matrix(groud_truth[,i])
            colnames(deconv_result) <- 'deconv_result'
            colnames(synthetic_comp) <- 'synthetic_comp'
            temp <- True_False_rate(deconv_result,synthetic_comp)
            True_rate[i] <- temp[1]+temp[2]
          }
          return(True_rate)
        })
        TR_ct_all <- matrix(unlist(TR_ct) ,nr=20) %>% as.data.frame()
        rownames(TR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
        colnames(TR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      }
    }
    # Fig. 2e summary overall
    {
      overall_summary <- df_Pearson_l_highres_avg[1:6,] %>% 
        mutate(algorithm = rownames(.),
               ext_low=rowSums(.<0.25),
               low=rowSums(.<0.5 & .>=0.25),
               high=rowSums(.<0.75 & .>=0.5),
               ext_high=rowSums(.>=0.75)) %>% select(ext_low,low,high,ext_high,algorithm)
      overall_summary[,1:4] <- overall_summary[,1:4]/20
      p_overall_summary <- pivot_longer(overall_summary,cols = colnames(overall_summary)[1:4],
                                        names_to = 'confidence_group',values_to = 'num')
      p_overall_summary$confidence_group <- factor(p_overall_summary$confidence_group,
                                                   levels = c('ext_high','high','low','ext_low'))
      p_overall_summary$algorithm <- factor(p_overall_summary$algorithm,
                                            levels = rev(c('UCASpatial','stereoscope','SPOTlight','cell2location',
                                                           'RCTD','CARD')))
      
      
      ggbarplot(p_overall_summary,x='algorithm',y='num',orientation = "horiz",color = 'confidence_group',
                fill = 'confidence_group')+xlab(NULL)+ylab(NULL)+
        theme_bw()+
        theme(legend.position = 'bottom',legend.justification = c(1.1,0))+
        scale_fill_manual(values = RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)])+
        scale_color_manual(values = rep(alpha('white',alpha = 0),4))
      ggsave('PCC_overall_summary_barplot.pdf',height = 4,width = 4.6)
    }

    
    
    
  }
  ## Fig SI. RMSE line chart
  {
    df_RMSE_l_highres <- RMSE_l_highres %>% na.omit() %>% dplyr::group_by(algorithm,index) %>% 
      dplyr::mutate(Mean = mean(value)) %>% dplyr::select(algorithm,index,Mean) %>% unique()
    df_RMSE_l_highres <- pivot_wider(df_RMSE_l_highres,names_from = 'index',values_from = 'Mean') %>% as.data.frame()
    
    rownames(df_RMSE_l_highres) <- as.character(df_RMSE_l_highres$algorithm)
    df_RMSE_l_highres <- df_RMSE_l_highres[,-1] 
    
    df_RMSE_l_highres[is.na(df_RMSE_l_highres)] <- 0
    # remove Tangram
    # df_RMSE_l_highres_avg <- df_RMSE_l_highres[-5,] %>% t() %>% as.data.frame()
    df_RMSE_l_highres_avg <- df_RMSE_l_highres %>% t() %>% as.data.frame()
    # rownames(df_RMSE_l_highres_avg) <- colnames(df_RMSE_l_highres)
    df_RMSE_l_highres_avg <- df_RMSE_l_highres_avg %>% 
      mutate(avg = rowMeans(.),minm= Rfast::rowMins(as.matrix(df_RMSE_l_highres_avg[,-1]),value=T)) %>%
      dplyr::arrange(minm) %>% t() %>% as.data.frame()
    
    df_RMSE_l_highres_normalized <- t(t(df_RMSE_l_highres_avg) - as.data.frame(t(df_RMSE_l_highres_avg))$UCASpatial) %>%
      as.data.frame()
    
    
    df_RMSE_l_highres_normalized <- df_RMSE_l_highres_normalized %>% 
      dplyr::select(c("Stromal.c1","Stromal.c2","Epithelial.c1","Epithelial.c2",
                      "CD8.Trm","CD8.Tem",'CD8.Tex',"NK","CD4.T.CCR7","CD4.Tfh","CD4.Treg",
                      "B.GC","B.memory.like","B.naive.like","Plasma",
                      "Macrophage","Monocyte","DC","pDC","Mast"))
    
    
    # example overall
    {
      example_overall <- df_RMSE_l_highres_avg %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = colnames(.)[1:20],names_to = 'Celltype',values_to = 'RMSE')
      # example_overall <- example_overall %>% filter(Algorithm%in%c('UCASpatial','minm','avg'))
      example_overall$Algorithm <- as.factor(example_overall$Algorithm)
      example_overall$Algorithm <- factor(example_overall$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                               'cell2location','CARD','stereoscope',
                                                                               'avg','minm'))
      example_overall$Celltype <- as.factor(example_overall$Celltype)
      example_overall$Celltype <- factor(example_overall$Celltype,levels = colnames(df_RMSE_l_highres_avg))
      
      
      
      
      p0 <- ggplot(example_overall,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                       linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.5,alpha=0.5)+ ylim(0,0.25)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,avg=2,CARD=3,stereoscope=3,dif = 2,minm = 2))+
        scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,avg=1,CARD=0.5,stereoscope=0.5,dif=1,minm = 1))+
        scale_color_manual(values = ggradar_col[c(1:4,6:7,16,9)])+
        theme_bw()+xlab(NULL)+ylab('RMSE')+
        # ggtitle(label = 'RMSE for')+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")
      # ggsave('RMSE_example_CD8T.pdf',height = 2.6,width = 4)
      p0
      ggsave('RMSE_overall_lineplot.pdf',height = 4,width = 5)
      
      
      
      example_CD8T <- df_RMSE_l_highres_normalized %>% dplyr::select("CD8.Trm","CD8.Tem",'CD8.Tex',"NK") %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = c("CD8.Trm","CD8.Tem",'CD8.Tex',"NK"),names_to = 'Celltype',values_to = 'RMSE')
      example_CD8T$Algorithm <- as.factor(example_CD8T$Algorithm)
      example_CD8T$Algorithm <- factor(example_CD8T$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                         'cell2location','Tangram','CARD','stereoscope'))
      
      p1 <- ggplot(example_CD8T,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                    linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.8)+ #ylim(0,0.8)+
        geom_line()+
        # scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
        #                                  cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
        scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
        scale_color_manual(values = ggradar_col)+
        theme_bw()+xlab(NULL)+ylab('Normalized RMSE')+
        # ggtitle(label = 'Gradient proportion')+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")
      # ggsave('RMSE_normalized_example_CD8T.pdf',height = 2.6,width = 4)
      
      }
    # example CD8T/NK
    {
      
      example_stro_mye <- df_RMSE_l_highres %>% dplyr::select("CD8.Tem","CD8.Tex","CD8.Trm","NK") %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = c("CD8.Tem","CD8.Tex","CD8.Trm","NK"),names_to = 'Celltype',values_to = 'RMSE')
      example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
      example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                 'cell2location','CARD','stereoscope'))
      
      p3 <- ggplot(example_stro_mye,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                        linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.8)+ ylim(0,0.1)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
        scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
        scale_color_manual(values = ggradar_col[-5])+
        theme_bw()+xlab(NULL)+ylab('RMSE')+labs(title = 'CD8T/NK')+
        # ggtitle(label = 'RMSE for')+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+NoLegend()
      p3
      ggsave('RMSE_example_2_CD8T.pdf',height = 4,width = 2)
      # ggsave('RMSE_example_CD8T.pdf',height = 2,width = 4)
      
      
      example_stro_mye <- df_RMSE_l_highres[-5,] %>% dplyr::select("Macrophage","Monocyte","DC","pDC") %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = c("Macrophage","Monocyte","DC","pDC"),names_to = 'Celltype',values_to = 'RMSE')
      example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
      example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                 'cell2location','CARD','stereoscope'))
      example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c("Macrophage","Monocyte","DC","pDC"))  
      
      p4 <- ggplot(example_stro_mye,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                        linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.8)+ ylim(0,0.1)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
        scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
        scale_color_manual(values = ggradar_col)+
        theme_bw()+xlab(NULL)+ylab('RMSE')+labs(title = 'Myeloid')+
        # ggtitle(label = 'Gradient proportion')+
        theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+NoLegend()
      # ggsave('RMSE_example_mye.pdf',height = 2,width = 4)
      
      
      example_stro_mye <- df_RMSE_l_highres %>% dplyr::select('CD4.Tfh','CD4.T.CCR7','CD4.Treg') %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = c('CD4.Tfh','CD4.T.CCR7','CD4.Treg'),names_to = 'Celltype',values_to = 'RMSE')
      example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
      example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                 'cell2location','CARD','stereoscope'))
      example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c('CD4.Tfh','CD4.T.CCR7','CD4.Treg'))  
      
      p5 <- ggplot(example_stro_mye,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                        linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.8)+ ylim(0,0.15)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
        scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
        scale_color_manual(values = ggradar_col)+
        theme_bw()+xlab(NULL)+ylab('RMSE')+labs(title = 'CD4T')+
        # ggtitle(label = 'Gradient proportion')+
        theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+NoLegend()
      # ggsave('RMSE_example_CD4T.pdf',height = 2,width = 4)
      
      
      
      example_stro_mye <- df_RMSE_l_highres %>% dplyr::select('B.GC','B.naive.like','B.memory.like','Plasma') %>% 
        mutate(Algorithm = rownames(.)) %>% 
        pivot_longer(cols = c('B.GC','B.naive.like','B.memory.like','Plasma'),names_to = 'Celltype',values_to = 'RMSE')
      example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
      example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                 'cell2location','CARD','stereoscope'))
      example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c('B.GC','B.naive.like','B.memory.like','Plasma'))  
      
      p5 <- ggplot(example_stro_mye,aes(Celltype,RMSE,group=Algorithm,color=Algorithm,
                                        linetype=Algorithm,linewidth=Algorithm))+
        geom_point(size=1.8)+ ylim(0,0.1)+
        geom_line()+
        scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                         cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
        scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                          cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
        scale_color_manual(values = ggradar_col)+
        theme_bw()+xlab(NULL)+ylab('RMSE')+labs(title = 'B')+
        # ggtitle(label = 'Gradient proportion')+
        theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "right")+NoLegend()
      # ggsave('RMSE_example_CD4T.pdf',height = 2,width = 4)
      
      
      
      
      
      
      
      
      p3+NoLegend()+p4
      ggsave('RMSE_example_2in1.pdf',height = 4,width = 5)
    }
    # False rate of CD8T/NK
    {
      source('/data/xy/scripts/Summary_True_False_rate.R')
      highres_list_TFR <- list(UCASpatial_highres_all,SPOTlight_highres_all,RCTD_highres_all,C2L_highres_all,
                               CARD_highres_all,stereoscope_highres_all)
      highres_list_TFR <- lapply(highres_list_TFR, function(x){
        a <- select(x,colnames(UCASpatial_highres_all)[1:20])
        return(a)
      })
      groud_truth <- select(groud_truth,colnames(UCASpatial_highres_all)[1:20])
      FR_ct <- lapply(highres_list_TFR, function(x){
        False_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
        temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
        for(i in 1:20)
        {
          deconv_result <- as.matrix(x[,i])
          deconv_result[deconv_result<normal_cutoff] <- 0
          synthetic_comp = as.matrix(groud_truth[,i])
          colnames(deconv_result) <- 'deconv_result'
          colnames(synthetic_comp) <- 'synthetic_comp'
          temp <- True_False_rate(deconv_result,synthetic_comp)
          False_rate[i] <- temp[3]+temp[4]
        }
        return(False_rate)
      })
      FR_ct_all <- matrix(unlist(FR_ct) ,nr=20) %>% as.data.frame()
      
      rownames(FR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
      colnames(FR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      
      p_FR_ct_all <- pivot_longer(as.data.frame(t(FR_ct_all)) %>% mutate(algorithm=rownames(.)),
                                  cols = rownames(FR_ct_all),
                                  names_to = 'cell_type',values_to = 'False_Counts')
      p_FR_ct_all_example <- filter(p_FR_ct_all,cell_type %in% c("CD8.Tem","CD8.Tex","CD8.Trm","NK"))
      p_FR_ct_all_example$algorithm <- factor(p_FR_ct_all_example$algorithm,
                                              levels = c('UCASpatial','stereoscope','cell2location','RCTD','CARD','SPOTlight'))
      ggbarplot(p_FR_ct_all_example,x='algorithm',y='False_Counts',
                fill = 'cell_type',color = 'cell_type')+
        xlab(NULL)+ylab('False Positive + False Negative')+
        theme_bw(base_size = 20)+
        theme(legend.position = 'top',legend.justification = c(1,0))+labs(fill = "Cell States",color= "Cell States")+
        theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values = RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],guide=guide_legend(ncol = 2))+
        scale_color_manual(values = alpha(RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],alpha = 0))+
        scale_y_continuous(position = "right")
      ggsave('PCC_example_FalseCounts_2_barplot.pdf',height = 8,width = 4)
      
      TR_ct <- lapply(highres_list_TFR, function(x){
        True_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
        temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
        for(i in 1:20)
        {
          deconv_result <- as.matrix(x[,i])
          deconv_result[deconv_result<normal_cutoff] <- 0
          synthetic_comp = as.matrix(groud_truth[,i])
          colnames(deconv_result) <- 'deconv_result'
          colnames(synthetic_comp) <- 'synthetic_comp'
          temp <- True_False_rate(deconv_result,synthetic_comp)
          True_rate[i] <- temp[1]+temp[2]
        }
        return(True_rate)
      })
      TR_ct_all <- matrix(unlist(TR_ct) ,nr=20) %>% as.data.frame()
      rownames(TR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
      colnames(TR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      
      
      
    }
    
  }
  ## Fig SI. F1 line chart
  {
    # Calculate F1 score celltype level
    {
      ## Accuracy level by cell types
      {
        source('/data/xy/scripts/Test_Accuracy.R')
        # summary Acur,F1,... 
        {
          cell_type_list <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
          All_auc_ct <- lapply(highres_list, function(x){
            auc_celltype <- lapply(cell_type_list, function(z){
              groud_truth_split_10 <- dplyr::filter(groud_truth_split,grepl(z, new_ct)) %>% arrange(new_ct)
              split_10 <-
                dplyr::filter(x, x$new_ct %in% groud_truth_split_10$new_ct) %>%
                dplyr::filter(grepl(z, new_ct)) %>% arrange(new_ct)
              Predict <- as.matrix(split_10$proportion_filter)
              Truth <- as.matrix(groud_truth_split_10$proportion)
              colnames(Predict) <- 'Value'
              colnames(Truth) <- 'Value'
              return(Test_Acurracy(Predict, Truth)[5])
            })
            return(auc_celltype)
          })
          
          auc_ct_md <- lapply(All_auc_ct,function(x){
            test_10 <- matrix(unlist(x) ,nr=1) %>% as.data.frame()
            return(test_10)
          })
          
          auc_ct_UCASpatial <- auc_ct_md[[1]]
          colnames(auc_ct_UCASpatial) <- cell_type_list
          
          auc_ct_SPOTlight <- auc_ct_md[[2]]
          colnames(auc_ct_SPOTlight) <- cell_type_list
          
          auc_ct_RCTD <- auc_ct_md[[3]]
          colnames(auc_ct_RCTD) <- cell_type_list
          
          auc_ct_C2L <- auc_ct_md[[4]]
          colnames(auc_ct_C2L) <- cell_type_list
          
          auc_ct_CARD <- auc_ct_md[[5]]
          colnames(auc_ct_CARD) <- cell_type_list
          
          auc_ct_stereoscope <- auc_ct_md[[6]]
          colnames(auc_ct_stereoscope) <- cell_type_list
          
          
          auc_ct_all <- rbind(auc_ct_UCASpatial,auc_ct_SPOTlight,auc_ct_RCTD,
                              auc_ct_C2L,auc_ct_CARD,auc_ct_stereoscope)
          rownames(auc_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
          # saveRDS(auc_ct_all,'auc_ct_all.rds')
          auc_ct_max <- matrix(unlist(auc_ct_md) ,nc=6) %>% as.data.frame()
          
          auc_ct_avg <- lapply(list(auc_ct_UCASpatial,auc_ct_SPOTlight,auc_ct_RCTD,
                                    auc_ct_C2L,auc_ct_CARD,auc_ct_stereoscope), function(x){
                                      avg <- rowSums(x,na.rm = T)/20
                                      return(avg)
                                    })
          auc_ct_avg2 <- matrix(unlist(auc_ct_avg) ,nc=6) %>% as.data.frame()
          
          # saveRDS(auc_ct_avg2,'auc_ct_avg.rds')
        }
        
        ## F1 line chart
        {
          auc_ct_all[is.na(auc_ct_all)] <- 0
          auc_ct_all_avg <- auc_ct_all %>% t() %>% as.data.frame()
          
          auc_ct_all_avg <- auc_ct_all_avg %>% 
            mutate(avg = rowMeans(.),maxm= Rfast::rowMaxs(as.matrix(auc_ct_all_avg[,-1]),value=T)) %>%
            dplyr::arrange(maxm) %>% t() %>% as.data.frame()
          
          
          auc_ct_all_normalized <- t(t(auc_ct_all_avg) - as.data.frame(t(auc_ct_all_avg))$UCASpatial) %>%
            as.data.frame()
          
          auc_ct_all_normalized <- auc_ct_all_normalized %>% 
            dplyr::select(c("Stromal.c1","Stromal.c2","Epithelial.c1","Epithelial.c2",
                            "CD8.Trm","CD8.Tem",'CD8.Tex',"NK","CD4.T.CCR7","CD4.Tfh","CD4.Treg",
                            "B.GC","B.memory.like","B.naive.like","Plasma",
                            "Macrophage","Monocyte","DC","pDC","Mast"))
          
          
          # example overall
          {
            example_overall <- auc_ct_all_avg %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = colnames(.)[1:20],names_to = 'Celltype',values_to = 'F1')
            # example_overall <- example_overall %>% filter(Algorithm%in%c('UCASpatial','minm','avg'))
            example_overall$Algorithm <- as.factor(example_overall$Algorithm)
            example_overall$Algorithm <- factor(example_overall$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                     'cell2location','CARD','stereoscope',
                                                                                     'avg','maxm'))
            example_overall$Celltype <- as.factor(example_overall$Celltype)
            example_overall$Celltype <- factor(example_overall$Celltype,levels = colnames(auc_ct_all_avg))
            
            
            
            
            p0 <- ggplot(example_overall,aes(Celltype,F1,group=Algorithm,color=Algorithm,
                                             linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.5,alpha=0.5)+ ylim(0,1)+
              geom_line()+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,CARD=3,stereoscope=3,avg=2,maxm = 2))+
              scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,CARD=0.5,stereoscope=0.5,avg=1,maxm = 1))+
              scale_color_manual(values = ggradar_col[c(1:4,6:7,16,9)])+
              theme_bw()+xlab(NULL)+ylab('F1 score')+
              # ggtitle(label = 'F1 for')+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")
            # ggsave('F1_example_CD8T.pdf',height = 2.6,width = 4)
            p0
            ggsave('F1_overall_lineplot.pdf',height = 4,width = 5)
            
            
            
            example_CD8T <- df_Pearson_l_highres_normalized %>% dplyr::select("CD8.Trm","CD8.Tem",'CD8.Tex',"NK") %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = c("CD8.Trm","CD8.Tem",'CD8.Tex',"NK"),names_to = 'Celltype',values_to = 'F1')
            example_CD8T$Algorithm <- as.factor(example_CD8T$Algorithm)
            example_CD8T$Algorithm <- factor(example_CD8T$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                               'cell2location','Tangram','CARD','stereoscope'))
            
            p1 <- ggplot(example_CD8T,aes(Celltype,F1,group=Algorithm,color=Algorithm,
                                          linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.8)+ #ylim(0,0.8)+
              geom_line()+
              # scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
              #                                  cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
              scale_color_manual(values = ggradar_col)+
              theme_bw()+xlab(NULL)+ylab('Normalized F1')+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")
            # ggsave('F1_normalized_example_CD8T.pdf',height = 2.6,width = 4)
            
            }
          # summary overall
          {
            overall_summary <- df_Pearson_l_highres_avg[1:6,] %>% 
              mutate(algorithm = rownames(.),
                     ext_low=rowSums(.<0.25),
                     low=rowSums(.<0.5 & .>=0.25),
                     high=rowSums(.<0.75 & .>=0.5),
                     ext_high=rowSums(.>=0.75)) %>% select(ext_low,low,high,ext_high,algorithm)
            overall_summary[,1:4] <- overall_summary[,1:4]/20
            p_overall_summary <- pivot_longer(overall_summary,cols = colnames(overall_summary)[1:4],
                                              names_to = 'confidence_group',values_to = 'num')
            p_overall_summary$confidence_group <- factor(p_overall_summary$confidence_group,
                                                         levels = c('ext_high','high','low','ext_low'))
            p_overall_summary$algorithm <- factor(p_overall_summary$algorithm,
                                                  levels = rev(c('UCASpatial','stereoscope','SPOTlight','cell2location',
                                                                 'RCTD','CARD')))
            
            
            ggbarplot(p_overall_summary,x='algorithm',y='num',orientation = "horiz",color = 'confidence_group',
                      fill = 'confidence_group')+xlab(NULL)+ylab(NULL)+
              theme_bw()+
              theme(legend.position = 'bottom',legend.justification = c(1.1,0))+
              scale_fill_manual(values = RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)])+
              scale_color_manual(values = rep(alpha('white',alpha = 0),4))
            ggsave('F1_overall_summary_barplot.pdf',height = 4,width = 4.6)
            
            ggradar_col
            
            
            1
          }
          # example CD8T/NK
          {
            
            example_stro_mye <- df_Pearson_l_highres[-5,] %>% dplyr::select("CD8.Tem","CD8.Tex","CD8.Trm","NK") %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = c("CD8.Tem","CD8.Tex","CD8.Trm","NK"),names_to = 'Celltype',values_to = 'F1')
            example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
            example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                       'cell2location','CARD','stereoscope'))
            
            p3 <- ggplot(example_stro_mye,aes(Celltype,F1,group=Algorithm,color=Algorithm,
                                              linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.8)+ ylim(0,0.65)+
              geom_line()+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
              scale_color_manual(values = ggradar_col[-5])+
              theme_bw()+xlab(NULL)+ylab('F1')+labs(title = 'CD8T/NK')+
              # ggtitle(label = 'F1 for')+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")+NoLegend()
            p3
            ggsave('F1_example_2_CD8T.pdf',height = 4,width = 2)
            # ggsave('F1_example_CD8T.pdf',height = 2,width = 4)
            
            
            example_stro_mye <- df_Pearson_l_highres[-5,] %>% dplyr::select("Macrophage","Monocyte","DC","pDC") %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = c("Macrophage","Monocyte","DC","pDC"),names_to = 'Celltype',values_to = 'F1')
            example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
            example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                       'cell2location','CARD','stereoscope'))
            example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c("Macrophage","Monocyte","DC","pDC"))  
            
            p4 <- ggplot(example_stro_mye,aes(Celltype,F1,group=Algorithm,color=Algorithm,
                                              linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.8)+ ylim(0,0.9)+
              geom_line()+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
              scale_color_manual(values = ggradar_col)+
              theme_bw()+xlab(NULL)+ylab('F1')+labs(title = 'Myeloid')+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")+NoLegend()
            # ggsave('PCC_example_mye.pdf',height = 2,width = 4)
            
            
            example_stro_mye <- df_Pearson_l_highres[-5,] %>% dplyr::select('CD4.Tfh','CD4.T.CCR7','CD4.Treg') %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = c('CD4.Tfh','CD4.T.CCR7','CD4.Treg'),names_to = 'Celltype',values_to = 'PCC')
            example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
            example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                       'cell2location','CARD','stereoscope'))
            example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c('CD4.Tfh','CD4.T.CCR7','CD4.Treg'))  
            
            p5 <- ggplot(example_stro_mye,aes(Celltype,PCC,group=Algorithm,color=Algorithm,
                                              linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.8)+ ylim(0,0.65)+
              geom_line()+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
              scale_color_manual(values = ggradar_col)+
              theme_bw()+xlab(NULL)+ylab('PCC')+labs(title = 'CD4T')+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")+NoLegend()
            # ggsave('PCC_example_CD4T.pdf',height = 2,width = 4)
            
            
            
            example_stro_mye <- df_Pearson_l_highres[-5,] %>% dplyr::select('B.GC','B.naive.like','B.memory.like','Plasma') %>% 
              mutate(Algorithm = rownames(.)) %>% 
              pivot_longer(cols = c('B.GC','B.naive.like','B.memory.like','Plasma'),names_to = 'Celltype',values_to = 'PCC')
            example_stro_mye$Algorithm <- as.factor(example_stro_mye$Algorithm)
            example_stro_mye$Algorithm <- factor(example_stro_mye$Algorithm,levels = c('UCASpatial','SPOTlight','RCTD',
                                                                                       'cell2location','CARD','stereoscope'))
            example_stro_mye$Celltype <- factor(example_stro_mye$Celltype,levels = c('B.GC','B.naive.like','B.memory.like','Plasma'))  
            
            p5 <- ggplot(example_stro_mye,aes(Celltype,PCC,group=Algorithm,color=Algorithm,
                                              linetype=Algorithm,linewidth=Algorithm))+
              geom_point(size=1.8)+ ylim(0,1)+
              geom_line()+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,Tangram=1,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 0.8, SPOTlight = 0.5, RCTD=0.5,
                                                cell2location = 0.5,Tangram=0.5,CARD=0.5,stereoscope=0.5))+
              scale_color_manual(values = ggradar_col)+
              theme_bw()+xlab(NULL)+ylab('PCC')+labs(title = 'B')+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 0,hjust = 0.5),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")+NoLegend()
            # ggsave('PCC_example_CD4T.pdf',height = 2,width = 4)
            
            
            
            
            
            
            
            
            p3+NoLegend()+p4
            ggsave('PCC_example_2in1.pdf',height = 4,width = 5)
          }
          # False rate of CD8T/NK
          {
            source('/data/xy/scripts/Summary_True_False_rate.R')
            highres_list_TFR <- list(UCASpatial_highres_all,SPOTlight_highres_all,RCTD_highres_all,C2L_highres_all,
                                     CARD_highres_all,stereoscope_highres_all)
            highres_list_TFR <- lapply(highres_list_TFR, function(x){
              a <- select(x,colnames(UCASpatial_highres_all)[1:20])
              return(a)
            })
            groud_truth <- select(groud_truth,colnames(UCASpatial_highres_all)[1:20])
            FR_ct <- lapply(highres_list_TFR, function(x){
              False_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
              temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
              for(i in 1:20)
              {
                deconv_result <- as.matrix(x[,i])
                deconv_result[deconv_result<normal_cutoff] <- 0
                synthetic_comp = as.matrix(groud_truth[,i])
                colnames(deconv_result) <- 'deconv_result'
                colnames(synthetic_comp) <- 'synthetic_comp'
                temp <- True_False_rate(deconv_result,synthetic_comp)
                False_rate[i] <- temp[3]+temp[4]
              }
              return(False_rate)
            })
            FR_ct_all <- matrix(unlist(FR_ct) ,nr=20) %>% as.data.frame()
            
            rownames(FR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
            colnames(FR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
            
            p_FR_ct_all <- pivot_longer(as.data.frame(t(FR_ct_all)) %>% mutate(algorithm=rownames(.)),
                                        cols = rownames(FR_ct_all),
                                        names_to = 'cell_type',values_to = 'False_Counts')
            p_FR_ct_all_example <- filter(p_FR_ct_all,cell_type %in% c("CD8.Tem","CD8.Tex","CD8.Trm","NK"))
            p_FR_ct_all_example$algorithm <- factor(p_FR_ct_all_example$algorithm,
                                                    levels = c('UCASpatial','stereoscope','cell2location','RCTD','CARD','SPOTlight'))
            ggbarplot(p_FR_ct_all_example,x='algorithm',y='False_Counts',
                      fill = 'cell_type',color = 'cell_type')+
              xlab(NULL)+ylab('False Positive + False Negative')+
              theme_bw(base_size = 20)+
              theme(legend.position = 'top',legend.justification = c(1,0))+labs(fill = "Cell States",color= "Cell States")+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              scale_fill_manual(values = RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],guide=guide_legend(ncol = 2))+
              scale_color_manual(values = alpha(RColorBrewer::brewer.pal(n=6,name = 'Paired')[c(6,5,1,2)],alpha = 0))+
              scale_y_continuous(position = "right")
            ggsave('PCC_example_FalseCounts_2_barplot.pdf',height = 8,width = 4)
            
            TR_ct <- lapply(highres_list_TFR, function(x){
              True_rate <- as.data.frame(matrix(nrow = 1,ncol = 20))
              temp <- as.data.frame(matrix(nrow = 1,ncol = 4))
              for(i in 1:20)
              {
                deconv_result <- as.matrix(x[,i])
                deconv_result[deconv_result<normal_cutoff] <- 0
                synthetic_comp = as.matrix(groud_truth[,i])
                colnames(deconv_result) <- 'deconv_result'
                colnames(synthetic_comp) <- 'synthetic_comp'
                temp <- True_False_rate(deconv_result,synthetic_comp)
                True_rate[i] <- temp[1]+temp[2]
              }
              return(True_rate)
            })
            TR_ct_all <- matrix(unlist(TR_ct) ,nr=20) %>% as.data.frame()
            rownames(TR_ct_all) <- colnames(UCASpatial_highres_all)[1:20]
            colnames(TR_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
            
            
            
          }
          
          
          
        }
        
        
      }
      
    }
  }
}

# Fig. 2g gradient proportion by cell type
{
  # UCASpatial results
  {
    UCASpatial_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(UCASpatial_highres_all) <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
    for(i in 1:15)
    {
      UCASpatial_highres_all <- rbind(UCASpatial_highres_all,UCASpatial_highres[[i]][[2]][,1:20])
    }
    UCASpatial_highres_split <- UCASpatial_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(UCASpatial_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    UCASpatial_highres_split$proportion_filter <- UCASpatial_highres_split$proportion
    UCASpatial_highres_split$proportion_filter[UCASpatial_highres_split$proportion_filter<0.05] <- 0
  }
  # SPOTlight results
  {
    SPOTlight_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(SPOTlight_highres_all) <- colnames(SPOTlight_highres[[1]][[2]])[1:20]
    for(i in 1:15)
    {
      SPOTlight_highres_all <- rbind(SPOTlight_highres_all,SPOTlight_highres[[i]][[2]][,1:20])
    }
    SPOTlight_highres_split <- SPOTlight_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(SPOTlight_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    SPOTlight_highres_split$proportion_filter <- SPOTlight_highres_split$proportion
    SPOTlight_highres_split$proportion_filter[SPOTlight_highres_split$proportion_filter<0.05] <- 0
  }
  # RCTD results
  {
    RCTD_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(RCTD_highres_all) <- colnames(RCTD_highres[[1]]@results[["weights"]])[1:20]
    for(i in 1:15)
    {
      RCTD_highres_all <- rbind(RCTD_highres_all,as.data.frame(RCTD_highres[[i]]@results[["weights"]]))
    }
    rownames(RCTD_highres_all) <- 1:750
    colnames(RCTD_highres_all) <- gsub("[_/ /+/-/(/)/-]",".",colnames(RCTD_highres_all))
    RCTD_highres_split <- RCTD_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(RCTD_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    RCTD_highres_split$proportion_filter <- RCTD_highres_split$proportion
    RCTD_highres_split$proportion_filter[RCTD_highres_split$proportion_filter<0.05] <- 0
  }
  # cell2location results
  {
    C2L_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(C2L_highres_all) <- colnames(C2L_highres[[1]][[2]])[1:20]
    for(i in 1:15)
    {
      C2L_highres_all <- rbind(C2L_highres_all,C2L_highres[[i]])
    }
    rownames(C2L_highres_all) <- 1:750
    colnames(C2L_highres_all) <- gsub("[_/ /+/-/(/)/-]",".",colnames(C2L_highres_all))
    C2L_highres_all[C2L_highres_all<0.1] <- 0
    C2L_highres_all <- as.data.frame(C2L_highres_all)/rowSums(as.matrix(C2L_highres_all))
    C2L_highres_split <- C2L_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(C2L_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    C2L_highres_split$proportion_filter <- C2L_highres_split$proportion
    C2L_highres_split$proportion_filter[C2L_highres_split$proportion_filter<0.1] <- 0
  }
  # CARD results
  {
    CARD_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(CARD_highres_all) <- colnames(CARD_highres[[1]])[1:20]
    for(i in 1:15)
    {
      CARD_highres_all <- rbind(CARD_highres_all,CARD_highres[[i]])
    }
    rownames(CARD_highres_all) <- 1:750
    colnames(CARD_highres_all) <- gsub("[_/ /+/-/(/)/-]",".",colnames(CARD_highres_all))
    CARD_highres_split <- CARD_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(CARD_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    CARD_highres_split$proportion_filter <- CARD_highres_split$proportion
    CARD_highres_split$proportion_filter[CARD_highres_split$proportion_filter<0.05] <- 0
  }
  # stereoscope results
  {
    stereoscope_highres_all <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(stereoscope_highres_all) <- colnames(stereoscope_highres[[1]])[1:20]
    for(i in 1:15)
    {
      stereoscope_highres_all <- rbind(stereoscope_highres_all,stereoscope_highres[[i]])
    }
    rownames(stereoscope_highres_all) <- 1:750
    colnames(stereoscope_highres_all) <- gsub("[_/ /+/-/(/)/-]",".",colnames(stereoscope_highres_all))
    stereoscope_highres_split <- stereoscope_highres_all %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(stereoscope_highres_all),
                   values_to = 'proportion',names_to = 'Celltype') %>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
    stereoscope_highres_split$proportion_filter <- stereoscope_highres_split$proportion
    stereoscope_highres_split$proportion_filter[stereoscope_highres_split$proportion_filter<0.05] <- 0
  }
  # groud truth results
  {
    groud_truth <- as.data.frame(matrix(nrow = 0,ncol=20))
    colnames(groud_truth) <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
    for(i in 1:15)
    {
      synthetic_comp <- as.matrix(simulatData_TME[[i]][["cell_composition"]][[3]] /
                                    rowSums(simulatData_TME[[i]][["cell_composition"]][[3]]))
      colnames(synthetic_comp) <- gsub('_','.',colnames(synthetic_comp))
      groud_truth <- rbind(groud_truth,synthetic_comp)
    }
    groud_truth_split <- groud_truth %>% mutate(spotid = rownames(.)) %>% 
      pivot_longer(cols = colnames(groud_truth),
                   values_to = 'proportion',names_to = 'Celltype')%>% 
      mutate(new_ct = paste(spotid,Celltype,sep = '_')) %>% 
      dplyr::select(new_ct,proportion) %>% arrange(new_ct)
  }
  highres_list <- list(UCASpatial_highres_split,SPOTlight_highres_split,RCTD_highres_split,C2L_highres_split,
                       CARD_highres_split,stereoscope_highres_split)
  
  # Fig. 2g RMSE
  {
    ## RMSE by cell types
    {
      cell_type_list <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
      
      RMSE_divided10_ct <- lapply(highres_list, function(x){
        RMSE_celltype <- lapply(cell_type_list, function(z){
          RMSE <- lapply(seq(0.1,1,0.1), function(y){
            groud_truth_split_10 <- filter(groud_truth_split,proportion <= y) %>% 
              filter(grepl(z,new_ct)) %>% arrange(new_ct)
            split_10 <- filter(x,x$new_ct %in% groud_truth_split_10$new_ct) %>% 
              filter(grepl(z,new_ct)) %>% arrange(new_ct)
            return(rmse(split_10$proportion,groud_truth_split_10$proportion))
          })
          return(RMSE)
        })
        return(RMSE_celltype)
      })
      
      RMSE_divided10_ct_md <- lapply(RMSE_divided10_ct,function(x){
        test_10 <- matrix(unlist(x) ,nr=10) %>% as.data.frame()
        return(test_10)
      })
      
      RMSE_10_ct_UCASpatial <- RMSE_divided10_ct_md[[1]]
      colnames(RMSE_10_ct_UCASpatial) <- cell_type_list
      rownames(RMSE_10_ct_UCASpatial) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_SPOTlight <- RMSE_divided10_ct_md[[2]]
      colnames(RMSE_10_ct_SPOTlight) <- cell_type_list
      rownames(RMSE_10_ct_SPOTlight) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_RCTD <- RMSE_divided10_ct_md[[3]]
      colnames(RMSE_10_ct_RCTD) <- cell_type_list
      rownames(RMSE_10_ct_RCTD) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_C2L <- RMSE_divided10_ct_md[[4]]
      colnames(RMSE_10_ct_C2L) <- cell_type_list
      rownames(RMSE_10_ct_C2L) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_CARD <- RMSE_divided10_ct_md[[5]]
      colnames(RMSE_10_ct_CARD) <- cell_type_list
      rownames(RMSE_10_ct_CARD) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_stereoscope <- RMSE_divided10_ct_md[[6]]
      colnames(RMSE_10_ct_stereoscope) <- cell_type_list
      rownames(RMSE_10_ct_stereoscope) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      RMSE_10_ct_max <- matrix(unlist(RMSE_divided10_ct_md) ,nc=6) %>% as.data.frame() %>% 
        filter(rownames(.) %in% c(seq(10,200,10)))
      
      RMSE_10_ct_avg <- lapply(list(RMSE_10_ct_UCASpatial,RMSE_10_ct_SPOTlight,RMSE_10_ct_RCTD,
                                    RMSE_10_ct_C2L,RMSE_10_ct_CARD,RMSE_10_ct_stereoscope), function(x){
                                      avg <- rowSums(x)/20
                                      return(avg)
                                    })
      RMSE_10_ct_avg2 <- matrix(unlist(RMSE_10_ct_avg) ,nc=6) %>% as.data.frame()
      rownames(RMSE_10_ct_avg2) <- paste('0-',seq(10,100,10),'%',sep = '')
      colnames(RMSE_10_ct_avg2) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      # saveRDS(RMSE_10_ct_avg2,'RMSE_10_ct_avg.rds')
    }
    # draw line chart of RMSE for all methods
    {
      library(Rmisc)
      RMSE_10_ct_all <- matrix(unlist(RMSE_divided10_ct_md) ,nc=6) %>% as.data.frame()
      colnames(RMSE_10_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      RMSE_10_ct_all_2 <- RMSE_10_ct_all %>% mutate(proportion = paste('0-',rep(seq(10,100,10),20),'%',sep = ''),
                                                    celltype=rep(colnames(UCASpatial_highres[[1]][[2]])[1:20],each = 10))
      RMSE_10_ct_all_p <- pivot_longer(RMSE_10_ct_all_2,cols = colnames(RMSE_10_ct_all_2)[1:6],names_to = 'algorithm',
                                       values_to = 'RMSE')
      RMSE_10_ct_all_p$proportion <- factor(RMSE_10_ct_all_p$proportion,
                                            levels = paste('0-',seq(10,100,10),'%',sep = ''))
      
      RMSE_10_ct_all_p$algorithm <- factor(RMSE_10_ct_all_p$algorithm,
                                           levels = c('UCASpatial','SPOTlight','RCTD',
                                                      'cell2location','CARD','stereoscope'))
      
      # overall summary
      {
        RMSE_10_ct_all_p_sum <- summarySE(RMSE_10_ct_all_p, measurevar = "RMSE",
                                          groupvars = c("algorithm","proportion"))
        RMSE_10_ct_all_p_sum$proportion <- factor(RMSE_10_ct_all_p_sum$proportion,
                                                  levels = paste('0-',seq(10,100,10),'%',sep = ''))
        
        RMSE_10_ct_all_p_sum$algorithm <- factor(RMSE_10_ct_all_p_sum$algorithm,
                                                 levels = c('UCASpatial','SPOTlight','RCTD',
                                                            'cell2location','CARD','stereoscope'))
        
        ggplot(RMSE_10_ct_all_p_sum,aes(proportion,RMSE,group=algorithm,color=algorithm,
                                        linetype=algorithm,linewidth=algorithm))+
          # geom_errorbar(aes(ymin = RMSE-se, ymax = RMSE+se), width = 0.1, size = 1.2)+
          geom_point(size=1.5)+ expand_limits(y=c(0,0.1))+
          geom_line(position = position_dodge(0.1),cex=0.5)+
          scale_color_manual(values = ggradar_col[-5])+
          theme_bw()+xlab(NULL)+
          scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                           cell2location = 3,CARD=3,stereoscope=3))+
          scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                            cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
          # ggtitle(label = 'Gradient proportion')+
          theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
          theme(legend.position = "right")
        ggsave('Gradient_proportion_avg_ct_line_chart_RMSE.pdf',height = 3,width = 4.5)
        
        RMSE_10_ct_all_p_sum_summ <- RMSE_10_ct_all_p_sum %>% dplyr::group_by(algorithm) %>% 
          dplyr::summarise(Mean =median(RMSE))
        RMSE_10_ct_all_p_sum_summ$Mean_div <- RMSE_10_ct_all_p_sum_summ$Mean-RMSE_10_ct_all_p_sum_summ$Mean[1]
        RMSE_10_ct_all_p_sum_summ$impr <- RMSE_10_ct_all_p_sum_summ$Mean_div/RMSE_10_ct_all_p_sum_summ$Mean*100
      }
      
      # for each cell type
      {
        dev.off()
        pdf('Gradient_proportion_avg_ct_line_chart_RMSE_each_celltype.pdf',height = 3,width = 4.5)
        p_RMSE_each_ct_grad_pro <- lapply(cell_type_list,function(x){
          RMSE_10_ct_all_p_temp <- RMSE_10_ct_all_p %>% filter(celltype == x)
          ggplot(RMSE_10_ct_all_p_temp,aes(proportion,RMSE,group=algorithm,color=algorithm,
                                           linetype=algorithm,linewidth=algorithm))+
            # geom_errorbar(aes(ymin = RMSE-se, ymax = RMSE+se), width = 0.1, size = 1.2)+
            geom_point(size=1.5)+ expand_limits(y=0)+
            geom_line(position = position_dodge(0.1),cex=0.5)+
            scale_color_manual(values = ggradar_col[-5])+
            theme_bw()+xlab(NULL)+labs(title=x)+
            scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                             cell2location = 3,CARD=3,stereoscope=3))+
            scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                              cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
            # ggtitle(label = 'Gradient proportion')+
            theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
            theme(legend.position = "right")
        })
        p_RMSE_each_ct_grad_pro
        dev.off()
        
        
        pdf('Gradient_proportion_avg_ct_line_chart_RMSE_Tex_Tem.pdf',height = 3,width = 4.5)
        p_RMSE_each_ct_grad_pro <- lapply(c('CD8.Tex','CD8.Tem'),function(x){
          RMSE_10_ct_all_p_temp <- RMSE_10_ct_all_p %>% filter(celltype == x)
          ggplot(RMSE_10_ct_all_p_temp,aes(proportion,RMSE,group=algorithm,color=algorithm,
                                           linetype=algorithm,linewidth=algorithm))+
            # geom_errorbar(aes(ymin = RMSE-se, ymax = RMSE+se), width = 0.1, size = 1.2)+
            geom_point(size=1.5)+ expand_limits(y=0)+
            geom_line(position = position_dodge(0.1),cex=0.5)+
            scale_color_manual(values = ggradar_col[-5])+
            theme_bw()+xlab(NULL)+#labs(title=x)+
            scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                             cell2location = 3,CARD=3,stereoscope=3))+
            scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                              cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
            # ggtitle(label = 'Gradient proportion')+
            theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
            theme(legend.position = "right")
        })
        p_RMSE_each_ct_grad_pro
        dev.off()
        
        dev.off()
        pdf('Gradient_proportion_avg_ct_line_chart_RMSE_each_celltype_in1.pdf',height = 15,width = 18)
        p_RMSE_each_ct_grad_pro <- lapply(cell_type_list,function(x){
          RMSE_10_ct_all_p_temp <- RMSE_10_ct_all_p %>% filter(celltype == x)
          ggplot(RMSE_10_ct_all_p_temp,aes(proportion,RMSE,group=algorithm,color=algorithm,
                                           linetype=algorithm,linewidth=algorithm))+
            # geom_errorbar(aes(ymin = RMSE-se, ymax = RMSE+se), width = 0.1, size = 1.2)+
            geom_point(size=1.5)+ expand_limits(y=0)+
            geom_line(position = position_dodge(0.1),cex=0.5)+
            scale_color_manual(values = ggradar_col[-5])+
            theme_bw()+xlab(NULL)+labs(title=x)+
            scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                             cell2location = 3,CARD=3,stereoscope=3))+
            scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                              cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
            # ggtitle(label = 'Gradient proportion')+
            theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
            theme(legend.position = "right")
        })
        cowplot::plot_grid(plotlist = p_RMSE_each_ct_grad_pro,ncol = 4)
        
        dev.off()
        
        RMSE_10_ct_all_p_summ <- RMSE_10_ct_all_p %>% dplyr::group_by(algorithm,celltype) %>% 
          dplyr::summarise(Mean=mean(RMSE))
        RMSE_10_ct_all_p_summ_Tex <- RMSE_10_ct_all_p_summ %>% 
          dplyr::filter(celltype %in% c('CD8.Tex')) %>% dplyr::group_by(algorithm) %>%  
          dplyr::summarise(Mean=mean(Mean,na.rm = T))
        RMSE_10_ct_all_p_summ_Tem <- RMSE_10_ct_all_p_summ %>% 
          dplyr::filter(celltype %in% c('CD8.Tem')) %>% dplyr::group_by(algorithm) %>%  
          dplyr::summarise(Mean=mean(Mean,na.rm = T))
        RMSE_10_ct_all_p_summ_Tex$div <- RMSE_10_ct_all_p_summ_Tex$Mean- RMSE_10_ct_all_p_summ_Tex$Mean[1]
        RMSE_10_ct_all_p_summ_Tex$impr <- RMSE_10_ct_all_p_summ_Tex$div/RMSE_10_ct_all_p_summ_Tex$Mean*100
        RMSE_10_ct_all_p_summ_Tem$div <- RMSE_10_ct_all_p_summ_Tem$Mean- RMSE_10_ct_all_p_summ_Tem$Mean[1]
        RMSE_10_ct_all_p_summ_Tem$impr <- RMSE_10_ct_all_p_summ_Tem$div/RMSE_10_ct_all_p_summ_Tem$Mean*100
      }
    }
  }
  
  # Fig. SI PCC / F1 score
  {
    
    # summary PCC with 10 division by cell types
    {
      cell_type_list <- colnames(UCASpatial_highres[[1]][[2]])[1:20]
      
      
      PCC_divided10_ct <- lapply(highres_list, function(x){
        PCC_celltype <- lapply(cell_type_list, function(z){
          PCC <- lapply(seq(0.1,1,0.1), function(y){
            groud_truth_split_10 <- filter(groud_truth_split,proportion <= y) %>% 
              filter(grepl(z,new_ct)) %>% arrange(new_ct)
            split_10 <- filter(x,x$new_ct %in% groud_truth_split_10$new_ct) %>% 
              filter(grepl(z,new_ct)) %>% arrange(new_ct)
            return(cor(split_10$proportion,groud_truth_split_10$proportion,method = c("pearson")))
          })
          return(PCC)
        })
        return(PCC_celltype)
      })
      
      PCC_divided10_ct_md <- lapply(PCC_divided10_ct,function(x){
        test_10 <- matrix(unlist(x) ,nr=10) %>% as.data.frame()
        return(test_10)
      })
      
      PCC_10_ct_UCASpatial <- PCC_divided10_ct_md[[1]]
      colnames(PCC_10_ct_UCASpatial) <- cell_type_list
      rownames(PCC_10_ct_UCASpatial) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_SPOTlight <- PCC_divided10_ct_md[[2]]
      colnames(PCC_10_ct_SPOTlight) <-  <- 
        rownames(PCC_10_ct_SPOTlight) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_RCTD <- PCC_divided10_ct_md[[3]]
      colnames(PCC_10_ct_RCTD) <- cell_type_list
      rownames(PCC_10_ct_RCTD) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_C2L <- PCC_divided10_ct_md[[4]]
      colnames(PCC_10_ct_C2L) <- cell_type_list
      rownames(PCC_10_ct_C2L) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_CARD <- PCC_divided10_ct_md[[5]]
      colnames(PCC_10_ct_CARD) <- cell_type_list
      rownames(PCC_10_ct_CARD) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_stereoscope <- PCC_divided10_ct_md[[6]]
      colnames(PCC_10_ct_stereoscope) <- cell_type_list
      rownames(PCC_10_ct_stereoscope) <- paste('0-',seq(10,100,10),'%',sep = '')
      
      PCC_10_ct_max <- matrix(unlist(PCC_divided10_ct_md) ,nc=6) %>% as.data.frame() %>% 
        filter(rownames(.) %in% c(seq(10,200,10)))
      
      PCC_10_ct_avg <- lapply(list(PCC_10_ct_UCASpatial,PCC_10_ct_SPOTlight,PCC_10_ct_RCTD,
                                   PCC_10_ct_C2L,PCC_10_ct_CARD,PCC_10_ct_stereoscope), function(x){
                                     avg <- rowSums(x)/20
                                     return(avg)
                                   })
      PCC_10_ct_avg2 <- matrix(unlist(PCC_10_ct_avg) ,nc=6) %>% as.data.frame()
      rownames(PCC_10_ct_avg2) <- paste('0-',seq(10,100,10),'%',sep = '')
      colnames(PCC_10_ct_avg2) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      # saveRDS(PCC_10_ct_avg2,'PCC_10_ct_avg.rds')
    }
    # draw line chart of PCC for all methods
    {
      library(Rmisc)
      
      PCC_10_ct_all <- matrix(unlist(PCC_divided10_ct_md) ,nc=6) %>% as.data.frame()
      colnames(PCC_10_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
      PCC_10_ct_all_2 <- PCC_10_ct_all %>% mutate(proportion = paste('0-',rep(seq(10,100,10),20),'%',sep = ''),
                                                  celltype=rep(colnames(UCASpatial_highres[[1]][[2]])[1:20],each = 10))
      PCC_10_ct_all_p <- pivot_longer(PCC_10_ct_all_2,cols = colnames(PCC_10_ct_all_2)[1:6],names_to = 'algorithm',
                                      values_to = 'PCC')
      PCC_10_ct_all_p$proportion <- factor(PCC_10_ct_all_p$proportion,
                                           levels = paste('0-',seq(10,100,10),'%',sep = ''))
      
      PCC_10_ct_all_p$algorithm <- factor(PCC_10_ct_all_p$algorithm,
                                          levels = c('UCASpatial','SPOTlight','RCTD',
                                                     'cell2location','CARD','stereoscope'))
      
      # overall summary
      {
        PCC_10_ct_all_p_sum <- summarySE(PCC_10_ct_all_p, measurevar = "PCC",
                                         groupvars = c("algorithm","proportion"))
        PCC_10_ct_all_p_sum$proportion <- factor(PCC_10_ct_all_p_sum$proportion,
                                                 levels = paste('0-',seq(10,100,10),'%',sep = ''))
        
        PCC_10_ct_all_p_sum$algorithm <- factor(PCC_10_ct_all_p_sum$algorithm,
                                                levels = c('UCASpatial','SPOTlight','RCTD',
                                                           'cell2location','CARD','stereoscope'))
        
        ggplot(PCC_10_ct_all_p_sum,aes(proportion,PCC,group=algorithm,color=algorithm,
                                       linetype=algorithm,linewidth=algorithm))+
          # geom_errorbar(aes(ymin = PCC-se, ymax = PCC+se), width = 0.1, size = 1.2)+
          geom_point(size=1.5)+ ylim(0,0.8)+
          geom_line(position = position_dodge(0.1),cex=0.5)+
          scale_color_manual(values = ggradar_col[-5])+
          theme_bw()+xlab(NULL)+
          scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                           cell2location = 3,CARD=3,stereoscope=3))+
          scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                            cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
          # ggtitle(label = 'Gradient proportion')+
          theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
          theme(legend.position = "right")
        ggsave('Gradient_proportion_avg_ct_line_chart_PCC.pdf',height = 3,width = 4.5)
        
      }
      
      # for each cell type
      {
        dev.off()
        pdf('Gradient_proportion_avg_ct_line_chart_PCC_each_celltype_in1.pdf',height = 15,width = 18)
        p_PCC_each_ct_grad_pro <- lapply(cell_type_list,function(x){
          PCC_10_ct_all_p_temp <- PCC_10_ct_all_p %>% filter(celltype == x)
          ggplot(PCC_10_ct_all_p_temp,aes(proportion,PCC,group=algorithm,color=algorithm,
                                          linetype=algorithm,linewidth=algorithm))+
            # geom_errorbar(aes(ymin = PCC-se, ymax = PCC+se), width = 0.1, size = 1.2)+
            geom_point(size=1.5)+ #ylim(0,0.8)+
            geom_line(position = position_dodge(0.1),cex=0.5)+
            scale_color_manual(values = ggradar_col[-5])+
            theme_bw()+xlab(NULL)+labs(title=x)+
            scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                             cell2location = 3,CARD=3,stereoscope=3))+
            scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                              cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
            # ggtitle(label = 'Gradient proportion')+
            theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
            theme(legend.position = "right")
        })
        cowplot::plot_grid(plotlist = p_PCC_each_ct_grad_pro,ncol = 4)
        dev.off()
      }
    }
    
    ## Accuracy level by cell types
    {
      source('/data/xy/scripts/Test_Accuracy.R')
      # summary Acur,F1,... (with 10 division or not)
      {
        
        All_auc_10_ct <- lapply(highres_list, function(x){
          auc_celltype <- lapply(cell_type_list, function(z){
            auc <- lapply(seq(0.1,1,0.1), function(y){
              groud_truth_split_10 <- filter(groud_truth_split,proportion <= y) %>% 
                filter(grepl(z,new_ct)) %>% arrange(new_ct)
              split_10 <- filter(x,x$new_ct %in% groud_truth_split_10$new_ct) %>% 
                filter(grepl(z,new_ct)) %>% arrange(new_ct)
              Predict <- as.matrix(split_10$proportion_filter)
              Truth <- as.matrix(groud_truth_split_10$proportion)
              colnames(Predict) <- 'Value'
              colnames(Truth) <- 'Value'
              return(Test_Acurracy(Predict,Truth)[5])
            })
            return(auc)
          })
          return(auc_celltype)
        })
        
        
        auc_divided10_ct_md <- lapply(All_auc_10_ct,function(x){
          test_10 <- matrix(unlist(x) ,nr=10) %>% as.data.frame()
          return(test_10)
        })
        
        auc_10_ct_UCASpatial <- auc_divided10_ct_md[[1]]
        colnames(auc_10_ct_UCASpatial) <- cell_type_list
        rownames(auc_10_ct_UCASpatial) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_SPOTlight <- auc_divided10_ct_md[[2]]
        colnames(auc_10_ct_SPOTlight) <- cell_type_list
        rownames(auc_10_ct_SPOTlight) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_RCTD <- auc_divided10_ct_md[[3]]
        colnames(auc_10_ct_RCTD) <- cell_type_list
        rownames(auc_10_ct_RCTD) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_C2L <- auc_divided10_ct_md[[4]]
        colnames(auc_10_ct_C2L) <- cell_type_list
        rownames(auc_10_ct_C2L) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_CARD <- auc_divided10_ct_md[[5]]
        colnames(auc_10_ct_CARD) <- cell_type_list
        rownames(auc_10_ct_CARD) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_stereoscope <- auc_divided10_ct_md[[6]]
        colnames(auc_10_ct_stereoscope) <- cell_type_list
        rownames(auc_10_ct_stereoscope) <- paste('0-',seq(10,100,10),'%',sep = '')
        
        auc_10_ct_max <- matrix(unlist(auc_divided10_ct_md) ,nc=6) %>% as.data.frame() %>% 
          filter(rownames(.) %in% c(seq(10,200,10)))
        
        auc_10_ct_avg <- lapply(list(auc_10_ct_UCASpatial,auc_10_ct_SPOTlight,auc_10_ct_RCTD,
                                     auc_10_ct_C2L,auc_10_ct_CARD,auc_10_ct_stereoscope), function(x){
                                       avg <- rowSums(x,na.rm = T)/20
                                       return(avg)
                                     })
        auc_10_ct_avg2 <- matrix(unlist(auc_10_ct_avg) ,nc=6) %>% as.data.frame()
        rownames(auc_10_ct_avg2) <- paste('0-',seq(10,100,10),'%',sep = '')
        colnames(auc_10_ct_avg2) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
        
        # saveRDS(auc_10_ct_avg2,'auc_10_ct_avg.rds')
      }
      # draw line chart of F1 for all methods
      {
        F1_10_ct_all <- matrix(unlist(auc_divided10_ct_md) ,nc=6) %>% as.data.frame()
        colnames(F1_10_ct_all) <- c('UCASpatial','SPOTlight','RCTD','cell2location','CARD','stereoscope')
        F1_10_ct_all_2 <- F1_10_ct_all %>% mutate(proportion = paste('0-',rep(seq(10,100,10),20),'%',sep = ''),
                                                  celltype=rep(colnames(UCASpatial_highres[[1]][[2]])[1:20],each = 10))
        F1_10_ct_all_p <- pivot_longer(F1_10_ct_all_2,cols = colnames(F1_10_ct_all_2)[1:6],names_to = 'algorithm',
                                       values_to = 'F1')
        F1_10_ct_all_p$proportion <- factor(F1_10_ct_all_p$proportion,
                                            levels = paste('0-',seq(10,100,10),'%',sep = ''))
        
        F1_10_ct_all_p$algorithm <- factor(F1_10_ct_all_p$algorithm,
                                           levels = c('UCASpatial','SPOTlight','RCTD',
                                                      'cell2location','CARD','stereoscope'))
        # overall summary
        {
          F1_10_ct_all_p_sum <- summarySE(F1_10_ct_all_p, measurevar = "F1",
                                          groupvars = c("algorithm","proportion"),na.rm = T)
          F1_10_ct_all_p_sum$proportion <- factor(F1_10_ct_all_p_sum$proportion,
                                                  levels = paste('0-',seq(10,100,10),'%',sep = ''))
          
          F1_10_ct_all_p_sum$algorithm <- factor(F1_10_ct_all_p_sum$algorithm,
                                                 levels = c('UCASpatial','SPOTlight','RCTD',
                                                            'cell2location','CARD','stereoscope'))
          
          ggplot(F1_10_ct_all_p_sum,aes(x=proportion,y=F1,group=algorithm,color=algorithm,
                                        linetype=algorithm,linewidth=algorithm))+
            geom_point(size=2)+ ylim(0,0.8)+
            geom_line(position = position_dodge(0.1),cex=0.5)+
            scale_color_manual(values = ggradar_col[-5])+
            theme_bw()+xlab(NULL)+ylab('F1 score')+
            scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                             cell2location = 3,CARD=3,stereoscope=3))+
            scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                              cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
            # ggtitle(label = 'Gradient proportion')+
            theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
            theme(legend.position = "right")
          ggsave('Gradient_proportion_avg_ct_line_chart_F1.pdf',height = 3,width = 4.5)
        }
        
        # for each cell type
        {
          dev.off()
          pdf('Gradient_proportion_avg_ct_line_chart_F1_each_celltype.pdf',height = 3,width = 4.5)
          p_F1_each_ct_grad_pro <- lapply(cell_type_list,function(x){
            F1_10_ct_all_p_temp <- F1_10_ct_all_p %>% filter(celltype == x)
            ggplot(F1_10_ct_all_p_temp,aes(x=proportion,y=F1,group=algorithm,color=algorithm,
                                           linetype=algorithm,linewidth=algorithm))+
              geom_point(size=2)+ #ylim(0,0.8)+
              geom_line(position = position_dodge(0.1),cex=0.5)+
              scale_color_manual(values = ggradar_col[-5])+
              theme_bw()+xlab(NULL)+ylab('F1 score')+labs(title = x)+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                                cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")
          })
          p_F1_each_ct_grad_pro
          dev.off()
          
          dev.off()
          pdf('Gradient_proportion_avg_ct_line_chart_F1_each_celltype_In1.pdf',height = 15,width = 18)
          p_F1_each_ct_grad_pro <- lapply(cell_type_list,function(x){
            F1_10_ct_all_p_temp <- F1_10_ct_all_p %>% filter(celltype == x)
            ggplot(F1_10_ct_all_p_temp,aes(x=proportion,y=F1,group=algorithm,color=algorithm,
                                           linetype=algorithm,linewidth=algorithm))+
              geom_point(size=2)+ #ylim(0,0.8)+
              geom_line(position = position_dodge(0.1),cex=0.5)+
              scale_color_manual(values = ggradar_col[-5])+
              theme_bw()+xlab(NULL)+ylab('F1 score')+labs(title = x)+
              scale_linetype_manual(values = c(UCASpatial = 1, SPOTlight = 3, RCTD=3,
                                               cell2location = 3,CARD=3,stereoscope=3))+
              scale_linewidth_manual(values = c(UCASpatial = 1, SPOTlight = 0.75, RCTD=0.75,
                                                cell2location = 0.75,CARD=0.75,stereoscope=0.75))+
              # ggtitle(label = 'Gradient proportion')+
              theme(axis.text.x = element_text(angle = 60,hjust = 1),plot.title = element_text(hjust = 0.5))+
              theme(legend.position = "right")
          })
          cowplot::plot_grid(plotlist = p_F1_each_ct_grad_pro,ncol = 4)
          dev.off()
        }
      }
    }
  }
}
