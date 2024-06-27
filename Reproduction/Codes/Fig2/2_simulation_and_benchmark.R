library(Seurat)
library(pacman)
p_unload(Seurat)
p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig2/')

source('../source/UCASpatial_final_v1.R')
col.sc <- readRDS('../../Data/Fig2/col.sc.re.rds')
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")

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

