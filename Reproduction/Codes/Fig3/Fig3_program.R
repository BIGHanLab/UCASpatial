### read st and sc data 
setwd('/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Codes/Fig3/')


######## step 1 read decov results #####
{
  ### UCASpatial
  ucaspatial.mat <- plot_data_espanno[,13:32] %>% as.matrix()
  colnames(ucaspatial.mat)
  all(row.names(ucaspatial.mat)==colnames(MRL.st))
  ### CARD
  card_res <- read.table('../../Data/Fig3/CARD/CRAD_regeneration_decon_result.txt',sep = '\t')
  head(card_res)
  card.mat <- card_res[row.names(card_res) %in% colnames(MRL.st),]
  colnames(card.mat)[c(1,2,3,6,13,10,11,14)] <- colnames(ucaspatial.mat)[c(2:5,12:15)]
  colnames(card.mat) <- gsub('\\.',' ',colnames(card.mat))
  card.mat <- card.mat[,colnames(ucaspatial.mat)]
  all(row.names(card.mat) == colnames(MRL.st))
  ### Cell2location
  cell2location_res <- read.csv('../../Data/Fig3/cell2location/result/Reg_C2L.csv',row.names = 1)
  head(cell2location_res)
  dim(cell2location_res)
  cell2location.mat <- cell2location_res[row.names(cell2location_res) %in% colnames(MRL.st),]
  all(row.names(cell2location.mat) == colnames(MRL.st))
  colnames(cell2location.mat)[c(1,6,7,9,11,12,17,19)] <- c('Cd163+ Mac',
                                                                    'Fmod+ fibroblast',
                                                                    'Lrrc15+ fibroblast',
                                                                    'Cd34+ fibroblast',
                                                                    'Mki67+ fibroblast',
                                                                    'MHCII+ Mac',
                                                                    'Il1b+ Mac',
                                                                    'Spp1+ Mac')
  colnames(cell2location.mat) <- gsub('\\.',' ',colnames(cell2location.mat))
  cell2location.mat <- cell2location.mat[,  colnames(ucaspatial.mat)]
  head(cell2location.mat)
  all(colnames(cell2location.mat) == levels(MRL.sc.data$Minor_class_lasted))
  
  ####RCTD
  rctd_red <- readRDS('../../Data/Fig3/RCTD/result_RCTD.rds')
  rctd_res_mat <- rctd_red@results[["weights"]] %>% data.frame()
  head(rctd_res_mat)
  RCTD_res_mat <- rctd_res_mat[colnames(MRL.st),]
  
  MRL.st.rctd <- subset(MRL.st, spot_id %in% row.names(RCTD_res_mat))
  RCTD_res_mat <- rctd_res_mat[colnames(MRL.st.rctd),]
  
  all(row.names(RCTD_res_mat) == colnames(MRL.st.rctd))
  RCTD_res_mat.norm <- RCTD_res_mat/rowSums(RCTD_res_mat)
  colnames(RCTD_res_mat)[c(2:5,12:15)] <- colnames(ucaspatial.mat)[c(2:5,12:15)]
  colnames(RCTD_res_mat) <- gsub('\\.',' ',colnames(RCTD_res_mat))
  RCTD_res_mat <- RCTD_res_mat[,colnames(ucaspatial.mat)]
  RCTD_res_mat <- RCTD_res_mat[row.names(RCTD_res_mat) %in% colnames(MRL.st.rctd),]
  all(row.names(RCTD_res_mat) == colnames(MRL.st.rctd))
  ####SPOTlight
  spotlight_res <- readRDS('../../Data/Fig3/SPOTlight/KRAS_SPOTlight_result.rds')
  spotlight_res[[2]]
  spotlight_res_mat <- spotlight_res[[2]] %>% data.frame()
  head(spotlight_res_mat)
  spotlight_res_mat <- spotlight_res_mat[,-ncol(spotlight_res_mat)]
  load('../../Data/Fig3/MRL.st.primary.Rdata')
  row.names(spotlight_res_mat) <- colnames(ear.integrated.filter)
  head(spotlight_res_mat)
  colnames(spotlight_res_mat)[c(1,19,12,17,9,6,11,7)] <- colnames(ucaspatial.mat)[c(2:5,12:15)]
  colnames(spotlight_res_mat) <- gsub('\\.',' ',colnames(spotlight_res_mat))
  spotlight_res_mat <- spotlight_res_mat[,match(colnames(ucaspatial.mat),colnames(spotlight_res_mat))]
  spotlight_res_mat <- spotlight_res_mat[row.names(spotlight_res_mat) %in% colnames(MRL.st),]
  all(row.names(spotlight_res_mat) == colnames(MRL.st))
  ##### stereoscope
  stereoscope_res <- read.table("../../Data/Fig3/stereoscope/stereo_decon_df.txt",sep = '\t')
  head(stereoscope_res)
  stereoscope_res_mat <- stereoscope_res[row.names(stereoscope_res) %in% colnames(MRL.st),]
  all(row.names(stereoscope_res_mat) == colnames(MRL.st))
  colnames(stereoscope_res_mat)[c(1,6,7,9,11,12,17,19)] <- c('Cd163+ Mac',
                                                             'Fmod+ fibroblast',
                                                             'Lrrc15+ fibroblast',
                                                             'Cd34+ fibroblast',
                                                             'Mki67+ fibroblast',
                                                             'MHCII+ Mac',
                                                             'Il1b+ Mac',
                                                             'Spp1+ Mac')
  
  colnames(stereoscope_res_mat) <- gsub('\\.',' ',colnames(stereoscope_res_mat))
  colnames(stereoscope_res_mat)
  stereoscope_res_mat <- stereoscope_res_mat[,colnames(ucaspatial.mat)]
  
  ct.name <- colnames(stereoscope_res_mat)
  ############
  deconv.list <- list(ucaspatial.mat,card.mat,cell2location.mat,
                      RCTD_res_mat,spotlight_res_mat,stereoscope_res_mat)
  names(deconv.list) <- c('UCASpatial','CARD','Cell2location','RCTD','SPOTlight','stereoscope')
  min.cutoff <- 0.03
  deconv.list <- mapply(function(x){
    if(x %in% c('Cell2location','RCTD')){
      decon_mtrx <- deconv.list[[x]]
      decon_mtrx <- decon_mtrx/rowSums(decon_mtrx)
      decon_mtrx[decon_mtrx < min.cutoff] <- 0
      decon_mtrx <- decon_mtrx/rowSums(decon_mtrx)
      decon_df <- decon_mtrx %>%
        data.frame()
      colnames(decon_df) <- ct.name
      decon_df$Fibroblast <- rowSums(decon_df[,12:15])
      decon_df$Macrophage <- rowSums(decon_df[,2:5])
      return(list(decon_df))
    }else{
      decon_mtrx <- deconv.list[[x]]
      decon_mtrx[decon_mtrx < min.cutoff] <- 0
      decon_mtrx <- decon_mtrx/rowSums(decon_mtrx)
      decon_df <- decon_mtrx %>%
        data.frame()
      colnames(decon_df) <- ct.name
      decon_df$Fibroblast <- rowSums(decon_df[,12:15])
      decon_df$Macrophage <- rowSums(decon_df[,2:5])
      return(list(decon_df))
    }
    
  },names(deconv.list))
  deconv.list[[1]] %>% head()
  save(deconv.list,file = 'deconv.list.Rdata')
  rowSums(deconv.list[[6]][,1:20]) %>% table
  MRL.st@assays$UCASpatial <- CreateAssayObject(counts = t(deconv.list$UCASpatial))
  MRL.st@assays$CARD <- CreateAssayObject(counts = t(deconv.list$CARD))
  MRL.st@assays$Cell2location <- CreateAssayObject(counts = t(deconv.list$Cell2location))
  MRL.st.rctd@assays$RCTD <- CreateAssayObject(counts = t(deconv.list$RCTD))
  MRL.st@assays$SPOTlight <- CreateAssayObject(counts = t(deconv.list$SPOTlight))
  MRL.st@assays$stereoscope <- CreateAssayObject(counts = t(deconv.list$stereoscope))
  save(MRL.st,MRL.st.rctd,file = 'MRL.st.add.decov.Rdata')
  all(row.names(deconv.list$stereoscope) == row.names(MRL.coord))
  rege.fro.spot <- row.names(plot_data_UCASpatial_sub)
  wound.center.spot <- intersect(rege.fro.spot,row.names(MRL.coord)[MRL.coord$Spot_group_v2 == 'wound'])
  plot_data_UCASpatial <- cbind(MRL.coord,deconv.list$UCASpatial)
  plot_data_UCASpatial$Method <- 'UCASpatial'
  head(plot_data_UCASpatial)
  plot_data_UCASpatial_sub <- plot_data_UCASpatial[wound.center.spot,]
  
  plot_data_CARD <- cbind(MRL.coord,deconv.list$CARD)
  plot_data_CARD$Method <- 'CARD'
  plot_data_CARD_sub<- plot_data_CARD[ wound.center.spot,]
  
  plot_data_Cell2location <- cbind(MRL.coord,deconv.list$Cell2location)
  plot_data_Cell2location$Method <- 'Cell2location'
  plot_data_Cell2location_sub<- plot_data_Cell2location[ wound.center.spot,]
  
  MRL.coord.rctd <- MRL.coord[row.names(deconv.list$RCTD),]
  plot_data_RCTD <- cbind(MRL.coord.rctd,deconv.list$RCTD)
  plot_data_RCTD$Method <- 'RCTD'
  plot_data_RCTD_sub <- plot_data_RCTD[ wound.center.spot,]
  
  plot_data_SPOTlight <- cbind(MRL.coord,deconv.list$SPOTlight)
  plot_data_SPOTlight$Method <- 'SPOTlight'
  plot_data_SPOTlight_sub<- plot_data_SPOTlight[wound.center.spot,]
  
  plot_data_stereoscope <- cbind(MRL.coord,deconv.list$stereoscope)
  plot_data_stereoscope$Method <- 'stereoscope'
  plot_data_stereoscope_sub <- plot_data_stereoscope[ wound.center.spot,]
  
  save(plot_data_CARD,plot_data_CARD_sub,plot_data_Cell2location,plot_data_Cell2location_sub,
       plot_data_RCTD,plot_data_RCTD_sub,plot_data_SPOTlight,plot_data_SPOTlight_sub,
       plot_data_stereoscope,plot_data_stereoscope_sub,plot_data_UCASpatial,plot_data_UCASpatial_sub,
       file = 'plot_data_6_methods.Rdata')
}
### ucaspatial without entropy 
{
  ### ucaspatial without entrophy 
  ucaspatial.mat.nw <- readRDS('/data/xy/Spatial_transcriptome/eWEIDE/20240517_UCASpatial_in_real_without_weight/without_weight_and_meta/regeneration/eSpanno_reg_noweight.rds')
  ucaspatial.mat.nw <- ucaspatial.mat.nw[[2]]
  row.names(ucaspatial.mat.nw) <- colnames(ear.integrated.filter)
  ucaspatial.mat.nw <- ucaspatial.mat.nw[row.names(ucaspatial.mat.nw) %in% colnames(MRL.st),]
  colnames(ucaspatial.mat.nw)[c(1,19,12,17,9,6,11,7)] <- colnames(ucaspatial.mat)[c(2:5,12:15)]
  colnames(ucaspatial.mat.nw) <- gsub('\\.',' ',colnames(ucaspatial.mat.nw))
  ucaspatial.mat.nw <- ucaspatial.mat.nw[,colnames(ucaspatial.mat)]
  all(row.names(ucaspatial.mat.nw) == colnames(MRL.st))
  {
    ## matrix
    decon_mtrx <- ucaspatial.mat.nw
    decon_mtrx[decon_mtrx < min.cutoff] <- 0
    decon_mtrx <- decon_mtrx/rowSums(decon_mtrx)
    decon_df <- decon_mtrx %>%
      data.frame()
    colnames(decon_df) <- ct.name
    decon_df$Fibroblast <- rowSums(decon_df[,12:15])
    decon_df$Macrophage <- rowSums(decon_df[,2:5])
  }
  deconv.list$UCASpatial_UW <- decon_df
  
  plot_data_ucaspatial_nw <- cbind(MRL.coord,deconv.list$UCASpatial_UW)
  plot_data_ucaspatial_nw$Method <- 'UCASpatial_UW'
  plot_data_UCASpatial_UW_sub<- plot_data_ucaspatial_nw[ wound.center.spot,]
  
  
  save(plot_data_CARD,plot_data_CARD_sub,plot_data_Cell2location,plot_data_Cell2location_sub,
       plot_data_RCTD,plot_data_RCTD_sub,plot_data_SPOTlight,plot_data_SPOTlight_sub,
       plot_data_stereoscope,plot_data_stereoscope_sub,plot_data_UCASpatial,plot_data_UCASpatial_sub,
       plot_data_ucaspatial_nw,plot_data_UCASpatial_UW_sub,
       file = 'plot_data_6_methods_plus_abolish.Rdata')
  
  MRL.st@assays$UCASpatial_UW <- CreateAssayObject(counts = t(deconv.list$UCASpatial_UW))
  
}
##### fig3c 
{
  p1 <- ggplot(data = plot_data_UCASpatial_sub)+
    theme_base(base_size = 15)+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(size = 15),
          plot.background =  element_blank()
    )+
    geom_point_rast( aes(x=row,y=col,fill=Fibroblast ,#alpha =Fibroblast
                         ),
                     size=4,shape=21)+
    guides(alpha = F)+
    coord_fixed(ratio = 0.5)+
    xlim(18,28)+ylim(-65,-47)+
    labs(fill='Proportion',title = 'Fibroblast')+
    theme(plot.title = element_text(face = 'plain',hjust = 0.5,size=15))+
    scale_fill_gradientn(labels = scales::label_percent(),
                         colors =colorRampPalette(c('white',brewer.pal(n = 9, name = "Oranges")))(n=100) )
  p1
  #ggsave('MRL_D7_proximal_Fibroblast.pdf',width = 4.51,height = 3.5)
  
  p2 <- ggplot(data = plot_data_UCASpatial_sub)+
    theme_base(base_size = 15)+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(size = 15),
          plot.background =  element_blank()
    )+
    geom_point_rast( aes(x=row,y=col,fill=Keratinocyte ,#alpha =Keratinocyte
                         ),
                     size=4,shape=21)+
    guides(alpha = F)+
    coord_fixed(ratio = 0.5)+
    xlim(18,28)+ylim(-65,-47)+
    labs(fill='Proportion',title = 'Keratinocyte')+
    theme(plot.title = element_text(face = 'plain',hjust = 0.5,size=15))+
    scale_fill_gradientn(labels = scales::label_percent(),
                         colors =colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100) )
  p2
  #ggsave('MRL_D7_proximal_Keratinocyte.pdf',width = 4.51,height = 3.5)
  
  p3 <- ggplot(data = plot_data_UCASpatial_sub)+
    theme_base(base_size = 15)+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(size = 15),
          plot.background =  element_blank()
    )+
    geom_point_rast( aes(x=row,y=col,fill=`Endothelial cell`  ,#alpha =`Endothelial cell` 
                         ),
                     size=4,shape=21)+
    guides(alpha = F)+
    coord_fixed(ratio = 0.5)+ 
    xlim(18,28)+ylim(-65,-47)+
    labs(fill='Proportion',title = 'Endothelial')+
    theme(plot.title = element_text(face = 'plain',hjust = 0.5,size = 15))+
    scale_fill_gradientn(labels = scales::label_percent(),
                         colors =colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100) )
  p3
  
  p_list <- list(p2,p3,p1)
  plot_grid(plotlist = p_list,nrow = 3,align = "hv")
  ggsave('fig3c_MRL_D7_Promixal_Fibroblast_endothelial_keratinocyte_v3.pdf',
         width = 3.98,height = 5.4)
  
}

### fig3d
{
  ### UCASpatial
  {
    p1 <- ggplot(data = plot_data_UCASpatial_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac` ,#alpha =`Cd163+ Mac`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           #limits = c(0, quantile(plot_data_UCASpatial_sub$Keratinocyte, 0.99)),
                           #limits = c(0,0.5),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    p2 <- ggplot(data = plot_data_UCASpatial_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac` ,#alpha =`Spp1+ Mac`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    p3 <- ggplot(data = plot_data_UCASpatial_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac` ,#alpha =`MHCII+ Mac`
      ),
      size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           #limits = c(0, quantile(plot_data_UCASpatial_sub$Keratinocyte, 0.99)),
                           #limits = c(0, 0.5),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    ##############################################################
    p4 <- ggplot(data = plot_data_UCASpatial_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac` ,#alpha =`Il1b+ Mac`
      ),
      size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           #limits = c(0, quantile(plot_data_UCASpatial_sub$`Endothelial cell`, 0.99)),
                           #limits = c(0, 0.5),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    p_list <- list(p1,p2,p3,p4)
    ucaspatial_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    ucaspatial_p
    ggsave('fig3d_MRL_D7_split_4_macrophage_espanno_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
    
  }
  ### cell2location
  {
    plot_data_Cell2location_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_Cell2location_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                             max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_RCTD_sub$`Cd163+ Mac`)
    plot_data_Cell2location_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_Cell2location_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                             max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_RCTD_sub$`MHCII+ Mac`)
    plot_data_Cell2location_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_Cell2location_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                            max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_RCTD_sub$`Il1b+ Mac`)
    plot_data_Cell2location_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_Cell2location_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                            max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_RCTD_sub$`Spp1+ Mac`)
    p1 <- ggplot(data = plot_data_Cell2location_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_Cell2location_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03` ,#alpha =`Spp1+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    
    p3 <- ggplot(data = plot_data_Cell2location_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_Cell2location_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    p_list <- list(p1,p2,p3,p4)
    cl_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    cl_p
    ggsave('fig3d_Cell2loaction_MRL_macophage_D7_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
  }
  ### RCTD
  {
    plot_data_RCTD_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_RCTD_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                    max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_RCTD_sub$`Cd163+ Mac`)
    plot_data_RCTD_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_RCTD_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                    max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_RCTD_sub$`MHCII+ Mac`)
    plot_data_RCTD_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_RCTD_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                   max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_RCTD_sub$`Il1b+ Mac`)
    plot_data_RCTD_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_RCTD_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                   max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_RCTD_sub$`Spp1+ Mac`)
    
    p1 <- ggplot(data = plot_data_RCTD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_RCTD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03`,#alpha =`Spp1+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    p3 <- ggplot(data = plot_data_RCTD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_RCTD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    
    p_list <- list(p1,p2,p3,p4)
    rctd_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    rctd_p
    ggsave('fig3d_RCTD_MRL_macophage_D7_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
  }
  ### SPOTlight
  {
    #### 
    plot_data_SPOTlight_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_SPOTlight_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                         max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_SPOTlight_sub$`Cd163+ Mac`)
    plot_data_SPOTlight_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_SPOTlight_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                         max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_SPOTlight_sub$`MHCII+ Mac`)
    plot_data_SPOTlight_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_SPOTlight_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                        max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_SPOTlight_sub$`Il1b+ Mac`)
    plot_data_SPOTlight_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_SPOTlight_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                        max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_SPOTlight_sub$`Spp1+ Mac`)
    
    p1 <- ggplot(data = plot_data_SPOTlight_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_SPOTlight_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03` ,#alpha =`Spp1+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    
    p2
    p3 <- ggplot(data = plot_data_SPOTlight_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_SPOTlight_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    
    p_list <- list(p1,p2,p3,p4)
    spotlight_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    spotlight_p
    ggsave('fig3d_SPOTlight_MRL_macophage_D7_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
  }
  ### stereoscope
  {
    #### 
    plot_data_stereoscope_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_stereoscope_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                           max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_stereoscope_sub$`Cd163+ Mac`)
    plot_data_stereoscope_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_stereoscope_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                           max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_stereoscope_sub$`MHCII+ Mac`)
    plot_data_stereoscope_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_stereoscope_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                          max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_stereoscope_sub$`Il1b+ Mac`)
    plot_data_stereoscope_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_stereoscope_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                          max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_stereoscope_sub$`Spp1+ Mac`)
    
    p1 <- ggplot(data = plot_data_stereoscope_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_stereoscope_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03` ,#alpha =`Spp1+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    p3 <- ggplot(data = plot_data_stereoscope_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_stereoscope_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    
    p_list <- list(p1,p2,p3,p4)
    stereo_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    stereo_p
    ggsave('fig3d_stereoscope_MRL_macophage_D7_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
  }
  ### CARD
  {
    ### CARD 
    plot_data_CARD_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_CARD_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                    max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_CARD_sub$`Cd163+ Mac`)
    plot_data_CARD_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_CARD_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                    max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_CARD_sub$`MHCII+ Mac`)
    plot_data_CARD_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_CARD_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                   max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_CARD_sub$`Il1b+ Mac`)
    plot_data_CARD_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_CARD_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                   max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_CARD_sub$`Spp1+ Mac`)
    
    head(plot_data_CARD_sub)
    p1 <- ggplot(data = plot_data_CARD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_CARD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03` ,#alpha =`Spp1+ Mac F0.03`
                      ),
                  size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    p3 <- ggplot(data = plot_data_CARD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_CARD_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
                      ),
                  size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    
    p_list <- list(p1,p2,p3,p4)
    card_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    card_p
    ggsave('fig3d_CARD_MRL_D7_mac_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
  }
  ### UCASpatial_UW
  {
    ##UCASpatial_UW
    plot_data_UCASpatial_UW_sub$`Cd163+ Mac F0.03` <- ifelse(plot_data_UCASpatial_UW_sub$`Cd163+ Mac` > max(plot_data_UCASpatial_sub$`Cd163+ Mac`),
                                                             max(plot_data_UCASpatial_sub$`Cd163+ Mac`),plot_data_RCTD_sub$`Cd163+ Mac`)
    plot_data_UCASpatial_UW_sub$`MHCII+ Mac F0.03` <- ifelse(plot_data_UCASpatial_UW_sub$`MHCII+ Mac` > max(plot_data_UCASpatial_sub$`MHCII+ Mac`),
                                                             max(plot_data_UCASpatial_sub$`MHCII+ Mac`),plot_data_RCTD_sub$`MHCII+ Mac`)
    plot_data_UCASpatial_UW_sub$`Il1b+ Mac F0.03` <- ifelse(plot_data_UCASpatial_UW_sub$`Il1b+ Mac` > max(plot_data_UCASpatial_sub$`Il1b+ Mac`),
                                                            max(plot_data_UCASpatial_sub$`Il1b+ Mac`),plot_data_RCTD_sub$`Il1b+ Mac`)
    plot_data_UCASpatial_UW_sub$`Spp1+ Mac F0.03` <- ifelse(plot_data_UCASpatial_UW_sub$`Spp1+ Mac` > max(plot_data_UCASpatial_sub$`Spp1+ Mac`),
                                                            max(plot_data_UCASpatial_sub$`Spp1+ Mac`),plot_data_RCTD_sub$`Spp1+ Mac`)
    p1 <- ggplot(data = plot_data_UCASpatial_UW_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Cd163+ Mac F0.03` ,#alpha =`Cd163+ Mac F0.03`
      ),
      size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Reds")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Cd163+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p1
    
    
    p2 <- ggplot(data = plot_data_UCASpatial_UW_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Spp1+ Mac F0.03` ,#alpha =`Spp1+ Mac F0.03`
      ),
      size=4,shape=21)+
      #scale_color_gradient2(low = "white", mid = "black", high =  "red")+
      scale_fill_gradient2(low = "white", high =  "gold2",
                           limits = c(0, max(plot_data_UCASpatial_sub$`Spp1+ Mac`)),
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p2
    
    
    p3 <- ggplot(data = plot_data_UCASpatial_UW_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`MHCII+ Mac F0.03` ,#alpha =`MHCII+ Mac F0.03`
      ),
      size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Blues")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`MHCII+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p3
    
    
    p4 <- ggplot(data = plot_data_UCASpatial_UW_sub)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15),
            plot.background =  element_blank()
      )+
      geom_point( aes(x=row,y=col,fill=`Il1b+ Mac F0.03` ,#alpha =`Il1b+ Mac F0.03`
      ),
      size=4,shape=21)+  
      scale_fill_gradientn(colours=colorRampPalette(c('white',brewer.pal(n = 9, name = "Purples")))(n=100),
                           limits = c(0, max(plot_data_UCASpatial_sub$`Il1b+ Mac`)),
                           na.value = "#9E0142",
                           labels = scales::label_percent())+
      guides(alpha = F)+
      xlim(18,28)+ylim(-65,-47)+
      coord_fixed(ratio = 0.5)+
      labs(fill ='Proportion')
    p4
    
    p_list <- list(p1,p2,p3,p4)
    UCASpatial_UW_p <- plot_grid(plotlist = p_list,ncol = 1,align = 'hv')
    UCASpatial_UW_p
    ggsave('fig3d_UCASpatial_UW_MRL_macophage_D7_cutoff0.03_v3.pdf',width = 4.98,height = 6.16)
    
  }
  plot_grid(ucaspatial_p,UCASpatial_UW_p,card_p,cl_p,rctd_p,spotlight_p,stereo_p,
            align = 'hv',nrow = 1)
  ggsave('fig3d_6_method_MRL_D7_mac_cutoff0.03_adding_ablation_v4.pdf',
         width = 21,height = 7)
  
}
##### fig 3c
{
  plot_data_UCASpatial_sub$Region <- gsub('Epithelium','Epidermis',plot_data_UCASpatial_sub$Region)
  plot_data_UCASpatial_sub$Region <- factor(plot_data_UCASpatial_sub$Region,
                                            levels = c('Epidermis','Dermis'))
  ggplot(plot_data_UCASpatial_sub,aes(x=row,y=col,fill=Region))+
    geom_point(shape=21,size=8)+
    theme_base()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(size = 15),
          plot.background =  element_blank()
    )+
    scale_fill_manual(values = c('#4DAF4A','#FF7F00'))+
    xlim(18,28)+ylim(-65,-47)+
    coord_fixed(ratio = 0.5)
  ggsave('fig3c_epidermis_dermis_region_in_wound_center.pdf',
         width = 4.55,height = 3.22)
}
## fig 3e
{
  UCASpatial_df <- espanno.df[,c(1:15,38:43)]
  UCASpatial_df <- cbind(UCASpatial_df,deconv.list$UCASpatial)
  UCASpatial.df.long.raw <- gather(data = UCASpatial_df,key = 'subset',value = 'Proportion',`Cd163+ Mac`:`Il1b+ Mac`)
  plot.data.raw <- aggregate(x=UCASpatial.df.long.raw$Proportion, by=list(UCASpatial.df.long.raw$Spot_group_v2,
                                                                       UCASpatial.df.long.raw$subset,
                                                                       UCASpatial.df.long.raw$Day),mean)
  head(plot.data.raw)
  colnames(plot.data.raw) <- c('Group',
                               'Subset','Day','Proportion')
  
  
  
  plot.data.ggplot.raw <- aggregate(x=UCASpatial.df.long.raw$Proportion, by=list(UCASpatial.df.long.raw$Spot_group_v2,
                                                                              UCASpatial.df.long.raw$subset,
                                                                              UCASpatial.df.long.raw$Day),mean)
  head(plot.data.ggplot.raw)
  colnames(plot.data.ggplot.raw) <- c('Group',
                                      'Subset','Day','Proportion')
  head(plot.data.ggplot.raw)
  ggplot(plot.data.ggplot.raw,aes(x=Day,y=Group,fill=Proportion))+
    #geom_tile(colour = 'black',size=0.5)+
    geom_raster()+
    facet_wrap(~Subset)+
    scale_fill_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"))+
    theme_base()+
    theme(plot.background = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())+
    labs(fill = 'Proportion')
  plot.data.ggplot.raw.fil <- plot.data.ggplot.raw[plot.data.ggplot.raw$Day !='D10' & plot.data.ggplot.raw$Group != 'Other',]
  p_list <- lapply(unique(plot.data.ggplot.raw.fil$Subset), function(x){
    tmp <- plot.data.ggplot.raw.fil[plot.data.ggplot.raw.fil$Subset == x ,]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      facet_wrap(~Subset)+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base()+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank())+
      labs(fill = 'Proportion')+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))
  })
  plot_grid(plotlist = p_list)
  plot.data.ggplot.raw.fil$Group <- sapply(as.character(plot.data.ggplot.raw.fil$Group),
  function(x){
    switch (x,
      'Uninjury' = 'Uninjury',
      'wound' = 'Wound center',
      'proximal' = 'Wound proximal',
      'distal' = 'Wound distal'
      
    )
  })
  plot.data.ggplot.raw.fil$Group <- factor(plot.data.ggplot.raw.fil$Group,
                                           levels = c('Uninjury','Wound center','Wound proximal','Wound distal'))
  plot.data.ggplot.raw.fil$Day <- gsub('D','Day ',as.character(plot.data.ggplot.raw.fil$Day))
  plot.data.ggplot.raw.fil$Day <- factor(plot.data.ggplot.raw.fil$Day,
                                         levels = c('Day 0','Day 3','Day 7','Day 15'))
  p_list <- lapply(c("Cd163+ Mac","Spp1+ Mac" ,"MHCII+ Mac" ,"Il1b+ Mac"), function(x){
    tmp <- plot.data.ggplot.raw.fil[plot.data.ggplot.raw.fil$Subset == x & plot.data.ggplot.raw.fil$Day != 'Day 0',]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 15)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12),
            legend.title = element_text(size = 12))+
      labs(fill = 'Average\nproportion',
           title = paste(x,'\n','UW:',round(plot.data.ggplot.raw.fil$Proportion[plot.data.ggplot.raw.fil$Day == 'Day 0' & plot.data.ggplot.raw.fil$Subset == x]*100,2),'%',sep=''))+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))
  })
  plot_grid(plotlist = p_list)
  ggsave('fig3e_mac_heatmap_diff_3_groups_v4_without normalization.pdf',width = 8.21,height = 3.88)
  
}

#### fig S3c
{
  ###### all cell types
  UCASpatial.df.long.all <- gather(data = UCASpatial_df,key = 'subset',value = 'Proportion',Monocyte:Macrophage)
  
  plot.data.all.raw <- aggregate(x=UCASpatial.df.long.all$Proportion, by=list(UCASpatial.df.long.all$Spot_group,
                                                                           UCASpatial.df.long.all$subset,
                                                                           UCASpatial.df.long.all$Day),mean)
  head(plot.data.all.raw)
  colnames(plot.data.all.raw) <- c('Group',
                                   'Subset','Day','Proportion')
  
  
  plot.data.all.raw$Group <- sapply(as.character(plot.data.all.raw$Group),
                                    function(x){
                                      switch (x,
                                              'Uninjury' = 'Unwounded',
                                              'Wound' = 'Wound center',
                                              'proximal' = 'Wound proximal',
                                              'distal' = 'Wound distal',
                                              'Other' = 'Other'
                                      )
                                    })
  plot.data.all.raw$Group <- factor(plot.data.all.raw$Group,
                                    levels = c('Unwounded','Wound center','Wound proximal','Wound distal','Other'))
  used.plot <- plot.data.all.raw[plot.data.all.raw$Subset %in% c("B cell" ,"Basophil" ,"Chondrocyte",
                                                                 "Dendritic cell","Endothelial cell",
                                                                 "Fibroblast","Keratinocyte","Macrophage","Mast cell",
                                                                 "Neutrophil","Perivascular cell","Satellite cell","T cell"),]
  
  used.plot$Subset <- factor(used.plot$Subset,
                             levels = c('Macrophage','Dendritic cell',"Neutrophil","T cell","B cell",
                                        "Basophil","Mast cell" , "Fibroblast"  ,"Chondrocyte",
                                        "Keratinocyte", "Perivascular cell","Endothelial cell","Satellite cell"))
  used.plot.fil <- used.plot[used.plot$Group != 'Other' & used.plot$Day != 'D10',]
  p_list <- lapply(levels(used.plot.fil$Subset), function(x){
    tmp <- used.plot.fil[used.plot.fil$Subset == x ,]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      facet_wrap(~Subset)+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 12)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12))+
      labs(fill = 'Proportion')+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))
    
  })
  plot_grid(plotlist = p_list,ncol = 3)
  #ggsave('heatmap_all_cell_types_v3.pdf',width = 10.97,height = 8)
  
  used.plot.fil$Day <- gsub('D','Day ',as.character(used.plot.fil$Day))
  used.plot.fil$Day <- factor(used.plot.fil$Day,
                                         levels = c('Day 0','Day 3','Day 7','Day 15'))
  
  p_list <- lapply(levels(used.plot.fil$Subset), function(x){
    tmp <- used.plot.fil[used.plot.fil$Subset == x & used.plot.fil$Day != 'Day 0',]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 12)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1))+
      labs(fill = 'Proportion',title = paste(x,'\n','UW:',round(used.plot.fil$Proportion[used.plot.fil$Day == 'Day 0' & used.plot.fil$Subset == x]*100,2),'%',sep=''))+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))+
      guides(fill = F)
    
  })
  plot_grid(plotlist = p_list,ncol = 6)
  ggsave('figS3c_heatmap_all_cell_types_v4.pdf',width = 15.58,height = 5.26)
}

### figS3c V2
{
  ###### all cell types
  UCASpatial.df.long.all <- gather(data = UCASpatial_df,key = 'subset',value = 'Proportion',colnames(UCASpatial_df)[c(22,27:41)])
  
  plot.data.all.raw <- aggregate(x=UCASpatial.df.long.all$Proportion, by=list(UCASpatial.df.long.all$Spot_group,
                                                                              UCASpatial.df.long.all$subset,
                                                                              UCASpatial.df.long.all$Day),mean)
  head(plot.data.all.raw)
  colnames(plot.data.all.raw) <- c('Group',
                                   'Subset','Day','Proportion')
  
  
  plot.data.all.raw$Group <- sapply(as.character(plot.data.all.raw$Group),
                                    function(x){
                                      switch (x,
                                              'Uninjury' = 'Unwounded',
                                              'Wound' = 'Wound center',
                                              'proximal' = 'Wound proximal',
                                              'distal' = 'Wound distal',
                                              'Other' = 'Other'
                                      )
                                    })
  plot.data.all.raw$Group <- factor(plot.data.all.raw$Group,
                                    levels = c('Unwounded','Wound center','Wound proximal','Wound distal','Other'))
  
  used.plot.fil <- plot.data.all.raw[plot.data.all.raw$Group != 'Other' & plot.data.all.raw$Day != 'D10',]
  used.plot.fil$Subset <- factor(used.plot.fil$Subset,
                                 levels = colnames(UCASpatial_df)[c(22,27:41)])
  p_list <- lapply(levels(used.plot.fil$Subset), function(x){
    tmp <- used.plot.fil[used.plot.fil$Subset == x ,]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      facet_wrap(~Subset)+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 12)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12))+
      labs(fill = 'Proportion')+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))
    
  })
  plot_grid(plotlist = p_list,ncol = 3)

  used.plot.fil$Day <- gsub('D','Day ',as.character(used.plot.fil$Day))
  used.plot.fil$Day <- factor(used.plot.fil$Day,
                              levels = c('Day 0','Day 3','Day 7','Day 15'))
  
  p_list <- lapply(levels(used.plot.fil$Subset), function(x){
    tmp <- used.plot.fil[used.plot.fil$Subset == x & used.plot.fil$Day != 'Day 0',]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 12)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1))+
      labs(fill = 'Proportion',title = paste(x,'\n','UW:',round(used.plot.fil$Proportion[used.plot.fil$Day == 'Day 0' & used.plot.fil$Subset == x]*100,2),'%',sep=''))+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))+
      guides(fill = F)
    
  })
  plot_grid(plotlist = p_list,ncol = 6)
  ggsave('figS3c_heatmap_all_cell_subpopulations.pdf',width = 15.58,height = 5.26)
}

#fig S3a
{
  col_all <- c(brewer.pal(n = 9, name = "Reds")[c(3,4,5,6,7,9)], #6
               brewer.pal(n = 9, name = "Greens")[c(2,3,5,7,8)] , # 4
               brewer.pal(n = 9, name = "Blues")[4:9] , #6
               brewer.pal(n = 9, name = "Purples")[2:9], #8
               brewer.pal(n = 9, name = "BrBG"))
  col_new <- col_all[c(1:16,27,19:21)]
  col_new <- c(col_all[c(1:16)],'#FFFF33','#FFCC33','#CC9933','#CC6633',col_all[19:21])

  UCASpatial.df_deld10 <- UCASpatial_df[UCASpatial_df$Day != 'D10',]
  UCASpatial.df_deld10$col[UCASpatial.df_deld10$Day == 'D15'] <- 
    UCASpatial.df_deld10$col[UCASpatial.df_deld10$Day == 'D15']+30
  ggplot2::ggplot() + 
    scatterpie::geom_scatterpie(data = UCASpatial.df_deld10,
                                mapping = ggplot2::aes(x = row,y = col),
                                cols = ct.name, color = NA, 
                                alpha = 1, pie_scale = 0.45,legend_name = 'Type') + 
    scale_fill_manual(values =col_new)+
    theme_base()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(size = 15),
          plot.background =  element_blank()
    )+
    coord_fixed(ratio = 0.5)+
    guides(fill = guide_legend(ncol = 2))

  ggsave('figS3a_UCASpatial_piechart_all_cellltype.pdf',width = 5.37,height = 5.16)
  
}

### figS3d
{
  plot_data_ucaspatial_nw$Method <- 'UCASpatial_ablation'
  ggradar_col <- c(brewer.pal(name = "Set3",n=12)[c(-1:-3)],brewer.pal(name = "Set1",n=9))
  ggradar_col <- ggradar_col[c(1,4,3,6,7,2,9)]
  ggradar_col <- ggradar_col[c(1,7,5,6,3,4,2)]
  plot.deconv.merge.neu <- rbind(plot_data_UCASpatial[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )],
                                 plot_data_ucaspatial_nw[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )],
                                 plot_data_CARD[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )],
                                 plot_data_Cell2location[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil")],
                                 plot_data_RCTD[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )],
                                 plot_data_SPOTlight[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )],
                                 plot_data_stereoscope[,c("row","col","sample","Day" ,"location", "Method"  ,"Neutrophil" )])
  table(plot.deconv.merge.neu$Method)
  head(plot.deconv.merge.neu)
  plot.deconv.merge.neu$Spot_id <- row.names(plot.deconv.merge.neu)
  plot.deconv.merge.neu <- plot.deconv.merge.neu[plot.deconv.merge.neu$Day != 'D10',]
  plot.deconv.merge.neu.long <-melt(plot.deconv.merge.neu,
                                    id.vars = c('Spot_id','row','col','Day','Method'),#需要保留不参与聚合的变量,
                                    measure.vars = c('Neutrophil' ),#用于聚合的变量,
                                    variable.name='Subset',
                                    value.name='Proportion')
  head(plot.deconv.merge.neu.long)
  
  plot.deconv.merge.neu.long$Method <- factor(plot.deconv.merge.neu.long$Method,
                                              levels = c("UCASpatial",'UCASpatial_ablation',"stereoscope","SPOTlight","RCTD",'CARD',"Cell2location" ))
  
  
  plot.deconv.merge.neu.long.SE <- Rmisc::summarySE(plot.deconv.merge.neu.long, measurevar="Proportion", 
                                                    groupvars=c("Method","Day","Subset"))
  
  
  plot.deconv.merge.neu.long$Day <- factor(plot.deconv.merge.neu.long$Day,
                                           levels = c('D0','D3','D7','D15'))
  
  
  plot.deconv.merge.neu.long.SE$Method <- factor(plot.deconv.merge.neu.long.SE$Method,
                                                 levels = c("UCASpatial",'UCASpatial_ablation',"stereoscope","SPOTlight","RCTD",'CARD',"Cell2location" ))
  plot.deconv.merge.neu.long.SE$Day <- gsub('D','Day ',plot.deconv.merge.neu.long.SE$Day)
  plot.deconv.merge.neu.long.SE$Day <- factor(plot.deconv.merge.neu.long.SE$Day,
                                              levels = c('Day 0','Day 3','Day 7','Day 15'))
  
  ggplot(plot.deconv.merge.neu.long.SE,aes(x=Day,y=Proportion*100,
                                           color=Method,group = Method))+
    geom_line(size=1)+
    theme_base(base_size = 16)+
    theme(axis.title.x = element_blank(),
          plot.background = element_blank(),
          plot.title = element_text(face = 'plain',hjust = 0.5,size=16))+
    scale_color_manual(values = ggradar_col)+
    ylab('Average Proportion')+
    labs(title = 'Neutrophil',color=NULL)
  ggsave('FigS3d_line_plot_neutrophil_in_all_methods_adding_ablation.pdf',width = 5.36,height = 2.7)
  
}

### figS3e
{
  plot.deconv.merge.mac <- rbind(plot_data_UCASpatial[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )],
                                 plot_data_ucaspatial_nw[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )],
                                 plot_data_CARD[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )],
                                 plot_data_Cell2location[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage")],
                                 plot_data_RCTD[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )],
                                 plot_data_SPOTlight[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )],
                                 plot_data_stereoscope[,c("row","col","sample","Day" ,"location", "Method"  ,"Macrophage" )])
  table(plot.deconv.merge.mac$Method)
  head(plot.deconv.merge.mac)
  plot.deconv.merge.mac$Spot_id <- row.names(plot.deconv.merge.mac)
  plot.deconv.merge.mac <- plot.deconv.merge.mac[plot.deconv.merge.mac$Day != 'D10',]
  plot.deconv.merge.mac.long <-melt(plot.deconv.merge.mac,
                                    id.vars = c('Spot_id','row','col','Day','Method'),#需要保留不参与聚合的变量,
                                    measure.vars = c('Macrophage' ),#用于聚合的变量,
                                    variable.name='Subset',
                                    value.name='Proportion')
  head(plot.deconv.merge.mac.long)
  
  plot.deconv.merge.mac.long$Method <- factor(plot.deconv.merge.mac.long$Method,
                                              levels = c("UCASpatial",'UCASpatial_ablation',"stereoscope","SPOTlight","RCTD",'CARD',"Cell2location" ))
  
  
  plot.deconv.merge.mac.long.SE <- Rmisc::summarySE(plot.deconv.merge.mac.long, measurevar="Proportion", 
                                                    groupvars=c("Method","Day","Subset"))
  
  
  plot.deconv.merge.mac.long$Day <- factor(plot.deconv.merge.mac.long$Day,
                                           levels = c('D0','D3','D7','D15'))
  
  ggradar_col <- c(brewer.pal(name = "Set3",n=12)[c(-1:-3)],brewer.pal(name = "Set1",n=9))
  ggradar_col <- ggradar_col[c(1,4,3,6,7,2,9)]
  ggradar_col <- ggradar_col[c(1,7,5,6,3,4,2)]
  
  plot.deconv.merge.mac.long.SE$Method <- factor(plot.deconv.merge.mac.long.SE$Method,
                                                 levels = c("UCASpatial",'UCASpatial_ablation',"stereoscope","SPOTlight","RCTD",'CARD',"Cell2location" ))
  plot.deconv.merge.mac.long.SE$Day <- gsub('D','Day ',plot.deconv.merge.mac.long.SE$Day)
  plot.deconv.merge.mac.long.SE$Day <- factor(plot.deconv.merge.mac.long.SE$Day,
                                              levels = c('Day 0','Day 3','Day 7','Day 15'))
  
  ggplot(plot.deconv.merge.mac.long.SE,aes(x=Day,y=Proportion*100,
                                           color=Method,group = Method))+
    geom_line(size=1)+
    theme_base(base_size = 16)+
    theme(axis.title.x = element_blank(),
          plot.background = element_blank(),
          plot.title = element_text(face = 'plain',hjust = 0.5,size=16))+
    scale_color_manual(values = ggradar_col)+
    ylab('Average Proportion')+
    labs(title = 'Macrophage',color=NULL)
  ggsave('FigS3e_line_plot_Macrophage_in_all_methods_adding_ablation.pdf',width = 5.36,height = 2.7)
  
}

### figSI plot all cell types
{
  {
    #UCASpatial
    plot_data_UCASpatial_deld10 <- plot_data_UCASpatial[plot_data_UCASpatial$Day != 'D10',]
    plot_data_UCASpatial_deld10$col[plot_data_UCASpatial_deld10$Day == 'D15'] <- 
      plot_data_UCASpatial_deld10$col[plot_data_UCASpatial_deld10$Day == 'D15']+30
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_UCASpatial_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('UCASpatial_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  
  ### UCASpatial NW
  {
    #UCASpatial NW
    plot_data_UCASpatial_UW_deld10 <- plot_data_ucaspatial_nw[plot_data_ucaspatial_nw$Day != 'D10',]
    plot_data_UCASpatial_UW_deld10$col[plot_data_UCASpatial_UW_deld10$Day == 'D15'] <- 
    plot_data_UCASpatial_UW_deld10$col[plot_data_UCASpatial_UW_deld10$Day == 'D15']+30
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_UCASpatial_UW_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('UCASpatial_UW_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  ###CARD
  {
    #CARD
    plot_data_CARD_deld10 <- plot_data_CARD[plot_data_CARD$Day != 'D10',]
    plot_data_CARD_deld10$col[plot_data_CARD_deld10$Day == 'D15'] <- 
      plot_data_CARD_deld10$col[plot_data_CARD_deld10$Day == 'D15']+30
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_CARD_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('CARD_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  
  ###RCTD
  {
    #RCTD
    
    plot_data_RCTD_deld10 <- plot_data_RCTD[plot_data_RCTD$Day != 'D10',]
    plot_data_RCTD_deld10$col[plot_data_RCTD_deld10$Day == 'D15'] <- 
      plot_data_RCTD_deld10$col[plot_data_RCTD_deld10$Day == 'D15']+30
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_RCTD_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('RCTD_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  
  ### cell2location
  {
    #cell2location
    plot_data_Cell2location_deld10 <- plot_data_Cell2location[plot_data_Cell2location$Day != 'D10',]
    plot_data_Cell2location_deld10$col[plot_data_Cell2location_deld10$Day == 'D15'] <- 
      plot_data_Cell2location_deld10$col[plot_data_Cell2location_deld10$Day == 'D15']+30
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_Cell2location_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('Cell2location_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  
  
  ### spotlight
  {
    #spotlight
    plot_data_SPOTlight_deld10 <- plot_data_SPOTlight[plot_data_SPOTlight$Day != 'D10',]
    plot_data_SPOTlight_deld10$col[plot_data_SPOTlight_deld10$Day == 'D15'] <- 
      plot_data_SPOTlight_deld10$col[plot_data_SPOTlight_deld10$Day == 'D15']+30
    
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_SPOTlight_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('Spotlight_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
  }
  
  ### stereoscope
  {
    #stereoscope
    plot_data_stereoscope_deld10 <- plot_data_stereoscope[plot_data_stereoscope$Day != 'D10',]
    plot_data_stereoscope_deld10$col[plot_data_stereoscope_deld10$Day == 'D15'] <- 
      plot_data_stereoscope_deld10$col[plot_data_stereoscope_deld10$Day == 'D15']+30
    
    p_list <- lapply(ct.name, function(x){
      plot_tmp <- plot_data_stereoscope_deld10[,c('row','col',x)]
      colnames(plot_tmp)[3] <- 'Type'
      ggplot(data = plot_tmp)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15),
              plot.background =  element_blank()
        )+
        geom_point_rast( aes(x=row,y=col,color=Type),
                         size=1,shape=19)+
        #scale_color_gradientn(colours = plasma(8),
        #                     guide = guide_colorbar(frame.colour = 'black',barwidth = 1))+
        scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                              guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                              #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                              na.value = '#DE1A37',
                              labels = scales::label_percent()
        )+
        labs(title = x,color = NULL)+
        coord_fixed(ratio = 0.5)
    })
    pdf('stereoscope_all_cell_tpye.pdf',width = 25,height = 20)
    plot_grid(plotlist = p_list,nrow = 5,align = 'hv')
    dev.off()
    
  }
  
  
}
plot_tmp <- plot_data_stereoscope_deld10[,c('row','col',x)]
colnames(plot_tmp)[3] <- 'Type'
ggplot(data = plot_data_stereoscope_deld10)+
  theme_base()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 15),
        plot.background =  element_blank()
  )+
  geom_point_rast( aes(x=row,y=col,color=Macrophage),
                   size=1,shape=19)+
  scale_color_gradientn(colors = c("#90BCD9","#C6E0ED","white","#FDBF6F","#DE1A37"),
                        guide = guide_colorbar(frame.colour = 'black',barwidth = 0.8),
                        #limits = c(0, quantile(plot_tmp$Type, 0.99)),
                        na.value = '#DE1A37',
                        labels = scales::label_percent()
  )
  scale_color_viridis_c(labels = scales::label_percent())+
  coord_fixed(ratio = 0.5)

### cluster amrkers
{
  source('/data/xy/scripts/eMark_v2.R')
  
  cluster_markers <- Seurat::FindAllMarkers(
    object = MRL.sc.data,
    only.pos = T,assay = 'RNA',slot = 'data',min.pct = 0.2,min.diff.pct = 0.1
  )
  cluster_markers_wetent <- CalculateWeightedEntropy(MRL.sc.data,cluster_markers,assay = 'RNA',
                                                     slot = 'data',unit = 'log2')
  saveRDS(cluster_markers_wetent,"w.ent_cluster_markers.rds")
  WriteXLS::WriteXLS(x = cluster_markers_wetent,ExcelFileName = 'Supplementary Table murine wound healing cell type markers in scRNA.xls',SheetNames = 'murine wound healing')
}

#### cell2location 
{
  ###### all cell types
  cell2location_df <- cbind(cell2location_df[,1:21],deconv.list$Cell2location)
  cell2location.df.long.all <- gather(data = cell2location_df,key = 'subset',value = 'Proportion',colnames(cell2location_df)[c(22:41)])
  
  plot.data.all.raw <- aggregate(x=cell2location.df.long.all$Proportion, by=list(cell2location.df.long.all$Spot_group,
                                                                                 cell2location.df.long.all$subset,
                                                                                 cell2location.df.long.all$Day),mean)
  head(plot.data.all.raw)
  colnames(plot.data.all.raw) <- c('Group',
                                   'Subset','Day','Proportion')
  
  
  plot.data.all.raw$Group <- sapply(as.character(plot.data.all.raw$Group),
                                    function(x){
                                      switch (x,
                                              'Uninjury' = 'Unwounded',
                                              'Wound' = 'Wound center',
                                              'proximal' = 'Wound proximal',
                                              'distal' = 'Wound distal',
                                              'Other' = 'Other'
                                      )
                                    })
  plot.data.all.raw$Group <- factor(plot.data.all.raw$Group,
                                    levels = c('Unwounded','Wound center','Wound proximal','Wound distal','Other'))
  
  used.plot.fil <- plot.data.all.raw[plot.data.all.raw$Group != 'Other' & plot.data.all.raw$Day != 'D10',]
  used.plot.fil$Subset <- factor(used.plot.fil$Subset,
                                 levels = colnames(cell2location_df)[c(22:41)])
  used.plot.fil$Day <- gsub('D','Day ',as.character(used.plot.fil$Day))
  used.plot.fil$Day <- factor(used.plot.fil$Day,
                              levels = c('Day 0','Day 3','Day 7','Day 15'))
  
  p_list <- lapply(levels(used.plot.fil$Subset), function(x){
    tmp <- used.plot.fil[used.plot.fil$Subset == x & used.plot.fil$Day != 'Day 0',]
    ggplot(tmp,aes(x=Day,y=Group,fill=round(Proportion*100,2)))+
      #geom_tile(colour = 'black',size=0.5)+
      geom_raster()+
      scale_fill_gradientn(colors = c("#90BCD9","white","#DE1A37"))+
      theme_base(base_size = 12)+
      theme(plot.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(face = "plain",size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1))+
      labs(fill = 'Proportion',title = paste(x,'\n','UW:',round(used.plot.fil$Proportion[used.plot.fil$Day == 'Day 0' & used.plot.fil$Subset == x]*100,2),'%',sep=''))+
      geom_text(data = tmp,
                aes(x=Day,y=Group,label= round(Proportion*100,2)))+
      guides(fill = F)
    
  })
  plot_grid(plotlist = p_list,ncol = 6)
  ggsave('cell2location_heatmap_all_cell_types.pdf',width = 15.58,height = 6.9)
}
  
### figS9e cell community
{
  head(plot_data_espanno_sub)
  ### wound center spot
  used.data <- plot_data_UCASpatial[row.names(MRL.coord)[MRL.coord$Spot_group_v2 == 'wound'],]
  used.data <- used.data[used.data$Day != 'D10',]
  used.data <- used.data[,17:36]
  pval_mat_2 <- matrix(NA, 20, 20)
  for(i in 1:20) {
    for(j in 1:20) {
      temp <- cor.test(used.data[,i], used.data[,j],method = "spearman")
      pval_mat_2[i,j] <- temp$p.value
    }
  }
  
  sig_mat_2 <- ifelse(pval_mat_2 < 0.001, "***", ifelse(pval_mat_2 < 0.01, "**", ifelse(pval_mat_2 < 0.05, "*", "")))
  rownames(sig_mat_2) = colnames(used.data)
  colnames(sig_mat_2) = colnames(used.data)
  PCC_matrix <- cor(used.data)
  sig_mat_2 <- sig_mat_2[rownames(PCC_matrix),colnames(PCC_matrix)]
  pdf('heatmap_cell_community_wound_center.pdf',width = 8.36,height = 5.36)
  pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                     cellwidth = 10,cellheight = 10,
                     fontsize = 10,
                     breaks = seq(-1, 1, length.out = 100),
                     border_color = "black",
                     display_numbers = sig_mat_2,
                     cutree_cols = 5,cutree_rows = 5,
                     clustering_method = "complete")
  dev.off()
  ##### D7 wound center
  used.data <- plot_data_UCASpatial[row.names(MRL.coord)[MRL.coord$Spot_group_v2 == 'wound' & MRL.coord$Day == 'D7'],]
  used.data <- used.data[used.data$Day != 'D10',]
  used.data <- used.data[,17:36]
  pval_mat_2 <- matrix(NA, 20, 20)
  for(i in 1:20) {
    for(j in 1:20) {
      temp <- cor.test(used.data[,i], used.data[,j],method = "spearman")
      pval_mat_2[i,j] <- temp$p.value
    }
  }
  
  sig_mat_2 <- ifelse(pval_mat_2 < 0.001, "***", ifelse(pval_mat_2 < 0.01, "**", ifelse(pval_mat_2 < 0.05, "*", "")))
  rownames(sig_mat_2) = colnames(used.data)
  colnames(sig_mat_2) = colnames(used.data)
  PCC_matrix <- cor(used.data)
  sig_mat_2 <- sig_mat_2[rownames(PCC_matrix),colnames(PCC_matrix)]
  PCC_matrix <- PCC_matrix[-9,-9]
  sig_mat_2 <- sig_mat_2[-9,-9]
  pdf('heatmap_cell_community_D7_wound_center.pdf',width = 8.36,height = 5.36)
  pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                     cellwidth = 10,cellheight = 10,
                     fontsize = 10,
                     breaks = seq(-1, 1, length.out = 100),
                     border_color = "black",
                     display_numbers = sig_mat_2,
                     cutree_cols = 5,cutree_rows = 5,
                     clustering_method = "complete")
  dev.off()
  
  ### total spot
  head(plot_data_UCASpatial)
  used.data <- plot_data_UCASpatial[plot_data_UCASpatial$Day != 'D10',]
  used.data <- used.data[,17:36]
  pval_mat_2 <- matrix(NA, 20, 20)
  for(i in 1:20) {
    for(j in 1:20) {
      temp <- cor.test(used.data[,i], used.data[,j],method = "spearman")
      pval_mat_2[i,j] <- temp$p.value
    }
  }
  
  sig_mat_2 <- ifelse(pval_mat_2 < 0.001, "***", ifelse(pval_mat_2 < 0.01, "**", ifelse(pval_mat_2 < 0.05, "*", "")))
  rownames(sig_mat_2) = colnames(used.data)
  colnames(sig_mat_2) = colnames(used.data)
  PCC_matrix <- cor(used.data)
  sig_mat_2 <- sig_mat_2[rownames(PCC_matrix),colnames(PCC_matrix)]
  pdf('heatmap_cell_community_total_spot.pdf',width = 8.36,height = 5.36)
  pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                     cellwidth = 10,cellheight = 10,
                     fontsize = 10,
                     breaks = seq(-1, 1, length.out = 100),
                     border_color = "black",
                     display_numbers = sig_mat_2,
                     cutree_cols = 5,cutree_rows = 5,
                     clustering_method = "complete")
  dev.off()
  
  ##### divide region: wound center, proximal, distal
  used.data <- plot_data_UCASpatial[plot_data_UCASpatial$Day != 'D10',]
  used.data$tmp_group <- paste(used.data$Day,used.data$Spot_group_v2,sep = ' ')
  used.data <- used.data[,c(17:36,40)]
  used.data %>%
    group_by(tmp_group) %>%
    summarize_all(mean) -> used.data
  used.data <- used.data[-c(2,9),] %>% data.frame()
  row.names(used.data) <- used.data$tmp_group
  used.data$tmp_group <- NULL
  colnames(used.data) <- colnames(plot_data_UCASpatial)[17:36]
  pval_mat_2 <- matrix(NA, 20, 20)
  for(i in 1:20) {
    for(j in 1:20) {
      temp <- cor.test(used.data[,i], used.data[,j],method = "spearman")
      pval_mat_2[i,j] <- temp$p.value
    }
  }
  
  sig_mat_2 <- ifelse(pval_mat_2 < 0.001, "***", ifelse(pval_mat_2 < 0.01, "**", ifelse(pval_mat_2 < 0.05, "*", "")))
  rownames(sig_mat_2) = colnames(used.data)
  colnames(sig_mat_2) = colnames(used.data)
  PCC_matrix <- cor(used.data)
  sig_mat_2 <- sig_mat_2[rownames(PCC_matrix),colnames(PCC_matrix)]
  pdf('heatmap_cell_community_divide_region_splot_day.pdf',width = 8.36,height = 5.36)
  pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                     cellwidth = 10,cellheight = 10,
                     fontsize = 10,
                     breaks = seq(-1, 1, length.out = 100),
                     border_color = "black",
                     display_numbers = sig_mat_2,
                     cutree_cols = 5,cutree_rows = 5,
                     clustering_method = "complete")
  dev.off()
  
  
  #######
  ##### divide region: wound center, proximal, distal
  used.data <- plot_data_UCASpatial[plot_data_UCASpatial$Day != 'D10',]
  used.data <- used.data[,c(16:36)]
  used.data %>%
    group_by(Spot_group_v2) %>%
    summarize_all(mean) -> used.data
  used.data <- used.data[-5,] %>% data.frame()
  row.names(used.data) <- used.data$Spot_group_v2
  used.data$Spot_group_v2 <- NULL
  colnames(used.data) <- colnames(plot_data_UCASpatial)[17:36]
  pval_mat_2 <- matrix(NA, 20, 20)
  for(i in 1:20) {
    for(j in 1:20) {
      temp <- cor.test(used.data[,i], used.data[,j],method = "spearman")
      pval_mat_2[i,j] <- temp$p.value
    }
  }
  
  sig_mat_2 <- ifelse(pval_mat_2 < 0.001, "***", ifelse(pval_mat_2 < 0.01, "**", ifelse(pval_mat_2 < 0.05, "*", "")))
  rownames(sig_mat_2) = colnames(used.data)
  colnames(sig_mat_2) = colnames(used.data)
  PCC_matrix <- cor(used.data)
  sig_mat_2 <- sig_mat_2[rownames(PCC_matrix),colnames(PCC_matrix)]
  pdf('heatmap_cell_community_divide_region_ward.D2.pdf',width = 8.36,height = 5.36)
  pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),
                     cellwidth = 10,cellheight = 10,
                     fontsize = 10,
                     breaks = seq(-1, 1, length.out = 100),
                     border_color = "black",
                     display_numbers = sig_mat_2,
                     cutree_cols = 4,cutree_rows = 4,
                     clustering_method = "ward.D2")
  dev.off()
}

### fig3b and figS3 macrophage
{
    p <- DimPlot(MRL_scRNA,reduction = 'umap')
    head(p$data)
    all_cluster.umap <- p$data
    head(all_cluster.umap)
    colnames(all_cluster.umap) <- c('UMAP_1','UMAP_2','Subset')
    
    col_all <- c(brewer.pal(n = 9, name = "Reds")[c(3,5,6,7,8)], #5
                 brewer.pal(n = 9, name = "Greens")[c(2,4,6,8)] , # 4
                 brewer.pal(n = 9, name = "Blues")[c(4,5,6,7,8)] , #5
                 brewer.pal(n = 9, name = "Purples")[c(2,4,6,8)], #4
                 brewer.pal(n = 9, name = "BrBG"))
    
    
    ggplot(all_cluster.umap,aes(x=UMAP_1,y=UMAP_2,fill=Subset))+
      geom_point_rast(shape=21,color='black',size=0.8,stroke = 0)+
      theme_base(base_size = 16)+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.background = element_blank(),
            panel.border = element_blank(),
            panel.background	 = element_blank())+
      guides(fill=guide_legend(ncol=2,override.aes = list(size = 3)))+
      scale_fill_manual(values = col_all)+
      coord_fixed(ratio = 0.8)+
      ggtitle(' ')+
      labs(x='UMAP_1',y='UMAP_2',fill=NULL)
    ggsave('/data/zhangyw/project/eSpanno/figure_final/MRL_umap_modify_all_clusters_v2.pdf',width = 10,height = 6)
    DimPlot(mac,cols = brewer.pal(9,'Set1')[c(1,2,4,5,7)],reduction = 'umap',raster = T,pt.size = 2)+
      theme_base(base_size = 16)+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.background = element_blank(),
            panel.border = element_blank(),
            panel.background	 = element_blank())+
      coord_fixed(ratio = 1.5)
    ggsave('/data/zhangyw/project/eSpanno/figure_final/umap_mac_raster_v2.pdf',width = 7,height = 2.94)
    
    FeaturePlot(mac,features = c('Ly6c2','Cd163','Spp1','H2-Aa','Il1b','Arg1'),)
    ##### macrophage features 
    FeaturePlot(mac,features = c('Ly6c2','Cd163','Spp1','H2-Aa','Il1b','Arg1'),ncol = 3,
                reduction = 'umap',min.cutoff = 'q5',max.cutoff = 'q95',raster = T,
                pt.size = 6,cols = rev(colorRampPalette(c("firebrick3", "white", "navy"))(100)))
    ggsave('/data/zhangyw/project/eSpanno/figure_final/featureplot_for_macrophages.pdf',width = 11,height = 5.74)
    
}
