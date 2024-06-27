my_spatial_featureplot <- function(spatial_data,data_coord,features,
                                   min.cutoff = NA,
                                   max.cutoff = NA,
                                   label.size = 5,
                                   ncol = 1){
  
  data <- FetchData(object = spatial_data, vars = features)
  features <- colnames(x = data)
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  
  data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
    data.feature <- as.vector(x = data[, index])
    min.use <- SetQuantile(cutoff = min.cutoff[index], 
                           data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[index], 
                           data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    return(data.feature)
  })
  colnames(x = data) <- features
  rownames(x = data) <- Cells(x = spatial_data)
  features <- colnames(x = data)
  plot_data <- cbind(data_coord, 
                     data)
  plot <- list()
  for (j in 1:length(x = features)) {
    tmp <- plot_data[,c('row','col','sample','sample_detail','mice',features[j])]
    colnames(tmp) <- c('row','col','sample','sample_detail','mice','gene')
    plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=gene))+
      geom_point(shape = 21)+
      facet_wrap(~mice)+
      theme_base()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(color='black',size=0.2),
            axis.ticks.y.right = element_line(color = 'black'),
            panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            strip.text = element_text(size = 15))+
      labs(fill=NULL) +
      scale_fill_gradientn(colours=colorRampPalette(rev(x = brewer.pal(n = 10, name = "Spectral")))(n=100),
                           #limits = c(0, quantile(tmp$gene, 0.99)+0.01),
                           na.value = "#9E0142")+
      coord_fixed(ratio = 0.5)+
      annotate('text',x=-15,y=12,label='D0',hjust=0,size=5)+
      annotate('text',x=-15,y=-16,label='D3',hjust=0,size=5)+
      annotate('text',x=-15,y=-45,label='D7',hjust=0,size=5)+
      annotate('text',x=-15,y=-75,label='D10',hjust=0,size=5)+
      annotate('text',x=-15,y=-106,label='D15',hjust=0,size=5)+
      geom_hline(yintercept = c(-10,-40,-70,-100),colour="grey", linetype="dashed")+
      ggtitle(features[j])+
      theme(plot.background = element_blank())
    
  }
  plot_grid(plotlist = plot,ncol = ncol)
}

SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}


my_spatialplot <- function(spatial.metadata,interset.cluster=NULL,
                           label.size=5,
                           ncol = 1,cols = NULL,
						               fill = cluster,
						               legend_title = NULL){

  if(is.null(interset.cluster)){
    plot <- list()
    if(is.null(cols)){
      plot[[1]] <- ggplot(data = spatial.metadata, aes(x=row,y=col,fill = fill))+
        geom_point(colour = "black",shape = 21)+
        facet_wrap(~mice)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line(color='black',size=0.2),
              axis.ticks.y.right = element_line(color = 'black'),
              panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              strip.text = element_text(size = 15))+
        coord_fixed(ratio = 0.5)+
        annotate('text',x=-15,y=12,label='D0',hjust=0,size=5)+
        annotate('text',x=-15,y=-16,label='D3',hjust=0,size=5)+
        annotate('text',x=-15,y=-45,label='D7',hjust=0,size=5)+
        annotate('text',x=-15,y=-75,label='D10',hjust=0,size=5)+
        annotate('text',x=-15,y=-106,label='D15',hjust=0,size=5)+
        geom_hline(yintercept = c(-10,-40,-70,-100),colour="grey", linetype="dashed")+
        theme(plot.background = element_blank())+
        scale_fill_manual(values = cols)+
        guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4))) 
    }else{
      plot[[1]] <- ggplot(data = spatial.metadata, aes(x=row,y=col,fill = fill))+
        geom_point(colour = "black",shape = 21)+
        facet_wrap(~mice)+
        theme_base()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line(color='black',size=0.2),
              axis.ticks.y.right = element_line(color = 'black'),
              panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              strip.text = element_text(size = 15))+
        coord_fixed(ratio = 0.5)+
        annotate('text',x=-15,y=12,label='D0',hjust=0,size=5)+
        annotate('text',x=-15,y=-16,label='D3',hjust=0,size=5)+
        annotate('text',x=-15,y=-45,label='D7',hjust=0,size=5)+
        annotate('text',x=-15,y=-75,label='D10',hjust=0,size=5)+
        annotate('text',x=-15,y=-106,label='D15',hjust=0,size=5)+
        geom_hline(yintercept = c(-10,-40,-70,-100),colour="grey", linetype="dashed")+
        theme(plot.background = element_blank())+
        scale_fill_manual(values = cols)+
        guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
    }
    
      
  }else{
    plot <- list()
    for (j in 1:length(x = interset.cluster)) {
      tmp <- spatial.metadata[,c('row','col','sample','sample_detail','mice',interset.cluster[j])]
      colnames(tmp) <- c('row','col','sample','sample_detail','mice','cluster')
      if(is.null(cols)){
        plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=fill))+
          geom_point(colour = "black",shape = 21)+
          facet_wrap(~mice)+
          theme_base()+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_line(color='black',size=0.2),
                axis.ticks.y.right = element_line(color = 'black'),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                strip.text = element_text(size = 15))+
          labs(fill=interset.cluster[j]) +
          coord_fixed(ratio = 0.5)+
          annotate('text',x=-15,y=12,label='D0',hjust=0,size=5)+
          annotate('text',x=-15,y=-16,label='D3',hjust=0,size=5)+
          annotate('text',x=-15,y=-45,label='D7',hjust=0,size=5)+
          annotate('text',x=-15,y=-75,label='D10',hjust=0,size=5)+
          annotate('text',x=-15,y=-106,label='D15',hjust=0,size=5)+
          geom_hline(yintercept = c(-10,-40,-70,-100),colour="grey", linetype="dashed")+
          ggtitle(interset.cluster[j])+
          theme(plot.background = element_blank())+
          scale_fill_manual(values = cols)+
          guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
      }else{
        plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=fill))+
          geom_point()+
          facet_wrap(~mice)+
          theme_base()+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_line(color='black',size=0.2),
                axis.ticks.y.right = element_line(color = 'black'),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                strip.text = element_text(size = 15))+
          labs(fill=interset.cluster[j]) +
          coord_fixed(ratio = 0.5)+
          annotate('text',x=-15,y=12,label='D0',hjust=0,size=5)+
          annotate('text',x=-15,y=-16,label='D3',hjust=0,size=5)+
          annotate('text',x=-15,y=-45,label='D7',hjust=0,size=5)+
          annotate('text',x=-15,y=-75,label='D10',hjust=0,size=5)+
          annotate('text',x=-15,y=-106,label='D15',hjust=0,size=5)+
          geom_hline(yintercept = c(-10,-40,-70,-100),colour="grey", linetype="dashed")+
          ggtitle(interset.cluster[j])+
          theme(plot.background = element_blank())+
          scale_fill_manual(values = cols)+
          guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
      }
      
    }
  }
  
  plot_grid(plotlist = plot,ncol = ncol)
}



my_spatial_eweideplot_scSample <- function(spatial_data,data_coord,features,
                                           min.cutoff = NA,
                                           max.cutoff = NA,
                                           label.size = 1,
                                           ncol = 1,
                                           pt.size = 2){
  
  data <- FetchData(object = spatial_data, vars = features)
  features <- colnames(x = data)
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  
  data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
    data.feature <- as.vector(x = data[, index])
    min.use <- SetQuantile(cutoff = min.cutoff[index], 
                           data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[index], 
                           data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    return(data.feature)
  })
  colnames(x = data) <- features
  rownames(x = data) <- Cells(x = spatial_data)
  features <- colnames(x = data)
  plot_data <- cbind(data_coord, 
                     data)
  plot <- list()
  for (j in 1:length(x = features)) {
    tmp <- plot_data[,c('row','col','sample','sample_detail','mice',features[j])]
    colnames(tmp) <- c('row','col','sample','sample_detail','mice','gene')
    plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=gene))+
      geom_point(shape = 21,size=pt.size)+
      theme_cowplot()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15))+
      labs(fill=NULL) +
      #scale_fill_gradientn(colours=colorRampPalette(rev(x = brewer.pal(n = 10, name = "Spectral")))(n=100),
      #limits = c(0, quantile(tmp$gene, 0.99)+0.01),
      #                     na.value = "#9E0142")+
      scale_fill_gradientn(colours = colorRampPalette(rev(x = brewer.pal(n = 10, name = "Spectral")))(n=100), 
                           limits = c(0, 1))+
      coord_fixed(ratio = 0.5)+
      ggtitle(features[j])+
      theme(plot.background = element_blank())
    
  }
  plot_grid(plotlist = plot,ncol = ncol)
}



my_spatialplot_scSample <- function(spatial.metadata,interset.cluster=NULL,
                                    label.size=5,
                                    ncol = 1,cols = NULL,
                                    fill = cluster,
                                    pt.size=1,
                                    legend_title = NULL){
  
  if(is.null(interset.cluster)){
    plot <- list()
    if(is.null(cols)){
      plot[[1]] <- ggplot(data = spatial.metadata, aes(x=row,y=col,fill = fill))+
        geom_point(colour = "black",shape = 21,size=pt.size)+
        theme_cowplot()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15))+
        coord_fixed(ratio = 0.5)+
        theme(plot.background = element_blank())+
        scale_fill_manual(values = cols)+
        guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4))) 
    }else{
      plot[[1]] <- ggplot(data = spatial.metadata, aes(x=row,y=col,fill = fill))+
        geom_point(colour = "black",shape = 21,size=pt.size)+
        theme_cowplot()+
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              strip.text = element_text(size = 15))+
        coord_fixed(ratio = 0.5)+
        theme(plot.background = element_blank())+
        scale_fill_manual(values = cols)+
        guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
    }
    
    
  }else{
    plot <- list()
    for (j in 1:length(x = interset.cluster)) {
      tmp <- spatial.metadata[,c('row','col','sample','sample_detail','mice',interset.cluster[j])]
      colnames(tmp) <- c('row','col','sample','sample_detail','mice','cluster')
      if(is.null(cols)){
        plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=fill))+
          geom_point(colour = "black",shape = 21,size=pt.size)+
          theme_cowplot()+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                strip.text = element_text(size = 15))+
          labs(fill=interset.cluster[j]) +
          coord_fixed(ratio = 0.5)+
          ggtitle(interset.cluster[j])+
          theme(plot.background = element_blank())+
          scale_fill_manual(values = cols)+
          guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
      }else{
        plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=fill))+
          geom_point(colour = "black",shape = 21,size=pt.size)+
          theme_cowplot()+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                strip.text = element_text(size = 15))+
          labs(fill=interset.cluster[j]) +
          coord_fixed(ratio = 0.5)+
          ggtitle(interset.cluster[j])+
          theme(plot.background = element_blank())+
          scale_fill_manual(values = cols)+
          guides(fill=guide_legend(title = legend_title,override.aes = list(size = 4)))
      }
      
    }
  }
  
  plot_grid(plotlist = plot,ncol = ncol)
}

my_spatial_featureplot_scSample <- function(spatial_data,data_coord,features,
                                            min.cutoff = NA,
                                            max.cutoff = NA,
                                            label.size = 1,
                                            ncol = 1,
                                            cols = colorRampPalette(rev(x = brewer.pal(n = 10, name = "Spectral")))(n=100),
                                            pt.size = 2){
  
  data <- FetchData(object = spatial_data, vars = features)
  features <- colnames(x = data)
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  
  data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
    data.feature <- as.vector(x = data[, index])
    min.use <- SetQuantile(cutoff = min.cutoff[index], 
                           data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[index], 
                           data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    return(data.feature)
  })
  colnames(x = data) <- features
  rownames(x = data) <- Cells(x = spatial_data)
  features <- colnames(x = data)
  plot_data <- cbind(data_coord, 
                     data)
  plot <- list()
  for (j in 1:length(x = features)) {
    tmp <- plot_data[,c(colnames(plot_data)[1:ncol(data_coord)],features[j])]
    colnames(tmp) <- c(colnames(plot_data)[1:ncol(data_coord)],'gene')
    plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=gene))+
      geom_point(shape = 21,size=pt.size)+
      theme_cowplot()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15))+
      labs(fill=NULL) +
      scale_fill_gradientn(colours=cols,
                           #limits = c(0, quantile(tmp$gene, 0.99)+0.01),
                           na.value = "#9E0142")+
      coord_fixed(ratio = 0.5)+
      ggtitle(features[j])+
      theme(plot.background = element_blank())
    
  }
  plot_grid(plotlist = plot,ncol = ncol)
}


my_spatial_devconplot_scSample <- function(spatial_data,data_coord,features,
                                            min.cutoff = NA,
                                            max.cutoff = NA,
                                            label.size = 1,
                                            assay = 'eSpanno',
                                            ncol = 1,
                                            pt.size = 2){
  ####
  DefaultAssay(spatial_data) <- assay
  if(length(features) > 1){
    data <- spatial_data@assays[[assay]]@counts[features,] %>% as.matrix() %>% data.frame() 
    rownames(data) <- features
  }else{
    data <- spatial_data@assays[[assay]]@counts[features,] %>% as.matrix()  %>% t() %>% data.frame( )
    rownames(data) <- features
  }
  colnames(data) <- colnames(spatial_data)
  data <- t(data)
  data <- data.frame(data)
  colnames(data) <- features
  #return(data)

  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  print(features)
  data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
    data.feature <- as.vector(x = data[, index])
    min.use <- SetQuantile(cutoff = min.cutoff[index], 
                           data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[index], 
                           data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    return(data.feature)
  })
  colnames(x = data) <- features
  rownames(x = data) <- Cells(x = spatial_data)
  features <- colnames(x = data)
  plot_data <- cbind(data_coord, 
                     data)
  plot <- list()
  for (j in 1:length(x = features)) {
    tmp <- plot_data[,c(colnames(plot_data)[1:ncol(data_coord)],features[j])]
    colnames(tmp) <- c(colnames(plot_data)[1:ncol(data_coord)],'gene')
    plot[[j]] <- ggplot(data = tmp, aes(x=row,y=col,fill=gene))+
      geom_point(shape = 21,size=pt.size)+
      theme_cowplot()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = 15))+
      labs(fill=NULL) +
      scale_fill_gradientn(colours=colorRampPalette(rev(x = brewer.pal(n = 10, name = "Spectral")))(n=100),
                           #limits = c(0, quantile(tmp$gene, 0.99)+0.01),
                           na.value = "#9E0142")+
      coord_fixed(ratio = 0.5)+
      ggtitle(features[j])+
      theme(plot.background = element_blank())
    
  }
  plot_grid(plotlist = plot,ncol = ncol)
}
