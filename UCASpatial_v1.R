## varible explaination:
# sc_ref: reference single cell RNA-seq data (Seurat object)
# st_vis: spatial transcriptomic data (Seurat object)
# clust_vr: cluster information of single cell reference data
# assay(optional): set the default assay for sc_ref. By default 'RNA'.
# slot(optional): set the default slot for sc_ref. By default 'counts'.
# cluster_markers(optional): give the marker genes for each cluster of sc_ref. If is NULL, it will be calculated use the default parameters. By default NULL.
# normalize: method for normalizing the count matrix: uv (unit variance), raw (no transformation applied). By default 'uv'.
# n_cluster(optional): set the maximum size for each cluster after downsample. By default '100'.
# n_top(optional): set the n_top marker gene used for each cluster. By default 'NULL'.
# min_cont(optional): set the cutoff threshold. By default '0.01'.


eWEIDE_deconv <- function (sc_ref, st_vis, clust_vr, assay = "RNA", slot = "data",output_path = NULL,
                           cluster_markers = NULL,min.pct = 0.2,logfc.threshold = 0.25,min.diff.pct = 0.1,
                           normalize = 'uv',downsample_n = 1 , n_cluster = 100,n_top = NULL,
                           remove.RPL=F,remove.MT=F,
                           cos.filter = T,cos.mu=1,cos.n_genes_user=900,marker.slot = 'data',
                           min_cont = 0.01,unit = "log2",random.seed = 10000,meta.filter = T,nmf.tol=1e-04,
                           meta.assay = 'integrated',meta.ndims = 30,meta.resolution = 100,meta.purity = 0.95,
                           ent.filter.threshold = 0.5,cos.filter.threshold = 0.05,weight.filter.threshold = 0.2) 
{
  # Step0: Preprocess
  cat("Load required packages...\n")
  suppressMessages(require(NMF))
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))
  suppressMessages(require(dplyr))
  suppressMessages(require(edgeR))
  suppressMessages(require(purrr))
  suppressMessages(require(tibble))
  suppressMessages(require(nnls))
  suppressMessages(require(SPOTlight))
  suppressMessages(require(reshape2))
  suppressMessages(require(RcppML))
  suppressMessages(require(COSG))
  
  cat("Step0    Check the variables............................\n")
  {
    if (is(sc_ref) != "Seurat") 
      stop("ERROR: sc_ref must be a Seurat object!\n")
    if (is.null(sc_ref@meta.data[[clust_vr]])) 
      stop("ERROR: clust_vr must be the cluster information of sc_ref!\n(path:sc_ref@meta.data[[clust_vr]])\n")
    if(is.null(output_path))
    {
      output_path <- getwd()
      cat('Warning: Using current path: \"',output_path,'\" as the \'output_path\'...\n',sep = '')
    }
    if (!is.numeric(n_cluster)) 
      stop("ERROR: n_cluster must be an integer!\n(n_cluster: maximum size for each cluster after downsample)\n")
    if (!is.character(normalize)) 
      stop("ERROR: normalize must be a character string!\n(method for normalizing the count matrix: uv (unit variance), raw (no transformation applied))\n")
    if(!is.null(n_top))
      if (!is.numeric(n_top)) 
        stop("ERROR: n_top must be a an integer!\n(set the n_top marker gene used for each cluster)\n")
    if (!is.numeric(min_cont)) 
    {
      cat("Warning message: min_cont should be numeric (instead with default 0.01).\n")
      min_cont <- 0.01
    }
  }
  cat("...........\nLoad the spatial expression matrix (row counts):st_vis@assays$Spatial@counts\n")
  st_vis_matr <- st_vis@assays$Spatial@counts
  if(remove.RPL)
  {
    cat("...........\nFilter Ribosome genes\n")
    ribosome_genes <- grep("^RPL|^RPS", rownames(st_vis), value = TRUE)
    st_vis_filtered <- st_vis[!rownames(st_vis) %in% ribosome_genes,]
    st_vis <- st_vis_filtered
    st_vis_matr <- st_vis@assays$Spatial@counts
  }
  if(remove.MT)
  {
    cat("...........\nFilter mitochondrial genes\n")
    mitochondrial_genes <- grep("^MT-", rownames(st_vis), value = TRUE)
    st_vis_filtered <- st_vis[!rownames(st_vis) %in% mitochondrial_genes, ]
    st_vis <- st_vis_filtered
    st_vis_matr <- st_vis@assays$Spatial@counts
  }
  
  DefaultAssay(sc_ref) <- assay
  #### Step0: Preprocess the sc_ref data. ####
  cat("Step0    Preprocess the sc_ref data.....................\n")
  ## Step0.1  Meta.cell filter
  cat("Step0.1  Meta.cell filter...............................\n")
  if(meta.filter == T)
  {
    sc_ref <- meta_filter(sc_ref = sc_ref,meta.assay = meta.assay,assay = assay,output_path = output_path,clust_vr,
                          meta.purity = meta.purity,meta.resolution = meta.resolution,meta.ndims = meta.ndims)
  }
  ## Step0.2  Downsample the sc_ref data by min.varience
  if(downsample_n != 0){
    cat("Step0.2  Downsample the sc_ref data by min.varience.....\n")
    sc_ref_down <- eWEIDE_downsample(sc_ref = sc_ref,clust_vr = clust_vr,random.seed = random.seed,
                                     downsample_n = downsample_n,n_cluster = n_cluster,assay = assay,slot = 'counts')
    cat(paste("Auto save the downsampled sc_ref under the path:\n",output_path,"/sc_ref_down.rds\n",sep = ""))
    saveRDS(sc_ref_down,paste(output_path,"/sc_ref_down.rds",sep = ""))
    sc_ref <- sc_ref_down
  }
  #### Step1: Calculate the markers and add weights. ####
  cat("Step1    Calculate the markers and add weights..........\n")
  if (is.null(cluster_markers))
  {
    # cat("Warrning message: Did not provide the marker gene sets.")
    cat("Step1.1  Calculate the marker genes.....................\n")
    cat(paste("Change the idents of sc_ref into 'sc_ref@meta.data$",clust_vr,"'...\n",sep = ""))
    Idents(sc_ref) <- sc_ref@meta.data[[clust_vr]]
    cluster_markers <- Seurat::FindAllMarkers(
      object = sc_ref,
      assay = assay,
      slot = marker.slot,
      min.pct = min.pct,
      only.pos = TRUE,
      logfc.threshold = logfc.threshold,
      min.diff.pct = min.diff.pct
    )
    cat(paste("Auto save the marker genes under the path:\n",output_path,"/cluster_markers.rds\n",sep = ""))
    saveRDS(cluster_markers,paste(output_path,"/cluster_markers.rds",sep = ""))
  }
  # if (is.null(cluster_markers$gene) | is.null(cluster_markers$cluster) | is.null(cluster_markers$p_val)) 
  #   stop("ERROR: Missing correspondence information for 'cluster_markers'\n")
  if (is.null(cluster_markers$weight)) 
  {
    cat("Step1.2  Calculate the entropy-based weight.............\n")
    if (is.null(cluster_markers$avg_logFC))
      cluster_markers <- cluster_markers %>% rename(avg_logFC = avg_log2FC)
    cluster_markers %>% dplyr::count(cluster)
    cluster_markers <-
      CalculateEntropy(
        sc_ref,
        cluster_markers,
        assay = assay,
        slot = marker.slot,
        unit = unit
      )
    # saveRDS(cluster_markers, paste(output_path, "/cluster_markers_ent.rds", sep = ""))
    if(cos.filter == T)
    {
      # Calculate the cosine score as the weight
      cat("Step1.3  Calculate the cosine-based weight..............\n")
      COSG_markers <- cosg(
        sc_ref,
        groups='all',
        assay=assay,
        slot=marker.slot,
        mu=cos.mu,
        n_genes_user=cos.n_genes_user)
      marker_list<-reshape2::melt(as.matrix(COSG_markers$names))[,2:3]
      score_list<-reshape2::melt(as.matrix(COSG_markers$scores))[,2:3]
      marker_list<-cbind(marker_list,score_list)
      marker_list<-marker_list[,c(1,2,4)]
      colnames(marker_list) <- c("cluster","gene","weight")
      
      # combined with raw cluster_markers
      cluster_markers_ent_cos_filtered <- 
        cluster_markers[cluster_markers$gene %in% marker_list$gene,]
      markerlist_2 <- (marker_list %>% group_by(gene) %>% filter(weight==max(weight)))
      colnames(markerlist_2)[3] <- "cos_weight"
      colnames(cluster_markers_ent_cos_filtered)[9] <- "ent_weight"
      cluster_markers_ent_COS_score <- merge(cluster_markers_ent_cos_filtered,markerlist_2[,2:3],by="gene")
      cluster_markers_ent_COS_score <- mutate(cluster_markers_ent_COS_score,weight = ent_weight)
      # save the cluster_markers_ent_COS
      # saveRDS(cluster_markers_ent_COS_score,file = paste(output_path, "/cluster_markers_ent_COS.rds", sep = ""))
      cluster_markers <- cluster_markers_ent_COS_score
    }
  }
  cat("Step1.4  Filter the markers.............................\n")
  if(cos.filter == T)
    cluster_markers <- cluster_markers %>% dplyr::filter(ent.adj > ent.filter.threshold & cos_weight > cos.filter.threshold & weight > weight.filter.threshold)
  else
    cluster_markers <- cluster_markers %>% dplyr::filter(ent.adj > ent.filter.threshold & weight > weight.filter.threshold)
  cluster_markers$cluster <- factor(cluster_markers$cluster,levels = levels(droplevels(cluster_markers$cluster)))
  # levels(cluster_markers$cluster) <- names(table(sc_ref@meta.data[[clust_vr]]))
  cat("Save the markers........................................\n")
  saveRDS(cluster_markers,file = paste(output_path, "/cluster_markers_final.rds", sep = ""))
  sc_genes <- unique(cluster_markers$gene)
  sc_ref <- subset(x = sc_ref, features = sc_genes)
  
  #### Step2: Train the nsNMF model ####
  cat('Step2  Train the nsNMF model............................\n')
  nsnmf_mod <- train_nsnmfmod(sc_ref = sc_ref, st_vis_matr = st_vis_matr,
                              cluster_markers = cluster_markers, clust_vr = clust_vr,
                              n_top = n_top, assay = assay, slot = slot,nmf.tol=nmf.tol)
  
  #### Step3: Calculate the cluster-topic profile ####
  cat('Step3  Calculate the cluster-topic profile..............\n')
  clus_topic_profile <- topic_profile_per_cluster_nmf(h = nsnmf_mod[[1]]@h, 
                                                                 train_cell_clust = nsnmf_mod[[2]])
  
  #### Step4: Use entropy-based weighted-NNLS to implement the deconvolution ####
  cat('Step4  Weighted-NNLS to implement the deconvolution.....\n')
  decon_matr <- eWEIDE_deconvolution_nmf(nmf_mod = nsnmf_mod[[1]], cluster_markers = cluster_markers,
                                         mixture_transcriptome = st_vis_matr, normalize = normalize, 
                                         reference_profiles = clus_topic_profile, min_cont = min_cont)
  return(list(nsnmf_mod, decon_matr))
}

eWEIDE_downsample <- function(sc_ref,clust_vr,random.seed,downsample_n,n_cluster,assay,slot)
{
  set.seed(random.seed)
  # 多轮抽样并选择方差最小的子集
  Idents(sc_ref) <- sc_ref@meta.data[[clust_vr]]
  sc_barcodes <- lapply(names(table(sc_ref@meta.data[[clust_vr]])),function(clus){
    sc_tmp <- subset(x = sc_ref, cells = WhichCells(sc_ref,idents = clus))
    if(ncol(sc_tmp) <= n_cluster)
      return(colnames(sc_tmp))
    samples <- replicate(downsample_n, sample(colnames(sc_tmp),if_else(ncol(sc_tmp) < n_cluster, as.numeric(ncol(sc_tmp)),as.numeric(n_cluster))), simplify = FALSE)
    variances <- lapply(samples,function(x){
      cell_expr <- GetAssayData(object = sc_tmp,assay = assay,slot = slot)
      cell_expr <- cell_expr[,x]
      sum(apply(cell_expr, 1, var))
    })
    min_variance_index <- which.min(variances)
    selected_sample <- samples[[min_variance_index]]
    min_variance <- variances[min_variance_index]
    return(selected_sample)
  })
  
  sc_ref_down <- subset(x = sc_ref,cells = unlist(sc_barcodes))
  return(sc_ref_down)
}


meta_filter <- function(sc_ref,meta.assay,meta.resolution,meta.purity,meta.ndims,assay,clust_vr,output_path)
{
  sc_ref_meta <- sc_ref
  DefaultAssay(sc_ref_meta) <- meta.assay
  
  # Judge the required information
  if(length(sc_ref_meta@reductions) == 0)
  {
    cat('Warning: Found no matched \'reductions\' in sc_ref! \nExecute FindVariableFeatures, ScaleData and RunPCA using default parameters......\n')
    sc_ref_meta <- FindVariableFeatures(sc_ref_meta, selection.method = "vst", nfeatures = 2000)
    sc_ref_meta <- ScaleData(sc_ref_meta)
    sc_ref_meta <- RunPCA(sc_ref_meta,assay = meta.assay,features = VariableFeatures(sc_ref_meta),verbose = F)
  }
  else if(sc_ref_meta@reductions$pca@assay.used!=meta.assay)
  {
    cat('Warning: Found no matched \'reductions\' in sc_ref! \nExecute FindVariableFeatures, ScaleData and RunPCA using default parameters......\n')
    sc_ref_meta <- FindVariableFeatures(sc_ref_meta, selection.method = "vst", nfeatures = 2000)
    sc_ref_meta <- ScaleData(sc_ref_meta)
    sc_ref_meta <- RunPCA(sc_ref_meta,assay = meta.assay,features = VariableFeatures(sc_ref_meta),verbose = F)
  }
  
  if(length(sc_ref_meta@graphs) == 0)
  {
    cat('Warning: Found no \'graph\' in sc_ref! \nExecute FindNeighbors using default parameters......\n')
    sc_ref_meta <- FindNeighbors(sc_ref_meta,dims = 1:meta.ndims,reduction = 'pca')
  }
  
  # Construct meta cells by Seurat::FindClusters, using provided meta.resolution.
  if(is.null(sc_ref_meta@meta.data[[paste(meta.assay,'_snn_res.',meta.resolution,sep = '')]]))
  {
    cat('Execute FindClusters using default parameters......\n')
    sc_ref_meta <- suppressMessages(FindClusters(object = sc_ref_meta, resolution = meta.resolution,verbose = T))
    sc_ref_meta[[paste('Meta_cell_',meta.resolution,sep = '')]] <- paste('Meta_cell',t(sc_ref_meta@meta.data[[paste(meta.assay,'_snn_res.',meta.resolution,sep = '')]]),sep = '_')
    saveRDS(sc_ref_meta,file = paste(output_path,'/sc_ref_meta_before_filter.rds',sep = ''))
  }
  # sc.mat <- sc_ref_meta@assays[[assay]]@counts %>% as.matrix() %>% t() %>% data.frame
  # sc.mat$Meta_cell <- sc_ref_meta[[paste(meta.assay,'_snn_res.',meta.resolution,sep = '')]]
  
  # Filter the meta cells with purity over meta.purity.
  sc_ref_meta[[paste('Meta_cell_',meta.resolution,sep = '')]] <- paste('Meta_cell',t(sc_ref_meta@meta.data[[paste(meta.assay,'_snn_res.',meta.resolution,sep = '')]]),sep = '_')
  meta_cell.count <- table(t(sc_ref_meta[[clust_vr]]),t(sc_ref_meta[[paste('Meta_cell_',meta.resolution,sep = '')]])) %>% data.frame()
  colnames(meta_cell.count) <- c('Cell_type','Meta_cell','Freq')
  meta_cell.count$Sum <- rep(table(sc_ref_meta[[paste('Meta_cell_',meta.resolution,sep = '')]]),each = length(unique(meta_cell.count$Cell_type)))
  meta_cell.count$Proportion <- meta_cell.count$Freq/meta_cell.count$Sum
  meta_cell.count.fil <- meta_cell.count[meta_cell.count$Proportion >= meta.purity,]
  
  # Return the filtered sc_ref
  cellid_filter <- colnames(sc_ref_meta)[which(t(sc_ref_meta[[paste('Meta_cell_',meta.resolution,sep = '')]]) %in% unique(meta_cell.count.fil$Meta_cell))]
  sc_ref <- subset(sc_ref_meta,cells = cellid_filter)
  DefaultAssay(sc_ref) <- assay
  saveRDS(sc_ref,file = paste(output_path,'/sc_ref_meta_filter.rds',sep = ''))
  return(sc_ref)
}




init_nmf_matr <- function (cluster_markers, sc_ref, clust_vr,assay = "RNA", slot = "counts", n_top = NULL) 
{
  sc_nmf_ready <-
    Matrix::Matrix(t(as.matrix(Seurat::GetAssayData(sc_ref, assay = assay,slot = slot))), sparse = T)
  k <- length(unique(sc_ref@meta.data[[clust_vr]]))
  # Take then_top marker genes for each cluster into consider
  if (is.null(n_top)) 
    n_top <- max(table(cluster_markers$cluster))
  cluster_markers_cut <- cluster_markers %>% dplyr::arrange(cluster, weight) %>% dplyr::group_by(cluster) %>% 
    dplyr::top_n(n_top) %>% dplyr::ungroup() %>% data.frame()
  
  ## uniq marker genes ————————————————————
  cluster_markers_uniq <- lapply(unique(cluster_markers_cut$cluster), 
                                 function(clust) {
                                   ls1 <- cluster_markers_cut[cluster_markers_cut$cluster == 
                                                                clust, "gene"]
                                   ls2 <- cluster_markers_cut[cluster_markers_cut$cluster != 
                                                                clust, "gene"]
                                   ls1_unique <- ls1[!ls1 %in% ls2]
                                   return(cluster_markers_cut[cluster_markers_cut$cluster == 
                                                                clust & cluster_markers_cut$gene %in% ls1_unique, 
                                                              ])
                                 }) %>% bind_rows()
  seedgenes <- matrix(nrow = k, ncol = ncol(sc_nmf_ready), 
                      data = 1e-10)
  colnames(seedgenes) <- colnames(sc_nmf_ready)
  for (i in seq_len(k)) {
    clust_row <- cluster_markers_uniq$cluster == as.character(unique(sc_ref@meta.data[, 
                                                                                      clust_vr])[[i]])
    seedgenes[i, as.character(cluster_markers_uniq[clust_row, 
                                                   "gene"])] = cluster_markers_uniq[clust_row, "weight"]
  }
  W <- t(seedgenes)
  H <- matrix(data = 1e-10, nrow = k, ncol = nrow(sc_nmf_ready))
  for (i in seq_len(nrow(sc_nmf_ready))) {
    h_row <- which(unique(sc_ref@meta.data[, clust_vr]) == 
                     sc_ref@meta.data[i, clust_vr])
    H[h_row, i] <- 1
  }
  rownames(W) <- rownames(sc_ref@assays$RNA@counts)
  colnames(H) <- colnames(sc_ref@assays$RNA@counts)
  return(list(W = W, H = H))
}


train_nsnmfmod <- function(sc_ref,st_vis_matr,cluster_markers,clust_vr,n_top = NULL,assay = "RNA",slot = "data",nmf.tol=1e-04) 
{
  library(RcppML)
  print("Preparing Gene set")
  sc_ref_matr <- as.matrix(Seurat::GetAssayData(sc_ref, assay = assay, 
                                            slot = slot))
  # fliter the genes with no expression in sc_ref and st_vis data
  sc_noexp_gene <- which(!rowSums(sc_ref_matr == 0) == ncol(sc_ref_matr))
  sc_ref <- sc_ref[sc_noexp_gene, ]
  
  st_noexp_gene <- which(!rowSums(as.matrix(st_vis_matr) == 
                                 0) == ncol(st_vis_matr))
  st_vis_matr <- st_vis_matr[st_noexp_gene, ]
  
  st_genes <- rownames(st_vis_matr)
  sc_genes <- rownames(Seurat::GetAssayData(sc_ref, assay = assay, 
                                            slot = slot))
  
  # select the overlaped genes between sc_ref and st_vis
  if (length(intersect(sc_genes, st_genes)) < 10) 
    cat("Warring: The overlaped genes between sc_ref and st_vis is less than 10.\n")
  sc_ref <- sc_ref[intersect(sc_genes, st_genes), ]
  sc_ref_matr <- as.matrix(Seurat::GetAssayData(sc_ref, assay = assay, 
                                            slot = slot))
  cluster_markers <- cluster_markers[cluster_markers$gene %in% 
                                       rownames(sc_ref), ]
  
  # Normalize the sc_ref_matr
  cat("Normalize the sc_ref matrix...\n")
  count_matr_t <- scale(t(sc_ref_matr), center = FALSE, scale = apply(sc_ref_matr,1, sd, na.rm = TRUE))
  count_matr <- t(count_matr_t)
  
  # initialize and train nsNMF model
  k <- length(unique(sc_ref@meta.data[, clust_vr]))
  start_t <- Sys.time()
  cat("Initialize the NMF matrices...\n")
  init_matr <-
    seed_init_mtrx_nmf(
      cluster_markers = cluster_markers,
      sc_ref = sc_ref,
      ntop = n_top,
      clust_vr = clust_vr
    )
  print("NMF Training...")
  nsnmf_mod <- RcppML::nmf(data = as.matrix(count_matr),k = k,seed = as.matrix(init_matr),tol = nmf.tol)
  total_t <- round(difftime(Sys.time(), start_t, units = "mins"),2)
  print(sprintf("Time to initialize and train NMF model was %smins", total_t))
  return(list(nsnmf_mod, as.vector(sc_ref@meta.data[[clust_vr]])))
}






predict_spatial_mixtures_nmf_weighted <- function (nmf_mod, cluster_markers, mixture_transcriptome, normalize) 
{
  keep_genes <- rownames(nmf_mod@w)[rownames(nmf_mod@w) %in% 
                                           rownames(mixture_transcriptome)]
  mixture_transcriptome_subs <- as.matrix(mixture_transcriptome[keep_genes,])
  {
    if (normalize == "cpm") {
      count_matr <- edgeR::cpm(mixture_transcriptome_subs, 
                               normalized.lib.sizes = FALSE)
    }
    else if (normalize == "uv") {
      count_matr <- scale(t(mixture_transcriptome_subs), center = FALSE, 
                          scale = apply(mixture_transcriptome_subs, 1, sd, 
                                        na.rm = TRUE))
      count_matr <- t(count_matr)
      pos_0 <- which(rowSums(is.na(count_matr)) == ncol(count_matr))
      count_matr[pos_0, ] <- 0
    }
    else if (normalize == "raw") {
      count_matr <- mixture_transcriptome_subs
    }
    else stop("Error non specified parameter passed for normalize!")
  }
  
  
  # 20211210 refilter
  cluster_markers2 <- cluster_markers %>% 
    group_by(gene) %>% 
    dplyr::top_n(n=1,wt=weight) %>% 
    dplyr::select(c(gene,weight)) %>% 
    unique()
  
  weight_matr <- subset(cluster_markers2,cluster_markers2$gene %in% keep_genes)
  weight_matr[which(weight_matr[,2] == Inf),2] = max(weight_matr[which(weight_matr[,2] != Inf),2])
  weight_matr <- diag(sqrt(weight_matr$weight))
  
  W <- weight_matr %*% nmf_mod@w
  count_matr <- weight_matr %*% count_matr
  coef_pred <- matrix(data = NA, nrow = ncol(W), ncol = ncol(count_matr))
  colnames(coef_pred) <- colnames(count_matr)
  for (i in seq_len(ncol(count_matr))) {
    nnls_pred <- nnls::nnls(A = W, b = count_matr[, i])
    coef_pred[, i] <- nnls_pred$x
  }
  return(coef_pred)
}


#' Run mixtures through the NMF model to get the cell type composition.
eWEIDE_deconvolution_nmf <- function(nmf_mod,
                                     cluster_markers,
                                     mixture_transcriptome,
                                     normalize,
                                     reference_profiles,
                                     min_cont = 0.01)
{
  profile_matr <- predict_spatial_mixtures_nmf_weighted(nmf_mod = nmf_mod,
                                                        cluster_markers = cluster_markers,
                                                        mixture_transcriptome = mixture_transcriptome,
                                                        normalize = normalize)
  
  # We add 1 extra column to add the residual error
  decon_matr <- matrix(data = NA,
                       nrow = ncol(profile_matr),
                       ncol = ncol(reference_profiles) + 1)
  colnames(decon_matr) <- c(colnames(reference_profiles), "res_ss")
  
  # create progress bar
  print("Deconvoluting spots")
  total <- ncol(profile_matr)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in seq_len(ncol(profile_matr))) {
    ## NNLS to get cell type composition
    nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_matr[, i])
    weights <- nnls_pred$x
    
    ## get proportions of each cell type
    comp <- weights / sum(weights)
    
    ## Remove cell types not contributing the minimum
    comp[comp < min_cont] <- 0
    weights[comp < min_cont] <- 0
    
    ### Updated proportions after filtering out minimum contributions
    comp_prop <- comp / sum(comp)
    comp_prop[is.na(comp_prop)] <- 0
    
    ## Get Total sum of squares
    fit_null <- 0
    tot_ss <- sum((profile_matr[, i] - fit_null) ^ 2)
    
    ## Get % of unexplained residuals
    unexpl_ss <- nnls_pred$deviance / tot_ss
    
    decon_matr[i, 1:(ncol(decon_matr) - 1)] <- comp_prop
    decon_matr[i, ncol(decon_matr)] <- unexpl_ss
    
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  # Close progress bar
  close(pb)
  
  return(decon_matr)
}





spatial_scatterpie <- function (st_vis, cell_types_all, img_path,cols = NULL, cell_types_interest = NULL, 
                                slice = NULL, scatterpie_alpha = 1, pie_scale = 1) 
{
  if (!is(st_vis, "Seurat")) 
    stop("ERROR: st_vis must be a Seurat object!")
  if (!is(cell_types_all, "vector")) 
    stop("ERROR: cell_types_all must be a vector/list object!")
  if (!is.character(img_path)) 
    stop("ERROR: must be a character string!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) 
    stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is.numeric(scatterpie_alpha)) 
    stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale)) 
    stop("ERROR: pie_scale must be numeric between 0 and 1!")
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))
  suppressMessages(require(png))
  suppressMessages(require(jpeg))
  suppressMessages(require(grid))
  metadata_ds <- data.frame(st_vis@meta.data)
  colnames(metadata_ds) <- colnames(st_vis@meta.data)
  if (is.null(cell_types_interest)) {
    cell_types_interest <- cell_types_all
  }
  if (!all(cell_types_all %in% cell_types_interest)) {
    metadata_ds <- metadata_ds %>% tibble::rownames_to_column("barcodeID") %>% 
      dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest, 
                                           drop = FALSE])) %>% dplyr::filter(rsum != 0) %>% 
      dplyr::select("barcodeID") %>% dplyr::left_join(metadata_ds %>% 
                                                        tibble::rownames_to_column("barcodeID"), by = "barcodeID") %>% 
      tibble::column_to_rownames("barcodeID")
  }
  if (is.null(slice) | (!is.null(slice) && !slice %in% names(st_vis@images))) {
    slice <- names(st_vis@images)[1]
    warning(sprintf("Using slice %s", slice))
  }
  spatial_coord <- data.frame(st_vis@images[[slice]]@coordinates) %>% 
    tibble::rownames_to_column("barcodeID") %>% 
    dplyr::mutate(imagerow_scaled = imagerow*st_vis@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol*st_vis@images[[slice]]@scale.factors$lowres) %>% 
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"), by = "barcodeID")
  
  img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
  if (img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  }
  else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }
  img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  if (is.null(cols)) 
  {
    scatterpie_plt <-
      suppressMessages(
        ggplot2::ggplot() + ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)
        ) +
          scatterpie::geom_scatterpie(
            data = spatial_coord,
            ggplot2::aes(x = imagecol_scaled,
                         y = imagerow_scaled),
            cols = cell_types_all,
            color = NA,
            alpha = scatterpie_alpha,
            pie_scale = pie_scale
          ) +
          ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img),
                                                     0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11,
                                                                                                                 rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(
                                                                                                                   ratio = 1,
                                                                                                                   xlim = NULL,
                                                                                                                   ylim = NULL,
                                                                                                                   expand = TRUE,
                                                                                                                   clip = "on"
                                                                                                                 )
      )
  }
  else
  {
    scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                         ggplot2::annotation_custom(grob = img_grob,xmin = 0,
                                                                    xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                                         scatterpie::geom_scatterpie(data = spatial_coord,
                                                                     mapping = ggplot2::aes(x = imagecol_scaled,y = imagerow_scaled),
                                                                     cols = cell_types_all, color = NA, 
                                                                     alpha = scatterpie_alpha, pie_scale = pie_scale) + 
                                         scale_fill_manual(values=cols)+
                                         ggplot2::scale_y_reverse() + 
                                         ggplot2::ylim(nrow(img), 0) + 
                                         ggplot2::xlim(0, ncol(img)) + 
                                         cowplot::theme_half_open(11, rel_small = 1) + 
                                         ggplot2::theme_void() + 
                                         ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))
    
  }
  return(scatterpie_plt)
}


### calculate the entropy for the markers _ Version 9 ###
## Input: 
# object:seurat object
# unit = "log2", method for calculate entropy
# exp = NULL / 2 / 'e', whether to exp the entropy
# unify = T, whether to replace the duplicant genes' weight with the max weight

# entropy : calculate by log2
# pct : max pct of gene in cluster
# entropy.adjust : entropy.max - entropy
# weight = entropy.adjust * pct * avg_logFC

CalculateEntropy <- function(object,features,assay = "RNA",slot = "data",unit = "log2")
{
  library(entropy)
  # 01 Get the candidate features
  features_raw <- features
  # library(SeuratObject)
  colnames(features_raw) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")
  
  data.use <- Seurat::GetAssayData(object = object ,  slot = slot,assay = assay)
  
  features <- features$gene
  # features <- features %||% rownames(x = data.use)
  thresh.min <- 0
  # library(pacman)
  # p_unload(Seurat)
  # # p_unload(SeuratObject)
  # p_load(Seurat) 
  
  # 02 Construct the features expression percentage of each cell type
  pct <- as.data.frame(features)
  tot = 1
  for(cells.t in levels(Idents(object)))
  {
    cells.t.name <- WhichCells(object,idents = cells.t)
    pct.1 <- round(x = rowSums(x = as.matrix(data.use[features, cells.t.name, drop = FALSE] > thresh.min)) / length(x = cells.t.name),
                   digits = 3)
    pct<- cbind(pct,pct.1)
    tot = tot + 1 
  }
  all(features %in% rownames(data.use))
  
  # 03 Calculate the entropy of each feature
  pct.uq <- unique(pct)
  rownames(pct.uq) <- pct.uq$features
  pct.data <- pct.uq[, 1:length(levels(Idents(object))) + 1]
  
  colnames(pct.data) <- levels(Idents(object))
  ent.data <- pct.data/rowSums(pct.data)
  colnames(ent.data) <- levels(Idents(object))
  
  # 03.1 filter the cluster corresponding to gene which has the max pct*log2(avg_logFC+1)
  features_max_avgLogFC <- features_raw %>%
    group_by(gene) %>%
    mutate(pre_weight = log2(pct.1+1) * log2(avg_logFC + 1)) %>%
    dplyr::arrange(-pre_weight, .by_group = T) %>%
    top_n(n = 1, wt = pre_weight) %>% 
    top_n(n = 1, wt = -p_val)
  if(sum(duplicated(features_max_avgLogFC$gene))!=0)
    features_max_avgLogFC <- features_max_avgLogFC[-which(duplicated(features_max_avgLogFC$gene)),]
  
  ent.result <- data.frame(gene = character(),entropy= numeric(),entropy.adjust= numeric(),pct= numeric(),cluster= character(), stringsAsFactors=F)
  entropy.max = entropy::entropy(c(rep(1,length(levels(Idents(object))))),unit = unit)
  pb <- txtProgressBar(style=3)
  for(i in 1:nrow(features_max_avgLogFC)){
    ent.result[i,1] <- rownames(ent.data)[i]
    ent.result[i,2] <- entropy::entropy(as.numeric(ent.data[i,]),unit = unit)
    ent.result[i,3] <- entropy.max - ent.result[i,2]
    ent.result[i,4] <- max(pct.data[i,])
    ent.result[i,5] <- as.character(features_max_avgLogFC$cluster[which(features_max_avgLogFC$gene == rownames(ent.data)[i])[1]])
    setTxtProgressBar(pb, i/length(rownames(pct.uq)))
  }
  close(pb)
  
  cat('\n')
  # 04 Calculate the weight of each feature based on the entropy, pct and avg_logFC
  ent.result.rank <- ent.result %>% 
    group_by(gene) %>% 
    dplyr::arrange(gene,.by_group = T)
  
  ent.result.rank$avg_logFC <- features_max_avgLogFC$avg_logFC
  
  # 20220623: Remove Inf in avg_logFC, replace with max avg_logFC
  max_avg_logFC <- max(ent.result.rank$avg_logFC[which(ent.result.rank$avg_logFC != Inf)])
  ent.result.rank$avg_logFC[which(ent.result.rank$avg_logFC == Inf)] <- max_avg_logFC
  
  ent.result.rank <- ent.result.rank %>% 
    mutate(weight = entropy.adjust * pct * avg_logFC)%>% 
    group_by(cluster) %>% 
    dplyr::arrange(-weight,.by_group = T)
  rownames(ent.result.rank) <- ent.result.rank$gene
  
  # 05 Rescue the duplicate features among all cell types
  features_new <- features_raw
  features_new$ent.adj <- 0
  features_new$weight <- 0
  features_new$avg_logFC[which(features_new$avg_logFC == Inf)] <- max_avg_logFC
  cat('\n')
  for(i in 1:nrow(features_new))
  {
    features_new$ent.adj[i] <- ent.result.rank$entropy.adjust[which(ent.result.rank$gene == features_new$gene[i])[1]]
    features_new$weight[i] <- 4* features_new$ent.adj[i]/entropy.max * features_new$pct.1[i] * features_new$avg_logFC[i]
    if(features_new$weight[i] > 1 )
      features_new$weight[i] <- 1
  }
  return(features_new)
  
}


# 20220102: This function can simulate a Spatial data set with (n = 1000) spots using single cell data. 
# (Default) First randomly select 1~(max.celltype = 4) celltype, then randomly select 1~(max.cell = 5) cells for each selected celltype.
new_simulate <- function(sc_ref,clust_vr,n = 1000,verbose = TRUE,max.celltype = 4,max.cell = 5, min.celltype = 1,min.cell = 1)
{
  # Check variables
  if (is(sc_ref) != "Seurat") stop("ERROR: sc_ref must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
  if (!is.numeric(n.celltype)) stop("ERROR: n.celltype must be an integer!")
  if (!is.numeric(n.cell)) stop("ERROR: n.cell must be an integer!")
  if (!is.logical(verbose)) stop("ERROR: verbose must be a logical object!")
  
  suppressMessages(require(DropletUtils)) # For the downsampling
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  
  sc_ref@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                       x = sc_ref@meta.data[, clust_vr],
                                       perl = TRUE)
  print("Generating synthetic test spots...")
  start_gen <- Sys.time()
  # create progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Save count matrix
  count_matr <- as.matrix(sc_ref@assays$RNA@counts)
  
  # Save celltype names
  celltype <- names(table(sc_ref[[clust_vr]]))
  
  ds_spots <- lapply(seq_len(n), function(i) {
    
    # Select between 1 and n.celltype randomly from the celltype
    selected_celltype <- sample(celltype, sample(x = min.celltype:max.celltype, size = 1))
    
    cell_pool <- NULL
    for (k in selected_celltype) {
      # construct selected count matrix
      selected_cells <- which(sc_ref[[clust_vr]]==k)
      
      # sc_ref[[clust_vr]][selected_cells,1]
      selected_count_matr <- count_matr[,selected_cells]
      
      # Select between 1 and n.cell randomly from the selected_count_matr
      cell_pool <- c(cell_pool,sample(colnames(selected_count_matr), sample(x = min.cell:max.cell, size = 1)))
    }
    
    
    
    # Determine the weight each cell will have on the synthetic spot
    # weigh <-runif(length(cell_pool))
    # weigh <- weigh/sum(weigh)
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_matr) %in% cell_pool)
    tmp_ds <- sc_ref@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- tmp_ds %>%
      dplyr::select(all_of(clust_vr), weight) %>%
      # dplyr::mutate(clust_vr = paste("clust_",
      #                                tmp_ds[, clust_vr], sep = "")) %>%
      dplyr::group_by(!! sym(clust_vr)) %>%
      dplyr::summarise(sum_weights = sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = all_of(clust_vr),
                         values_from = sum_weights) %>%
      dplyr::mutate(name = name_simp)
    
    # Generate synthetic spot
    
    ## Here we multiply each vector by its weight
    # weighted <- lapply(1:length(cell_pool), function(ii) expr <- as.integer(round(weigh[[ii]]*count_matr[hvg,cell_pool[[ii]]],0)))
    
    ## Next step we add all the vectors by position
    # syn_spot <- Reduce(`+`,weighted)
    # ret_ds <- data.frame(gene=hvg, tmp=syn_spot)
    
    ## Here we add up the counts of each cell
    syn_spot <- rowSums(as.matrix(count_matr[, cell_pool])); sum(syn_spot)
    names_genes <- names(syn_spot)
    
    ## Downsample
    ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k
    if (sum(syn_spot) > 25000) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, sparse = T),
                                                        prop = 20000 / sum(syn_spot))
    } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    return(list(syn_spot_sparse, spot_ds))
  })
  
  ds_syn_spots <- purrr::map(ds_spots, 1) %>%
    base::Reduce(function(m1, m2) cbind(unlist(m1), unlist(m2)), .)
  
  # Generate dataframe of spot characteristic
  ds_spots_metadata <- purrr::map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    data.frame()
  
  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
  
  # change column order so that its progressive
  lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(sc_ref@meta.data[, clust_vr]))
  colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
  # all_cn <- c(paste("clust_", lev_mod, sep = ""), "name") # This was set to deal when cluster names are numeric
  
  # Check if there are missing columns (Cell types not selected) and add them as all 0s
  if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(sc_ref@meta.data[, clust_vr])) + 1)) {
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  } else {
    
    missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}



# 20230222: Generate in silico data sets: 
# 1. cell type complexity
# 2. low percentage cluster
# 3. cell type sensitivity
Simu_cell_type_complexity <- function(sc_ref,clust_vr,n = 1000,verbose = TRUE,
                                      cell_type_complexity = 4,max.cell = 5, 
                                      min.celltype = 1,min.cell = 1)
{
  # Check variables
  if (is(sc_ref) != "Seurat") stop("ERROR: sc_ref must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
  if (!is.numeric(n.celltype)) stop("ERROR: n.celltype must be an integer!")
  if (!is.numeric(n.cell)) stop("ERROR: n.cell must be an integer!")
  if (!is.logical(verbose)) stop("ERROR: verbose must be a logical object!")
  
  suppressMessages(require(DropletUtils)) # For the downsampling
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tidyr))
  
  sc_ref@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                       x = sc_ref@meta.data[, clust_vr],
                                       perl = TRUE)
  print("Generating synthetic test spots...")
  start_gen <- Sys.time()
  # create progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Save count matrix
  count_matr <- as.matrix(sc_ref@assays$RNA@counts)
  
  # Save celltype names
  celltype <- names(table(sc_ref[[clust_vr]]))
  
  ds_spots <- lapply(seq_len(n), function(i) {
    
    # Select between 1 and n.celltype randomly from the celltype
    selected_celltype <- sample(celltype, sample(x = min.celltype:max.celltype, size = 1))
    
    cell_pool <- NULL
    for (k in selected_celltype) {
      # construct selected count matrix
      selected_cells <- which(sc_ref[[clust_vr]]==k)
      
      # sc_ref[[clust_vr]][selected_cells,1]
      selected_count_matr <- count_matr[,selected_cells]
      
      # Select between 1 and n.cell randomly from the selected_count_matr
      cell_pool <- c(cell_pool,sample(colnames(selected_count_matr), sample(x = min.cell:max.cell, size = 1)))
    }
    
    
    
    # Determine the weight each cell will have on the synthetic spot
    # weigh <-runif(length(cell_pool))
    # weigh <- weigh/sum(weigh)
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_matr) %in% cell_pool)
    tmp_ds <- sc_ref@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- tmp_ds %>%
      dplyr::select(all_of(clust_vr), weight) %>%
      # dplyr::mutate(clust_vr = paste("clust_",
      #                                tmp_ds[, clust_vr], sep = "")) %>%
      dplyr::group_by(!! sym(clust_vr)) %>%
      dplyr::summarise(sum_weights = sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = all_of(clust_vr),
                         values_from = sum_weights) %>%
      dplyr::mutate(name = name_simp)
    
    # Generate synthetic spot
    
    ## Here we multiply each vector by its weight
    # weighted <- lapply(1:length(cell_pool), function(ii) expr <- as.integer(round(weigh[[ii]]*count_matr[hvg,cell_pool[[ii]]],0)))
    
    ## Next step we add all the vectors by position
    # syn_spot <- Reduce(`+`,weighted)
    # ret_ds <- data.frame(gene=hvg, tmp=syn_spot)
    
    ## Here we add up the counts of each cell
    syn_spot <- rowSums(as.matrix(count_matr[, cell_pool])); sum(syn_spot)
    names_genes <- names(syn_spot)
    
    ## Downsample
    ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k
    if (sum(syn_spot) > 25000) {
      syn_spot_sparse <- DropletUtils::downsampleMatrix(Matrix::Matrix(syn_spot, sparse = T),
                                                        prop = 20000 / sum(syn_spot))
    } else {
      syn_spot_sparse <- Matrix::Matrix(syn_spot, sparse = T)
    }
    
    rownames(syn_spot_sparse) <- names_genes
    colnames(syn_spot_sparse) <- name_simp
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    return(list(syn_spot_sparse, spot_ds))
  })
  
  ds_syn_spots <- purrr::map(ds_spots, 1) %>%
    base::Reduce(function(m1, m2) cbind(unlist(m1), unlist(m2)), .)
  
  # Generate dataframe of spot characteristic
  ds_spots_metadata <- purrr::map(ds_spots, 2) %>%
    dplyr::bind_rows() %>%
    data.frame()
  
  ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
  
  # change column order so that its progressive
  lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(sc_ref@meta.data[, clust_vr]))
  colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
  # all_cn <- c(paste("clust_", lev_mod, sep = ""), "name") # This was set to deal when cluster names are numeric
  
  # Check if there are missing columns (Cell types not selected) and add them as all 0s
  if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(sc_ref@meta.data[, clust_vr])) + 1)) {
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  } else {
    
    missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
    ds_spots_metadata[missing_cols] <- 0
    ds_spots_metadata <- ds_spots_metadata[, lev_mod]
  }
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}




# 20220110: visualization of simulate method
# see it on the KRAS_Pt8_benchmark.Rmd


# SPOTlight functions
topic_profile_per_cluster_nmf <- function (h, train_cell_clust) 
{
  if (!is(h, "matrix")) 
    stop("ERROR: h must be a matrix object!")
  if (!is(train_cell_clust, "vector")) 
    stop("ERROR: train_cell_clust must be a vector/list object!")
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))
  h_ds <- data.frame(t(h))
  h_ds[, "clust_vr"] <- train_cell_clust
  ct_topic_profiles <- h_ds %>% dplyr::group_by(clust_vr) %>% 
    dplyr::summarise_all(list(median)) %>% tibble::remove_rownames() %>% 
    tibble::column_to_rownames(var = "clust_vr") %>% as.matrix()
  ct_topic_profiles_t <- t(ct_topic_profiles)
  colnames(ct_topic_profiles_t) <- gsub("[[:punct:]]|[[:blank:]]", 
                                        ".", colnames(ct_topic_profiles_t))
  return(ct_topic_profiles_t)
}


seed_init_mtrx_nmf <- function (cluster_markers, sc_ref, clust_vr, ntop = NULL) 
{
  sc_nmf_ready <- prep_seobj_topic_fun(sc_ref = sc_ref)
  k <- length(unique(sc_ref@meta.data[, clust_vr]))
  if (is.null(ntop)) 
    ntop <- max(table(cluster_markers$cluster))
  cluster_markers_cut <- suppressMessages(cut_markers2(markers = cluster_markers, 
                                                       ntop = ntop))
  cluster_markers_uniq <- lapply(unique(cluster_markers_cut$cluster), 
                                 function(clust) {
                                   ls1 <- cluster_markers_cut[cluster_markers_cut$cluster == 
                                                                clust, "gene"]
                                   ls2 <- cluster_markers_cut[cluster_markers_cut$cluster != 
                                                                clust, "gene"]
                                   ls1_unique <- ls1[!ls1 %in% ls2]
                                   return(cluster_markers_cut[cluster_markers_cut$cluster == 
                                                                clust & cluster_markers_cut$gene %in% ls1_unique, 
                                                              ])
                                 }) %>% bind_rows()
  seedgenes <- matrix(nrow = k, ncol = ncol(sc_nmf_ready), 
                      data = 1e-10)
  colnames(seedgenes) <- colnames(sc_nmf_ready)
  for (i in seq_len(k)) {
    clust_row <- cluster_markers_uniq$cluster == as.character(unique(sc_ref@meta.data[, 
                                                                                      clust_vr])[[i]])
    seedgenes[i, as.character(cluster_markers_uniq[clust_row, 
                                                   "gene"])] = cluster_markers_uniq[clust_row, "weight"]
  }
  W <- t(seedgenes)
  rownames(W) <- rownames(sc_ref@assays$RNA@counts)
  return(W)
}

prep_seobj_topic_fun <- function (sc_ref) 
{
  if (is(sc_ref) != "Seurat") 
    stop("ERROR: sc_ref must be a Seurat object!")
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))
  count_mtrx <- t(as.matrix(sc_ref@assays$RNA@counts))
  count_mtrx <- Matrix::Matrix(count_mtrx, sparse = T)
  return(count_mtrx)
}

cut_markers2 <- function (markers, ntop) 
{
  if (!is.data.frame(markers)) 
    stop("ERROR: markers must be a data.frame object!")
  if (!is.numeric(ntop)& !is.null(ntop)) 
    stop("ERROR: ntop must be a an integer!")
  if (!"cluster" %in% colnames(markers)) 
    stop("ERROR: cluster information needs to be provided in markers.")
  if (!"gene" %in% colnames(markers)) 
    stop("ERROR: gene information needs to be provided in markers.")
  suppressMessages(require(dplyr))
  suppressMessages(require(entropy))
  tmp_markers <- markers %>% dplyr::arrange(cluster, weight) %>% 
    dplyr::group_by(cluster) %>% dplyr::top_n(ntop) %>% 
    dplyr::ungroup() %>% dplyr::select(gene,weight,  cluster) %>% data.frame()
  return(tmp_markers)
}


# 20230320: Test accuracy, recall, F1 score and so on.
test_MultiIndex <- function (deconv_result, synthetic_comp) 
{
  if (!is.matrix(deconv_result)) 
    stop("ERROR: deconv_result must be a matrix object!")
  if (!is.matrix(synthetic_comp)) 
    stop("ERROR: synthetic_comp must be a matrix object!")
  colnames(synthetic_comp) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
                                          ".", x = colnames(synthetic_comp), perl = TRUE)
  colnames(deconv_result) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", 
                                             ".", x = colnames(deconv_result), perl = TRUE)
  suppressMessages(require(philentropy))
  true_jsd_mtrx <- matrix(nrow = nrow(deconv_result), 
                          ncol = 1)
  tp <- 0
  tn <- 0
  fp <- 0
  fn <- 0
  for (i in seq_len(nrow(synthetic_comp))) {
    x <- rbind(synthetic_comp[i, ], deconv_result[i,])
    if (sum(synthetic_comp[i, ]) > 0) {
      true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, unit = "log2", est.prob = "empirical"))
    }
    else {
      true_jsd_mtrx[i, 1] <- 1
    }
    for (index in colnames(synthetic_comp)) {
      if (x[1, index] > 0 & x[2, index] > 0) {
        tp <- tp + 1
      }
      else if (x[1, index] == 0 & x[2, index] == 0) {
        tn <- tn + 1
      }
      else if (x[1, index] > 0 & x[2, index] == 0) {
        fn <- fn + 1
      }
      else if (x[1, index] == 0 & x[2, index] > 0) {
        fp <- fp + 1
      }
    }
    rm(index)
  }
  rm(i)
  accuracy <- round((tp + tn)/(tp + tn + fp + fn), 2)
  specificity <- round(tn/(tn + fp), 2)
  precision <- round(tp/(tp + fp), 2)
  recall <- round(tp/(tp + fn), 2)
  F1 <- round(2 * ((precision * recall)/(precision + recall)), 2)
  quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx, na.rm = TRUE), c(0.25, 0.5, 0.75)), 5)
  # cat(sprintf("The following summary statistics are obtained:
  #             Accuracy: %s,
  #             Specificity: %s,
  #             precision: %s,
  #             recall: %s,
  #             F1 score: %s,
  #             JSD quantiles: %s[%s-%s]", 
  #             accuracy, specificity, precision, recall, F1,
  #             quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), 
  #     sep = "\n")
  result <- data.frame(matrix(nrow = 1,ncol = 8))
  result <- cbind(accuracy, specificity, precision, recall, F1,quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]])
  colnames(result) <- c('Accuracy','Specificity','precision','recall','F1_score','JSD_quantiles','JSD_quantiles_1','JSD_quantiles_2')
  return(result)
}


# 20220609: Test Pearson correlation for each cell type
test_Pearson <- function(deconv_result,synthetic_comp){
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = ncol(deconv_result2))
  for(i in 1:ncol(deconv_result2))
  {
    t <- cbind(deconv_result2[,i],synthetic_comp[,i])
    if(sum(synthetic_comp[,i])==0)
      result_list[i] <- NA
    else
    {
      dataPearson <- round(cor(t,method = c("pearson")) , 2)
      if(is.na(dataPearson[1,2]))
        result_list[i] <- 0
      else
        result_list[i] <- dataPearson[1,2]
    }
    
  }
  names(result_list) <- colnames(synthetic_comp)
  return(result_list)
}

test_Pearson_spotlevel <- function(deconv_result,synthetic_comp){
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = nrow(deconv_result2))
  for(i in 1:nrow(deconv_result2))
  {
    t <- cbind(deconv_result2[i,],synthetic_comp[i,])
    dataPearson <- round(cor(t,method = c("pearson")) , 2)
    if(is.na(dataPearson[1,2]))
      result_list[i] <- 0
    else
      result_list[i] <- dataPearson[1,2]
  }
  names(result_list) <- rownames(synthetic_comp)
  return(result_list)
}

test_RMSE_spotlevel <- function(deconv_result,synthetic_comp){
  library(Metrics)
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = nrow(deconv_result2))
  for(i in 1:nrow(deconv_result2))
  {
    result_list[i] <- rmse(deconv_result2[i,],synthetic_comp[i,])
  }
  names(result_list) <- rownames(synthetic_comp)
  return(result_list)
}

test_RMSE <- function(deconv_result,synthetic_comp){
  library(Metrics)
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = ncol(deconv_result2))
  for(i in 1:ncol(deconv_result2))
  {
    result_list[i] <- rmse(deconv_result2[,i],synthetic_comp[,i])
  }
  names(result_list) <- rownames(synthetic_comp)
  return(result_list)
}



dot_plot_profiles_fun <- function (h, train_cell_clust)
{
  suppressMessages(require(stringr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))
  suppressMessages(require(tidyr))
  suppressMessages(require(ggplot2))
  h_df <- data.frame(t(h))
  colnames(h_df) <- gsub(".", " ", colnames(h_df), fixed = TRUE)
  h_ds <- round(h_df / rowSums(h_df), 4)
  h_ds[, "clust_vr"] <- train_cell_clust
  train_cells_plt <- h_ds %>% tibble::rowid_to_column("id") %>%
    tidyr::pivot_longer(
      cols = -c(clust_vr, id),
      names_to = "topics",
      values_to = "weights"
    ) %>% dplyr::mutate(weights_txt = dplyr::if_else(weights >0.1, round(weights, 2), 0)) %>%
    ggplot2::ggplot(aes(x = id,y = topics)) + 
    ggplot2::geom_point(aes(size = weights,colour = weights)) + 
    ggplot2::facet_wrap(clust_vr ~., scales = "free") + 
    ggplot2::scale_color_continuous(low = "grey",high = "#59b371") + 
    ggplot2::theme_classic() + 
    ggplot2::labs(title = "NMF: Topic proportion within cell types") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 20),
      axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5),
      axis.text = ggplot2::element_text(size = 15)) +
    ggplot2::scale_size(range = c(0, 5)) + ggplot2::guides(
      colour = ggplot2::guide_legend("Proportion"),
      size = ggplot2::guide_legend("Proportion")
    )
  ct_topic_profiles <- h_ds %>% dplyr::group_by(clust_vr) %>%
    dplyr::summarise_all(list(median)) %>% tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("clust_vr")
  ct_topic_profiles <- ct_topic_profiles / rowSums(ct_topic_profiles)
  ct_topic_profiles[is.na(ct_topic_profiles)] <- 0
  
  topic_order <- apply(ct_topic_profiles, 1, function(x) order(x, decreasing = TRUE)[1])
  cell_type_order <- apply(ct_topic_profiles, 2, function(x) order(x, decreasing = TRUE)[1])
  ct_topic_profiles <- ct_topic_profiles[,topic_order]

  cell_type_pf <- round(ct_topic_profiles, 2) %>% tibble::rownames_to_column("Cell type") %>%
    tidyr::pivot_longer(cols = -`Cell type`, names_to = "Topics") %>%
    dplyr::mutate(value_txt = dplyr::if_else(value > 0.1,round(value, 2), 0),
                  Topics = factor(x = Topics,levels = colnames(ct_topic_profiles)),
                  `Cell type` = factor(x = `Cell type`,levels = rownames(ct_topic_profiles))) 
  
  
  cell_type_plt <- cell_type_pf %>% 
    ggplot2::ggplot(ggplot2::aes(x = `Cell type`, y = Topics)) + 
    ggplot2::geom_point(ggplot2::aes(size = value,colour = value)) + 
    ggplot2::scale_color_continuous(low = "grey",high = "#59b371") + 
    ggplot2::theme_classic() + ggplot2::labs(title = "NMF: Topic profiles by cell type") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5,size = 20),
      axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5),
      axis.text = ggplot2::element_text(size = 15)
    ) +
    ggplot2::scale_size(range = c(0, 10)) + ggplot2::guides(
      colour = ggplot2::guide_legend("Proportion"),
      size = ggplot2::guide_legend("Proportion")
    )
  
  return(list(train_cells_plt, cell_type_plt))
}


