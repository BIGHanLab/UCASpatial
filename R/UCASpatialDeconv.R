#' @title Ultra-high resolution Cellular annotation of Spatial transcriptomics - UCASpatial
#'
#' @description This package could deconvolute the spatial transcriptomics data into cell-type-composition based on a reference scRNA-seq data with appropriate cell-type labels.
#'
#' @param sc_ref reference single cell RNA-seq data: Seurat object
#' @param st_vis spatial transcriptomic data: Seurat object
#' @param clust_vr colname of the meta.data in sc_ref which include the cluster information
#' @param assay optional: set the default assay for sc_ref. By default is 'RNA'.
#' @param slot optional: set the default slot for sc_ref. By default is 'counts'.
#' @param output_path optional: set the output path for the intermediate files. If is NULL, it will use current path by 'getwd()' function. By default is 'NULL'.
#' @param cluster_markers optional: give the marker genes for each cluster of sc_ref. If is NULL, it will be calculated use the default parameters. By default is 'NULL'.
#' @param min.pct optional: parameter for Seurat::FindAllMarkers
#' @param logfc.threshold optional: parameter for Seurat::FindAllMarkers
#' @param min.diff.pct optional: parameter for Seurat::FindAllMarkers
#' @param normalize optional: method for normalizing the count matrix, including uv: unit variance, raw: no transformation applied. By default is 'uv'.
#' @param downsample_n optional: set the downsample number of times. By default is '1'.
#' @param n_cluster optional: set the maximum size for each cluster after downsample. By default is '100'.
#' @param n_top optional: set the n_top marker gene used for each cluster. By default is 'NULL'.
#' @param min_cont optional: set the cutoff threshold. By default is '0.01'.
#' @param remove.RPL whether remove the RPL-related genes. By default is 'FALSE'.
#' @param remove.MT whether remove the MT-related genes. By default is 'FALSE'.
#' @param cos.filter optional: whether use the cosg methods to filter some features. By default is 'TRUE'.
#' @param cos.mu optional: parameter for cosg.
#' @param cos.n_genes_user optional: parameter for cosg.
#' @param marker.slot optional: set the slot of sc_ref for marker finding. By default is 'data'.
#' @param min_cont optional: parameter for calculating entropy. By default is '0.01'.
#' @param unit optional: parameter for calculating entropy. By default is 'log2'.
#' @param random.seed optional: set the random seed for randomly downsampling. By default is '10000'.
#' @param meta.filter optional: whether execute meta-purity filtering process. By default is 'TRUE'.
#' @param nmf.tol optional: parameter for the nmf iteration. By default is '1e-04'.
#' @param meta.assay optional: assay parameter for meta-purity filtering. By default is 'integrated'.
#' @param meta.ndims optional: ndims parameter for meta-purity filtering. By default is '30'.
#' @param meta.resolution optional: resolution parameter for meta-purity filtering. By default is '100'.
#' @param meta.purity optional: purity parameter for meta-purity filtering. By default is '0.95'.
#' @param ent.filter.threshold optional: parameter for filter low quality features based on entropy.
#' @param cos.filter.threshold optional: parameter for filter low quality features based on cosine similarity.
#' @param weight.filter.threshold optional: parameter for filter low quality features based on weight.
#'
#' @importFrom dplyr %>%
#' @return a list including NMF result and deconvolution result
#' @export


UCASpatial_deconv <- function (sc_ref, st_vis, clust_vr, assay = "RNA", slot = "data",output_path = NULL,
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
  # suppressMessages(require(SPOTlight))
  suppressMessages(require(reshape2))
  suppressMessages(require(RcppML))

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
    sc_ref_down <- UCASpatial_downsample(sc_ref = sc_ref,clust_vr = clust_vr,random.seed = random.seed,
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
  decon_matr <- UCASpatial_deconvolution_nmf(nmf_mod = nsnmf_mod[[1]], cluster_markers = cluster_markers,
                                         mixture_transcriptome = st_vis_matr, normalize = normalize,
                                         reference_profiles = clus_topic_profile, min_cont = min_cont)
  return(list(nsnmf_mod, decon_matr))
}

UCASpatial_downsample <- function(sc_ref,clust_vr,random.seed,downsample_n,n_cluster,assay,slot)
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
  suppressMessages(require(RcppML))
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


# Run mixtures through the NMF model to get the cell type composition.
UCASpatial_deconvolution_nmf <- function(nmf_mod,
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


cosg <- function (object, groups = "all", assay = "RNA", slot = "data",
          mu = 1, remove_lowly_expressed = TRUE, expressed_pct = 0.1,
          n_genes_user = 100)
{
  genexcell <- Seurat::GetAssayData(object = object[[assay]],
                                    slot = slot)
  if (length(groups) > 1) {
    object <- subset(x = object, idents = groups)
    group_info <- Seurat::Idents(object = object)
  }
  else {
    if (groups == "all") {
      group_info <- Seurat::Idents(object = object)
    }
    else {
      stop("Cannot perform marker gene identification on a single cluster. Please reset the groups variable.")
    }
  }
  groups_order = sort(unique(group_info))
  n_cluster = length(groups_order)
  if (n_cluster == 1) {
    stop("Cannot perform marker gene identification on a single cluster.")
  }
  n_cell = ncol(genexcell)
  n_gene = nrow(genexcell)
  gene_name = rownames(genexcell)
  if (n_genes_user > n_gene) {
    n_genes_user = n_gene
  }
  cluster_mat = matrix(0, nrow = n_cluster, ncol = n_cell)
  order_i = 1
  for (group_i in groups_order) {
    idx_i = group_info == group_i
    cluster_mat[order_i, idx_i] = 1
    order_i = order_i + 1
  }
  cluster_mat_sparse = as(cluster_mat, "dgCMatrix")
  cosine_sim = proxyC::simil(genexcell, cluster_mat_sparse,
                             method = "cosine", drop0 = TRUE)
  pos_nonzero = cosine_sim != 0
  pos_nonzero = which(as.matrix(pos_nonzero), arr.ind = TRUE)
  genexlambda = cosine_sim * cosine_sim
  e_power2_sum = Matrix::rowSums(genexlambda)
  if (mu == 1) {
    genexlambda[pos_nonzero] = genexlambda[pos_nonzero]/(replicate(ncol(genexlambda),
                                                                   e_power2_sum)[as.matrix(pos_nonzero)])
  }
  else {
    genexlambda[pos_nonzero] = genexlambda[pos_nonzero]/(((1 -
                                                             mu) * genexlambda[pos_nonzero] + mu * (replicate(ncol(genexlambda),
                                                                                                              e_power2_sum)[as.matrix(pos_nonzero)])))
  }
  genexlambda = genexlambda * cosine_sim
  rank_stats_names = data.frame(matrix(matrix(), n_genes_user,
                                       length(groups_order), dimnames = list(seq(1, n_genes_user),
                                                                             groups_order)), stringsAsFactors = F)
  rank_stats_scores = data.frame(matrix(matrix(), n_genes_user,
                                        length(groups_order), dimnames = list(seq(1, n_genes_user),
                                                                              groups_order)), stringsAsFactors = F)
  order_i = 1
  for (group_i in groups_order) {
    idx_i = group_info == group_i
    scores = genexlambda[, order_i]
    if (remove_lowly_expressed) {
      n_cells_expressed = tabulate(genexcell[, idx_i]@i +
                                     1)
      n_cells_i = sum(idx_i)
      scores[n_cells_expressed < n_cells_i * expressed_pct] = -1
    }
    global_indices = select_top_n(scores, n_genes_user)
    rank_stats_names[, order_i] = gene_name[global_indices]
    rank_stats_scores[, order_i] = scores[global_indices]
    order_i = order_i + 1
  }
  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order
  ranks_stats = list(names = rank_stats_names, scores = rank_stats_scores)
  return(ranks_stats)
}


select_top_n <- function (scores, n_top)
{
  d <- data.frame(x = data.table::copy(scores), indice = seq(1,
                                                             length(scores)))
  data.table::setDT(d)
  data.table::setorder(d, -x)
  n_top_indice <- d$indice[1:n_top]
  return(n_top_indice)
}





