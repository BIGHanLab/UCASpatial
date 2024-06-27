library(Seurat)
library(pacman)
p_unload(Seurat)
p_load(Seurat)
setwd('/data/xy/Spatial_transcriptome/eWEIDE/20240511_Codes_for_upload/Fig2/')

source('../source/UCASpatial_final_v1.R')
source('../source/color_ref.R')
sc.silico <- readRDS("../../Data/Fig2/sc.silico.recluster.rds")

#### generate gradient resolution of sc.silico ####
# high res (recluster)
{
  sc.silico$high_res_ident <- sc.silico$new_ident
  DimPlot(sc.silico,label = T,repel = T,cols = col.sc$colors)
  sc.silico$high_res_ident <- sc.silico$new_ident
  DimPlot(sc.silico,group.by = 'high_res_ident',label = T,repel = T,label.size = 5,cols = col.sc$Zurui_color) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom') +
    ggtitle('High resolution')
}

# low res
{
  sc.silico$low_res_ident <- sc.silico$cluster_last_version
  
  table(sc.silico$low_res_ident)
  sc.silico$low_res_ident <- as.factor(sc.silico$low_res_ident)
  levels(sc.silico$low_res_ident)
  levels(sc.silico$low_res_ident)[1:3] <- 'B_cells'
  levels(sc.silico$low_res_ident)[2:4] <- 'CD4_T'
  levels(sc.silico$low_res_ident)[3:4] <- 'CD8_T'
  levels(sc.silico$low_res_ident)[c(4,7,8,9,11)] <- 'Myeloid_cells'
  levels(sc.silico$low_res_ident)[5:6] <- 'Epithelial'
  levels(sc.silico$low_res_ident)[8:9] <- 'Stromal'
  levels(sc.silico$low_res_ident)[6] <- 'CD8_T'
  levels(sc.silico$low_res_ident)
  DimPlot(sc.silico,group.by = 'low_res_ident',label = T,repel = T,label.size = 5) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom') +
    ggtitle('Low resolution')
  # ggsave('Low_resolution.pdf',height = 6,width = 6)
}

# medium res
{
  sc.silico$med_res_ident <- sc.silico$cluster_last_version
  
  table(sc.silico$med_res_ident)
  sc.silico$med_res_ident <- as.factor(sc.silico$med_res_ident)
  levels(sc.silico$med_res_ident)
  levels(sc.silico$med_res_ident)[2:3] <- 'memory_naive_B_cells'
  levels(sc.silico$med_res_ident)[3:4] <- 'CD4_T_CCR7'
  levels(sc.silico$med_res_ident)[5:6] <- 'CD8_T'
  levels(sc.silico$med_res_ident)[c(6,9,11,13)] <- 'Myeloid_cells'
  levels(sc.silico$med_res_ident)[7:8] <- 'Epithelial'
  levels(sc.silico$med_res_ident)[11:12] <- 'Stromal'
  levels(sc.silico$med_res_ident)
  DimPlot(sc.silico,group.by = 'med_res_ident')
  DimPlot(sc.silico,group.by = 'med_res_ident',label = T,repel = T,label.size = 5) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom') +
    ggtitle('Medium resolution')
  # ggsave('Medium_resolution.pdf',height = 6,width = 6)
}

# saveRDS(sc.silico,'sc.silico.RDS')

#### generate TME ####

random.seed <- c(19980205,19990911,19960330,20000502,20230509)

## tumor core:  0.5X ~ 1.5X
# 50% - 80% epithelial
# 0% - 20% myeloid
# 0% - 10% stromal
# 0% - 20% CD4 T
# 0% - 20% CD8 T
# 0% B, 0% Plasma
TME_simu_TC <- function(sc_ref,clust_vr,clust_vr_list,n = 50,verbose = TRUE,n.cell = 10,seeds = 20230508)
{
  # Check variables
  if (is(sc_ref) != "Seurat") stop("ERROR: sc_ref must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
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
  
  set.seed(seeds)
  # Save count matrix
  count_matr <- as.matrix(sc_ref@assays$RNA@counts)
  
  # Save celltype names
  celltype <- names(table(sc_ref[[clust_vr]]))
  
  ds_spots <- lapply(seq_len(n), function(i) {
    # epithelial
    n.epi <- sample(x = as.integer(0.5*n.cell):as.integer(0.8*n.cell), size = 1)
    n.mye <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.str <- sample(x = as.integer(0 * n.cell):as.integer(0.1 * n.cell), size = 1)
    n.CD4 <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.CD8 <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    cell_pool <- NULL
    n.select.ct <- list(n.epi,n.mye,n.str,n.CD4,n.CD8)
    celltype.list <- list('Epithelial','Myeloid.cells','Stromal', 'CD4.T' ,'CD8.T')
    cell_pool <- unlist(mapply(function(x,y){
      selected_cells <- which(sc_ref[[clust_vr]]==y)
      selected_count_matr <- count_matr[,selected_cells]
      cell_pool <- sample(colnames(selected_count_matr),x)
    },n.select.ct,celltype.list))
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_matr) %in% cell_pool)
    tmp_ds <- sc_ref@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    
    spot_ds <- lapply(clust_vr_list,function(clust_vr){
      tmp_ds %>%
        dplyr::select(all_of(clust_vr), weight) %>%
        # dplyr::mutate(clust_vr = paste("clust_",
        #                                tmp_ds[, clust_vr], sep = "")) %>%
        dplyr::group_by(!! sym(clust_vr)) %>%
        dplyr::summarise(sum_weights = sum(weight)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = all_of(clust_vr),
                           values_from = sum_weights) %>%
        dplyr::mutate(name = name_simp)
    })
    
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
  ds_spots_metadata <- lapply(1:length(clust_vr_list), function(k){
    ds_spots_metadata <- purrr::map(purrr::map(ds_spots, 2),k) %>%
      dplyr::bind_rows() %>%
      data.frame()
    ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
    return(ds_spots_metadata)
  })
  
  
  # change column order so that its progressive
  ds_spots_metadata <- mapply(function(ds_spots_metadata,clust_vr){
    lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(sc_ref@meta.data[, clust_vr]))
    colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
    if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(sc_ref@meta.data[, clust_vr])) + 1)) {
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    } else {
      missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
      ds_spots_metadata[missing_cols] <- 0
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    }
    return(ds_spots_metadata)
  },ds_spots_metadata,clust_vr_list)
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}

## Tumor margin:  0.5X ~ 2X
# 30% - 60% epithelial
# 10% - 30% myeloid
# 10% - 30% stromal
# 0% - 20% CD4 T
# 0% - 20% CD8 T
# 0% - 20% B cell
# 0% - 20% Plasma
TME_simu_TM <- function(sc_ref,clust_vr,clust_vr_list,n = 50,verbose = TRUE,n.cell = 10,seeds = 20230508)
{
  # Check variables
  if (is(sc_ref) != "Seurat") stop("ERROR: sc_ref must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
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
  
  set.seed(seeds)
  # Save count matrix
  count_matr <- as.matrix(sc_ref@assays$RNA@counts)
  
  # Save celltype names
  celltype <- names(table(sc_ref[[clust_vr]]))
  
  ds_spots <- lapply(seq_len(n), function(i) {
    # epithelial
    n.epi <- sample(x = as.integer(0.3*n.cell):as.integer(0.6*n.cell), size = 1)
    n.mye <- sample(x = as.integer(0.1 * n.cell):as.integer(0.3 * n.cell), size = 1)
    n.str <- sample(x = as.integer(0.1 * n.cell):as.integer(0.3 * n.cell), size = 1)
    n.CD4 <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.CD8 <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.B <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.Pla <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    cell_pool <- NULL
    n.select.ct <- list(n.epi,n.mye,n.str,n.CD4,n.CD8,n.B,n.Pla)
    celltype.list <- list('Epithelial','Myeloid.cells','Stromal', 'CD4.T' ,'CD8.T','B.cells','Plasma')
    cell_pool <- unlist(mapply(function(x,y){
      selected_cells <- which(sc_ref[[clust_vr]]==y)
      selected_count_matr <- count_matr[,selected_cells]
      cell_pool <- sample(colnames(selected_count_matr),x)
    },n.select.ct,celltype.list))
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_matr) %in% cell_pool)
    tmp_ds <- sc_ref@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- lapply(clust_vr_list,function(clust_vr){
      tmp_ds %>%
        dplyr::select(all_of(clust_vr), weight) %>%
        # dplyr::mutate(clust_vr = paste("clust_",
        #                                tmp_ds[, clust_vr], sep = "")) %>%
        dplyr::group_by(!! sym(clust_vr)) %>%
        dplyr::summarise(sum_weights = sum(weight)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = all_of(clust_vr),
                           values_from = sum_weights) %>%
        dplyr::mutate(name = name_simp)
    })
    
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
  ds_spots_metadata <- lapply(1:length(clust_vr_list), function(k){
    ds_spots_metadata <- purrr::map(purrr::map(ds_spots, 2),k) %>%
      dplyr::bind_rows() %>%
      data.frame()
    ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
    return(ds_spots_metadata)
  })
  
  
  # change column order so that its progressive
  ds_spots_metadata <- mapply(function(ds_spots_metadata,clust_vr){
    lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(sc_ref@meta.data[, clust_vr]))
    colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
    if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(sc_ref@meta.data[, clust_vr])) + 1)) {
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    } else {
      missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
      ds_spots_metadata[missing_cols] <- 0
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    }
    return(ds_spots_metadata)
  },ds_spots_metadata,clust_vr_list)
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}


## Tumor stroma:  0.5X ~ 1.9X
# 0% - 20% epithelial
# 10% - 20% myeloid
# 20% - 50% stromal
# 10% - 30% CD4 T
# 10% - 30% CD8 T
# 0% - 20% B cell
# 0% - 20% Plasma
TME_simu_TS <- function(sc_ref,clust_vr,clust_vr_list,n = 50,verbose = TRUE,n.cell = 10,seeds = 20230508)
{
  # Check variables
  if (is(sc_ref) != "Seurat") stop("ERROR: sc_ref must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.numeric(n)) stop("ERROR: n must be an integer!")
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
  
  set.seed(seeds)
  # Save count matrix
  count_matr <- as.matrix(sc_ref@assays$RNA@counts)
  
  # Save celltype names
  celltype <- names(table(sc_ref[[clust_vr]]))
  
  ds_spots <- lapply(seq_len(n), function(i) {
    # epithelial
    n.epi <- sample(x = as.integer(0*n.cell):as.integer(0.2*n.cell), size = 1)
    n.mye <- sample(x = as.integer(0.1 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.str <- sample(x = as.integer(0.2 * n.cell):as.integer(0.5 * n.cell), size = 1)
    n.CD4 <- sample(x = as.integer(0.1 * n.cell):as.integer(0.3 * n.cell), size = 1)
    n.CD8 <- sample(x = as.integer(0.1 * n.cell):as.integer(0.3 * n.cell), size = 1)
    n.B <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    n.Pla <- sample(x = as.integer(0 * n.cell):as.integer(0.2 * n.cell), size = 1)
    cell_pool <- NULL
    n.select.ct <- list(n.epi,n.mye,n.str,n.CD4,n.CD8,n.B,n.Pla)
    celltype.list <- list('Epithelial','Myeloid.cells','Stromal', 'CD4.T' ,'CD8.T','B.cells','Plasma')
    cell_pool <- unlist(mapply(function(x,y){
      selected_cells <- which(sc_ref[[clust_vr]]==y)
      selected_count_matr <- count_matr[,selected_cells]
      cell_pool <- sample(colnames(selected_count_matr),x)
    },n.select.ct,celltype.list))
    
    # We're not going to sum the reads as is bc spots are **enriched**
    # so we'll add up the counts and downsample to the ~depth of a typical spot.
    
    # Create a name for the spot with the necessary info to deconvolute it
    pos <- which(colnames(count_matr) %in% cell_pool)
    tmp_ds <- sc_ref@meta.data[pos, ] %>% mutate(weight = 1)
    # tmp_ds[, "weight"] <- weigh
    name_simp <- paste("spot_", i, sep = "")
    
    spot_ds <- lapply(clust_vr_list,function(clust_vr){
      tmp_ds %>%
        dplyr::select(all_of(clust_vr), weight) %>%
        # dplyr::mutate(clust_vr = paste("clust_",
        #                                tmp_ds[, clust_vr], sep = "")) %>%
        dplyr::group_by(!! sym(clust_vr)) %>%
        dplyr::summarise(sum_weights = sum(weight)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = all_of(clust_vr),
                           values_from = sum_weights) %>%
        dplyr::mutate(name = name_simp)
    })
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
  ds_spots_metadata <- lapply(1:length(clust_vr_list), function(k){
    ds_spots_metadata <- purrr::map(purrr::map(ds_spots, 2),k) %>%
      dplyr::bind_rows() %>%
      data.frame()
    ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
    return(ds_spots_metadata)
  })
  
  
  # change column order so that its progressive
  ds_spots_metadata <- mapply(function(ds_spots_metadata,clust_vr){
    lev_mod <- gsub("[\\+|\\ |\\/]", ".", unique(sc_ref@meta.data[, clust_vr]))
    colnames(ds_spots_metadata) <- gsub("[\\+|\\ |\\/]", ".", colnames(ds_spots_metadata))
    if (sum(lev_mod %in% colnames(ds_spots_metadata)) == (length(unique(sc_ref@meta.data[, clust_vr])) + 1)) {
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    } else {
      missing_cols <- lev_mod[which(!lev_mod %in% colnames(ds_spots_metadata))]
      ds_spots_metadata[missing_cols] <- 0
      ds_spots_metadata <- ds_spots_metadata[, lev_mod]
    }
    return(ds_spots_metadata)
  },ds_spots_metadata,clust_vr_list)
  
  # Close progress bar
  close(pb)
  
  print(sprintf("Generation of %s test spots took %s mins", n,
                round(difftime(Sys.time(), start_gen, units = "mins"), 2)))
  print("output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot")
  return(list(topic_profiles = ds_syn_spots, cell_composition = ds_spots_metadata))
}

clust_vr_list <- c('low_res_ident','med_res_ident','high_res_ident')

simu_TC <- lapply(random.seed,function(seeds){
  silico <- TME_simu_TC(sc_ref = sc.silico,clust_vr = 'low_res_ident',clust_vr_list = clust_vr_list,seeds = seeds)
  return(silico)
})
saveRDS(simu_TC,'simu_TC.rds')
simu_TM <- lapply(random.seed,function(seeds){
  silico <- TME_simu_TM(sc_ref = sc.silico,clust_vr = 'low_res_ident',clust_vr_list = clust_vr_list,seeds = seeds)
  return(silico)
})
saveRDS(simu_TM,'simu_TM.rds')
simu_TS <- lapply(random.seed,function(seeds){
  silico <- TME_simu_TS(sc_ref = sc.silico,clust_vr = 'low_res_ident',clust_vr_list = clust_vr_list,seeds = seeds)
  return(silico)
})
saveRDS(simu_TS,'simu_TS.rds')


#### Transfer to scanpy data ####

sc.silico <- readRDS("sc.silico.RDS")
simu_TC <- readRDS("simu_TC.rds")
simu_TM <- readRDS("simu_TM.rds")
simu_TS <- readRDS("simu_TS.rds")

# sc_ref porocess
{
  sc.silico.scanpy <- readRDS("sc.silico.RDS")
  for(i in colnames(sc.silico.scanpy@meta.data))
  {
    if(length(sc.silico.scanpy@meta.data[[i]][is.na(sc.silico.scanpy@meta.data[[i]])]) != 0)
      sc.silico.scanpy@meta.data[[i]] <- NULL
  }
  
  DefaultAssay(sc.silico.scanpy) <- "RNA"
  sc.silico.scanpy@assays$RNA@meta.features$'GeneID-2' <- rownames(sc.silico.scanpy)
  sc.silico.scanpy@assays$RNA@meta.features$SYMBOL <- rownames(sc.silico.scanpy)
  sc.silico.scanpy@assays$RNA@data <- sc.silico.scanpy@assays$RNA@counts 
  
  sc.silico.scanpy.loom <- as.loom(x = sc.silico.scanpy, filename = "sc.silico.loom", verbose = T)
  sc.silico.scanpy.loom
  sc.silico.scanpy.loom$close_all()
  write.csv(sc.silico.scanpy@meta.data,'sc.silico_meta.data.csv')
  write.csv(sc.silico.scanpy@assays$RNA@meta.features,'sc.silico_meta.features.csv')
}

# st_vis process
{
  
  for(i in 1:5)
  {
    saveRDS(as.data.frame(simu_TC[[i]][[1]]),paste('simu_TC_',i,'.rds',sep = ''))
    saveRDS(as.data.frame(simu_TM[[i]][[1]]),paste('simu_TM_',i,'.rds',sep = ''))
    saveRDS(as.data.frame(simu_TS[[i]][[1]]),paste('simu_TS_',i,'.rds',sep = ''))
  }
  
}
