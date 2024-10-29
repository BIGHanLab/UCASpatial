#' @title Generate simulated spatial transcriptomics data - UCASpatial
#'
#' @description This function can simulate a Spatial data set with n spots using single cell data.
#'
#' @param sc_ref reference single cell RNA-seq data: Seurat object
#' @param clust_vr colname of the meta.data in sc_ref which include the cluster information
#' @param n optional: set the spots number for generation. By default is '1000'.
#' @param max.celltype optional: set the max cell type included for a spot (randomly select from min~max). By default is '4'.
#' @param max.cell optional: set the max cell number included for a spot (randomly select from min~max). By default is '5'.
#' @param min.celltype optional: set the min cell type included for a spot (randomly select from min~max). By default is '1'.
#' @param min.cell optional: set the min cell number included for a spot (randomly select from min~max). By default is '1'.
#'
#' @importFrom dplyr %>%
#' @importFrom DropletUtils downsampleMatrix
#' @return a list including 1: simulated spatial transcriptomics expression matrix; 2: spot-cell-type composition
#' @export


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





