#' @title Calculate the entropy weight for the cGEPs - Version 9 - UCASpatial
#'
#' @description This function can simulate a Spatial data set with n spots using single cell data.
#'
#' @param object reference single cell RNA-seq data: Seurat object
#' @param features pre-calculated marker gene dataframe by Seurat FindAllMarkers, which should include p_val, avg_logFC/avg_log2FC
#' @param assay optional: set the default assay of object for calculating marker weights. By default is 'RNA'.
#' @param slot optional: set the default slot of object for calculating marker weights. By default is 'data'.
#' @param unit optional: method for calculate entropy. By default is 'log2'.
#'
#' @importFrom dplyr %>%
#' @return a list including 1: simulated spatial transcriptomics expression matrix; 2: spot-cell-type composition
#' @export

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
  suppressMessages(require(entropy))
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
