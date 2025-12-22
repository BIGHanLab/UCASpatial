#' @title Test Performance In-Silico with Multi-metrics
#'
#' @description This function can test the performance of deconvlution methods using: Pearson correlation, Spearman correlation, F1 score, RMSE, and so on.
#' @importFrom dplyr %>%
#'
#' @param deconv_result A matrix represents the deconvolution result of a method
#' @param synthetic_comp The ground truth matrix
#' @param metrics set the evaluation metrics using 'Pearson', 'RMSE', 'Spearman',or 'MultiIndex'.
#' @param evalu_level set the evaluation level using 'spot' or 'celltype'. By default is 'spot'.
#'
#' @return The performance of deconvolution result using given metrics.
#' @export


test_Performance <- function(deconv_result, synthetic_comp,metrics,evalu_level = 'spot')
{
  # metrics == 'Pearson'
  if(metrics == 'Pearson')
  {
    if(evalu_level == 'spot')
    {
      return(test_Pearson_spotlevel(deconv_result,synthetic_comp))
    }
    else
      return(test_Pearson(deconv_result,synthetic_comp))
  }
  # metrics == 'Spearman'
  if(metrics == 'Spearman')
  {
    if(evalu_level == 'spot')
    {
      return(test_Spearman_spotlevel(deconv_result,synthetic_comp))
    }
    else
      return(test_Spearman(deconv_result,synthetic_comp))
  }
  # metrics == 'RMSE'
  if(metrics == 'RMSE')
  {
    if(evalu_level == 'spot')
    {
      return(test_RMSE_spotlevel(deconv_result,synthetic_comp))
    }
    else
      return(test_RMSE(deconv_result,synthetic_comp))
  }
  # metrics == 'MultiIndex'
  if(metrics == 'MultiIndex')
  {
    if(evalu_level == 'spot')
    {
      stop("ERROR: 'evalu_level' for MultiIndex was only available for 'celltype'!\n")
    }
    # return(test_MultiIndex(deconv_result,synthetic_comp))
    return(Test_Acurracy(deconv_result,synthetic_comp))
  }
}



test_MultiIndex <- function(deconv_result, synthetic_comp)
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

Test_Acurracy<- function (deconv_result, synthetic_comp)
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
  tp <- 0
  tn <- 0
  fp <- 0
  fn <- 0
  for (i in seq_len(nrow(synthetic_comp))) {
    x <- rbind(synthetic_comp[i, ], deconv_result[i,])
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
  result <- cbind(accuracy, specificity, precision, recall,F1)
  colnames(result) <- c('Accuracy','Specificity','precision','recall','F1_score')
  return(result)
}

test_Spearman <- function(deconv_result,synthetic_comp){
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = ncol(deconv_result2))
  for(i in 1:ncol(deconv_result2))
  {
    t <- cbind(deconv_result2[,i],synthetic_comp[,i])
    if(sum(synthetic_comp[,i])==0)
      result_list[i] <- NA
    else
    {
      dataSpearman <- round(cor(t,method = c("spearman")) , 2)
      if(is.na(dataSpearman[1,2]))
        result_list[i] <- 0
      else
        result_list[i] <- dataSpearman[1,2]
    }

  }
  names(result_list) <- colnames(synthetic_comp)
  return(result_list)
}

test_Spearman_spotlevel <- function(deconv_result,synthetic_comp){
  deconv_result2 <- deconv_result[, colnames(synthetic_comp)]
  result_list <- vector("list", length = nrow(deconv_result2))
  for(i in 1:nrow(deconv_result2))
  {
    t <- cbind(deconv_result2[i,],synthetic_comp[i,])
    dataSpearman <- round(cor(t,method = c("spearman")) , 2)
    if(is.na(dataSpearman[1,2]))
      result_list[i] <- 0
    else
      result_list[i] <- dataSpearman[1,2]
  }
  names(result_list) <- rownames(synthetic_comp)
  return(result_list)
}
