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
