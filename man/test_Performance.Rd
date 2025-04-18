% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Test_Metrics_InSilico_v2.R
\name{test_Performance}
\alias{test_Performance}
\title{Test Performance In-Silico with Multi-metrics}
\usage{
test_Performance(deconv_result, synthetic_comp, metrics, evalu_level = "spot")
}
\arguments{
\item{deconv_result}{A matrix represents the deconvolution result of a method}

\item{synthetic_comp}{The ground truth matrix}

\item{metrics}{set the evaluation metrics using 'Pearson', 'RMSE', 'Spearman',or 'MultiIndex'.}

\item{evalu_level}{set the evaluation level using 'spot' or 'celltype'. By default is 'spot'.}
}
\value{
The performance of deconvolution result using given metrics.
}
\description{
This function can test the performance of deconvlution methods using: Pearson correlation, Spearman correlation, F1 score, RMSE, and so on.
}
