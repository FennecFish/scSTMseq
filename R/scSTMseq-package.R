#' Single Cell Structural Topic Model Sequencing
#'
#' This package implements the Structural Topic Model, adapted to single cell RNAseq data,
#' a general approach to including cell-level metadata within mixed-membership topic
#' models. This package is built upon package stm, and borrowed functions from Splatter
#'
#' Functions to simulate single cell data: \code{\link{SimPairedData}}
#'
#' Functions to fit the model: \code{\link{scSTM}} \code{\link{selectModel}}
#' \code{\link{selectModel_parallel}}
#'
#' Functions to summarize a model: \code{\link{labelTopics}}
#' \code{\link{summary.scSTM}}
#'
#' Functions for Post-Estimation: \code{\link{estimateEffect}}
#' \code{\link{topicCorr}}
#'
#' Plotting Functions: \code{\link{plot.scSTM}} \code{\link{structure_plot}}
#' \code{\link{plot.estimateEffect}}
#'
#'
#' @name scSTMseq-package
#' @docType _PACKAGE
#' @author Author: Euphy Wu, Didong Li, Naim Rashid
#' @keywords _PACKAGE
#'
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline axis hist legend lines par plot points segments smoothScatter text title
#' @importFrom stats aggregate as.formula coef cor cov lm loess median model.frame model.response model.matrix na.omit optim optimize pchisq predict quantile rbinom rgamma rnorm runif terms
#' @useDynLib scSTMseq, .registration = TRUE
NULL
