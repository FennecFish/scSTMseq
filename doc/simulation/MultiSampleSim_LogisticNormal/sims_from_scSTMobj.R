# Goal: This script is to generate a single patient simulation
# From a previous run scSTMobj
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(SingleCellExperiment)
library(MASS)
library(VariantAnnotation)
library(checkmate)
library(MCMCpack)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")
gc <- as.integer(args[2])

select_top_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values, na.rm = T)
    if(length(which(all_values == max_value)) > 1){
      max_position_in_vector <- which(all_values == max_value)[1]
    } else{
      max_position_in_vector <- which(all_values == max_value)
    }
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  return(scSTMobj)
}
set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SinglePatient/nSample1_nCellType5_noBatch_StromalCell/"
scSTMobj <- readRDS(paste0(dir, "scSTM_LinearRegression_noContent_Prevalence_TimeandResponse/", file_name))
scSTMobj <- select_top_scSTM(scSTMobj)
theta <- scSTMobj$theta * 1e5
beta <- exp(scSTMobj$beta$logbeta[[1]])
counts <- theta %*% beta
counts <- round(counts/min(counts))
