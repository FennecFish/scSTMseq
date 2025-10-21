# Goal: This script is to generate a single patient simulation
# both pre and post timepoints proportion is generated from a Logistic Normal
# Gamma is drawn from a multivariate normal
setwd("/proj/milovelab/wu/scSTMseq")
library(checkmate)
library(Rcpp)
library(SingleCellExperiment)
library(makeCluster)
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(MASS)
library(VariantAnnotation)
library(MCMCpack)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sims <- SimPairedData(nSample = 6, nGenes = 300,
                      nCell = c(100))

sourceCpp("src/STMCfuns.cpp")
ngroup = 10
scSTM.mod <- selectModel(sce = sims, sample = "Sample",
                         K = ngroup, prevalence = ~Time, content = NULL,
                         gamma.prior = "Pooled",
                         N = 1, ts_runs = 1, random_run = 1,
                         max.em.its = 1, net.max.em.its = 2)

scSTM.mod <- selectModel_parallel(sce = sims, sample = "Sample",
                                  K = ngroup, prevalence = ~Time, content = NULL,
                                  gamma.prior = "Pooled",
                                  N = 5, ts_runs = 30, random_run = 30,
                                  max.em.its = 100, net.max.em.its = 15, gc = 1)


