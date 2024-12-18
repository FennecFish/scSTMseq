# this script is to run scSTMseq with filtered gene
# without content
setwd("~/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stm)
library(scater)
library(scran)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")
gc <- as.integer(args[2])

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  "sims_1723778736_pos_L1_c9_nsample12.rds")

sims <- readRDS("data/sims_1723778736_pos_L1_c9_nsample12.rds")
sims$Batch_ID <- paste0(sims$Sample, "_", sims$Time)
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=800)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
# scSTM.mod <- stm(sce = sims, K = ngroup, prevalence = ~time, init.type = "Spectral")
# scSTM.mod <- selectModel(sce = sims,
#                          K = ngroup, prevalence = ~Batch, content = ~Batch,
#                          N = 1, ts_runs = 1, random_run = 1,
#                          max.em.its = 100, net.max.em.its = 5) 
scSTM.mod <- selectModel_parallel(sce = sims,sample = "Sample",
                                  K = ngroup, prevalence = ~Time, content = ~Batch_ID,
                                  N = 2, ts_runs = 0, random_run = 1,
                                  max.em.its = 2, net.max.em.its = 1)

scSTM.mod <- selectModel(sce = sims,sample = "Sample",
                         K = ngroup, prevalence = ~Time, content = ~Batch_ID,
                         N = 2, ts_runs = 0, random_run = 1,
                         max.em.its = 2, net.max.em.its = 1)

scSTM.mod <- scSTMseq(sce = sims, sample = "Sample", #init.type = "TopicScore",
                                  K = ngroup, prevalence = ~Time, content = ~Batch_ID,
                                  max.em.its = 100, seed = 1034)

# msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
# saveRDS(scSTM.mod, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/scSTM_Content_Batch_Prevalence_Time/", 
#                                  "scSTM_", set_level, ".rds"))

