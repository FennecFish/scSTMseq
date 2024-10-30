setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(fastTopics)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample12_nCellType8_Batch_CancerCell/"
sims <- readRDS(paste0(dir, "sims/", file_name))

sims <- sims[rowSums(counts(sims)) != 0,]
counts <- t(counts(sims))
ngroup <- length(unique(sims$Group))
fit <- fit_topic_model(counts,k = ngroup)
saveRDS(fit, file = paste0(dir, "fastTopics/", 
                                   "fastTopics_", set_level, ".rds"))
