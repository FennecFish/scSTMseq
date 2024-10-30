# Goal: This script is to generate multiple patient simulation
# the sample variation is from batch effect using Splatter
# each patient has two batches, representing the paired timepoints
# there is also batch effect among patients
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

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

################################################################################
# first start with small sample difference
# while the majority of the variation should come from cell type and treatment
# we use regular spat
sims_scRNA <- function(nsample, sampleSize, nGenes, nCellType, de.prob, 
                       de.facLoc, sample.facLoc, trt.facLoc, 
                       seed){
  
  first.prob <- rep(1/nCellType, nCellType)
  neg.prob <- rep(1/nCellType, nCellType)
  
  if(nCellType == 5){
    pos.prob <- first.prob + c(-0.15, 0.1, 0.2, -0.15, 0)
  }
  if(nCellType == 9){
    pos.prob <- first.prob + c(0.1, -0.05, -0.08, 0.03, -0.05, 0.02, 0.02, 0.01, 0)
  }
  if(nCellType == 13){
    pos.prob <- first.prob + c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
                               -0.04, 0)
  }
  
  nGenes <- nGenes
  de.facScale = 0.1
  
  # we use batch as samples
  batchCells <- rep(sampleSize, nsample*2)
  batch.facLoc = rep(c(sample.facLoc, trt.facLoc), each = nsample)
  
  params <- newSplatParams()
  params <- setParams(params, nGenes = nGenes,
                      group.prob = list(first.prob, pos.prob),
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      batch.facLoc = batch.facLoc, batchCells=batchCells,
                      seed = seed)
  sims <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = FALSE)
  
  sims$Time <- ifelse(sims$Batch %in% paste0("Batch", 1:nsample),
                      "Time1", "Time2")
  sim.sample.number <- as.numeric(sub("Batch", "", sims$Batch))
  sim.sample.number <- ifelse(sim.sample.number > nsample,
                              sim.sample.number - nsample, sim.sample.number)
  sims$Sample <- paste0("Sample", sim.sample.number)
  
  return(sims)
  # df <- colData(s) %>% 
  #   as.data.frame() %>% 
  #   count(Batch, Group, time) %>%
  #   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0))
  # total_1 <- sum(df$`1`)
  # total_2 <- sum(df$`2`)
  # result <- df %>%
  #   mutate(
  #     Proportion_1 = `1` / total_1,
  #     Proportion_2 = `2` / total_2,
  #     Ratio = Proportion_1 / Proportion_2
  #   )
}

# sims <- logNormCounts(sims)
# sims <- runPCA(sims)
# plotPCA(sims, colour_by = "Group")
# plotPCA(sims, colour_by = "Sample")
# plotPCA(sims, colour_by = "Time")

nsample <- c(3, 6, 12)
nCellType <- c(5, 9, 13)

nparam <- expand.grid(nsample, nCellType)
colnames(nparam) <- c("nsample", "nCellType")

for(i in 1:dim(nparam)[1]){
  
  nsample <- nparam[i,1]
  nCellType <- nparam[i,2]
  nGenes <- 1000
  sampleSize <- 300

  #################################### level 1 ########################################
  ## cell Type  variation dominates
  de.prob <- 0.3
  de.facLoc <- 1
  sample.facLoc <- 0.1
  trt.facLoc <- 0.1
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)

  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/",
                                   "sims_", seed, "_pos_L1_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L1\n")
  rm(sims)
  
  #################################### level 2 ########################################
  ## large cell type variation, and increased sample batch effect
  de.prob <- 0.3
  de.facLoc <- 1
  sample.facLoc <- 0.5
  trt.facLoc <- 0.1
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)
  
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/",
                              "sims_", seed, "_pos_L2_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L2\n")
  rm(sims)
  
  #################################### level 3 ########################################
  ## large cell type variation, with increased treatment and sample effect
  de.prob <- 0.3
  de.facLoc <- 1
  sample.facLoc <- 0.5
  trt.facLoc <- 0.3
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)
  
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/",
                              "sims_", seed, "_pos_L3_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L3 \n")
  rm(sims)
  
  #################################### level 4 ########################################
  ## smaller cell type variation, with increased treatment and sample effect
  de.prob <- 0.3
  de.facLoc <- 0.5
  sample.facLoc <- 0.5
  trt.facLoc <- 0.3
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)
  
  saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/",
                              "sims_", seed, "_pos_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4 \n")
  rm(sims)
}



############### code to fix issue in the current simulation ######
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/", pattern = "sims*")

safe_readRDS <- function(file_path) {
  tryCatch({
    # Attempt to read the RDS file
    data <- readRDS(file_path)
    return(data)
  }, error = function(e) {
    # Handle the error
    message(paste("Error reading RDS file:", file_path))
    message("Skipping to the next file.")
    return(NULL)  # Return NULL if an error occurs
  })
}

for(file_name in files){
  dat <- data.frame()
  file_path <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/", file_name)
  sims <- safe_readRDS(file_path)
  
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  nsample <- as.numeric(str_extract(set_level, "(?<=nsample)\\d+"))
  sim.sample.number <- as.numeric(sub("Batch", "", sims$Batch))
  sim.sample.number <- ifelse(sim.sample.number > nsample,
                              sim.sample.number - nsample, sim.sample.number)
  sims$Sample <- paste0("Sample", sim.sample.number)
  saveRDS(sims, file_path)
  cat(file_name, "\n")
}
