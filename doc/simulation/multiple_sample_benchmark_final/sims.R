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
  
  if(nCellType == 17){
    pos.prob <- first.prob + c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
                               -0.04, 0.01, 0.02, -0.03, 0, 0)
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
  sims.pos <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = FALSE)
  
  sims.pos$Time <- ifelse(sims.pos$Batch %in% paste0("Batch", 1:nsample),
                      "Time1", "Time2")
  sim.pos.sample.number <- as.numeric(sub("Batch", "", sims.pos$Batch))
  sim.pos.sample.number <- ifelse(sim.pos.sample.number > nsample,
                                  sim.pos.sample.number - nsample, sim.pos.sample.number)
  sims.pos$Sample <- paste0("Sample", sim.pos.sample.number)
  
  params <- newSplatParams()
  params <- setParams(params, nGenes = nGenes,
                      group.prob = list(first.prob, neg.prob),
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      batch.facLoc = batch.facLoc, batchCells=batchCells,
                      seed = seed)
  sims.neg <- splatSimulate(params, method = "groups",
                            verbose = TRUE, batch.rmEffect = FALSE)
  
  sims.neg$Time <- ifelse(sims.neg$Batch %in% paste0("Batch", 1:nsample),
                          "Time1", "Time2")
  sim.neg.sample.number <- as.numeric(sub("Batch", "", sims.neg$Batch))
  sim.neg.sample.number <- ifelse(sim.neg.sample.number > nsample,
                                  sim.neg.sample.number - nsample, sim.neg.sample.number)
  sims.neg$Sample <- paste0("Sample", sim.neg.sample.number)
  
  return(list(sims.pos, sims.neg))
  
}

# sims <- logNormCounts(sims)
# sims <- runPCA(sims)
# plotPCA(sims, colour_by = "Group")
# plotPCA(sims, colour_by = "Sample")
# plotPCA(sims, colour_by = "Time")

nsample <- c(6, 12, 24, 36, 48)
nCellType <- c(5, 9, 13, 17)

nparam <- expand.grid(nsample, nCellType)
colnames(nparam) <- c("nsample", "nCellType")

for(i in 1:dim(nparam)[1]){
  
  nsample <- nparam[i,1]
  nCellType <- nparam[i,2]
  nGenes <- 2000
  sampleSize <- 500

  #################################### level 1 ########################################
  ## cell Type  variation dominates
  de.prob <- 0.3
  de.facLoc <- 1
  sample.facLoc <- 0.1
  trt.facLoc <- 0.1
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)

  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                                   "sims_", seed, "_pos_L1_c", nCellType, "_nsample", nsample, ".rds"))
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                                   "sims_", seed, "_neg_L1_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L1\n")
  rm(sims)
  
  # proportion <- colData(sims[[2]]) %>%
  #   as.data.frame() %>%
  #   group_by(Sample, Time, Group) %>%
  #   summarise(count = n()) %>%
  #   mutate(proportion = count / sum(count)) %>%
  #   arrange(Sample, Group, Time) %>%
  #   dplyr::select(-count) %>%
  #   pivot_wider(names_from = Time, values_from = proportion) %>%
  #   mutate(ratio = Time1 / Time2)

  #################################### level 2 ########################################
  ## large cell type variation, and increased sample batch effect
  de.prob <- 0.3
  de.facLoc <- 1
  sample.facLoc <- 0.5
  trt.facLoc <- 0.1
  sims <- sims_scRNA(nsample, sampleSize, nGenes, nCellType, de.prob, 
                     de.facLoc, sample.facLoc, trt.facLoc, 
                     seed)
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                              "sims_", seed, "_pos_L2_c", nCellType, "_nsample", nsample, ".rds"))
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                                   "sims_", seed, "_neg_L2_c", nCellType, "_nsample", nsample, ".rds"))
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
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                              "sims_", seed, "_pos_L3_c", nCellType, "_nsample", nsample, ".rds"))
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                                   "sims_", seed, "_neg_L3_c", nCellType, "_nsample", nsample, ".rds"))
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
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                              "sims_", seed, "_pos_L4_c", nCellType, "_nsample", nsample, ".rds"))
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                                   "sims_", seed, "_neg_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4 \n")
  rm(sims)
}

# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
#                     pattern = "^sims.*pos_L4_c5_nsample6\\.rds$")
# file_name <- files[1]
# sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/", file_name))
# sims <- logNormCounts(sims)
# sims <- runPCA(sims)
# # plotPCA(sims, colour_by = "Sample", shape_by = "Group", size = "Time")
# # plotPCA(sims, colour_by = "Sample")
# # plotPCA(sims, colour_by = "Time")
# 
# pca_data <- sims@int_colData@listData$reducedDims@listData$PCA %>% as.data.frame()
# pca_data$Sample <- sims$Sample
# pca_data$Group <- sims$Group
# pca_data$Time <- sims$Time
# 
# png("res/PCA_splatter_Sample_Group.png", width = 2500, height = 2000, res = 300)
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample, shape = Group)) +
#   geom_point() +
#   labs(title = "PCA Plot Simulating Technical Effects Using Splatter", 
#        subtitle = "With Sample and Group Effect",
#        x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1") +
#   theme(
#     plot.title = element_text(size = 18, face = "bold"),  # Title font size
#     axis.title = element_text(size = 14),  # Axis titles font size
#     axis.text = element_text(size = 14),   # Axis text font size
#     legend.title = element_text(size = 14),  # Legend title font size
#     legend.text = element_text(size = 12)    # Legend text font size
#   )
# dev.off()
# 
# png("res/PCA_splatter_Time_Group.png", width = 2500, height = 2000, res = 300)
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Time, shape = Group)) +
#   geom_point() +
#   labs(title = "PCA Plot Simulating Technical Effects Using Splatter", 
#        subtitle = "With Group and Time Effect",
#        x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1") +
#   theme(
#     plot.title = element_text(size = 18, face = "bold"),  # Title font size
#     axis.title = element_text(size = 14),  # Axis titles font size
#     axis.text = element_text(size = 14),   # Axis text font size
#     legend.title = element_text(size = 14),  # Legend title font size
#     legend.text = element_text(size = 12)    # Legend text font size
#   )
# dev.off()
