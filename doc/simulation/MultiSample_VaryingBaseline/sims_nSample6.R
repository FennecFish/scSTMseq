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
library(MCMCpack)

# Each patient has a different starting baseline
# for simplicity, only two cellType changed, one + one -
# within the simulation, there is Respond and nonRespond patients
# R patients have change, the NE patients do not have changes
# the change depends on the patient's baseline

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

########################## some functions ######################################
# the following function is to generate post treatment proportion
post_proportion_generation <- function(baseline_proportion, effectSize){
  # order group 1 from max to min
  sorted_indices <- order(baseline_proportion[,1], decreasing = T)
  
  # Split the dataset into the top 3 and bottom 3 rows based on group 1
  top3_indices <- sorted_indices[1:(nrow(baseline_proportion)/2)]
  bottom3_indices <- sorted_indices[(nrow(baseline_proportion)/2 + 1):nrow(baseline_proportion)]
  
  
  post_proportion <- baseline_proportion
  for (i in top3_indices) {
    # calculate group 2 decrease. Think it as cancer cell decrease due to trt
    col2_decrease <- post_proportion[i, 2] * effectSize
    # Increase column 1 to keep the row sum to 1
    col1_increase <- post_proportion[i, 2] * (1 - effectSize)
    
    
    post_proportion[i, 2] <- col2_decrease
    post_proportion[i, 1] <- post_proportion[i, 1] + col1_increase
  }
  response_labels <- rep("NR", nrow(baseline_proportion)) 
  response_labels[top3_indices] <- "R"
  return(list(post_proportion = post_proportion, response_label = response_labels))
}

# Write a function to check sims proportion change
proportion_check <- function(sims){
  sims_df <- as.data.frame(colData(sims))
  
  # Filter for Time1 and Time2
  sims_time1 <- sims_df[sims_df$Time == "Time1", ]
  sims_time2 <- sims_df[sims_df$Time == "Time2", ]
  
  # Count the number of cells in each group for Time1 and Time2
  group_counts_time1 <- table(sims_time1$Sample, sims_time1$Group)
  group_counts_time2 <- table(sims_time2$Sample, sims_time2$Group)
  
  # Convert the counts to proportions (relative frequency)
  prop_time1 <- prop.table(group_counts_time1, margin = 1)
  prop_time2 <- prop.table(group_counts_time2, margin = 1)
  
  # Calculate the proportion change for each group between Time1 and Time2
  proportion_change <- prop_time2 - prop_time1
  
  return(proportion_change)
  # Convert to data frame for better readability
  # proportion_change_df <- as.data.frame(as.table(proportion_change))
}

# take the same baseline_proportion but varying effectSize to simulate data
delta_sim <- function(baseline_proportion, post_proportion,
                      label = label, effectSize, seed = seed, nSample = nSample, 
                      nGenes = nGenes, nCellType = nCellType, 
                      de.prob = de.prob, de.facLoc = de.facLoc, 
                      batchCells = batchCells, batch.facLoc = batch.facLoc,
                      cancerCellGroup = NULL, cancerCell.facLoc = NULL, cancerCell.facScale = NULL, 
                      batch.rmEffect = TRUE){
  params <- newSplatParams()
  params <- setParams(params, nGenes = nGenes,
                      group.prob = list(baseline_proportion, post_proportion),
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      batchCells=batchCells, batch.facLoc = batch.facLoc,
                      seed = seed)
  sims <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = batch.rmEffect,
                        sampleLabel = label, cancerCellGroup = cancerCellGroup, 
                        cancerCell.facLoc = cancerCell.facLoc, cancerCell.facScale = cancerCell.facScale)
  return(sims)
}

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType5_Batch_CancerCell/sims/"
batch.rmEffect = FALSE
nSample = 6
nGenes = 3000
nCellType = 5
de.prob <- 0.3
de.facLoc <- 0.5
cancerCellGroup = 2
cancerCell.facLoc = de.facLoc
cancerCell.facScale = 0.2
batchCells <- rep(c(500,500), each = nSample)
batch.facLoc <- runif(nSample, min = 0, max = 0.5)
batch.facLoc <- rep(batch.facLoc, times =2)
# generate baseline for every sample
baseline_proportion <- rdirichlet(nSample, alpha = rep(1, nCellType))

if(batch.rmEffect){save_batch <- "noBatch"} else{save_batch <- "Batch"}
if(is.null(cancerCellGroup)){save_cancer <- "StromalCell"} else{save_cancer <- "CancerCell"}
# No treatment effect
effectSize = 1
post_trt <- post_proportion_generation(baseline_proportion, effectSize = effectSize)
post_proportion <- post_trt$post_proportion
label <- post_trt$response_label

# r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
# sapply(r.file, source)
sims <- delta_sim(baseline_proportion, post_proportion,
                  label = label, effectSize = effectSize, 
                  seed = seed, nSample = nSample, 
                  nGenes = nGenes, nCellType = nCellType, 
                  de.prob = de.prob, de.facLoc = de.facLoc,
                  batchCells = batchCells, batch.facLoc = batch.facLoc,
                  cancerCellGroup = cancerCellGroup, cancerCell.facLoc = cancerCell.facLoc, cancerCell.facScale = cancerCell.facScale, 
                  batch.rmEffect = batch.rmEffect)
saveRDS(sims, file = paste0(dir, "MultiSample_VaryingBaseline_nSample", nSample, 
                            "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/",
                                 "sims_", seed, "_effectSize", effectSize, ".rds"))
cat("Generation Completed for effectSize", effectSize, "\n")
rm(sims)

# small treatment effect, cancer decreased 0.9
effectSize = 0.7
post_trt <- post_proportion_generation(baseline_proportion, effectSize = effectSize)
post_proportion <- post_trt$post_proportion
label <- post_trt$response_label

sims <- delta_sim(baseline_proportion, post_proportion,
                  label = label, effectSize = effectSize, 
                  seed = seed, nSample = nSample, 
                  nGenes = nGenes, nCellType = nCellType, 
                  de.prob = de.prob, de.facLoc = de.facLoc,
                  batchCells = batchCells, batch.facLoc = batch.facLoc,
                  cancerCellGroup = cancerCellGroup, cancerCell.facLoc = cancerCell.facLoc, cancerCell.facScale = cancerCell.facScale, 
                  batch.rmEffect = batch.rmEffect)
saveRDS(sims, file = paste0(dir, "MultiSample_VaryingBaseline_nSample", nSample, 
                            "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/",
                            "sims_", seed, "_effectSize", effectSize, ".rds"))
cat("Generation Completed for effectSize", effectSize, "\n")
rm(sims)

# medium treatment effect, cancer decreased 0.6
effectSize = 0.4
post_trt <- post_proportion_generation(baseline_proportion, effectSize = effectSize)
post_proportion <- post_trt$post_proportion
label <- post_trt$response_label

sims <- delta_sim(baseline_proportion, post_proportion,
                  label = label, effectSize = effectSize, 
                  seed = seed, nSample = nSample, 
                  nGenes = nGenes, nCellType = nCellType, 
                  de.prob = de.prob, de.facLoc = de.facLoc,
                  batchCells = batchCells, batch.facLoc = batch.facLoc,
                  cancerCellGroup = cancerCellGroup, cancerCell.facLoc = cancerCell.facLoc, cancerCell.facScale = cancerCell.facScale, 
                  batch.rmEffect = batch.rmEffect)
saveRDS(sims, file = paste0(dir, "MultiSample_VaryingBaseline_nSample", nSample, 
                            "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/",
                            "sims_", seed, "_effectSize", effectSize, ".rds"))
cat("Generation Completed for effectSize", effectSize, "\n")
rm(sims)

# large treatment effect, cancer decreased by 0.1
effectSize = 0.1
post_trt <- post_proportion_generation(baseline_proportion, effectSize = effectSize)
post_proportion <- post_trt$post_proportion
label <- post_trt$response_label

sims <- delta_sim(baseline_proportion, post_proportion,
                  label = label, effectSize = effectSize, 
                  seed = seed, nSample = nSample, 
                  nGenes = nGenes, nCellType = nCellType, 
                  de.prob = de.prob, de.facLoc = de.facLoc,
                  batchCells = batchCells, batch.facLoc = batch.facLoc,
                  cancerCellGroup = cancerCellGroup, cancerCell.facLoc = cancerCell.facLoc, cancerCell.facScale = cancerCell.facScale, 
                  batch.rmEffect = batch.rmEffect)
saveRDS(sims, file = paste0(dir, "MultiSample_VaryingBaseline_nSample", nSample, 
                            "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/",
                            "sims_", seed, "_effectSize", effectSize, ".rds"))
cat("Generation Completed for effectSize", effectSize, "\n")
rm(sims)

# sims <- sims[,sims$Group == "Group2"]
# sims <- logNormCounts(sims)
# sims <- runPCA(sims)
# plotPCA(sims, colour_by = "Sample", shape_by = "Time")
# plotPCA(sims, colour_by = "Group", shape_by = "Batch")
# plotPCA(sims, colour_by = "Group")



# sim_change <- proportion_check(sims)
# true_change <- post_proportion - baseline_proportion
