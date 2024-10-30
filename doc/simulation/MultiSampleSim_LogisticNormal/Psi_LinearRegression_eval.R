# Goal: This script is trying assess compositional change using gamma
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
require(pals)
library(MASS)
library(tibble)
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
# first read in scSTM data
files <- list.files(path = paste0(dir, "sims/"))
res <- vector(mode = "list")

for(file_name in files){
  # read in sims data
  set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
  
  ######## Pooled Version ####################
  # read in scSTMseq data
  scSTM_name <- paste0(dir, "scSTM_LinearRegression_noContent_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }
  
  psi <- scSTMobj$psi$alpha
  rownames(psi) <- unique(scSTMobj$sampleID)
  colnames(psi) <- paste0("topic_", 1:ncol(psi))
  # map group to Topics
  likelihood <- OneToOne_Mapping_Topics(scSTMobj)
  
  psi <- psi[,match(likelihood$Topic, colnames(psi))]
  
  trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  truePsi <- trueParam$psi
  colnames(truePsi) <- paste0("Group", 1:ncol(truePsi))
  rownames(truePsi) <- paste0("Sample", 1:nrow(truePsi))
  
  # truePsi <- truePsi %>% 
  #   as.data.frame() %>%
  #   rownames_to_column("Sample") 
  
  
  res[[set_level]] <- list(InferredPsi = psi, TrueTheta = truePsi)
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/NormalGamma_SingleResponse_LinearRegression_Null_nSample3_nCellType5_noBatch_StromalCell_Psi.rds")

# ############################# analysis #########################################
res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearRegression_Null_nSample3_nCellType5_noBatch_StromalCell_Psi.rds")

dat <- lapply(res, function(mat) {
  colnames(mat$InferredPsi) <- paste0("Group", 1:ncol(mat$InferredPsi)) 
  common_groups <- intersect(colnames(mat$InferredPsi), colnames(mat$TrueTheta))
  
  # Subtract InferredPsi from TrueTheta for common groups
  difference <- mat$InferredPsi[, common_groups] - mat$TrueTheta[, common_groups]
  difference <- difference %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  return(difference)
})

dat <- do.call(rbind, dat)
dat <- dat %>%
  pivot_longer(cols = starts_with("Group"), names_to = "Group", values_to = "Value")

ggplot(dat, aes(x = Group, y = Value)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_boxplot() +
  labs(title = "Difference Between Simulated Sample Random Variaiton Versus Estimated Psi Using Simple Linear Regression", x = "K", y = "Proportion Difference") +
  theme_minimal() +
  facet_grid(~Sample)

hist(dat$Value, 
     main = "Histogram of Difference between True Sample Variation versus Inferred Psi",
     xlab = "Difference")
