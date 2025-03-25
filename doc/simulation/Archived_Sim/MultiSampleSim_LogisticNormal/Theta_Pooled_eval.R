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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SinglePatient/nSample1_nCellType5_noBatch_StromalCell/"
# first read in scSTM data
files <- list.files(path = paste0(dir, "sims/"), pattern = "Null")
res <- vector(mode = "list")
nsims <- 1e6

for(file_name in files){
  # read in sims data
  set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
  
  ######## Pooled Version ####################
  # read in scSTMseq data
  # scSTM_name <- paste0(dir, "scSTM_Pooled_noSample_noContent_Prevalence_TimeandResponse/scSTM_",
  #                      set_level,".rds")
  scSTM_name <- paste0(dir, "scSTM_Pooled_noContent_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }
  
  # map group to Topics
  likelihood <- OneToOne_Mapping_Topics(scSTMobj)
  
  K <- scSTMobj$settings$dim$K
  sce <- scSTMobj$settings$sce
  theta <- scSTMobj$theta
  rownames(theta) <- scSTMobj$DocName
  theta_t1 <- theta[match(sce[,sce$Time == "Time1"]$Cell, rownames(theta)),]
  theta_t2 <- theta[match(sce[,sce$Time == "Time2"]$Cell, rownames(theta)),]
  metadata <- data.frame(
    Cell = sce$Cell,
    Sample = sce$Sample,
    Group = sce$Group,
    Timepoint = sce$Time
  )

  proportion_results <- data.frame()
  for (timepoint in unique(metadata$Timepoint)) {
    for (sample in unique(metadata$Sample)) {
      # Filter for the current Sample and Timepoint
      subset_data <- metadata %>%
        filter(Timepoint == timepoint & Sample == sample)
      
      if (nrow(subset_data) > 0) {
        theta_subset <- theta[rownames(theta) %in% subset_data$Cell, ]
        
        # Calculate proportions for each group
        proportions <- colSums(theta_subset) / sum(colSums(theta_subset))
        proportion_results <- rbind(proportion_results, data.frame(Timepoint = timepoint, Sample = sample, t(proportions)))
      }
    }
  }
  colnames(proportion_results) <- c("Timepoint", "Sample", paste0("topic_", 1:K))
  proportion_topics <- proportion_results[,3:ncol(proportion_results)]
  proportion_topics <- proportion_topics[,match(likelihood$Topic, colnames(proportion_topics))]
  proportion_results <- cbind(proportion_results[,1:2], proportion_topics)
  # proportion_results <- proportion_results %>%
  #   column_to_rownames("Timepoint")
  # # Loop through each unique Sample and Timepoint combination
  # for (sample in unique(scSTMobj$sampleID)) {
  #   for (timepoint in unique(metadata$Timepoint)) {
  #     # Filter for the current Sample and Timepoint
  #     subset_data <- metadata %>%
  #       filter(Sample == sample, Timepoint == timepoint)
  #     theta_subset <- theta[rownames(theta) %in% subset_data$Cell, ]
  # 
  #     # Calculate proportions for each group
  #     proportions <- colSums(theta_subset) / sum(colSums(theta_subset))
  #     proportion_results <- rbind(proportion_results, data.frame(Sample = sample, Timepoint = timepoint, t(proportions)))
  #   }
  # }
  # 
  # # Set appropriate column names for the proportion results
  # colnames(proportion_results)[3:ncol(proportion_results)] <- paste0("Group", 1:(ncol(proportion_results) - 2))
  
  trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  trueTheta <- trueParam$theta
  trueTheta <- do.call(rbind, trueTheta)
  colnames(trueTheta) <- paste0("Group", 1:K)
  trueTheta <- trueTheta %>% 
    as.data.frame() %>%
    rownames_to_column("RowName") %>%
    mutate(Sample = sub("_.*", "", RowName),
           Time = gsub("t", "Time", sub(".*_", "", RowName))) %>%
    dplyr::select(Time, Sample, Group1:Group5)
  
  res[[set_level]] <- list(InferredTheta = proportion_results, TrueTheta = trueTheta)
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/NormalGamma_SinglePatient/NormalGamma_nSample1_nCellType5_noBatch_StromalCell_Theta.rds")

# # ############################# analysis #########################################
# # # # res <- readRDS("res/composition_change/NormalGamma_ZeroPsi_nSample3_nCellType5_noBatch_StromalCell.rds")
# res <- readRDS("res/composition_change/NormalGamma_SinglePatient_Pooled_Null_nSample1_nCellType5_noBatch_StromalCell_Theta.rds")
# 
# dat <- lapply(res, function(mat) {
#   n <- ncol(mat$InferredTheta)
#   diff <- mat$TrueTheta[,3:n] - mat$InferredTheta[,3:n]
#   diff <- cbind(mat$TrueTheta[,1:2], diff)
#   # diff <- diff %>% rownames_to_column("Time")
#   return(diff)
# })
# 
# dat <- do.call(rbind, dat)
# dat <- dat %>%
#   pivot_longer(cols = starts_with("Group"), names_to = "Group", values_to = "Value")
# 
# ggplot(dat, aes(x = Group, y = Value, fill = Time)) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#   geom_boxplot() +
#   labs(title = "Difference Between Simulated Proportion Versus Estimated Proportion Using Normal-Gamma", x = "K", y = "Proportion Difference") +
#   theme_minimal() +
#   facet_grid(~Sample)
