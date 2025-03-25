# This script is to evalulate different configuration of scSTM
# To find out why it has decreased performance compared to fastTopics
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(cluster)
library(monocle3)
library(sctransform)
library(fastTopics)
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

process_scSTM <- function(scSTMobj) {
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
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

cluster_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

################################################################################################################

dir = "MultiSample_VaryingBaseline_Batch0.5_CancerCell"
files <- list.files(path = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_Content_Prevalence_Time/"), pattern = "scSTM*")

res <- data.frame()

for(file_name in files){
  # create a new dataframe to calculate adjRandIndex
  dat <- data.frame()
  
  set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)
  # cat(set_level, "\n")
  
  file_path <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/sims/sims_", set_level, ".rds")
  sims <- safe_readRDS(file_path)
  if (is.null(sims)) {
    next  # Skip to the next iteration if reading the file failed
  }
  dat <- colData(sims) %>% data.frame() 
  
  ##### scSTM_noContent_noPrevalanence
  scSTM_nC_nP_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_noContent_noPrevalence/scSTM_",set_level,".rds")
  if(file.exists(scSTM_nC_nP_name)){
    scSTM_nC_nP <- readRDS(scSTM_nC_nP_name)
    scSTM_nC_nP_cluster <- process_scSTM(scSTM_nC_nP)
    dat$scSTM_nC_nP_cluster <- scSTM_nC_nP_cluster[match(rownames(dat), names(scSTM_nC_nP_cluster))]
  } else {dat$scSTM_nC_nP_cluster <- NA}
  
  ##### scSTM_noContent_Prevalanence_Time
  scSTM_nC_P_Time_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_noContent_Prevalence_Time/scSTM_",set_level,".rds")
  if(file.exists(scSTM_nC_P_Time_name)){
    scSTM_nC_P_Time <- readRDS(scSTM_nC_P_Time_name)
    scSTM_nC_P_Time_cluster <- process_scSTM(scSTM_nC_P_Time)
    dat$scSTM_nC_P_Time_cluster <- scSTM_nC_P_Time_cluster[match(rownames(dat), names(scSTM_nC_P_Time_cluster))]
  } else {dat$scSTM_nC_P_Time_cluster <- NA}
  
  ##### scSTM_noContent_Prevalanence_TimeandResponse
  scSTM_nC_P_TandR_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_noContent_Prevalence_TimeandResponse/scSTM_",set_level,".rds")
  if(file.exists(scSTM_nC_P_TandR_name)){
    scSTM_nC_P_TandR <- readRDS(scSTM_nC_P_TandR_name)
    scSTM_nC_P_TandR_cluster <- process_scSTM(scSTM_nC_P_TandR)
    dat$scSTM_nC_P_TandR_cluster <- scSTM_nC_P_TandR_cluster[match(rownames(dat), names(scSTM_nC_P_TandR_cluster))]
  } else {dat$scSTM_nC_P_TandR_cluster <- NA}
  
  ##### scSTM_Content_noPrevalanence
  scSTM_C_nP_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_Content_noPrevalence/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_nP_name)){
    scSTM_C_nP <- readRDS(scSTM_C_nP_name)
    scSTM_C_nP_cluster <- process_scSTM(scSTM_C_nP)
    dat$scSTM_C_nP_cluster <- scSTM_C_nP_cluster[match(rownames(dat), names(scSTM_C_nP_cluster))]
  } else {dat$scSTM_C_nP_cluster <- NA}
  
  ##### scSTM_Content_Prevalanence_Time
  scSTM_C_P_Time_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_Content_Prevalence_Time/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_P_Time_name)){
    scSTM_C_P_Time <- readRDS(scSTM_C_P_Time_name)
    scSTM_C_P_Time_cluster <- process_scSTM(scSTM_C_P_Time)
    dat$scSTM_C_P_Time_cluster <- scSTM_C_P_Time_cluster[match(rownames(dat), names(scSTM_C_P_Time_cluster))]
  } else {dat$scSTM_C_P_Time_cluster <- NA}
  
  ##### scSTM_Content_Prevalanence_TimeandResponse
  scSTM_C_P_TandR_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/scSTM_Content_Prevalence_TimeandResponse/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_P_TandR_name)){
    scSTM_C_P_TandR <- readRDS(scSTM_C_P_TandR_name)
    scSTM_C_P_TandR_cluster <- process_scSTM(scSTM_C_P_TandR)
    dat$scSTM_C_P_TandR_cluster <- scSTM_C_P_TandR_cluster[match(rownames(dat), names(scSTM_C_P_TandR_cluster))]
  } else {dat$scSTM_C_P_TandR_cluster <- NA}
  
  ###### fastTopics
  fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/fastTopics/fastTopics_", set_level, ".rds")
  if(file.exists(fasttopic_name)){
    nmf.sims <- readRDS(fasttopic_name)
    max_indices <- apply(nmf.sims$L, 1, which.max)
    fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    names(fastTopics_cluster) <- rownames(nmf.sims$L)
    dat$fastTopics_cluster <- fastTopics_cluster[match(rownames(dat), names(fastTopics_cluster))]
  } else{ dat$fastTopics_cluster <- NA }
  
  
  # monocle3
  monocle_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/", dir, "/monocle3/monocle3_", set_level, ".rds")
  if(file.exists(monocle_name)){
    monocle3 <- readRDS(monocle_name)
    dat$monocle3_cluster <- partitions(monocle3)[match(rownames(dat), names(partitions(monocle3)))]
  } else {dat$monocle3_cluster <- NA}
  
  adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("cluster")), function(x) {
    adjustedRandIndex(x, sims$Group)
  })
  # adjusted_rand_indices <- adjustedRandIndex(dat[, "scSTM_Content_Sample_cluster"], sims$Group)
  
  # Create a data frame to store results
  res.temp <- data.frame(
    sim = set_level,
    seed = unlist(strsplit(set_level, "_"))[1],
    effectSize = unlist(strsplit(set_level, "_"))[2],
    # scSTM_Content_Sample_cluster = adjusted_rand_indices
    t(adjusted_rand_indices)
    # t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  )
  
  res <- bind_rows(res, res.temp)
  cat(file_name, "\n")
  rm(sims)
}
# 
# pre_res <- read.csv("res/adjRandIndex_multiple_sample_varyingBaseline.csv")
# # res <- read.csv("res/adjRandIndex_multiple_sample_varyingBaseline_update.csv")
# # # # 
# # # res.update <- rbind(pre_res, res)
# res.update <- res %>%
#  # dplyr::select(sim, -matches("cluster$"), matches("cluster$")) %>%
#   dplyr::full_join(pre_res, by = c("sim", "seed", "effectSize")) %>%
#  dplyr::select(sim, seed, effectSize,
#                -matches("cluster$"), matches("cluster$"))

write.csv(res, file = paste0("res/scSTM_adjRandIndex_", dir, ".csv"))

