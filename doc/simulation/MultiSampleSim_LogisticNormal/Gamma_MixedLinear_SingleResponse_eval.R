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
  
  # read in scSTMseq data
  scSTM_name <- paste0(dir, "scSTM_LinearMixed_Content_Sample_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }

  formula = ~Time + (1|Sample) # mixed linear model formula
  res <- GammaError(model, formula = formula, nsims = 100)
  saveRDS(res, file = paste0(dir, "GammaError_MixedModel/GammaError_",set_level,".rds"))
  
  # 
  # # map group to Topics
  # likelihood <- OneToOne_Mapping_Topics(scSTMobj)
  # 
  # Gamma.model <- scSTMobj$mu$lm.model
  # Gamma.est <- scSTMobj$mu$gamma
  # Gamma.conf <- lapply(Gamma.model, function(res){
  #   conf.int <- confint(res)
  #   conf.int <- conf.int[4:nrow(conf.int),]
  #   return(conf.int)
  # })
  # 
  # gamma.conf <- lapply(seq_along(Gamma.conf), function(i){
  #   mat <- Gamma.conf[[i]]
  #   mat <- t(mat) %>%
  #     as.data.frame() %>% 
  #     mutate(Estimate = Gamma.est[-1,i]) %>%
  #     mutate(Topic = paste0("topic_", i ))
  #   colnames(mat) <- c("lowBound", "HighBound", "Estimate", "Group")
  #   return(mat)
  # })
  # 
  # names(gamma.conf) <- paste0("topic_", 1:length(gamma.conf))
  # gamma.conf <- gamma.conf[match(likelihood$Topic, names(gamma.conf))]
  # names(gamma.conf) <- paste0("Group", 1:length(gamma.conf))
  # 
  # trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  # true_gamma <- as.data.frame(scSTMobj$settings$sce@metadata$TrueParams$gamma)
  # colnames(true_gamma) <- gsub("K", "Group", colnames(true_gamma))
  # 
  # Type = sapply(strsplit(set_level, "_"), `[`, 2)
  # 
  # if(Type == "Null"){
  #   sim_df <- lapply(names(gamma.conf), function(name){
  #     mat <- gamma.conf[[name]] %>% as.data.frame() 
  #     mat <- mat %>% mutate(True_Gamma = 0)
  #     return(mat)
  #   })
  # }else{
  #   sim_df <- lapply(names(gamma.conf), function(name){
  #     mat <- gamma.conf[[name]] %>% as.data.frame() 
  #     if(nrow(mat) != 0 & name %in% colnames(true_gamma)) {mat <- mat %>% mutate(True_Gamma = true_gamma[,name])}
  #     return(mat)
  #   })
  # }
  # names(sim_df) <- paste0("Group", 1:length(sim_df))
  # 
  # res[[set_level]] <- sim_df
  cat(file_name, "\n")
  # rm(scSTMobj)
  # rm(Gamma.conf)
  # rm(gamma.conf)
  # rm(sim_df)
}

saveRDS(res, file = "res/composition_change/NormalGamma_SingleResponse_LinearMixed_nSample3_nCellType5_noBatch_CancerCell_Gamma.rds")

# # ############################# analysis #########################################
# res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_nSample3_nCellType5_Batch_StromalCell_Gamma.rds")
# # res_null <- res[!grepl("Null", names(res))]
# 
# dat <- lapply(res, function(mat) {
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
#   mat <- do.call(rbind, mat)
#   # colnames(mat) <- c("lowBound", "HighBound")
#   mat <- mat %>%
#     as.data.frame() %>%
#     mutate(exist_change = ifelse(lowBound < 0 & HighBound > 0, FALSE, TRUE))
#   if(sum(mat$exist_change, na.rm = T) == 0){exist_change = FALSE} else{exist_change = TRUE}
#   return(exist_change)
# })
# 
# dat <-do.call(rbind, dat)
# dat <- dat %>% as.data.frame() %>% mutate(Type = sapply(strsplit(rownames(dat), "_"), `[`, 2))
# colnames(dat) <- c("exist_change", "Type")
# 
# ### Type I error ####
# dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# 
# 
# #### if we look at each topic individually #####
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
#   mat <- do.call(rbind, mat)
#   # colnames(mat) <- c("lowBound", "HighBound")
#   mat <- mat %>%
#     as.data.frame() %>%
#     mutate(exist_change = ifelse(lowBound < 0 & HighBound > 0, FALSE, TRUE))
#   # Add Type from the provided name
#   mat$Type <- strsplit(name, "_")[[1]][2]
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# 
# ### Type I error ####
# dat  %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE, na.rm = T))
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE, na.rm = T))
# 
# #### Check if the True Gamma is Within 95% Interval #####
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat %>%
#     as.data.frame() %>%
#     mutate(Recovered = ifelse(lowBound < True_Gamma & HighBound > True_Gamma, TRUE, FALSE))
#   # Add Type from the provided name
#   mat$Type <- strsplit(name, "_")[[1]][2]
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(Recovered == TRUE, na.rm = T))
# 
# #### Check if the True Gamma is in the same direction of the estimation #####
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat %>%
#     as.data.frame() %>%
#     mutate(Type = strsplit(name, "_")[[1]][2]) %>%
#     filter(Type == "HighVar") %>%
#     mutate(Recovered = ifelse(sign(Estimate) == sign(True_Gamma), TRUE, FALSE))
#   # Add Type from the provided name
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(Recovered == TRUE, na.rm = T))
