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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_MultiResponse/nSample6_nCellType5_Batch_StromalCell/"
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

  # map group to Topics
  likelihood <- OneToOne_Mapping_Topics(scSTMobj)
  
  Gamma.model <- scSTMobj$mu$lm.model
  Gamma.est <- scSTMobj$mu$gamma
  Gamma.conf <- lapply(Gamma.model, function(res){
    conf.int <- confint(res)
    conf.int <- conf.int[4:nrow(conf.int),]
    return(conf.int)
  })
  
  gamma.conf <- lapply(seq_along(Gamma.conf), function(i){
    mat <- Gamma.conf[[i]]
    mat <- mat %>%
      as.data.frame() %>% 
      mutate(Gamma.est[-1,i]) %>%
      mutate(Topic = paste0("topic_", i ))
    colnames(mat) <- c("lowBound", "HighBound", "Estimate", "Group")
    return(mat)
  })
  
  names(gamma.conf) <- paste0("topic_", 1:length(gamma.conf))
  gamma.conf <- gamma.conf[match(likelihood$Topic, names(gamma.conf))]
  names(gamma.conf) <- paste0("Group", 1:length(gamma.conf))
  
  trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  true_gamma <- as.data.frame(scSTMobj$settings$sce@metadata$TrueParams$gamma)
  colnames(true_gamma) <- gsub("K", "Group", colnames(true_gamma))
  
  Type = sapply(strsplit(set_level, "_"), `[`, 2)
  
  if(Type == "Null"){
    sim_df <- lapply(names(gamma.conf), function(name){
      mat <- gamma.conf[[name]] %>% as.data.frame() 
      mat <- mat %>% mutate(True_Gamma = 0)
      return(mat)
    })
  }else{
    sim_df <- lapply(names(gamma.conf), function(name){
      mat <- gamma.conf[[name]] %>% as.data.frame() 
      if(nrow(mat) != 0 & name %in% colnames(true_gamma)) {mat <- mat %>% mutate(True_Gamma = true_gamma[,name])}
      return(mat)
    })
  }
  names(sim_df) <- paste0("Group", 1:length(sim_df))
  
  res[[set_level]] <- sim_df
  rm(scSTMobj)
  rm(Gamma.conf)
  rm(gamma.conf)
  rm(sim_df)
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/NormalGamma_MultiResponse_LinearMixed_nSample6_nCellType5_Batch_StromalCell_Gamma.rds")

# # # ############################# analysis #########################################
# res <- readRDS("res/composition_change/NormalGamma_MultiResponse_LinearMixed_nSample6_nCellType5_noBatch_StromalCell_Gamma.rds")
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# # res_null <- res[!grepl("Null", names(res))]
# 
# dat <- lapply(res, function(mat) {
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
#   mat <- lapply(mat, function(df) {
#     df <- df %>%
#       as.data.frame() %>%
#       mutate(exist_change = ifelse(lowBound < 0 & HighBound > 0, FALSE, TRUE)) %>%
#       rownames_to_column(var = "cov")
#     return(df)
#   })
# 
#   mat <- do.call(rbind, mat)
#   
#   unique_covariates <- unique(mat$cov)
#   res_cov <- data.frame()
#   for (covariate in unique_covariates) {
#     cov_data <- mat[mat$cov == covariate, ]
#     exist_change <- ifelse(all(cov_data$exist_change == FALSE), FALSE, TRUE)
#     temp <- c(covariate, exist_change)
#     res_cov <- rbind(res_cov, temp)
#   }
#   
#   colnames(res_cov) <- c("covariate", "exist_change") 
#   return(res_cov)
# })
# 
# dat <- lapply(names(dat), function(name){
#   mat <- dat[[name]]
#   mat <- mat %>% as.data.frame() %>% mutate(Type = sapply(strsplit(name, "_"), `[`, 2))
#   return(mat)
# })
# dat <-do.call(rbind, dat)
# 
# ### Type I error ####
# dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & covariate == "TimeTime2") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & covariate == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & covariate == "ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# 
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & covariate == "TimeTime2") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & covariate == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & covariate == "ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# 
# # ggplot(proportion_data, aes(x = Condition, y = Proportion, color = Type)) +
# #   geom_point(size = 4) +
# #   facet_wrap(~Type, scales = "free_x") + # Separate plots for Type I error and Power
# #   theme_minimal() +
# #   labs(title = "Type I error and Power For Sample Level Testing",
# #        x = "Condition", y = "Proportion") +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# #   scale_color_manual(values = c("Type I error" = "red",
# #                                 "Power" = "blue"))+ # Assign different colors to dots
# #   scale_y_continuous(breaks = seq(0, 1, by = 0.1), # Define y-axis ticks at 0.1 intervals
# #                      limits = c(0, 1))
# ### if we look at each topic individually #####
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
# 
#   mat <- lapply(mat, function(df) {
#     df <- df %>%
#       as.data.frame() %>%
#       mutate(exist_change = ifelse(lowBound < 0 & HighBound > 0, FALSE, TRUE)) %>%
#       rownames_to_column(var = "cov")
#     return(df)
#   })
#   
#   mat <- do.call(rbind, mat) %>%
#     as.data.frame() %>%
#     mutate(Type = strsplit(name, "_")[[1]][2])
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# 
# ### Type I error ####
# dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & cov == "TimeTime2") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & cov == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "NullModel" & cov == "ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# 
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "ResponseResponse") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
# 
# 
# #### Check if the True Gamma is Within 95% Interval #####
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
#   
#   mat <- lapply(mat, function(df) {
#     df <- df %>%
#       as.data.frame() %>%
#       mutate(Recovered = ifelse(lowBound < True_Gamma & HighBound > True_Gamma, TRUE, FALSE)) %>%
#       rownames_to_column(var = "cov")
#     return(df)
#   })
#   
#   mat <- do.call(rbind, mat) %>%
#     as.data.frame() %>%
#     mutate(Type = strsplit(name, "_")[[1]][2])
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(Recovered == TRUE, na.rm = T))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "ResponseResponse") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
# 
# #### Check if the True Gamma is in the same direction of the estimation #####
# # res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# dat <- res[grep("HighVar", names(res))]
# 
# dat <- lapply(names(res), function(name) {
#   mat <- res[[name]]
#   mat <- mat[sapply(mat, function(df) "True_Gamma" %in% colnames(df))]
#   mat <- lapply(mat, function(df) {
#     df <- df %>%
#       as.data.frame() %>%
#       mutate(Type = strsplit(name, "_")[[1]][2]) %>%
#       mutate(Recovered = ifelse(sign(Estimate) == sign(True_Gamma), TRUE, FALSE))%>%
#       rownames_to_column(var = "cov")
#     return(df)
#   })
#   mat <- do.call(rbind, mat)
#   return(mat)
# })
# 
# dat <-do.call(rbind, dat)
# ### Power ####
# dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(Recovered == TRUE, na.rm = T))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "TimeTime2:ResponseResponse") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
# dat %>% dplyr::filter(Type == "HighVar" & cov == "ResponseResponse") %>% dplyr::summarise(proportion = mean(Recovered == TRUE))
