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
  scSTM_name <- paste0(dir, "scSTM_LinearMixed_noSample_noContent_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }

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


  trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  true_gamma <- as.data.frame(scSTMobj$settings$sce@metadata$TrueParams$gamma)
  # true_psi <- scSTMobj$settings$sce@metadata$TrueParams$psi
  # true_theta <- sce@metadata$TrueParams$theta
  true_gamma_long <- true_gamma %>%
    rownames_to_column(var = "Cov") %>%
    pivot_longer(cols = -Cov, names_to = "K", values_to = "True_Gamma")
  true_gamma_long$K <- as.numeric(gsub("K", "", true_gamma_long$K))
  
  
  lm_res <- scSTMobj$mu$lm.model
  lm_res <- lapply(lm_res, function(x) summary(x))
  tvalue_list <- lapply(lm_res, function(x) x$coefficients["X", "t value"])
  tvalue_list <- do.call(c, tvalue_list)
  pvalue_list <- lapply(lm_res, function(x) x$coefficients["X", "Pr(>|t|)"])
  pvalue_list <- do.call(c, pvalue_list)
  sim_df <- data.frame(K = 1:ncol(scSTMobj$mu$gamma),
                       Mean = scSTMobj$mu$gamma[2,],
                       SD = scSTMobj$mu$std.gamma[2,],
                       tValue = tvalue_list,
                       pValue = pvalue_list) %>%
    full_join(true_gamma_long, by = "K")
  
  
  res[[set_level]] <- sim_df
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/NormalGamma_SingleResponse_LinearMixed_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")

# ############################# analysis #########################################
# # # res <- readRDS("res/composition_change/NormalGamma_ZeroPsi_nSample3_nCellType5_noBatch_StromalCell.rds")
res <- readRDS("res/composition_change/NormalGamma_SingleResponse_LinearRegression_Null_nSample3_nCellType5_noBatch_StromalCell_pValue.rds")
# res_null <- res[!grepl("Null", names(res))]

dat <- lapply(res, function(mat) {
  mat$adj.pValue <- p.adjust(mat$pValue, method = "BH")
  mat$Significant <- ifelse(mat$adj.pValue < 0.05,
                        TRUE,
                        FALSE)
  if(sum(mat$Significant) == 0){exist_change = FALSE} else{exist_change = TRUE}
  return(exist_change)
})

dat <-do.call(rbind, dat)
dat <- dat %>% as.data.frame() %>% mutate(Type = sapply(strsplit(rownames(dat), "_"), `[`, 2))
colnames(dat) <- c("exist_change", "Type")

### Type I error ####
dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
### Power ####
dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))


#### if we look at each topic individually #####
res <- readRDS("res/composition_change/NormalGamma_SinglePatient_LinearRegression_nSample1_nCellType5_noBatch_StromalCell_pValue.rds")
dat <- lapply(names(res), function(name) {
  mat <- res[[name]]
  mat$adj.pValue <- p.adjust(mat$pValue, method = "BH")
  mat$Significant <- ifelse(mat$adj.pValue < 0.05, TRUE, FALSE)

  # Determine if there's any significant change
  exist_change <- ifelse(sum(mat$Significant) == 0, FALSE, TRUE)

  # Add Type from the provided name
  mat$Type <- strsplit(name, "_")[[1]][2]
  return(mat)
})

dat <-do.call(rbind, dat)

### Type I error ####
dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(Significant == TRUE))
### Power ####
dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(Significant == TRUE))
