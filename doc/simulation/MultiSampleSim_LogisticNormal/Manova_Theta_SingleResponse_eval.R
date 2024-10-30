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
library(MANOVA.RM)

# args <- commandArgs(trailingOnly = TRUE)
# file_name <- args[1]
# file_name <- basename(file_name)
# cat(file_name, "\n")
# 
# set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)
# dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample20_nCellType5_noBatch_CancerCell/"
# 
# r.file <- paste0("R/",list.files("R/"))
# sapply(r.file, source)
# 
# scSTM_name <- paste0(dir, "scSTM_LinearMixed_noContent_Prevalence_Time/scSTM_",
#                      set_level,".rds")
# if(file.exists(scSTM_name)){
#   scSTMobj <- readRDS(scSTM_name)
#   scSTMobj <- select_top_scSTM(scSTMobj)
# }else{
#   next
# }
# 
# ThetaManova <- ThetaManova(model = scSTMobj, nsims = 100)
# saveRDS(ThetaManova, file = paste0(dir, "ThetaManova/ThetaManova_",set_level,".rds"))


##### analysis ####'
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample20_nCellType5_noBatch_CancerCell/"
files <- list.files(paste0(dir, "ThetaManova/"))
# compDTU.mixed <- vector(mode = "list")
# compDTU <- vector(mode = "list")
res <- vector(mode = "list")
# manova.dm <- vector(mode = "list")
for(file in files){
  set_level <- sub("^ThetaManova_(.*)\\.rds$", "\\1", file)
  ThetaManova <- readRDS(paste0(dir, "ThetaManova/ThetaManova_", set_level, ".rds"))
  # compDTU.mixed[[set_level]] <- ThetaManova$compDTU.mixed
  # compDTU[[set_level]] <- ThetaManova$compDTU
  rm.wts <- lapply(ThetaManova$manova.fit, function(x){
    # return(list(statistics = x$MATS, p.value = x$resampling[2]))
    dat <- x$resampling[1]
    return(dat)
  })
  
  rm.mats <- lapply(ThetaManova$manova.fit, function(x){
    dat <- x$resampling[2]
    return(dat)
  })
  mcmc <- lapply(ThetaManova$mcmcglmm.fit, function(x){
    # return(list(statistics = x$MATS, p.value = x$resampling[2]))
    res <- MCMCglmm::summary.MCMCglmm(x)
    dat <- res$solutions[2,"pMCMC"]
    return(dat)
  })
  
  dat.clean <- function(given_list){
    given_list <- do.call(rbind, given_list)
    rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
    colnames(given_list) <- "pValue"
    return(given_list)
  }
  res[[set_level]] <- list(rm.wts = dat.clean(rm.wts),
                           rm.mats = dat.clean(rm.mats),
                           mcmc.glmm = dat.clean(mcmc))
  cat(file, "\n")
}
name <- basename(dir)
# saveRDS(compDTU.mixed, file = paste0("res/compDTUmixed_MixedLinearModelComparison_",name,".rds"))
# saveRDS(compDTU, file = paste0("res/compDTU_MixedLinearModelComparison_",name,".rds"))
# saveRDS(manova.dm, file = paste0("res/manova.dm_MixedLinearModelComparison_",name,".rds"))
# saveRDS(mcmcglmm, file = paste0("res/mcmcGLMM_MixedLinearModelComparison_",name,".rds"))
saveRDS(res, file = paste0("res/MixedManovaComparison_",name,".rds"))

#################################################################################
################################## plot #########################################
#################################################################################
num_elements <- function(obj) {
  if (is.data.frame(obj)) {
    return(nrow(obj) * ncol(obj))  # For data frames, return the total number of elements (rows * columns)
  } else {
    return(length(obj))  # For vectors, return the length (number of elements)
  }
}
plot_data_generate <- function(thresholds, pValueList){

  exist.change.list <- lapply(thresholds, function(thr) {
    lapply(pValueList, function(method) {
      ifelse(method[,1] < thr, TRUE, FALSE)
    })
  })
  TypeIError_data <- data.frame()

  for (i in seq_along(threshold)) {
    threshold <- thresholds[i]

    for (method in names(exist.change.list[[i]])) {
      positive <- sum(exist.change.list[[i]][[method]], na.rm = T)
      total <- num_elements(exist.change.list[[i]][[method]])

      TypeIError <- positive/total

      # Store the results in a dataframe
      TypeIError_data <- rbind(TypeIError_data, data.frame(
        Method = method,
        Threshold = threshold,
        TypeIError = TypeIError
      ))
    }
  }
  return(TypeIError_data)
}

# compDTU.mixed <- readRDS("res/compDTUmixed_MixedLinearModelComparison_nSample20_nCellType5_noBatch_CancerCell.rds")
# compDTU <- readRDS("res/compDTU_MixedLinearModelComparison_nSample20_nCellType5_noBatch_CancerCell.rds")
# manova.dm <- readRDS("res/manova.dm_MixedLinearModelComparison_nSample20_nCellType5_noBatch_CancerCell.rds")
# mcmc.dm <- readRDS("res/mcmcGLMM_MixedLinearModelComparison_nSample20_nCellType5_noBatch_CancerCell.rds")
dat <- readRDS("res/MixedManovaComparison_nSample20_nCellType5_noBatch_CancerCell.rds")
null.dat <- dat[grep("NullModel",names(dat))]
alt.dat <- dat[grep("HighVar",names(dat))]
# compDTU.mixed.null <- compDTU.mixed[grep("NullModel",names(compDTU.mixed))]
# compDTU.mixed.null <- do.call(rbind, compDTU.mixed.null)
# colnames(compDTU.mixed.null)[1] <- "p.value"
# compDTU.null <- compDTU[grep("NullModel",names(compDTU))]
# compDTU.null <- do.call(rbind, compDTU.null)
# colnames(compDTU.null)[1] <- "p.value"
manova.wts.pValue <- lapply(null.dat, function(x) {
  mat <- x$rm.wts
  # mean.chisq <- mean(mat)
  # pchisq(mean.chisq, df = mean(x$df), lower.tail = FALSE)
  quantile(mat, probs = 0.50)
})
manova.wts.pValue <- do.call(rbind, manova.wts.pValue)
colnames(manova.wts.pValue) <- "p.value"

manova.mats.pValue <- lapply(null.dat, function(x) {
  mat <- x$rm.mats
  # mean.chisq <- mean(mat)
  # pchisq(mean.chisq, df = mean(x$df), lower.tail = FALSE)
  quantile(mat, probs = 0.50)
})
manova.mats.pValue <- do.call(rbind, manova.mats.pValue)
colnames(manova.mats.pValue) <- "p.value"

mcmc.glmm.pValue <- lapply(null.dat, function(x) {
  mat <- x$mcmc.glmm
  quantile(mat, probs = 0.50)
})
mcmc.glmm.pValue <- do.call(rbind, mcmc.glmm.pValue)
colnames(mcmc.glmm.pValue) <- "p.value"

# res.pValue <- list(compDTU.reduced = compDTU.mixed.null, compDTU = compDTU.null, manova.dm = manova.dm.null)
res.pValue <- list(manova.rm.wts = manova.wts.pValue, 
                   manova.rm.mats = manova.mats.pValue,
                   mcmc.glmm = mcmc.glmm.pValue)
threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = res.pValue)

ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  # geom_point(size = 3, position = position_jitter(height = 1e-10)) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Type I Error by Methods and Thresholds",
       x = "Threshold",
       y = "Type I Error") +
  theme_minimal() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14)   # Increase Y axis labels size
  )

################## power #############################
# compDTU.mixed.alt <- compDTU.mixed[grep("HighVar",names(compDTU.mixed))]
# compDTU.mixed.alt <- do.call(rbind, compDTU.mixed.alt)
# colnames(compDTU.mixed.alt)[1] <- "p.value"
# compDTU.alt <- compDTU[grep("HighVar",names(compDTU))]
# compDTU.alt <- do.call(rbind, compDTU.alt)
# colnames(compDTU.alt)[1] <- "p.value"
# manova.dm.alt <- manova.dm[grep("HighVar",names(manova.dm))]
# manova.dm.alt <- do.call(rbind, manova.dm.alt)
# colnames(manova.dm.alt) <- "p.value"

manova.wts.pValue <- lapply(alt.dat, function(x) {
  mat <- x$rm.wts
  # mean.chisq <- mean(mat)
  # pchisq(mean.chisq, df = mean(x$df), lower.tail = FALSE)
  quantile(mat, probs = 0.50)
})
manova.wts.pValue <- do.call(rbind, manova.wts.pValue)
colnames(manova.wts.pValue) <- "p.value"

manova.mats.pValue <- lapply(alt.dat, function(x) {
  mat <- x$rm.mats
  # mean.chisq <- mean(mat)
  # pchisq(mean.chisq, df = mean(x$df), lower.tail = FALSE)
  quantile(mat, probs = 0.50)
})
manova.mats.pValue <- do.call(rbind, manova.mats.pValue)
colnames(manova.mats.pValue) <- "p.value"

mcmc.glmm.pValue <- lapply(alt.dat, function(x) {
  mat <- x$mcmc.glmm
  quantile(mat, probs = 0.50)
})
mcmc.glmm.pValue <- do.call(rbind, mcmc.glmm.pValue)
colnames(mcmc.glmm.pValue) <- "p.value"

# res.pValue <- list(compDTU.reduced = compDTU.mixed.null, compDTU = compDTU.null, manova.dm = manova.dm.null)
res.pValue <- list(manova.rm.wts = manova.wts.pValue, 
                   manova.rm.mats = manova.mats.pValue,
                   mcmc.glmm = mcmc.glmm.pValue)
threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = res.pValue)

colnames(topic_plot)[3] <- "Power"
ggplot(topic_plot, aes(x = Threshold, y = Power, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  geom_point(size = 3, position = position_jitter(height = 1e-10)) +
  # geom_point(size = 3) +
  labs(title = "Power by Methods and Thresholds",
       x = "Threshold",
       y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14)   # Increase Y axis labels size
  )
