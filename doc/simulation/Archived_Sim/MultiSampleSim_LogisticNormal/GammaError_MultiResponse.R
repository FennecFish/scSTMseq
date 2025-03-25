# Goal: This script is trying to run GammaError for all scSTMseq model
# specifically those on the 
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

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_MultiResponse/nSample6_nCellType5_noBatch_StromalCell/"

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

scSTM_name <- paste0(dir, "scSTM_LinearMixed_noContent_Prevalence_TimeandResponse/scSTM_",
                     set_level,".rds")
if(file.exists(scSTM_name)){
  scSTMobj <- readRDS(scSTM_name)
  scSTMobj <- select_top_scSTM(scSTMobj)
}else{
  next
}
formula = ~Time*Response + (1|Sample) # mixed linear model formula
null.formula  = ~Time + Response + (1|Sample)
res <- GammaError(model = scSTMobj, formula = formula, null_formula = null.formula, nsims = 100)
saveRDS(res, file = paste0(dir, "GammaError_MixedModel/GammaError_",set_level,".rds"))

rm(res)
formula = ~Time*Response + (1|Sample)
null.formula  = ~Time + Response + (1|Sample)
res <- GammaError(model = scSTMobj, formula = formula, null_formula = null.formula, nsims = 100)
saveRDS(res, file = paste0(dir, "GammaError_MixedModel/GammaError_",set_level,".rds"))


# 
# ##### analysis ####'
# dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
# files <- list.files(paste0(dir, "GammaError_MixedModel/"))
# res <- vector("list")
# for(file in files){
#   set_level <- sub("^GammaError_(.*)\\.rds$", "\\1", file)
#   model <- readRDS(paste0(dir, "GammaError_MixedModel/", file))
#   simplex <- dim(model$SimEta[[1]])[2]
#   pValue_output_lists <- vector("list", simplex)
#   Chisq_output_lists <- vector("list", simplex)
#   # Loop over the K lists and extract the respective nested list from 'model'
#   for (i in 1:simplex) {
#     pValue.list <- lapply(model$model.diff, function(x) {
#       k.model <- x[[i]]
#       pValue <- k.model$`Pr(>Chisq)`[2]
#       pValue
#     })
#     pValue.list <- do.call(c, pValue.list)
#     
#     Chisq.list <- lapply(model$model.diff, function(x) {
#       k.model <- x[[i]]
#       Chisq <- k.model$Chisq[2]
#       Chisq
#     })
#     Chisq.list <- do.call(c, Chisq.list)
#     
#     pValue_output_lists[[i]] <- pValue.list
#     Chisq_output_lists[[i]] <- Chisq.list
#   }
#   pValue <- do.call(cbind, pValue_output_lists) %>% as.data.frame()
#   colnames(pValue) <- paste0("topic_", 1:ncol(pValue))
#   rownames(pValue) <- paste0("Replicate_", 1:nrow(pValue))
#   
#   Chisq <- do.call(cbind, Chisq_output_lists) %>% as.data.frame()
#   colnames(Chisq) <- paste0("topic_", 1:ncol(Chisq))
#   rownames(Chisq) <- paste0("Replicate_", 1:nrow(Chisq))
#   
#   res[[set_level]] <- list(pValue = pValue, Chisq = Chisq)
# }
# name <- basename(dir)
# saveRDS(res, file = paste0("res/GammaError_MixedLinearModelComparison_",name,".rds"))
# 
# # # #### plot ####
# res <- readRDS("res/GammaError_MixedLinearModelComparison_nSample3_nCellType5_noBatch_StromalCell.rds")
# null.res <- res[grep("NullModel",names(res))]
# alt.res <- res[grep("HighVar",names(res))]
# 
# ### approach 1, take the mean of chisq among replicates
# null.chisq <- lapply(null.res, function(x) {
#   mat <- x$Chisq
#   mean.chisq <- colMeans(mat)
#   mean.chisq
# })
# null.chisq <- do.call(rbind, null.chisq)
# null.chisq.pvalues <- apply(null.chisq, c(1, 2), function(x) pchisq(x, df = 1, lower.tail = FALSE))
# # adjust pvalue
# null.exist_change <- apply(null.chisq.pvalues, 1, function(row) {
#   p.adj <- p.adjust(row, method = "fdr")
#   exist_change <- ifelse(any(p.adj<0.05), TRUE, FALSE)
#   exist_change
#   })
# 
# ##### approach 2 Pval50Perc
# null.res.pValue <- lapply(null.res, function(x){
#   mat <- x$pValue
#   browser()
#   Pval50Perc = quantile(x, probs = 0.50)
# })
# 
# 
# null <- do.call(rbind, null.res)
# alt <- do.call(rbind, alt.res)
