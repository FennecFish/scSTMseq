setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
library(MANOVA.RM)
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType5_noBatch_StromalCell/Manualsims"
# dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType10_noBatch_StromalCell/sims"
files <- list.files(path = dir, pattern = "Null")
pValue <- vector(mode = "list")
for (file_name in files){
  #file_name <- files[i]
  sims <- readRDS(paste0(dir, "/", file_name))
  
  theta.collapsed <- colData(sims) %>%
    as.data.frame() %>%
    group_by(Time, Sample, Group) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(Time, Sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    dplyr::select(-count) %>%
    pivot_wider(names_from = Group, values_from = proportion, values_fill = 0)
  
  
  ##### Fit Manova.RM #######
  response_vars <- grep("^Group", colnames(theta.collapsed), value = TRUE)
  fit <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
                data = theta.collapsed, subject = "Sample", within = "Time", iter = 1000)
  pValue[[file_name]] <- summary(fit)[2] # extract MATS pValue
  rm(sims)
  cat(file_name, "\n")
}
saveRDS(pValue, "res/simsManual_Manova_pValue_Null_nSample20_nCellType5_noBatch_StromalCell.rds")

# 
# res <- readRDS("res/simsManual_Manova_pValue_Null_nSample20_nCellType5_noBatch_StromalCell.rds")
# res <- as.numeric(do.call(rbind, res))
# hist(res)
