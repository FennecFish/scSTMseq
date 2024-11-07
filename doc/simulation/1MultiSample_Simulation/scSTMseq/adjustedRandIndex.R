setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
design <- "nSample20_nCellType5_noBatch_StromalCell/"
pooled <- "1000scSTM_Pooled_noContent_Prevalence_Time/"
# lm <-"scSTM_LinearMixed_noContent_Prevalence_Time/"
path = paste0(dir, design, pooled)
files <- list.files(path)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

res <- data.frame()

for(file_name in files){
  
  set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file_name)
  
  # first read in pooled data
  file_path <- paste0(path, file_name)
  pooled_scSTMobj <- safe_readRDS(file_path)
  pooled_scSTMobj <- select_top_scSTM(pooled_scSTMobj)
  pooled_cluster <- cluster_scSTM(pooled_scSTMobj)
  
  dat <- colData(pooled_scSTMobj$settings$sce) %>% data.frame() 
  dat$pooled_scSTMseq <- pooled_cluster[match(rownames(dat), names(pooled_cluster))]
  
  # then read in linear mixed model
  # lm_scSTMobj_path <- paste0(dir, design, lm, "scSTM_", set_level, ".rds")
  # if(file.exists(lm_scSTMobj_path)){
  #   lm_scSTMobj <- readRDS(lm_scSTMobj_path)
  #   lm_cluster <- process_scSTM(lm_scSTMobj)
  #   dat$lm_scSTMseq <- lm_cluster[match(rownames(dat), names(lm_cluster))]
  # } else {dat$lm_scSTMseq <- NA}
  
  adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("scSTMseq")), function(x) {
    adjustedRandIndex(x, dat$Group)
  })

  # Create a data frame to store results
  res.temp <- data.frame(
    modelType = unlist(strsplit(set_level, "_"))[2],
    seed = unlist(strsplit(set_level, "_"))[1],
    nCellType = as.numeric(gsub("nCellType", "", unlist(strsplit(design, "_"))[2])),
    nSample = as.numeric(gsub("nSample", "", unlist(strsplit(design, "_"))[1])),
    Batch = ifelse(unlist(strsplit(design, "_"))[3]=="noBatch", FALSE, TRUE),
    CancerType = ifelse(unlist(strsplit(design, "_"))[4]=="StromalCell/", FALSE, TRUE),
    t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  )
  
  res <- bind_rows(res, res.temp)
  cat(file_name, "\n")
  rm(pooled_scSTMobj)
  rm(dat)
}

# pre_res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final.csv")
# res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final_update.csv")
# 
# res.update <- rbind(pre_res, res)
# res.update <- res %>%
#  # dplyr::select(sim, -matches("cluster$"), matches("cluster$")) %>%
#   dplyr::full_join(pre_res, by = c("sim", "seed", "control", "level", "nCellType", "nSample")) %>%
#  dplyr::select(sim, seed, control, level, nCellType, nSample,
#                -matches("cluster$"), matches("cluster$"))

save_path <- paste0("res/adjRandIndex_",basename(dir), "_", basename(design), "_", basename(pooled),".csv")
write.csv(res, file = save_path)


