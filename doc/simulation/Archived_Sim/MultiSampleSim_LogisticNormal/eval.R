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

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

################################################################################################################

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_MultiResponse/nSample6_nCellType5_noBatch_StromalCell/"
files <- list.files(path = paste0(dir, "sims/"), pattern = "sims*")

res <- data.frame()

for(file_name in files){
  # create a new dataframe to calculate adjRandIndex
  dat <- data.frame()
  
  set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
  # cat(set_level, "\n")
  
  file_path <- paste0(dir, "sims/sims_", set_level, ".rds")
  sims <- safe_readRDS(file_path)
  if (is.null(sims)) {
    next  # Skip to the next iteration if reading the file failed
  }
  dat <- colData(sims) %>% data.frame() 

  # ##### scSTM_Pooled
  # scSTM_pooled_name <- paste0(dir, "scSTM_Pooled_Sample_noContent_Prevalence_TimeandResponse/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_pooled_name)){
  #   scSTM_pooled <- readRDS(scSTM_pooled_name)
  #   scSTM_pooled_cluster <- process_scSTM(scSTM_pooled)
  #   dat$scSTM_pooled_cluster <- scSTM_pooled_cluster[match(rownames(dat), names(scSTM_pooled_cluster))]
  # } else {dat$scSTM_pooled_cluster <- NA}

  #### scSTM_LinearMixed
  scSTM_Mixed_name <- paste0(dir, "scSTM_LinearMixed_noSample_noContent_Prevalence_TimeandResponse/scSTM_",set_level,".rds")
  if(file.exists(scSTM_Mixed_name)){
    scSTM_Mixed <- readRDS(scSTM_Mixed_name)
    scSTM_Mixed_cluster <- process_scSTM(scSTM_Mixed)
    dat$scSTM_Mixed_cluster <- scSTM_Mixed_cluster[match(rownames(dat), names(scSTM_Mixed_cluster))]
  } else {dat$scSTM_Mixed_cluster <- NA}
  # 
  # #### scSTM_LinearRegression
  # scSTM_LR_name <- paste0(dir, "scSTM_LinearRegression_noContent_Prevalence_TimeandResponse/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_LR_name)){
  #   scSTM_LR <- readRDS(scSTM_LR_name)
  #   scSTM_LR_cluster <- process_scSTM(scSTM_LR)
  #   dat$scSTM_LR_cluster <- scSTM_LR_cluster[match(rownames(dat), names(scSTM_LR_cluster))]
  # } else {dat$scSTM_LR_cluster <- NA}


  # ###### fastTopics
  # fasttopic_name <- paste0(dir, "fastTopics/fastTopics_", set_level, ".rds")
  # if(file.exists(fasttopic_name)){
  #   nmf.sims <- readRDS(fasttopic_name)
  #   max_indices <- apply(nmf.sims$L, 1, which.max)
  #   fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
  #   names(fastTopics_cluster) <- rownames(nmf.sims$L)
  #   dat$fastTopics_cluster <- fastTopics_cluster[match(rownames(dat), names(fastTopics_cluster))]
  # } else{ dat$fastTopics_cluster <- NA }
  # 
  # # Seurat
  # seurat_name <- paste0(dir, "seurat/seurat_", set_level, ".rds")
  # if(file.exists(seurat_name)){
  #   seurat.sims <- readRDS(seurat_name)
  #   smeta <- seurat.sims@meta.data %>% as.data.frame()
  #   sub_sims <- sims[,rownames(smeta)] # filter by the rows
  #   seurat.adj <- sapply(smeta[,4:7], function(x) {
  #     adjustedRandIndex(x, sub_sims$Group)
  #   })
  #   # select the resolution that has the highest ARI. When there are multiple, select the first one
  #   best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
  #   seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
  #   dat$seurat_cluster <- seurat_cluster[match(rownames(dat), rownames(seurat_cluster)),]
  #   rm(sub_sims)
  # } else{dat$seurat_cluster <- NA}
  # 
  # 
  # # sctransform
  # sctf_name <- paste0(dir, "sctransform/sctransform_", set_level, ".rds")
  # if(file.exists(sctf_name)){
  #   sctf <- readRDS(sctf_name)
  #   sct_meta <- sctf@meta.data %>% as.data.frame()
  #   sub_sims <- sims[,rownames(sct_meta)] # filter by the rows
  #   sctf.adj <- sapply(sct_meta[,6:9], function(x) {
  #     adjustedRandIndex(x, sub_sims$Group)
  #   })
  #   # select the resolution that has the highest ARI. When there are multiple, select the first one
  #   best_res <- names(sctf.adj)[sctf.adj == max(sctf.adj)][1]
  #   sctf_cluster <- sctf@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
  #   dat$sctransform_cluster <- sctf_cluster[match(rownames(dat), rownames(sctf_cluster)),]
  #   rm(sub_sims)
  # }else{dat$sctransform_cluster <- NA}
  # 
  # # monocle3
  # monocle_name <- paste0(dir, "monocle3/monocle3_", set_level, ".rds")
  # if(file.exists(monocle_name)){
  #   monocle3 <- readRDS(monocle_name)
  #   dat$monocle3_cluster <- partitions(monocle3)[match(rownames(dat), names(partitions(monocle3)))]
  # } else {dat$monocle3_cluster <- NA}
  # 
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

sim_setup <- basename(dir)
write.csv(res, file = paste0("res/adjRandIndex_", "NormalGamma_MultiResponse", sim_setup, ".csv"))

