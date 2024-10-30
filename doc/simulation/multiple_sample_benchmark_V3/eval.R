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
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/", pattern = "sims*")
files <- list.files(
  path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  pattern = "^sims.*(nsample6\\.rds|nsample12\\.rds)$")
# files <- list.files(
#   path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
#   pattern = "^sims.*(nsample12\\.rds)$")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_Time_Prevalence_Time/")

res <- data.frame()

for(file_name in files){
  # create a new dataframe to calculate adjRandIndex
  dat <- data.frame()
  
  # set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  # setup <- sub("sims_(?:multiCellType_)?([^.]*)\\.rds", "\\1", file_name)
  file_path <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/sims_", set_level, ".rds")
  sims <- safe_readRDS(file_path)
  if (is.null(sims)) {
    next  # Skip to the next iteration if reading the file failed
  }
  dat <- colData(sims) %>% data.frame() 
  
  # scSTM_C_Time_P # without interaction
  # scSTM_C_TimeNI_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_TimeNoInteraction_Prevalence_Time/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_C_TimeNI_P_name)){
  #   scSTM_C_TimeNI_P <- readRDS(scSTM_C_TimeNI_P_name)
  #   scSTM_C_TimeNI_P_cluster <- process_scSTM(scSTM_C_TimeNI_P)
  #   dat$scSTM_C_TimeNI_P_cluster <- scSTM_C_TimeNI_P_cluster[match(rownames(dat), names(scSTM_C_TimeNI_P_cluster))]
  # } else {dat$scSTM_C_TimeNI_P_cluster <- NA}

  # 
  # # scSTM_C_Sample_P
  # scSTM_C_Sample_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_Sample_Prevalence_Time/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_C_Sample_P_name)){
  #   scSTM_C_Sample_P <- readRDS(scSTM_C_Sample_P_name)
  #   scSTM_C_Sample_P_cluster <- process_scSTM(scSTM_C_Sample_P)
  #   dat$scSTM_C_Sample_P_cluster <- scSTM_C_Sample_P_cluster[match(rownames(dat), names(scSTM_C_Sample_P_cluster))]
  # } else {dat$scSTM_C_Sample_P_cluster <- NA}

  # scSTM_C_NIBatch_P
  scSTM_C_BatchNI_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_BatchNoInteraction_Prevalence_Time/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_BatchNI_P_name)){
    scSTM_C_BatchNI_P <- readRDS(scSTM_C_BatchNI_P_name)
    scSTM_C_BatchNI_P_cluster <- process_scSTM(scSTM_C_BatchNI_P)
    dat$scSTM_C_BatchNI_P_cluster <- scSTM_C_BatchNI_P_cluster[match(rownames(dat), names(scSTM_C_BatchNI_P_cluster))]
  } else {dat$scSTM_C_BatchNI_P_cluster <- NA}

  
  # scSTM_nC_P
  # scSTM_nC_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_noContent_Prevalance/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_nC_P_name)){
  #   scSTM_nC_P <- readRDS(scSTM_nC_P_name)
  #   scSTM_nC_P_cluster <- process_scSTM(scSTM_nC_P)
  #   dat$scSTM_nC_P_cluster <- scSTM_nC_P_cluster[match(rownames(dat), names(scSTM_nC_P_cluster))]
  # } else {dat$scSTM_nC_P_cluster <- NA}

  # # scSTM_C_nP
  # scSTM_C_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_noPrevalance/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_C_P_name)){
  #   scSTM_C_nP <- readRDS(scSTM_C_P_name)
  #   scSTM_C_nP_cluster <- process_scSTM(scSTM_C_nP)
  #   dat$scSTM_C_nP_cluster <- scSTM_C_nP_cluster[match(rownames(dat), names(scSTM_C_nP_cluster))]
  # } else {dat$scSTM_C_nP_cluster <- NA}
  # 
  # # scSTM_C_P
  # scSTM_C_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_Prevalance/scSTM_",set_level,".rds")
  # if(file.exists(scSTM_C_P_name)){
  #   scSTM_C_P <- readRDS(scSTM_C_P_name)
  #   # scSTM_C_P_cluster <- cluster_scSTM(scSTM_C_P$runout[[10]])
  #   # scSTM_C_P_cluster <- scSTM_C_P_cluster[match(rownames(dat), names(scSTM_C_P_cluster))]
  #   # adjustedRandIndex(scSTM_C_P_cluster, dat$Group)
  #   scSTM_C_P_cluster <- process_scSTM(scSTM_C_P)
  #   dat$scSTM_C_P_cluster <- scSTM_C_P_cluster[match(rownames(dat), names(scSTM_C_P_cluster))]
  # } else {dat$scSTM_C_P_cluster <- NA}
  # 
  # fastTopics
  fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/fastTopics/fastTopics_", set_level, ".rds")
  if(file.exists(fasttopic_name)){
    nmf.sims <- readRDS(fasttopic_name)
    max_indices <- apply(nmf.sims$L, 1, which.max)
    fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    names(fastTopics_cluster) <- rownames(nmf.sims$L)
    dat$fastTopics_cluster <- fastTopics_cluster[match(rownames(dat), names(fastTopics_cluster))]
  } else{ dat$fastTopics_cluster <- NA }

  # Seurat
  seurat_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/seurat/seurat_", set_level, ".rds")
  if(file.exists(seurat_name)){
    seurat.sims <- readRDS(seurat_name)
    smeta <- seurat.sims@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(smeta)] # filter by the rows
    seurat.adj <- sapply(smeta[,4:7], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
    seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
    dat$seurat_cluster <- seurat_cluster[match(rownames(dat), rownames(seurat_cluster)),]
    rm(sub_sims)
  } else{dat$seurat_cluster <- NA}


  # # sctransform
  # sctf_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sctransform/sctransform_", set_level, ".rds")
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

  # monocle3
  monocle_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/monocle3/monocle3_", set_level, ".rds")
  if(file.exists(monocle_name)){
    monocle3 <- readRDS(monocle_name)
    dat$monocle3_cluster <- partitions(monocle3)[match(rownames(dat), names(partitions(monocle3)))]
  } else {dat$monocle3_cluster <- NA}

  # # monocle3_time
  # monocle_time_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/monocle3_time/monocle3_", set_level, ".rds")
  # if(file.exists(monocle_time_name)){
  #   monocle3_time <- readRDS(monocle_time_name)
  #   dat$monocle3_time_cluster <- partitions(monocle3_time)[match(rownames(dat), names(partitions(monocle3_time)))]
  # } else {dat$monocle3_time_cluster <- NA}
  
  adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("cluster")), function(x) {
    adjustedRandIndex(x, sims$Group)
  })
  # adjusted_rand_indices <- adjustedRandIndex(dat[, "monocle3_time_cluster"], sims$Group)

  # Create a data frame to store results
  res.temp <- data.frame(
    sim = set_level,
    seed = unlist(strsplit(set_level, "_"))[1],
    control = unlist(strsplit(set_level, "_"))[2],
    level = unlist(strsplit(set_level, "_"))[3],
    nCellType = unlist(strsplit(set_level, "_"))[4],
    nSample = unlist(strsplit(set_level, "_"))[5],
    t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  )
  
  res <- bind_rows(res, res.temp)
  cat(file_name, "\n")
  rm(sims)
}
# 
# pre_res <- read.csv("res/sim_cluster_benchmark/adjRandIndex_multiple_sample_benchmark_V3.csv")
# res <- read.csv("res/sim_cluster_benchmark/adjRandIndex_multiple_sample_benchmark_V3_update.csv")
# res.update <- res %>%
#   dplyr::select(sim, -matches("cluster$"), matches("cluster$")) %>%
#   dplyr::right_join(pre_res, by = "sim") %>%
#   dplyr::select(sim, seed, control, level, nCellType, nSample, -matches("cluster$"), matches("cluster$"))

write.csv(res, file = "res/sim_cluster_benchmark/adjRandIndex_multiple_sample_benchmark_V3.csv")
