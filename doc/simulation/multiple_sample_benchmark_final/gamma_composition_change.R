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
# functions to read simulated data
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

select_top_scSTM <- function(scSTMobj) {
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
  return(scSTMobj)
}

cluster_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

true_sim_RR <- list(c5 = c(-0.15, 0.1, 0.2, -0.15, 0)/(1/5),
                c9 = c(0.1, -0.05, -0.08, 0.03, -0.05, 0.02, 0.02, 0.01, 0)/(1/9),
                c13 = c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
                        -0.04, 0)/(1/13),
                c17 = c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
                        -0.04, 0.01, 0.02, -0.03, 0, 0)/(1/17))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

# first read in scSTM data
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/", 
                    pattern = "^sims.*(nsample6\\.rds|nsample12\\.rds)$")
# files <- files[80:110]
res <- vector(mode = "list")
nsims <- 100
for(file_name in files){
  # read in sims data
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  
  # read in scSTMseq data
  scSTM_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/scSTM_Content_Batch_Prevalence_Time/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }
  K <- scSTMobj$settings$dim$K
  relative_change <- gammaPosterior(model = scSTMobj, nsims = nsims)
  est_cluster <- cluster_scSTM(scSTMobj)
  
  # calculate true proportion change
  sims <- scSTMobj$settings$sce

  # The following code is to map topic back to original cell types
  max_topic_per_cell <- cluster_scSTM(scSTMobj) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Cell") %>%
    dplyr::rename(Topic = ".") 
  matched_data <- colData(sims) %>%
    as.data.frame %>%
    dplyr::right_join(max_topic_per_cell, by = "Cell") 
  topic_group_counts <- matched_data %>%
    group_by(Group, Topic) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Step 2: Calculate the total number of cells in each group
  group_totals <- matched_data %>%
    group_by(Group) %>%
    summarise(total = n(), .groups = 'drop')
  
  # Step 3: Calculate the likelihood of each topic belong to the given group
  likelihoods <- topic_group_counts %>%
    dplyr::left_join(group_totals, by = "Group") %>%
    dplyr::mutate(likelihood = count / total) %>%
    dplyr::select(Group, Topic, likelihood) %>%
    group_by(Group) %>%
    slice_max(likelihood, with_ties = FALSE) %>%
    ungroup()
  
  
  proportion_df <- colData(sims) %>%
    as.data.frame() %>%
    # group_by(Sample, Time, Group) %>%
    group_by(Time, Group) %>%
    summarise(Count = n()) %>%
    # group_by(Sample, Time) %>%
    group_by(Time) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  proportion_change_df <- proportion_df %>%
    # select(Sample, Time, Group, Proportion) %>%
    dplyr::select(Time, Group, Proportion) %>%
    tidyr::spread(key = Time, value = Proportion) %>%
    dplyr::mutate(True_RR = (Time2 - Time1)/Time1) %>%
    dplyr::right_join(likelihoods, by = "Group")
  
  # this is the RR from simulated data
  colnames(relative_change) <- paste0("topic_", 1:K)
  relative_change <- relative_change[,match(proportion_change_df$Topic, colnames(relative_change))]
  RR_from_simDat <- proportion_change_df$True_RR# [match(colnames(relative_change), proportion_change_df$Topic)]
  
  # this is the RR calculated from the parameter
  name_split <- unlist(strsplit(set_level, "_"))
  if(name_split[2] == "pos"){
    nCellType <- name_split[4]
    RR_from_param <- data.frame(true_sim_RR = true_sim_RR[[nCellType]],
                                Group = levels(sims$Group)) %>%
      dplyr::right_join(likelihoods, by = "Group")
    RR_from_param <- RR_from_param$true_sim_RR[match(proportion_change_df$Group, RR_from_param$Group)]
  }else{
    RR_from_param <- rep(0, K)
  }
  
  relative_change <- rbind(RR_from_param, RR_from_simDat, relative_change)
  rownames(relative_change) <- c("RR_from_param", "RR_from_simDat",paste0("Gamma", 1:nsims))
  res[[set_level]] <- relative_change
  
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/gammaPosterior_multiple_patient_final_orderByGroup.rds")

############################# analysis #########################################
res <- readRDS("res/composition_change/gammaPosterior_multiple_patient_final_orderByGroup.rds")

res <- lapply(res, function(x) {
  colnames(x) <- paste0("Group", 1:ncol(x))
  return(x)
})
calc_mean_sd <- function(x, name, true_RR){
  # truth <- x[nrow(x),]
  sim <- x[grep("Gamma", rownames(x)),]
  # mean <- apply(sim, 2, mean)
  # sd <- apply(sim, 2, sd)
  name_split <- unlist(strsplit(name, "_"))
  seed <- name_split[1]
  control <- name_split[2]
  level <- name_split[3]
  nCellType <- name_split[4]
  nSample <- name_split[5]
  # RR_from_simDat <- x["RR_from_simDat",]
  RR_from_param <- x["RR_from_param",]
  # bias <- abs(sim - RR_from_simDat)
  mean <- apply(sim, 2, mean)
  sd <- apply(sim, 2, sd)
  # dat <- data.frame(name = name,
  #                   seed = seed,
  #                   control = control,
  #                   level = level,
  #                   nCellType = nCellType,
  #                   nSample = nSample,
  #                   mean = mean,
  #                   sd = sd)
  dat <- data.frame(name = name,
                    seed = seed,
                    control = control,
                    level = level,
                    nCellType = nCellType,
                    nSample = nSample,
                    mean = mean,
                    sd = sd,
                    RR_from_param = RR_from_param,
                    Group = names(mean))
  dat$Topic <- rownames(dat)
  return(dat)
}

mean_sd <- do.call(rbind, lapply(names(res), function(name) {
  calc_mean_sd(res[[name]], name)
}))

RR_dat <- mean_sd %>% dplyr::select("control", "level", "nCellType", "nSample", "Group", "RR_from_param")

aggregated_mean_sd <- mean_sd %>%
  group_by(control, nCellType, level, nSample, Group) %>%
  summarize(mean = mean(mean),
            sd = mean(sd),
            .groups = 'drop') %>%
  left_join(RR_dat, by = c("control", "level", "nCellType", "nSample", "Group"))

aggregated_mean_sd$Group <- factor(aggregated_mean_sd$Group, levels = paste0("Group", 1:17))
aggregated_mean_sd$nCellType <- factor(aggregated_mean_sd$nCellType, levels= paste0("c", c(5, 9, 13, 17)))
# mean_sd <- mean_sd %>%
#   mutate(diff_mean_RR = mean - RR_from_param)

png("res/composition_change/composition_change_mean_plot.png", res= 300, width = 6000, height = 3000)
ggplot(aggregated_mean_sd, aes(x = Group, y = mean, fill = control)) +
  geom_point(aes(y = mean), size = 1, color = "blue") + # Mean points
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) + # Error bars for mean ± sd
  geom_point(aes(y = RR_from_param), size = 1, color = "red", shape = 17) + # Points for true RR
  facet_grid(level + nSample ~ nCellType + control,  scales = "free_x", space = "free_x", drop = TRUE) + 
  theme_minimal() + # Use a minimal theme for clarity
  labs(title = "Comparison of Mean ± SD and True RR",
       x = "Group", y = "Mean and True RR",
       subtitle = "Blue points: mean, Red triangles: true RR, Error bars: mean ± SD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()
