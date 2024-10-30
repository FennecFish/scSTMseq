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

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType5_noBatch_StromalCell/"
# first read in scSTM data
files <- list.files(path = paste0(dir, "sims/"))
res <- vector(mode = "list")
nsims <- 100

for(file_name in files){
  # read in sims data
  set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
  
  ######## Pooled Version ####################
  # read in scSTMseq data
  scSTM_name <- paste0(dir, "scSTM_LinearRegression_noContent_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }
  
  K <- scSTMobj$settings$dim$K
  sim_gamma <- gammaPosterior(model = scSTMobj, nsims = nsims, category = "ResponseR")
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
  
  
  # colData(sims) %>%
  #   as.data.frame() %>%
  #   dplyr::filter(Response == "R" & Group == "Group2" & Sample == "Sample2") %>%
  #   group_by(Time) %>%
  #   summarise(Count = n()) 
  
  proportion_df <- colData(sims) %>%
    as.data.frame() %>%
    group_by(Time, Sample, Group, Response) %>%
    summarise(Count = n()) %>%
    group_by(Time, Sample) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  proportion_change_df <- proportion_df %>%
    # select(Sample, Time, Group, Proportion) %>%
    dplyr::select(Time, Group, Response, Sample, Proportion) %>%
    tidyr::spread(key = Time, value = Proportion) %>%
    dplyr::mutate(True_RR = (Time2 - Time1)/Time1) %>%
    group_by(Group, Response) %>%
    summarize(mean_True_RR = mean(True_RR, na.rm = TRUE))%>%
    dplyr::right_join(likelihoods, by = "Group")
  
  R_from_simDat <- proportion_change_df %>% 
    dplyr::filter(Response == "R") 
  R_relative_change_from_simDat <- R_from_simDat$mean_True_RR
  NR_from_simDat <- proportion_change_df %>% 
    dplyr::filter(Response == "NR") 
  NR_relative_change_from_simDat <- NR_from_simDat$mean_True_RR
  
  # Calculate Relative Change from sim_gamma
  R_relative_change <- (sim_gamma$t2_alt-sim_gamma$t1_alt)/sim_gamma$t1_alt
  NR_relative_change <- (sim_gamma$t2_ref-sim_gamma$t1_ref)/sim_gamma$t1_ref
  colnames(R_relative_change) <- paste0("topic_", 1:K)
  colnames(NR_relative_change) <- paste0("topic_", 1:K)
  R_relative_change <- R_relative_change[,match(R_from_simDat$Topic, colnames(R_relative_change))]
  NR_relative_change <- NR_relative_change[,match(NR_from_simDat$Topic, colnames(NR_relative_change))]
  

  
  # this is the RR calculated from the parameter
  name_split <- unlist(strsplit(set_level, "_"))
  # if(name_split[2] == "pos"){
  #   nCellType <- name_split[4]
  #   RR_from_param <- data.frame(true_sim_RR = true_sim_RR[[nCellType]],
  #                               Group = levels(sims$Group)) %>%
  #     dplyr::right_join(likelihoods, by = "Group")
  #   RR_from_param <- RR_from_param$true_sim_RR[match(proportion_change_df$Group, RR_from_param$Group)]
  # }else{
  #   RR_from_param <- rep(0, K)
  # }
  
  R_relative_change <- rbind(R_relative_change_from_simDat, R_relative_change)
  NR_relative_change <- rbind(NR_relative_change_from_simDat, NR_relative_change)
  rownames(R_relative_change) <- c("relative_change_simDat",paste0("Gamma", 1:nsims))
  rownames(NR_relative_change) <- c("relative_change_simDat",paste0("Gamma", 1:nsims))
  res[[set_level]] <- list(Response = R_relative_change, NonResponse = NR_relative_change)
  
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/MultiSample_VaryingBaseline_noBatch_StromalCell_LinearRegression_RelativeChange.rds")

############################# analysis #########################################
res <- readRDS("res/composition_change/MultiSample_VaryingBaseline_noBatch_StromalCell_LinearRegression_RelativeChange.rds")

res <- lapply(res, function(x) {
  lapply(x, function(y) {
      colnames(y) <- paste0("Group", 1:ncol(y))
    return(y)
  })
})

calc_mean_sd <- function(dat_list, name, true_RR){
  name_split <- unlist(strsplit(name, "_"))
  seed <- name_split[1]
  effectSize <- name_split[2]
  dat <- data.frame()

  result_list <- lapply(names(dat_list), function(r) {
    x <- dat_list[[r]]
    sim <- x[grep("Gamma", rownames(x)), ]

    RR_from_simDat <- x["relative_change_simDat", ]
    mean <- apply(sim, 2, mean)
    sd <- apply(sim, 2, sd)

    # Create a temporary data frame
    dat.temp <- data.frame(
      name = name,
      seed = seed,
      Response = r,
      effectSize = effectSize,
      mean = mean,
      sd = sd,
      relative_change_simDat = RR_from_simDat,
      Group = names(mean)
    )

    return(dat.temp)
  })

  dat <- do.call(rbind, result_list)
  return(dat)
}

mean_sd <- do.call(rbind, lapply(names(res), function(name) {
  calc_mean_sd(res[[name]], name)
}))
mean_sd$Group <- as.factor(mean_sd$Group)
mean_sd$effectSize <- factor(mean_sd$effectSize, levels = paste0("effectSize", c(1, 0.7, 0.4, 0.1)))

mean_sd$within_range <- with(mean_sd,
                             relative_change_simDat >= (mean - 2*sd) &
                               relative_change_simDat <= (mean + 2*sd))

power_data <- aggregate(within_range ~ Group + Response, data = mean_sd, FUN = mean)

# Plot the power for each Group and Response
ggplot(power_data, aes(x = Group, y = within_range, fill = Response)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Power by Group and Response",
       x = "Group",
       y = "Proportion within Mean ± SD (Power)") +
  theme_minimal()

unique_seeds <- unique(mean_sd$seed)
plots_list <- list()
for (s in unique_seeds) {
  sub_mean_sd <- mean_sd %>% filter(seed == s)

  p <- ggplot(sub_mean_sd, aes(x = Group, y = mean, color = Group)) +
    geom_point(size = 3) +  # Plot the means
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
    geom_point(aes(y = relative_change_simDat), shape = 4, size = 3, color = "black") +  # Add relative_change_simDat
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add a horizontal line at y = 0
    facet_grid(effectSize ~ Response, scales = "free_y") +  # Facet by effectSize and Response
    theme_minimal() +
    labs(title = paste("Mean and SD by Group for Seed", s),  # Include seed in the title
         y = "Mean with SD",
         x = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Store the plot in the list with the seed as the key
  plots_list[[as.character(s)]] <- p
}

sub_mean_sd <- mean_sd %>% dplyr::filter(seed == "1726189900")
ggplot(sub_mean_sd, aes(x = Group, y = mean, color = Group)) +
  geom_point(size = 3) +  # Plot the means
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
  geom_point(aes(y = relative_change_simDat), shape = 4, size = 3, color = "black") +  # Add relative_change_simDat
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(effectSize ~ Response, scales = "free_y") +  # Facet by name and effectSize
  theme_minimal() +
  labs(title = "Mean and SD by Group for each name and effectSize",
       y = "Mean with SD",
       x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
