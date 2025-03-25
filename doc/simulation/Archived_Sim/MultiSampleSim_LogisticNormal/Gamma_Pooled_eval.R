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
files <- list.files(path = paste0(dir, "sims/"), pattern = "Null")
res <- vector(mode = "list")
nsims <- 1e6

for(file_name in files){
  # read in sims data
  set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
  
  ######## Pooled Version ####################
  # read in scSTMseq data
  # scSTM_name <- paste0(dir, "scSTM_Pooled_noSample_noContent_Prevalence_TimeandResponse/scSTM_",
  #                      set_level,".rds")
  scSTM_name <- paste0(dir, "scSTM_Pooled_Sample_noContent_Prevalence_TimeandResponse/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  }else{
    next
  }
  
  K <- scSTMobj$settings$dim$K
  sce <- scSTMobj$settings$sce
  # theta <- scSTMobj$theta
  # rownames(theta) <- scSTMobj$DocName 
  # theta_t1 <- theta[match(sce[,sce$Time == "Time1"]$Cell, rownames(theta)),]
  # theta_t2 <- theta[match(sce[,sce$Time == "Time2"]$Cell, rownames(theta)),]
  # metadata <- data.frame(
  #   Cell = sce$Cell,
  #   Sample = sce$Sample,
  #   Group = sce$Group,
  #   Timepoint = sce$Time
  # )
  # 
  # proportion_results <- data.frame()
  # for (timepoint in unique(metadata$Timepoint)) {
  #   # Filter for the current Sample and Timepoint
  #   subset_data <- metadata %>%
  #     filter(Timepoint == timepoint)
  #   theta_subset <- theta[rownames(theta) %in% subset_data$Cell, ]
  #   
  #   # Calculate proportions for each group
  #   proportions <- colSums(theta_subset) / sum(colSums(theta_subset))
  #   proportion_results <- rbind(proportion_results, data.frame(Timepoint = timepoint, t(proportions)))
  # }
  # # Loop through each unique Sample and Timepoint combination
  # for (sample in unique(scSTMobj$sampleID)) {
  #   for (timepoint in unique(metadata$Timepoint)) {
  #     # Filter for the current Sample and Timepoint
  #     subset_data <- metadata %>%
  #       filter(Sample == sample, Timepoint == timepoint)
  #     theta_subset <- theta[rownames(theta) %in% subset_data$Cell, ]
  # 
  #     # Calculate proportions for each group
  #     proportions <- colSums(theta_subset) / sum(colSums(theta_subset))
  #     proportion_results <- rbind(proportion_results, data.frame(Sample = sample, Timepoint = timepoint, t(proportions)))
  #   }
  # }
  # 
  # # Set appropriate column names for the proportion results
  # colnames(proportion_results)[3:ncol(proportion_results)] <- paste0("Group", 1:(ncol(proportion_results) - 2))
  
  trueParam <- scSTMobj$settings$sce@metadata$TrueParams
  true_gamma <- as.data.frame(scSTMobj$settings$sce@metadata$TrueParams$gamma)
  # true_psi <- scSTMobj$settings$sce@metadata$TrueParams$psi
  # true_theta <- sce@metadata$TrueParams$theta
  true_gamma_long <- true_gamma %>%
    rownames_to_column(var = "Cov") %>%
    pivot_longer(cols = -Cov, names_to = "K", values_to = "True_Gamma")
  true_gamma_long$K <- as.numeric(gsub("K", "", true_gamma_long$K))
  
  sim_gamma <- gammaPosterior(model = scSTMobj, nsims = nsims)
  mean_list <- list()
  sd_list <- list()

  # Loop over each list in sim_gamma
  for (i in seq_along(sim_gamma)) {
    # Calculate the column-wise mean for each matrix
    mean_list[[i]] <- apply(sim_gamma[[i]], 2, mean)

    # Calculate the column-wise standard deviation for each matrix
    sd_list[[i]] <- apply(sim_gamma[[i]], 2, sd)
  }
  sim_df <- data.frame(K = rep(seq_along(sim_gamma), each = length(mean_list[[1]])),
                       Cov = names(unlist(mean_list)),
                       Mean = unlist(mean_list),
                       SD = unlist(sd_list)) %>%
    dplyr::filter(Cov != "(Intercept)") %>%
    dplyr::mutate(Cov = dplyr::case_when(
      Cov == "TimeTime2" ~ "Timepoint",
      Cov == "ResponseResponse" ~ "Response",
      Cov == "TimeTime2:ResponseResponse" ~ "Timepoint_Response",
      TRUE ~ Cov
    )) %>%
    full_join(true_gamma_long, by = c("Cov", "K"))
  
  # lm_res <- scSTMobj$mu$lm.model
  # lm_res <- lapply(lm_res, function(x) summary(x))
  # tvalue_list <- lapply(lm_res, function(x) x$coefficients["X", "t value"])
  # tvalue_list <- do.call(c, tvalue_list)
  # pvalue_list <- lapply(lm_res, function(x) x$coefficients["X", "Pr(>|t|)"])
  # pvalue_list <- do.call(c, pvalue_list)
  # sim_df <- data.frame(K = 1:ncol(scSTMobj$mu$gamma),
  #                      Mean = scSTMobj$mu$gamma[2,],
  #                      SD = scSTMobj$mu$std.gamma[2,],
  #                      tValue = tvalue_list,
  #                      pValue = pvalue_list) %>%
  #   full_join(true_gamma_long, by = "K")
  
  
  res[[set_level]] <- sim_df
  cat(file_name, "\n")
}

saveRDS(res, file = "res/composition_change/NormalGamma_SingleResponse_Pooled_Null_nSample3_nCellType5_noBatch_StromalCell.rds")

# ############################# analysis #########################################
# # # res <- readRDS("res/composition_change/NormalGamma_ZeroPsi_nSample3_nCellType5_noBatch_StromalCell.rds")
res <- readRDS("res/composition_change/NormalGamma_SingleResponse_Pooled_Null_nSample3_nCellType5_noBatch_StromalCell.rds")
# res_null <- res[!grepl("Null", names(res))]

dat <- lapply(names(res), function(name) {
  mat <- res[[name]]
  mat$exist_change <- ifelse(0 > mat$Mean - 2*mat$SD &
                          0 < mat$Mean + 2*mat$SD,
                        FALSE,
                        TRUE)
  mat$Type <- strsplit(name, "_")[[1]][2]
  if(sum(mat$exist_change) == 0){exist_change = FALSE} else{exist_change = TRUE}
  return(exist_change)
})

dat <- do.call(rbind, dat)
dat <- dat %>% as.data.frame() %>% mutate(Type = sapply(strsplit(names(res), "_"), `[`, 2))
colnames(dat) <- c("exist_change", "Type")
### Type I error ####
dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
### Power ####
dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))

#### if we look at each topic individually #####
# res <- readRDS("res/composition_change/NormalGamma_SingleResponse_Pooled_Null_nSample3_nCellType5_noBatch_StromalCell.rds")
dat <- lapply(names(res), function(name) {
  mat <- res[[name]]
  mat$exist_change <- ifelse(0 > mat$Mean - 2*mat$SD &
                               0 < mat$Mean + 2*mat$SD,
                             FALSE,
                             TRUE)
  mat$Type <- strsplit(name, "_")[[1]][2]
  if(sum(mat$exist_change) == 0){exist_change = FALSE} else{exist_change = TRUE}
  return(mat)
})

dat <-do.call(rbind, dat)

### Type I error ####
dat %>% dplyr::filter(Type == "NullModel") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))
### Power ####
dat %>% dplyr::filter(Type == "HighVar") %>% dplyr::summarise(proportion = mean(exist_change == TRUE))

ggplot(res_null[["1727199560_NullModel"]], aes(x = K, y = Mean)) +
  geom_point(color = "blue") +            # Points for the mean values
  geom_errorbar(aes(ymin = Mean - 2*SD, ymax = Mean + 2*SD), width = 0.2) + # Error bars for the SD
  geom_point(aes(y = True_Gamma), color = "red") + # Points for true_gamma in red
  labs(x = "K", y = "Mean Gamma", title = "Mean and SD of Inferred Gamma (Blue), and True Gamma (Red) for Each Topic") +
  theme_minimal() +
  facet_grid(~Cov)


# #
# # res <- lapply(res, function(x) {
# #   lapply(x, function(y) {
# #       colnames(y) <- paste0("Group", 1:ncol(y))
# #     return(y)
# #   })
# # })
# #
# # calc_mean_sd <- function(dat_list, name, true_RR){
# #   name_split <- unlist(strsplit(name, "_"))
# #   seed <- name_split[1]
# #   effectSize <- name_split[2]
# #   dat <- data.frame()
# #
# #   result_list <- lapply(names(dat_list), function(r) {
# #     x <- dat_list[[r]]
# #     sim <- x[grep("Sim", rownames(x)), ]
# #
# #     RR_from_simDat <- x["relative_change_simDat", ]
# #     mean <- apply(sim, 2, mean)
# #     sd <- apply(sim, 2, sd)
# #
# #     # Create a temporary data frame
# #     dat.temp <- data.frame(
# #       name = name,
# #       seed = seed,
# #       Response = r,
# #       effectSize = effectSize,
# #       mean = mean,
# #       sd = sd,
# #       relative_change_simDat = RR_from_simDat,
# #       Group = names(mean)
# #     )
# #
# #     return(dat.temp)
# #   })
# #
# #   dat <- do.call(rbind, result_list)
# #   return(dat)
# # }
# #
# # mean_sd <- do.call(rbind, lapply(names(res), function(name) {
# #   calc_mean_sd(res[[name]], name)
# # }))
# # mean_sd$Group <- as.factor(mean_sd$Group)
# # mean_sd$effectSize <- factor(mean_sd$effectSize, levels = paste0("effectSize", c(1, 0.7, 0.4, 0.1)))
# #
# # mean_sd$within_range <- with(mean_sd,
# #                              relative_change_simDat >= (mean - 2*sd) &
# #                                relative_change_simDat <= (mean + 2*sd))
# #
# # mean_sd_G2 <- mean_sd %>% dplyr::filter(Group == "Group2" & Response == "Response")
# # power_data <- aggregate(within_range ~ effectSize, data = mean_sd_G2, FUN = mean)
# # ggplot(power_data, aes(x = effectSize, y = within_range)) +
# #   geom_bar(stat = "identity", position = "dodge") +
# #   labs(title = "Proportion of True Value In 2 SD of Mean for Group 2 Response Only",
# #        x = "Group 2",
# #        y = "Proportion") +
# #   theme_minimal()
# #
# #
# # power_data <- aggregate(within_range ~ Group + Response, data = mean_sd, FUN = mean)
# #
# # # Plot the power for each Group and Response
# # ggplot(power_data, aes(x = Group, y = within_range, fill = Response)) +
# #   geom_bar(stat = "identity", position = "dodge") +
# #   labs(title = "Proportion of True Value In 2 SD of Mean",
# #        x = "Group",
# #        y = "Proportion") +
# #   theme_minimal()
# #
# # unique_seeds <- unique(mean_sd$seed)
# # plots_list <- list()
# # for (s in unique_seeds) {
# #   sub_mean_sd <- mean_sd %>% filter(seed == s)
# #
# #   p <- ggplot(sub_mean_sd, aes(x = Group, y = mean, color = Group)) +
# #     geom_point(size = 3) +  # Plot the means
# #     geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
# #     geom_point(aes(y = relative_change_simDat), shape = 4, size = 3, color = "black") +  # Add relative_change_simDat
# #     geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add a horizontal line at y = 0
# #     facet_grid(effectSize ~ Response, scales = "free_y") +  # Facet by effectSize and Response
# #     theme_minimal() +
# #     labs(title = paste("Mean and SD by Group for Seed", s),  # Include seed in the title
# #          y = "Mean with SD",
# #          x = "Group") +
# #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #
# #   # Store the plot in the list with the seed as the key
# #   plots_list[[as.character(s)]] <- p
# # }
# #
# # sub_mean_sd <- mean_sd %>% dplyr::filter(seed == "1726189900")
# # ggplot(sub_mean_sd, aes(x = Group, y = mean, color = Group)) +
# #   geom_point(size = 3) +  # Plot the means
# #   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
# #   geom_point(aes(y = relative_change_simDat), shape = 4, size = 3, color = "black") +  # Add relative_change_simDat
# #   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
# #   facet_grid(effectSize ~ Response, scales = "free_y") +  # Facet by name and effectSize
# #   theme_minimal() +
# #   labs(title = "Mean and SD by Group for each name and effectSize",
# #        y = "Mean with SD",
# #        x = "Group") +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #
# # RR_dat <- mean_sd %>% dplyr::select("control", "level", "nCellType", "nSample", "Group", "RR_from_param")
# #
# # aggregated_mean_sd <- mean_sd %>%
# #   group_by(control, nCellType, level, nSample, Group) %>%
# #   summarize(mean = mean(mean),
# #             sd = mean(sd),
# #             .groups = 'drop') %>%
# #   left_join(RR_dat, by = c("control", "level", "nCellType", "nSample", "Group"))
# #
# # aggregated_mean_sd$Group <- factor(aggregated_mean_sd$Group, levels = paste0("Group", 1:17))
# # aggregated_mean_sd$nCellType <- factor(aggregated_mean_sd$nCellType, levels= paste0("c", c(5, 9, 13, 17)))
# # # mean_sd <- mean_sd %>%
# # #   mutate(diff_mean_RR = mean - RR_from_param)
# #
# # png("res/composition_change/composition_change_mean_plot.png", res= 300, width = 6000, height = 3000)
# # ggplot(aggregated_mean_sd, aes(x = Group, y = mean, fill = control)) +
# #   geom_point(aes(y = mean), size = 1, color = "blue") + # Mean points
# #   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) + # Error bars for mean ± sd
# #   geom_point(aes(y = RR_from_param), size = 1, color = "red", shape = 17) + # Points for true RR
# #   facet_grid(level + nSample ~ nCellType + control,  scales = "free_x", space = "free_x", drop = TRUE) +
# #   theme_minimal() + # Use a minimal theme for clarity
# #   labs(title = "Comparison of Mean ± SD and True RR",
# #        x = "Group", y = "Mean and True RR",
# #        subtitle = "Blue points: mean, Red triangles: true RR, Error bars: mean ± SD") +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # dev.off()
