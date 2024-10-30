# Goal: To test change between time in cell type proportion 
#       using estimateEffect from STM
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
require(pals)
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

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

# first read in scSTM data
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/", 
                    pattern = "^sims.*(nsample6\\.rds|nsample12\\.rds)$")
res_all_file <- vector(mode = "list")
for(file_name in files){
  # read in sims data
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  
  # read in scSTMseq data
  scSTM_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/scSTM_Content_Batch_Prevalence_Time/scSTM_",
                       set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  } 
  K <- scSTMobj$settings$dim$K
  ref <- "Time1"
  names(ref) <- "Time"
  time_effect <- estimateEffect(1:K ~ Time, ref.vec = ref,
                                stmobj = scSTMobj, 
                                uncertainty = "Global", nsims=30)
  
  res <- summary(time_effect)
  results_list <- lapply(res$tables, function(tbl) {
    tbl["TimeTime2", ]  # Select the row corresponding to "TimeTime2"
  })
  res <- do.call(cbind, results_list)
  colnames(res) <- paste0("Topic ", seq_along(results_list))
  rownames(res) <- c("Estimate", "Std", "t_value", "p_value")
  res_all_file[[set_level]] <- res
  cat(file_name, "\n")
}

saveRDS(res_all_file, file = "res/composition_change/estimateEffect_multiple_patient_minusSampleVariation_final.rds")

# ################# plot an example ################################
files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims/",
                    pattern = "^sims.*c13_(nsample6\\.rds|nsample12\\.rds)$")
file_name <- files[10]
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

# read in scSTMseq data
scSTM_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/scSTM_Content_Batch_Prevalence_Time/scSTM_",
                     set_level,".rds")
if(file.exists(scSTM_name)){
  scSTMobj <- readRDS(scSTM_name)
  scSTMobj <- select_top_scSTM(scSTMobj)
}
K <- scSTMobj$settings$dim$K
sims <- scSTMobj$settings$sce
data <- as.data.frame(colData(sims))
time1_data <- subset(data, Time == "Time1")
time2_data <- subset(data, Time == "Time2")
time1_props <- table(time1_data$Group) / nrow(time1_data)
time2_props <- table(time2_data$Group) / nrow(time2_data)
prop_diff <- time2_props - time1_props
prop_diff

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
structure_plot(scSTMobj, topics = 1:K, grouping = sims$Group, n = 2000, gap = 5)
ref <- "Time1"
names(ref) <- "Time"
time_effect <- estimateEffect(1:K ~ Time, ref.vec = ref,
                              stmobj = scSTMobj,
                              uncertainty = "Global", nsims=2)
plot(time_effect, covariate = "Time", model = scSTMobj,
     #method = "difference",
     method = "difference", # ref="Time1",alt="Time2",
     printlegend = T, labeltype = "custom", custom.labels = paste0("Topic ", 1:K),
     linecol = "black")


     # xlab = "proportion change between timepoints", moderator = "expansion",
     # moderator.value = "E", linecol = "skyblue", printlegend = FALSE,
     # main = "Expected Cell Group Proportion Change Before and After Treatment")
#### test #####
# I want to test the following
# Four out of five topics have significant p-values.
# Among the four significant topics, two should have negative estimates and two should have positive estimates.
# The last topic should have an insignificant p-value.
library(stringr)
res_all_file <- readRDS("res/composition_change/estimateEffect_multiple_patient_final.rds")
res_all_file_neg <- res_all_file[grep("neg", names(res_all_file))]
res_all_file_pos <- res_all_file[grep("pos", names(res_all_file))]

############ first check for positive condition ####################
check_pos_conditions <- function(df) {
  nCellType <- dim(df)[2]
  significance_threshold <- 0.05
  
  # Count significant p-values
  significant_indices <- which(df["p_value", ] < significance_threshold)
  insignificant_indices <- which(df["p_value", ] >= significance_threshold)
  
  # Extract the significant estimates
  significant_estimates <- df["Estimate", significant_indices]
  negative_estimates <- sum(significant_estimates < 0)
  
  # Define conditions based on nCellType
  if (nCellType == 5 && length(insignificant_indices) == 1 && negative_estimates == 2) {
    return(TRUE)
  } else if (nCellType == 9 && length(insignificant_indices) == 1 && negative_estimates == 5) {
    return(TRUE)
  } else if (nCellType == 13 && length(insignificant_indices) == 1 && negative_estimates == 4) {
    return(TRUE)
  } else if (nCellType == 17 && length(insignificant_indices) == 2 && negative_estimates == 6) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
check_neg_conditions <- function(df) {
  significance_threshold <- 0.05
  
  # Count significant p-values
  significant_indices <- which(df["p_value", ] < significance_threshold)

  # Define conditions based on nCellType
  if (length(significant_indices) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# Apply the function to each list in your dataset
results <- lapply(res_all_file_pos, check_pos_conditions)
results <- lapply(res_all_file_neg, check_neg_conditions)
                  
c5_list <- list()
c9_list <- list()
c13_list <- list()
c17_list <- list()
# Loop through each element in the original list
for (name in names(results)) {
  if (grepl("c5", name)) {
    c5_list[[name]] <- results[[name]]
  } else if (grepl("c9", name)) {
    c9_list[[name]] <- results[[name]]
  } else if (grepl("c13", name)) {
    c13_list[[name]] <- results[[name]]
  }else if (grepl("c17", name)) {
    c17_list[[name]] <- results[[name]]
  }
}

c5_list <- do.call(rbind, c5_list)
c9_list <- do.call(rbind, c9_list)
c13_list <- do.call(rbind, c13_list)
c17_list <- do.call(rbind, c17_list)

table(c17_list)

check_half_positive <- function(df) {
  estimates <- df["Estimate",]

  # Count the number of positive and negative values
  positive_count <- sum(estimates > 0)
  negative_count <- sum(estimates < 0)

  # Calculate the absolute difference
  absolute_difference <- abs(positive_count - negative_count)
  absolute_difference
}

results <- lapply(res_all_file_pos, check_half_positive)


c5_list <- list()
c9_list <- list()
c13_list <- list()
c17_list <- list()
# Loop through each element in the original list
for (name in names(results)) {
  if (grepl("c5", name)) {
    c5_list[[name]] <- results[[name]]
  } else if (grepl("c9", name)) {
    c9_list[[name]] <- results[[name]]
  } else if (grepl("c13", name)) {
    c13_list[[name]] <- results[[name]]
  }else if (grepl("c17", name)) {
    c17_list[[name]] <- results[[name]]
  }
}

c5_list <- do.call(rbind, c5_list)
c9_list <- do.call(rbind, c9_list)
c13_list <- do.call(rbind, c13_list)
c17_list <- do.call(rbind, c17_list)

table(c5_list ==1)
table(c9_list ==2)
table(c13_list ==4)
table(c17_list)
