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
    max_value <- max(all_values)
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

files <- list.files(
  path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims/",
  pattern = "^sims.*(nsample6\\.rds|nsample12\\.rds)$")

res_all_file <- vector(mode = "list")
for(file_name in files){
  # read in sims data
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  
  scSTM_C_Time_P_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_Time_Prevalence_Time/scSTM_",set_level,".rds")
  if(file.exists(scSTM_C_Time_P_name)){
    scSTMobj <- readRDS(scSTM_C_Time_P_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  } 
  K <- scSTMobj$settings$dim$K
  time_effect <- estimateEffect(1:K ~ Time, 
                                stmobj = scSTMobj, 
                                uncertainty = "Global", nsims=100)
  

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

saveRDS(res_all_file, file = "res/composition_change/estimateEffect_multiple_patient_V3.rds")


#### test #####
# I want to test the following 
# Four out of five topics have significant p-values.
# Among the four significant topics, two should have negative estimates and two should have positive estimates.
# The last topic should have an insignificant p-value.
res_all_file <- readRDS("res/composition_change/estimateEffect_multiple_patient_V3.rds")

check_conditions <- function(df) {
  # Define the threshold for significance
  significance_threshold <- 0.05
  # Count significant p-values
  significant_indices <- which(df["p_value", ] < significance_threshold)
  insignificant_indices <- which(df["p_value", ] >= significance_threshold)
  num_sig <- dim(df)[2]-1
  # Check if four out of five p-values are significant
  if (length(significant_indices) == num_sig && length(insignificant_indices) == 1) {
    # Check if two of the significant estimates are negative and two are positive
    significant_estimates <- df["Estimate", significant_indices]
    negative_estimates <- sum(significant_estimates < 0)
    positive_estimates <- sum(significant_estimates > 0)

    if (negative_estimates == num_sig/2 && positive_estimates == num_sig/2) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Apply the function to each list in your dataset
results <- lapply(res_all_file, check_conditions)

c5_list <- list()
c9_list <- list()
c13_list <- list()

# Loop through each element in the original list
for (name in names(results)) {
  if (grepl("c5", name)) {
    c5_list[[name]] <- results[[name]]
  } else if (grepl("c9", name)) {
    c9_list[[name]] <- results[[name]]
  } else if (grepl("c13", name)) {
    c13_list[[name]] <- results[[name]]
  }
}

c5_list <- do.call(rbind, c5_list)
c9_list <- do.call(rbind, c9_list)
c13_list <- do.call(rbind, c13_list)

cat(table(c13_list))


check_half_positive <- function(df) {
  estimates <- x["Estimate",]

  # Count the number of positive and negative values
  positive_count <- sum(estimates > 0)
  negative_count <- sum(estimates < 0)

  # Calculate the absolute difference
  absolute_difference <- abs(positive_count - negative_count)
  absolute_difference
}
results <- lapply(res_all_file, check_half_positive)
c5_list <- list()
c9_list <- list()
c13_list <- list()

# Loop through each element in the original list
for (name in names(results)) {
  if (grepl("c5", name)) {
    c5_list[[name]] <- results[[name]]
  } else if (grepl("c9", name)) {
    c9_list[[name]] <- results[[name]]
  } else if (grepl("c13", name)) {
    c13_list[[name]] <- results[[name]]
  }
}

c5_list <- do.call(rbind, c5_list)
c9_list <- do.call(rbind, c9_list)
c13_list <- do.call(rbind, c13_list)

table(c13_list ==1)
