# Goal: To test change between time in cell type proportion 
#       using estimateEffect from STM
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mcluster)
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

cluster_scSTM <- function(scSTMobj, sims) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  res_cluster <- res_cluster[match(colnames(sims), names(res_cluster))]
  return(adjustedRandIndex(res_cluster, sims$Group))
}


r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
n <- 50
# first read in scSTM data
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/", pattern = "sims*")
files <- list.files(
  path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/",
  pattern = "^sims.*(nsample6\\.rds|nsample12\\.rds)$")

res <- vector(mode = "list")

for(file_name in files){
  # read in sims data
  set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
  
  # read in scSTM_nC_P data
  scSTM_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/scSTM_Content_Batch_Prevalence_Time/scSTM_",set_level,".rds")
  if(file.exists(scSTM_name)){
    scSTMobj <- readRDS(scSTM_name)
    scSTMobj <- select_top_scSTM(scSTMobj)
  } 
  # 
  # file_path <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims/sims_", set_level, ".rds")
  # sims <- safe_readRDS(file_path)
  # # meta <- colData(sims) %>% 
  # #   data.frame()
  # # meta <- meta[rownames(meta) %in% scSTM_nC_P$DocName,]
  
  # structure plot
  sims <- scSTMobj$settings$sce
  K <- length(unique(sims$Group))
  
  r.file <- paste0("R/",list.files("R/"))
  sapply(r.file, source)
  
  x <- structure_plot(scSTMobj, topics = 1:K, grouping = sims$Group, n = 2000, gap = 5)
  max_topic_per_cell <- x$data %>%
    group_by(source_Doc) %>%
    slice_max(prop, with_ties = FALSE) %>%
    ungroup()
  matched_data <- colData(sims) %>%
    as.data.frame %>%
    dplyr::right_join(max_topic_per_cell, by = c("Cell" = "source_Doc")) 
  
  
  topic_group_counts <- matched_data %>%
    group_by(Group, topic) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Step 2: Calculate the total number of cells in each group
  group_totals <- matched_data %>%
    group_by(Group) %>%
    summarise(total = n(), .groups = 'drop')
  
  # Step 3: Calculate the likelihood
  likelihoods <- topic_group_counts %>%
    dplyr::left_join(group_totals, by = "Group") %>%
    dplyr::mutate(likelihood = count / total) %>%
    dplyr::select(Group, topic, likelihood) %>%
    group_by(Group) %>%
    slice_max(likelihood, with_ties = FALSE) %>%
    ungroup()

  ############## ############## ############## ############## ##############
  ############## extract DE genes in truth ##########################
  ############## ############## ############## ############## ############## ##############
  truth <- rowData(sims) %>% as.data.frame()
  # truth <- truth[rownames(truth) %in% scSTM_C_P$vocab,]
  groupDE_cols <- grep("DEFacGroup", colnames(truth))
  # subset to only upregulated genes
  # rows_to_keep <- apply(truth[, groupDE_cols], 1, function(row) any(row > 1))
  # up_gene <- truth[rows_to_keep, groupDE_cols]

  est_scSTMobj <- labelTopics(scSTMobj, n = n)
  rownames(est_scSTMobj$topics) <- paste0("topic_", rownames(est_scSTMobj$topics))
  scSTM_de_gene <- data.frame(
    sim = set_level,
    seed = unlist(strsplit(set_level, "_"))[1],
    control = unlist(strsplit(set_level, "_"))[2],
    level = unlist(strsplit(set_level, "_"))[3],
    nCellType = unlist(strsplit(set_level, "_"))[4],
    nSample = unlist(strsplit(set_level, "_"))[5],
    adjustedRandIndex = cluster_scSTM(scSTMobj, sims)
  )
  for(i in 1:nrow(likelihoods)){
    group <- likelihoods$Group[i]
    topics <- likelihoods$topic[i]
    up_gene_topic <- truth %>% 
      dplyr::select(paste0("DEFac", group)) %>%
      dplyr::filter(.[[1]] > 1) # select the genes that are up-regulated 

    top_DEgene_ratio = paste0(
      sum(est_scSTMobj$topics[topics,] %in% rownames(up_gene_topic)), 
      "/", n)
    colname <- paste0("DERatio", group)
    scSTM_de_gene[[colname]] <- top_DEgene_ratio
  }
  list_name <- paste0("nCellType",nrow(likelihoods))
  res[[list_name]] <- rbind(res[[list_name]], scSTM_de_gene)
  # # for batch
  # batchDE_cols <- grep("BatchFac", colnames(truth))
  # batch_up_gene <- truth[apply(truth[, batchDE_cols], 1, function(row) any(row > 1)), ]
  cat(file_name, "\n")
  rm(scSTMseq)
}
saveRDS(res, file = "res/de_genes_eval/multi_sample_benchmark_V4.rds")

################################################################################
########################### generate plot ######################################
################################################################################

res <- readRDS("res/de_genes_eval/multi_sample_benchmark_V4.rds")

# Function to convert proportion to numeric
convert_proportion <- function(x) {
  sapply(strsplit(x, "/"), function(y) as.numeric(y[1]) / as.numeric(y[2]))
}

# Convert proportions to numeric
df_list <- vector("list", length(res))

# Loop through each data frame in the list `res`
for (i in seq_along(res)) {
  df_list[[i]] <- res[[i]] %>%
    mutate(across(starts_with("DERatio"), convert_proportion)) %>%
    rename_with(~ gsub("DERatio", "", .), starts_with("DERatio")) %>%
    mutate(nCellType = recode(nCellType, "c5" = "Num_of_CellGroup = 5", "c9" = "Num_of_CellGroup = 9", 
                              "c13" = "Num_of_CellGroup = 13")) %>%
    pivot_longer(cols = starts_with("Group"),
                 names_to = "DEG_by_Group",
                 values_to = "value")
}

# Combine all the transformed data frames into one data frame
df_combined <- bind_rows(df_list)
df_combined$nCellType <- factor(df_combined$nCellType, levels = c("Num_of_CellGroup = 5", "Num_of_CellGroup = 9", "Num_of_CellGroup = 13"))
df_combined$DEG_by_Group <- factor(df_combined$DEG_by_Group, levels = paste0("Group", 1:13))
# Create the box plot
png("res/de_genes_eval/multi_sample_benchmark_V4_DEGs_proportion_fig.png", width = 1300, height = 600, res =120)
ggplot(df_combined, aes(x = DEG_by_Group, y = value)) +
  geom_boxplot() +
  facet_wrap(~ nCellType, scales = "free_x") +
  labs(title = "Overlap of Topic-Defined and Upregulated Genes by Cell Group",
       x = "Cell Group",
       y = "Proportion of Overlapping Genes") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 13)
  )
dev.off()


################################################################################
########## plot C_P different proportion, separate by nSample, nCellType ######
###############################################################################
# there are a few C_nP had very bad performance. Not sure why but remove them
C_P_df <- C_P_df %>% filter(adjustedRandIndex > 0)

# Convert proportions to numeric
C_P_df <- C_P_df %>%
  mutate(across(starts_with("num_of_top_"), convert_proportion))

# Transform data to long format
C_P_df <- C_P_df %>%
  rename_with(~ gsub("num_of_", "", .), starts_with("num_of_")) %>%
  mutate(nCellType = recode(nCellType, "c2" = "Number_of_CellType = 2", "c4" = "Number_of_CellType = 4"),
         nSample = recode(nSample, "nsample1" = "Number_of_Sample = 1", "nsample5" = "Number_of_Sample = 5")) %>%
  rename("top_DEgene" = "CellType_Genes",
         "top_Timegene" = "Content_Genes") %>%
  pivot_longer(cols = c(CellType_Genes, Content_Genes),
               names_to = "Top_Genes",
               values_to = "value")

# Create the box plot
png("res/de_genes_eval/C_P_DEGs_proportion_fig.png", width = 2500, height = 2000, res =300)
ggplot(C_P_df, aes(x = Top_Genes, y = value)) +
  geom_boxplot() +
  facet_wrap(~ nCellType + nSample) +
  labs(title = "Proportion of Top Genes Defined by Different Parts of the Model",
       subtitle = "With the Content, With the Prevalence Model",
       x = "Models",
       y = "Proportion of Overlapping Genes") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 13)
  )
dev.off()

