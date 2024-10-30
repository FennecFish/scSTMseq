# Goal: This script is trying to run GammaError for all scSTMseq model
# specifically those on the 
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
library(iCOBRA)
# 
# args <- commandArgs(trailingOnly = TRUE)
# file_name <- args[1]
# file_name <- basename(file_name)
# cat(file_name, "\n")
# 
# set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)
# dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
# 
# r.file <- paste0("R/",list.files("R/"))
# sapply(r.file, source)
# 
# scSTM_name <- paste0(dir, "scSTM_LinearMixed_noSample_noContent_Prevalence_TimeandResponse/scSTM_",
#                      set_level,".rds")
# if(file.exists(scSTM_name)){
#   scSTMobj <- readRDS(scSTM_name)
#   scSTMobj <- select_top_scSTM(scSTMobj)
# }else{
#   next
# }
# formula = ~Time + (1|Sample) # mixed linear model formula
# res <- GammaError(model = scSTMobj, formula = formula, null_formula = ~(1|Sample), nsims = 100)
# saveRDS(res, file = paste0(dir, "GammaError_MixedModel/GammaError_",set_level,".rds"))
# 

# ##### analysis ####'
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
files <- list.files(paste0(dir, "GammaError_MixedModel/"))
res <- vector("list")
for(file in files){
  set_level <- sub("^GammaError_(.*)\\.rds$", "\\1", file)
  model <- readRDS(paste0(dir, "GammaError_MixedModel/", file))
  simplex <- dim(model$SimEta[[1]])[2]
  pValue_output_lists <- vector("list", simplex)
  Chisq_output_lists <- vector("list", simplex)
  # Loop over the K lists and extract the respective nested list from 'model'
  for (i in 1:simplex) {
    pValue.list <- lapply(model$model.diff, function(x) {
      k.model <- x[[i]]
      pValue <- k.model$`Pr(>Chisq)`[2]
      pValue
    })
    pValue.list <- do.call(c, pValue.list)

    Chisq.list <- lapply(model$model.diff, function(x) {
      k.model <- x[[i]]
      Chisq <- k.model$Chisq[2]
      Chisq
    })
    Chisq.list <- do.call(c, Chisq.list)

    pValue_output_lists[[i]] <- pValue.list
    Chisq_output_lists[[i]] <- Chisq.list
  }
  pValue <- do.call(cbind, pValue_output_lists) %>% as.data.frame()
  colnames(pValue) <- paste0("topic_", 1:ncol(pValue))
  rownames(pValue) <- paste0("Replicate_", 1:nrow(pValue))

  Chisq <- do.call(cbind, Chisq_output_lists) %>% as.data.frame()
  colnames(Chisq) <- paste0("topic_", 1:ncol(Chisq))
  rownames(Chisq) <- paste0("Replicate_", 1:nrow(Chisq))

  res[[set_level]] <- list(pValue = pValue, Chisq = Chisq)
}
name <- basename(dir)
saveRDS(res, file = paste0("res/GammaError_MixedLinearModelComparison_",name,".rds"))


##### plotting estimated gamma distribution ###
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
files <- list.files(paste0(dir, "GammaError_MixedModel/"), pattern = "HighVar")
files <- files[1:5]
plot_list <- list()
for(file in files){
  set_level <- sub("^GammaError_(.*)\\.rds$", "\\1", file)
  model <- readRDS(paste0(dir, "GammaError_MixedModel/", file))

  gamma_topics <- lapply(1:4, function(topic_num) {
    topic_name <- paste0("topic_", topic_num)
    lapply(gamma, function(x) x[topic_name, "TimeTime2"])
  })
  
  # Combine the values for all topics into a single vector and label each with its topic
  gamma_topics_combined <- do.call(rbind, lapply(1:4, function(i) {
    data.frame(Gamma = do.call(c, gamma_topics[[i]]), Topic = paste0("Topic_", i))
  }))
  
  p <- ggplot(gamma_topics_combined, aes(x = Gamma)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +  # Adjust binwidth as needed
    facet_wrap(~ Topic, ncol = 2) +  # Facet by topic, arranging in 2 columns
    labs(title = "Distribution of Estimated Gamma Fitted with Posterior Eta in Null Model", 
         subtitle = "set_level",
         x = "Estimated Gamma", 
         y = "Frequency") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),  # Centered title with larger font
      axis.title.x = element_text(size = 14),  # X axis title size
      axis.title.y = element_text(size = 14),  # Y axis title size
      strip.text = element_text(size = 14)  # Facet label text size
    )
  plot_list[[length(plot_list) + 1]] <- p
}

saveRDS(plot_list, "res/alt_plot.rds")


#################################################################################
#################################  Functions ####################################
#################################################################################
num_elements <- function(obj) {
  if (is.data.frame(obj)) {
    return(nrow(obj) * ncol(obj))  # For data frames, return the total number of elements (rows * columns)
  } else {
    return(length(obj))  # For vectors, return the length (number of elements)
  }
}

plot_data_generate <- function(threshold, pValueList, topic_level = TRUE){
  
  if(topic_level){
    exist.change.list <- lapply(thresholds, function(thr) {
      lapply(pValueList, function(method) {
        t(apply(method, 1, function(row) ifelse(row < thr, TRUE, FALSE)))
      })
    })
  }else{
    exist.change.list <- lapply(thresholds, function(thr) {
      lapply(pValueList, function(method) {
        apply(method, 1, function(row) ifelse(any(row < thr), TRUE, FALSE))
      })
    })
  }
  
  
  TypeIError_data <- data.frame()
  
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    
    for (method in names(exist.change.list[[i]])) {
      positive <- sum(exist.change.list[[i]][[method]])
      total <- num_elements(exist.change.list[[i]][[method]])
      
      TypeIError <- positive/total
      
      # Store the results in a dataframe
      TypeIError_data <- rbind(TypeIError_data, data.frame(
        Method = method,
        Threshold = threshold,
        TypeIError = TypeIError
      ))
    }
  }
  return(TypeIError_data)
}

#################################################################################
#################################  Type I Error Plots  ##########################
#################################################################################

res <- readRDS("res/GammaError_MixedLinearModelComparison_nSample3_nCellType5_noBatch_StromalCell.rds")
null.res <- res[grep("NullModel",names(res))]




### approach 1, take the mean of chisq among replicates
# appraoch 2 is 50 percentile of pvalue
# appraoch 3 is 75 percentile of pvalue
null.chisq <- lapply(null.res, function(x) {
  mat <- x$Chisq
  mean.chisq <- colMeans(mat)
  mean.chisq
})
null.chisq <- do.call(rbind, null.chisq)
null.chisq.pvalues <- apply(null.chisq, c(1, 2), function(x) pchisq(x, df = 1, lower.tail = FALSE))

Pval50Perc <- lapply(null.res, function(x){
  mat <- x$pValue
  apply(mat, 2, function(col) quantile(col, probs = 0.50)) 
})
Pval50Perc <- do.call(rbind, Pval50Perc)

Pval75Perc <- lapply(null.res, function(x){
  mat <- x$pValue
  apply(mat, 2, function(col) quantile(col, probs = 0.75)) 
})
Pval75Perc <- do.call(rbind, Pval75Perc)

res.pValue <- list(MeanChisqPvalue = null.chisq.pvalues, Pval50Perc = Pval50Perc, Pval75Perc = Pval75Perc)

##### adjust pvalue ####

res.adjpValue <- lapply(res.pValue, function(rawPvalue){
  # rawPvalue <- res.pValue[[name]]
  adj.pValue <- apply(rawPvalue, 1, function(row) p.adjust(row, method = "fdr"))
  t(adj.pValue)
})

#### prepring for plots ######
threshold <- c(0.01, 0.05, 0.1)

topic_plot <- plot_data_generate(threshold, pValueList = res.adjpValue, topic_level = TRUE)
p1 <- ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  geom_point(size = 3, position = position_jitter(height = 1e-10)) + 
  labs(title = "Type I Error on Topic Level", 
       x = "Threshold", 
       y = "Type I Error") +
  theme_minimal() +
  scale_x_continuous(breaks = thresholds) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    legend.position = "none"  # Remove the legend
  )

sample_plot <- plot_data_generate(threshold, pValueList = res.adjpValue, topic_level = FALSE)
p2 <- ggplot(sample_plot, aes(x = Threshold, y = TypeIError, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  geom_point(size = 3, position = position_jitter(height = 1e-10)) + 
  labs(title = "Type I Error on Sample Level", 
       x = "Threshold", 
       y = "Type I Error") +
  theme_minimal() +
  scale_x_continuous(breaks = thresholds) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14), 
    legend.title = element_blank()  # Remove legend title
  )

require(gridExtra)
grid.arrange(p1, p2, ncol=2)



#################################################################################
##################################  Power Plots  ###############################
#################################################################################

res <- readRDS("res/GammaError_MixedLinearModelComparison_nSample3_nCellType5_noBatch_StromalCell.rds")
alt.res <- res[grep("HighVar",names(res))]

chisq <- lapply(alt.res, function(x) {
  mat <- x$Chisq
  mean.chisq <- colMeans(mat)
  mean.chisq
})
chisq <- do.call(rbind, chisq)
chisq.pvalues <- apply(chisq, c(1, 2), function(x) pchisq(x, df = 1, lower.tail = FALSE))

Pval50Perc <- lapply(alt.res, function(x){
  mat <- x$pValue
  apply(mat, 2, function(col) quantile(col, probs = 0.50)) 
})
Pval50Perc <- do.call(rbind, Pval50Perc)

Pval75Perc <- lapply(alt.res, function(x){
  mat <- x$pValue
  apply(mat, 2, function(col) quantile(col, probs = 0.75)) 
})
Pval75Perc <- do.call(rbind, Pval75Perc)

res.pValue <- list(MeanChisqPvalue = chisq.pvalues, Pval50Perc = Pval50Perc, Pval75Perc = Pval75Perc)

##### adjust pvalue ####

res.adjpValue <- lapply(res.pValue, function(rawPvalue){
  # rawPvalue <- res.pValue[[name]]
  adj.pValue <- apply(rawPvalue, 1, function(row) p.adjust(row, method = "fdr"))
  t(adj.pValue)
})

#### prepring for plots ######
threshold <- c(0.01, 0.05, 0.1)

topic_plot <- plot_data_generate(threshold, pValueList = res.adjpValue, topic_level = TRUE)
colnames(topic_plot)[3] <- "Power"
p1 <- ggplot(topic_plot, aes(x = Threshold, y = Power, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  geom_point(size = 3, position = position_jitter(height = 1e-10)) + 
  labs(title = "Power on Topic Level", 
       x = "Threshold", 
       y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = thresholds) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    legend.position = "none"  # Remove the legend
  )

sample_plot <- plot_data_generate(threshold, pValueList = res.adjpValue, topic_level = FALSE)
colnames(sample_plot)[3] <- "Power"
p2 <- ggplot(sample_plot, aes(x = Threshold, y = Power, color = Method, group = Method)) +
  # geom_line(size = 1) +  # Line plot for each method
  geom_point(size = 3) + 
  labs(title = "Power on Sample Level", 
       x = "Threshold", 
       y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = thresholds) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14), 
    legend.title = element_blank()  # Remove legend title
  )

require(gridExtra)
grid.arrange(p1, p2, ncol=2)

