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
library(MANOVA.RM)

##### analysis ####'
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
# manova.dir <- "1000ManovaTheta_Pooled_noContent_Prevalence_Time"
manova.dir <- "ManovaTheta_Pooled_noContent_Prevalence_Time"
paths <- Sys.glob(file.path(dir, "*", manova.dir)) 
files <- unlist(lapply(paths, function(path) list.files(path, pattern = "Manova_pValue", full.names = TRUE)))

res <- vector(mode = "list")

for(file_name in files){
  set_level <- sub("^Manova_pValue_(.*)\\.rds$", "\\1", basename(file_name))
  design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
  nCellType = unlist(strsplit(design, "_"))[2]
  ThetaManova <- readRDS(file_name)
  # rm.wts <- lapply(ThetaManova$manova.fit, function(x){
  #   dat <- x$resampling[1]
  #   return(dat)
  # })
  # 
  # rm.mats <- lapply(ThetaManova$manova.fit, function(x){
  #   dat <- x$resampling[2]
  #   return(dat)
  # })
  # 
  # dat.clean <- function(given_list){
  #   given_list <- do.call(rbind, given_list)
  #   rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
  #   colnames(given_list) <- "pValue"
  #   return(given_list)
  # }
  temp <- data.frame(pValue.wts = ThetaManova$rm.wts,
                     pValue.mats = ThetaManova$rm.mats)
  colnames(temp) <- c("pValue.wts", "pValue.mats")
  list.name <- paste0(nCellType, "_", set_level)
  res[[list.name]] <- temp
  cat(file_name, "\n")
}
name <- sub("^ManovaTheta_", "",basename(manova.dir))
saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_",basename(dir), "_", name,".rds"))

#################################################################################
################################## plot #########################################
#################################################################################
plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
  exist.change.list <- lapply(thresholds, function(thr) {
    res <- pValueList %>%
      mutate(wts.sig = wts < thr,
             mats.sig = mats < thr)
    return(res)
  })
  ProportionOfTrue <- data.frame()

  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    if(is.null(by_group)){
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        summarize(
          wts_sig_proportion = mean(wts.sig),
          mats_sig_proportion = mean(mats.sig)
        ) %>%
        mutate(Threshold = threshold)
    }else{
      # if we want to group the power/ type I error by effect size
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        group_by(EffectSize, nCellType) %>%
        summarize(
          wts_sig_proportion = mean(wts.sig),
          mats_sig_proportion = mean(mats.sig),
          .groups = "drop"
        ) %>%
        mutate(Threshold = threshold)
    }
    # Store the results in a dataframe
    ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
  }
  return(ProportionOfTrue)
}

# For power only
dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_SingleResponse_noBatch_StromalCell_Pooled_noContent_Prevalence_Time.rds")
# For Type I Error of 1000 sims
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_SingleResponse_1000ManovaTheta_Pooled_noContent_Prevalence_Time.rds")
null.dat <- dat[grep("NullModel",names(dat))]

check.na <- lapply(null.dat, function(x){
  y <- is.na(x[,1])
  return(y)
})
check.na <- do.call(rbind, check.na)
table(check.na)
null.dat.dist <- do.call(rbind, null.dat)
manova_long <- null.dat.dist %>%
  pivot_longer(cols = starts_with("pValue"),
               names_to = "Type",
               values_to = "pValue") %>%
  mutate(Type = sub("pValue\\.", "",Type))
# manova_long <- data.frame(
#   Value = c(null.dat.dist[,1], null.dat.dist[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(null.dat.dist))
# )
png("res/1MultiSample_SingleResponse_Simulation/Histogram_pValue_AllReplicates_SingleResponse_noBatch_StromalCell_1000scSTM_Pooled_noContent_Prevalence_Time.png", width = 1500, height = 1500, res = 220)
ggplot(manova_long, aes(x = pValue, fill = Type)) +
  geom_histogram(bins = 20, color = "black", aes(fill = Type)) +
  labs(x = "p-value", y = "Frequency", title = "Histogram of p-values under the null across all replicates") +
  scale_fill_manual(values = c("wts" = "skyblue", "mats" = "coral")) + 
  facet_wrap(~ Type) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()

manova.pValue <- lapply(null.dat, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
manova.pValue <- do.call(rbind, manova.pValue)

manova_long <- manova.pValue %>%
  pivot_longer(cols = everything(),
               names_to = "Type",
               values_to = "pValue") %>%
  mutate(Type = sub("pValue\\.", "",Type))
# manova_long <- data.frame(
#   Value = c(null.dat.dist[,1], null.dat.dist[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(null.dat.dist))
# )
png("res/1MultiSample_SingleResponse_Simulation/Histogram_pValue_50thQuantile_SingleResponse_noBatch_StromalCell_1000scSTM_Pooled_noContent_Prevalence_Time.png", width = 1500, height = 1500, res = 220)
ggplot(manova_long, aes(x = pValue, fill = Type)) +
  geom_histogram(bins = 20, color = "black", aes(fill = Type)) +
  labs(x = "p-value", y = "Frequency", title = "Histogram of p-values under the null using medium of replicates") +
  scale_fill_manual(values = c("wts" = "skyblue", "mats" = "coral")) + 
  facet_wrap(~ Type) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()


threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = manova.pValue)
topic_plot <- topic_plot %>%
  rename(WTS = wts_sig_proportion,
         MATS = mats_sig_proportion) %>%
  pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "TypeIError")

png("res/1MultiSample_SingleResponse_Simulation/TypeIError_SingleResponse_noBatch_StromalCell_1000scSTM_Pooled_noContent_Prevalence_Time.png", 
    width = 1500, height = 1200, res = 220)
ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Type I Error for Single Response Simulation Under Various Thresholds",
       subtitle = "Time Covariate Included, No Batch Effect and Stromal Cells Only",
       x = "Threshold",
       y = "Type I Error") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),  # Increase title text size
    axis.title.x = element_text(size = 14),  # Increase X axis title size
    axis.title.y = element_text(size = 14),  # Increase Y axis title size
    axis.text.x = element_text(size = 14),   # Increase X axis labels size
    axis.text.y = element_text(size = 14)   # Increase Y axis labels size
  )
dev.off()
################## power #############################
alt.dat <- dat[grep("HighVar",names(dat))]
# check if any wts are NA
check.na <- lapply(alt.dat, function(x){
  y <- is.na(x[,1])
  return(y)
})
check.na <- do.call(rbind, check.na)
table(check.na)

manova.pValue <- lapply(alt.dat, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
manova.pValue <- do.call(rbind, manova.pValue)
manova.pValue$EffectSize <- as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", rownames(manova.pValue)))
manova.pValue$nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", rownames(manova.pValue)))
# threshold <- c(0.01, 0.05, 0.1)
threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = manova.pValue, by_group = "EffectSize")
topic_plot <- topic_plot %>%
  rename(WTS = wts_sig_proportion,
         MATS = mats_sig_proportion) %>%
  pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "Power")

png("res/1MultiSample_SingleResponse_Simulation/Power_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png", 
    width = 2500, height = 1800, res = 220)
ggplot(topic_plot, aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize, 
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Response Simulation with No Batch Effect and Stromal Cells only",
       subtitle = "Time Variable Included with Pooled Gamma Model",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
dev.off()
