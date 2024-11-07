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
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType5_noBatch_StromalCell/1000ManovaTheta_Pooled_noContent_Prevalence_Time/"
files <- list.files(dir, pattern = "Manova_pValue")
res <- vector(mode = "list")

for(file in files){
  set_level <- sub("^Manova_pValue_(.*)\\.rds$", "\\1", file)
  ThetaManova <- readRDS(paste0(dir, "Manova_pValue_", set_level, ".rds"))
  rm.wts <- lapply(ThetaManova$manova.fit, function(x){
    # return(list(statistics = x$MATS, p.value = x$resampling[2]))
    dat <- x$resampling[1]
    return(dat)
  })

  rm.mats <- lapply(ThetaManova$manova.fit, function(x){
    dat <- x$resampling[2]
    return(dat)
  })

  dat.clean <- function(given_list){
    given_list <- do.call(rbind, given_list)
    rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
    colnames(given_list) <- "pValue"
    return(given_list)
  }
  temp <- data.frame(pValue.wts = dat.clean(rm.wts),
                     pValue.mats = dat.clean(rm.mats))
  colnames(temp) <- c("pValue.wts", "pValue.mats")
  res[[set_level]] <- temp
  cat(file, "\n")
}
name <- basename(dir)
saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/",name,".rds"))

#################################################################################
################################## plot #########################################
#################################################################################
# plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
#   exist.change.list <- lapply(thresholds, function(thr) {
#     res <- pValueList %>%
#       mutate(wts.sig = wts < thr,
#              mats.sig = mats < thr)
#     return(res)
#   })
#   ProportionOfTrue <- data.frame()
# 
#   for (i in seq_along(thresholds)) {
#     threshold <- thresholds[i]
#     if(is.null(by_group)){
#       prop_truth <- exist.change.list[[i]] %>%
#         as.data.frame() %>%
#         summarize(
#           wts_sig_proportion = mean(wts.sig),
#           mats_sig_proportion = mean(mats.sig)
#         ) %>%
#         mutate(Threshold = threshold)
#     }else{
#       # if we want to group the power/ type I error by effect size
#       prop_truth <- exist.change.list[[i]] %>%
#         as.data.frame() %>%
#         group_by(EffectSize) %>%
#         summarize(
#           wts_sig_proportion = mean(wts.sig),
#           mats_sig_proportion = mean(mats.sig)
#         ) %>%
#         mutate(Threshold = threshold)
#     }
#     # Store the results in a dataframe
#     ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
#   }
#   return(ProportionOfTrue)
# }
# 
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/ManovaRM_ManovaTheta_Pooled_noContent_Prevalence_Time.rds")
# null.dat <- dat[grep("NullModel",names(dat))]
# 
# null.dat.dist <- do.call(rbind, null.dat)
# manova_long <- data.frame(
#   Value = c(null.dat.dist[,1], null.dat.dist[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(null.dat.dist))
# )
# ggplot(manova_long, aes(x = Value, fill = Type)) +
#   geom_histogram(bins = 10, color = "black") +
#   labs(title = "Distribution of WTS and MATS pValues Including All Replicates",
#        x = "p-value",
#        y = "Frequency") +
#   facet_wrap(~ Type) +
#   theme_minimal()
# 
# manova.pValue <- lapply(null.dat, function(x) {
#   data.frame(wts =quantile(x$pValue.wts, probs = 0.50), 
#              mats = quantile(x$pValue.mats, probs = 0.50))
# })
# manova.pValue <- do.call(rbind, manova.pValue)
# 
# manova_long <- data.frame(
#   Value = c(manova.pValue[,1], manova.pValue[,2]),
#   Type = rep(c("WTS", "MATS"), each = nrow(manova.pValue))
# )
# 
# ### first check the distirbution of pvalue (should be uniform) ##
# ggplot(manova_long, aes(x = Value, fill = Type)) +
#   geom_histogram(bins = 10, color = "black") +
#   labs(title = "Distribution of WTS and MATS pValues by Taking the Median",
#        x = "Medium p-value",
#        y = "Frequency") +
#   facet_wrap(~ Type) +
#   theme_minimal()
# 
# 
# threshold <- c(0.01, 0.05, 0.1)
# topic_plot <- plot_data_generate(threshold, pValueList = manova.pValue)
# topic_plot <- topic_plot %>% 
#   rename(WTS = wts_sig_proportion,
#          MATS = mats_sig_proportion) %>%
#   pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "TypeIError")
# 
# ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
#   geom_point(size = 3) +
#   geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
#   labs(title = "Type I Error for Single Response Simulation with No Batch Effect and Stromal Cells only",
#        subtitle = "Time Variable Included with Pooled Gamma Model",
#        x = "Threshold",
#        y = "Type I Error") +
#   theme_bw() +
#   scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
#   theme(
#     plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
#     axis.title.x = element_text(size = 16),  # Increase X axis title size
#     axis.title.y = element_text(size = 16),  # Increase Y axis title size
#     axis.text.x = element_text(size = 14),   # Increase X axis labels size
#     axis.text.y = element_text(size = 14)   # Increase Y axis labels size
#   )
# 
# ################## power #############################
# alt.dat <- dat[grep("HighVar",names(dat))]
# manova.pValue <- lapply(alt.dat, function(x) {
#   data.frame(wts =quantile(x$pValue.wts, probs = 0.50), 
#              mats = quantile(x$pValue.mats, probs = 0.50))
# })
# manova.pValue <- do.call(rbind, manova.pValue)
# manova.pValue$EffectSize <- as.numeric(gsub(".*HighVar", "", rownames(manova.pValue)))
# threshold <- c(0.01, 0.05, 0.1)
# topic_plot <- plot_data_generate(threshold, pValueList = manova.pValue, by_group = "EffectSize")
# topic_plot <- topic_plot %>% 
#   rename(WTS = wts_sig_proportion,
#          MATS = mats_sig_proportion) %>%
#   pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "Power")
# 
# ggplot(topic_plot, aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
#   geom_point(size = 3) +
#   facet_grid(~EffectSize, labeller = labeller(EffectSize = function(x) paste("Effect Size =", x))) + 
#   labs(title = "Power for Single Response Simulation with No Batch Effect and Stromal Cells only",
#        subtitle = "Time Variable Included with Pooled Gamma Model",
#        x = "Thresholds",
#        y = "Power") +
#   theme_bw() +
#   scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
#   theme(
#     plot.title = element_text(size = 18, hjust = 0.5),  # Increase title text size
#     axis.title.x = element_text(size = 16),  # Increase X axis title size
#     axis.title.y = element_text(size = 16),  # Increase Y axis title size
#     axis.text.x = element_text(size = 14, angle = 30, hjust = 1),   
#     axis.text.y = element_text(size = 14)   # Increase Y axis labels size
#   )
