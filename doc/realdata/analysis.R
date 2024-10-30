# Goal: This script is used to analyze scSTM from real data PD1
# 1) compare assigned cluster from paper cluster
# 2) changes between groups (time/expansion)
setwd("/proj/milovelab/wu/scLDAseq")
# setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)

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

cluster_scSTM <- function(scSTMobj) {
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- scSTMobj$DocName
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

# scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")
scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_BatchNoInteraction_Prevalence_TimeResponseInteraction.rds")
scSTMobj <- select_top_scSTM(scSTMobj)
sims <- scSTMobj$settings$sce
# calculate adjusted rand index
scSTM_cluster <- cluster_scSTM(scSTMobj)
scSTM_cluster <- scSTM_cluster[match(rownames(colData(sims)), names(scSTM_cluster))]
adjustedRandIndex(scSTM_cluster, sims$cellType) # 0.6372073 for with batch interaction, 0.7622637 for without batch interaction

K <- length(unique(sims$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

png("res/PD1/structure_plot_BatchNoInteraction.png", height = 800, width = 1200, res = 250)
structure_plot(scSTMobj, topics = 1:K, grouping = sims$cellType, n = 2000, gap = 20)
dev.off()

# time effect
ref <- c("Pre", "NE")
names(ref) <- c("timepoint", "expansion")
time_effect <- estimateEffect(1:K ~ timepoint*expansion,
                              stmobj = scSTMobj, ref.vec = ref,
                              uncertainty = "Global",
                              nsims = 30)
saveRDS(time_effect, file = "res/PD1/effect_BatchNoInteraction_nsim30.rds")

# 
# # time_effect <- readRDS("res/PD1/effect_BatchNoInteraction_nsim30.rds")
summary(time_effect)
# plot.estimateEffect(time_effect, covariate = "timepoint", model=scSTMobj,
#                     method="difference",cov.value1="Pre",cov.value2="On", linecol = "black")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

png("res/PD1/estimateEffect_rep30.png", width = 1200, height = 800, res = 150)
plot(time_effect, covariate = "timepoint", model = scSTMobj,
     #method = "difference",
     method = "difference", # cov.value1="Pre",cov.value2="On",
     xlab = "proportion change between timepoints", moderator = "expansion",
     moderator.value = "E", linecol = "skyblue", printlegend = FALSE,
     main = "Expected Cell Group Proportion Change Before and After Treatment")

plot(time_effect, covariate = "timepoint", model = scSTMobj,
     #method = "difference",
     method = "difference", # cov.value1="Pre",cov.value2="On",
     xlab = "proportion change between timepoints", moderator = "expansion",
     moderator.value = "NE", linecol = "orange", add = T,
     printlegend = T, labeltype = "custom", custom.labels = paste0("Topic ", 1:8))
legend("bottomright", legend = c("expansion", "non-expansion"), lwd = 2, col = c("skyblue", "orange"))
dev.off()

dat <- as.data.frame(t(labelTopics(scSTMobj, topics = 1:8, n = 50)$topics))
colnames(dat) <- paste0("Topic_", 1:8)
write.csv(dat, file = "res/PD1/top_50genes_in_topics.csv", row.names = F)

library(xtable)
latex_code <- xtable(dat)

labelTopics(scSTMobj, topics = 7, n = 5)$topics
# # gene_list <- labelTopics(scSTMobj, topics = 5, n = 200)
# # gene_list <- as.vector(gene_list$topics)
# #
# # labelTopics(scSTMobj, topics = 1:8, n = 5)$topics
# # gene_universe <- labelTopics(scSTMobj, topics = 5, n = length(scSTMobj$vocab))
# # gene_universe <- as.vector(gene_universe$topics)
# #
# # # gene_list <- mapIds(org.Hs.eg.db,
# # #                      keys = gene_list,
# # #                      column = "ENTREZID",
# # #                      keytype = "SYMBOL",
# # #                      multiVals = "first")
# # library(clusterProfiler)
# # ggo <- enrichGO(gene          = gene_list,
# #                 universe      = gene_universe,
# #                 OrgDb         = org.Hs.eg.db,
# #                 keyType       = "SYMBOL",
# #                 ont           = "BP",
# #                 pAdjustMethod = "BH",
# #                 pvalueCutoff  = 0.05,
# #                 qvalueCutoff  = 1,
# #                 readable      = TRUE)
# #
# # goplot(ggo)
