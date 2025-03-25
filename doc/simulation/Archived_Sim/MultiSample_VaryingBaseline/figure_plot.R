setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)
library(tidyr)

dat <- read.csv("res/adjRandIndex_MultiSample_VaryingBaseline_Batch0.5_StromalCell.csv") %>%
  dplyr::select(-matches("cluster$"), matches("cluster$"))
dat$effectSize <- factor(dat$effectSize, levels = paste0("effectSize", c(1, 0.7, 0.4, 0.1)))
dat_long <- dat %>% 
  # rename("scSTMseq" = "scSTM_C_Batch_P_cluster") %>%
  rename_with(~ gsub("_cluster", "", .), .cols = ends_with("_cluster")) %>%
  gather(method, adjR, scSTM_noContent_Sample:monocle3, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(effectSize, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  mutate(line_type = ifelse(grepl("scSTM", method), "dashed", "solid"))

# png("res/sim_cluster_benchmark/adjustedRandIndex_lineplot_multiple_sample_benchmark_final_figure.png", width = 1500, height = 800, res = 150)
ggplot(dat_long, aes(x = effectSize, y = mean_adjR, color = method, group = method)) +
  geom_line(aes(linetype = line_type), linewidth = 1.2) +  # Increased line width
  # facet_grid(nSample ~ nCellType) +
  labs(
    x = "Noise Level",
    y = "Mean Adjusted Rand Index",
    title = "Benchmarking Clustering Accuracy Using Mean Adjusted Rand Index",
    subtitle = "Across Different Sample Sizes and Cell Typeswith Splatter Based Simulations"
  ) +
  theme_minimal(base_size = 12) +  # Base size adjusted for consistency
  theme_bw() + 
  theme(
    text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),  # Subtitle styling
    # plot.caption = element_text(size = 10, face = "italic", hjust = 0),  # Caption styling
    legend.position = "bottom",  # Move legend to the bottom for better readability
    legend.box = "horizontal",  # Arrange legend items horizontally
    panel.grid.major = element_line(color = "grey80"),  # Subtle major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    strip.text = element_text(size = 12, face = "bold")  # Facet label styling
  ) +
  guides(
    color = guide_legend(title = "Method"),  # Keep only color legend
    linetype = "none"  # Remove linetype from the legend
  )
# dev.off()

dat <- read.csv("res/adjRandIndex_MultiSample_VaryingBaseline_nSample12_nCellType8_Batch_CancerCell.csv") %>%
  dplyr::select(-matches("cluster$"), matches("cluster$"))
dat$effectSize <- factor(dat$effectSize, levels = paste0("effectSize", c(1, 0.7, 0.4, 0.1)))
dat_long <- dat %>% 
  dplyr::select(-c("scSTM_L1_cluster", "scSTM_LR_cluster")) %>%
  # dplyr::select(effectSize, matches("^scSTM")) %>%
  # rename("scSTMseq" = "scSTM_C_Batch_P_cluster") %>%
  rename_with(~ gsub("_cluster", "", .), .cols = ends_with("_cluster")) %>%
  # gather(method, adjR, scSTM_pooled:monocle3, factor_key=TRUE) %>%
  gather(method, adjR, scSTM_pooled:monocle3, factor_key=TRUE) %>%
  dplyr::filter(!is.na(adjR))
  
  # group_by(effectSize, method) %>%
  # summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  # mutate(line_type = ifelse(grepl("scSTM", method), "dashed", "solid"))

# png("res/sim_cluster_benchmark/adjustedRandIndex_boxplot_multiple_sample_benchmark_final.png", width = 1500, height = 800, res = 160)
ggplot(dat_long, aes(x = method, y = adjR, fill = method)) +
  geom_boxplot(width = 0.7) +  # Remove outliers from the boxplot, adjust box width
  # geom_jitter(width = 0.2, size = 1, alpha = 0.6, aes(color = method)) +  # Add jittered points for data distribution
  facet_grid(~effectSize) +  
  theme_bw() +
  theme_minimal(base_size = 14) +  # Set a base font size for better readability
  labs(
    title = "Boxplot of Clustering Accuracy Comparing scSTM with different Prevalence Model",
    subtitle = "Batch Effect with One Cancer Cell Group, with total of 12 Samples and 8 Cell Types",
    x = "Method",
    y = "Adjusted Rand Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14),  # Center and bold the title
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Adjust x-axis text for better readability
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  
  )
# dev.off()
