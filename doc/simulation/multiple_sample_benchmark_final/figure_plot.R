setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)
library(tidyr)
########### adjusted Rand Index Plot for Single Patient Simulation ###############
dat <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final.csv") %>%
  dplyr::select(-matches("cluster$"), matches("cluster$"))

dat_long <- dat %>% 
  rename("scSTMseq" = "scSTM_C_Batch_P_cluster") %>%
  rename_with(~ gsub("_cluster", "", .), .cols = ends_with("_cluster")) %>%
  mutate(nCellType = recode(nCellType, "c5" = "Num_of_CellGroup = 5", "c9" = "Num_of_CellGroup = 9", 
                            "c13" = "Num_of_CellGroup = 13",  "c17" = "Num_of_CellGroup = 17")) %>%
  gather(method, adjR, scSTMseq:monocle3, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(nCellType, nSample, level, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  mutate(line_type = ifelse(grepl("scSTM", method), "dashed", "solid"))
  
dat_long$nCellType <- factor(dat_long$nCellType, levels = c("Num_of_CellGroup = 5", "Num_of_CellGroup = 9", "Num_of_CellGroup = 13", "Num_of_CellGroup = 17"))
dat_long$nSample <- factor(dat_long$nSample, levels = c("nsample6", "nsample12", "nsample24", "nsample36", "nsample48"))

dat_long <- dat_long %>% dplyr::filter(nSample %in% c("nsample6", "nsample12"))

png("res/sim_cluster_benchmark/adjustedRandIndex_lineplot_multiple_sample_benchmark_final_figure.png", width = 1500, height = 800, res = 150)
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method)) +
  geom_line(aes(linetype = line_type), linewidth = 1.2) +  # Increased line width
  facet_grid(nSample ~ nCellType) +
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
dev.off()

png("res/sim_cluster_benchmark/adjustedRandIndex_boxplot_multiple_sample_benchmark_final.png", width = 1500, height = 800, res = 160)
ggplot(dat_long, aes(x = method, y = mean_adjR)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # Remove outliers from the boxplot, adjust box width
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = method)) +  # Add jittered points for data distribution
  facet_grid(nSample ~ nCellType) +  # Use facet_grid for a more structured layout
  theme_minimal(base_size = 14) +  # Set a base font size for better readability
  labs(
    title = "Boxplot of Clustering Accuracy by Methods with Splatter Based Simulations",
    x = "Method",
    y = "Mean Adjusted Rand Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),  # Adjust x-axis text for better readability
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 12),  # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()
