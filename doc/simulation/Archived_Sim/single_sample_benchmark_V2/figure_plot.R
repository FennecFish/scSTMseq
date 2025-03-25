setwd("/proj/milovelab/wu/scLDAseq")
library(ggplot2)
library(dplyr)

########### adjusted Rand Index Plot for Single Patient Simulation ###############
dat <- read.csv("res/single_patient_simulation_V2/adjRandIndex_single_sample_benchmark_V2.csv")
# dat <- read.csv("res/sim_cluster_benchmark/adjRandIndex_single_sample_multiCellType_benchmark_V2.csv")

dat_long <- dat %>% 
  dplyr::filter(level != "L5") %>%
  dplyr::select(sim, seed, control, level, nCellType, nSample, monocle3_cluster, scSTM_C_P_cluster, seurat_cluster, fastTopics_cluster) %>%
  rename( "scSTM_C_P_cluster" = "scSTMseq") %>%
  rename_with(~ gsub("_cluster", "", .), .cols = ends_with("_cluster")) %>%
  mutate(nCellType = recode(nCellType, "c5" = "Num_of_CellGroup = 5", "c9" = "Num_of_CellGroup = 9", 
                            "c13" = "Num_of_CellGroup = 13")) %>%
  gather(method, adjR, monocle3:fastTopics, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(nCellType, nSample, level, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
  mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))


# dat_long <- dat %>% 
#   dplyr::select(-sctransform_cluster) %>%
#   gather(method, adjR, scSTM_nC_P_cluster:monocle3_cluster, factor_key=TRUE) %>%
#   drop_na()  %>%
#   group_by(nCellType, level, control, method) %>%
#   summarise(mean_adjR = mean(adjR, na.rm = TRUE)) %>%
#   mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))
#   
dat_long$nCellType <- factor(dat_long$nCellType, levels = c("Num_of_CellGroup = 5", "Num_of_CellGroup = 9", "Num_of_CellGroup = 13"))

png("res/single_patient_simulation_V2/adjustedRandIndex_single_sample_benchmark_V2_figure.png", width = 1200, height = 900, res = 200)
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method)) +
  geom_line(aes(linetype = line_type)) +
  facet_wrap(~ nCellType) +
  labs(
    x = "Noise Level",
    y = "Mean Adjusted Rand Index",
    title = "Benchmarking Clustering Accuracy Using Mean Adjusted Rand Index",
    subtitle = "Single Sample Simulation"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic"), 
    legend.position = "bottom",  # Move legend to the bottom for better readability
    legend.box = "horizontal",  # Arrange legend items horizontally
    panel.grid.major = element_line(color = "grey80"),  # Subtle major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    strip.text = element_text(size = 10, face = "bold")  # Facet label styling
  ) +
  guides(linetype = "none") # Remove linetype from the legend
dev.off()
# first demonstrate different version of scSTM
res_stats <- res_stats %>% #filter(grepl("^scSTM", method)) %>%
  mutate(line_type = ifelse(grepl("stm_cluster", method), "solid", "dashed")) %>%
  mutate(control = sapply(strsplit(level, "_"), `[`, 1),
         level = sapply(strsplit(level, "_"), `[`, 2))

# res_stats$control <- factor(res_stats$control)
png("res/figure1_scSTM_adjR.png", width = 2500, height = 1500, res = 300)
ggplot(dat_long, aes(x = level, y = mean_adjR, color = method, group = method, linetype = line_type)) +
  geom_line(linewidth = 1) +
  labs(title = "Mean Adjusted Rand Index by scSTM Methods and Level of Noise",
       x = "Level of Noise",
       y = "Mean adjRandIndex") +
  facet_grid(~control) +
  #  geom_errorbar(aes(ymin = mean_adjR - sd, ymax = mean_adjR + sd), width = 0.2) +
  theme_minimal() +
  # scale_color_brewer(palette = "Set5") +
  # scale_color_viridis_d() + 
  scale_color_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
  theme(axis.text.x = element_text(hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis text
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase size of x-axis title
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase size of y-axis title
        plot.title = element_text(size = 12, face = "bold"),  # Increase size and bold the plot title
        strip.text = element_text(size = 12, face = "bold"))
dev.off()
