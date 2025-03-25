setwd("/proj/milovelab/wu/scLDAseq")
library(monocle3)
library(DirichletReg)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/"
files <- list.files(path = paste0(dir, "sims/"), pattern = "Null")

file_name <- files[1]
set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)

file_path <- paste0(dir, "sims/sims_", set_level, ".rds")
sims <- safe_readRDS(file_path)
if (is.null(sims)) {
  next  # Skip to the next iteration if reading the file failed
}
dat <- colData(sims) %>% data.frame() 

obj <- readRDS(paste0(dir, "monocle3/monocle3_", set_level, ".rds"))
dat$monocle3_cluster <- partitions(obj)[match(rownames(dat), names(partitions(obj)))]

res <- table(dat$monocle3_cluster, dat$Time) %>% as.data.frame()
colnames(res) <- c("Group", "Time", "Freq")

res <- dat %>% 
  group_by(Time, monocle3_cluster) %>%
  summarise(count = n()) %>%                     
  mutate(total_by_time = sum(count),              
         proportion = count / total_by_time) %>%  
  ungroup() %>%
  select(Time, monocle3_cluster, proportion) %>%
  pivot_wider(names_from = monocle3_cluster, values_from = proportion)

res$AL <- DR_data(res[,2:6])
mod1 <- DirichReg(AL ~ Time, res)
