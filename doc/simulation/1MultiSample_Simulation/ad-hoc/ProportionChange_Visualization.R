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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
design <- "nSample20_nCellType5_noBatch_StromalCell/"
scSTM <- "scSTM_Pooled_noContent_Prevalence_Time/"
path = paste0(dir, design, scSTM)
files <- list.files(path)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

scSTMobj <- readRDS(paste0(path, files[1]))
model <- select_top_scSTM(scSTMobj)
topics <- model$settings$dim$K

SimEta <- etaPosterior(model, nsims = 100)
metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, Sample, Time, Response)

collaposed.diff <- lapply(SimEta, function(x){
  x <- cbind(x,0)
  theta.old <- exp(x - log(rowSums(exp(x))))
  rownames(theta.old) <- model$DocName
  theta.old <- theta.old %>%
    as.data.frame() %>%
    rownames_to_column("Cell") %>%
    left_join(metadata, by = "Cell")

  theta.collapsed <- theta.old %>%
    group_by(Sample, Time, Response) %>%
    summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
    ungroup()
  
  theta.new <- theta.collapsed %>%
    rowwise() %>%
    mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
    ungroup()
  
  # take difference between timepoints
  diff <- theta.new %>%
    group_by(Sample, Response) %>%
    summarise(across(starts_with("V"), ~ .[Time == "Time2"] - .[Time == "Time1"]), .groups = 'drop') %>%
    group_by(Response) %>%
    summarise(across(starts_with("V"), mean), .groups = 'drop')
  return(diff)
})

collaposed.diff <- do.call(rbind, collaposed.diff)
ci.level <- 0.95
offset <- (1-ci.level)/2
#For each topic,  Find means and cis
means = list()
cis = list()
dat <- collaposed.diff %>% filter(Response == "t1") %>% dplyr::select(starts_with("V"))
for(i in 1:ncol(dat)){
  value <- unlist(dat[,i])
  means[[i]] <- mean(value)
  cis[[i]] = quantile(value, c(offset,1-offset))
}

##### plot ####
main = "Example Visualization for Proportion Change"
ylim <- c(0, topics+1)
xlim <- c(min(unlist(cis) , na.rm=T), max(unlist(cis) + 0.005, na.rm=T))
xlab <- "Estimated Topic Proportions"
ylab <- ""
plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main,
     xlab=xlab, ylab=ylab,yaxt="n")
#Add a line for the 0 x-axis
lines(c(0,0), c(0, topics+2), lty=2)
# #Create labels
# labels = createLabels(covariate=covariate, method="difference",
#                       cdata=cdata, ref=ref, alt=alt,model=model,n=n,
#                       topics=topics,custom.labels=custom.labels, frexw=frexw, verbose.labels=verbose.labels)
add = T
linecol = grDevices::rainbow(length(topics))
printlegend = F
#Plot everything
it <- topics
for(i in 1:topics){
  if(add){
    # Shift the line for better visualization
    shift_amount <- 0.1  
    shifted_y_values <- it - shift_amount
    points(means[[i]], shifted_y_values, pch = 16, col = linecol)
    lines(c(cis[[i]][1],cis[[i]][2]),c(shifted_y_values,shifted_y_values), col = linecol)
    # Add label for topic number on the right side of the line
    text(as.numeric(cis[[i]][2]), shifted_y_values, labels = paste("Topic", i), pos = 4, col = linecol)
  }else{
    points(means[[i]], it, pch=16, col = linecol)
    lines(c(cis[[i]][1],cis[[i]][2]),c(it,it), col = linecol)
    # Add label for topic number on the right side of the line
    text(cis[[i]][2], it, labels = paste("Topic", i), pos = 4, col = linecol)
  }
  it = it-1
}

# time_effect <- estimateEffect(1:topics ~ Time,
#                               stmobj = model, 
#                               uncertainty = "Global",
#                               nsims = 30)
# summary(time_effect)
# plot(time_effect, covariate = "Time", model = model,
#      method = "difference", ref="Time1",alt="Time2",
#      xlab = "proportion change between timepoints",
#      main = "Expected Cell Group Proportion Change Before and After Treatment",
#      verbose.labels=F)
