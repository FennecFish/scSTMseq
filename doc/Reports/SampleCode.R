# Goal: This script is to generate multiple patient simulation
# both pre and post timepoints proportion is generated from a Logistic Normal
# Gamma is drawn from a multivariate normal
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(SingleCellExperiment)
library(MASS)
library(VariantAnnotation)
library(checkmate)
library(MCMCpack)

###############################################################################
########################### Simulation ########################################
###############################################################################
# extract seed from the cluster ID and add time for more randomness
args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))
set.seed(seed)

# As a test, one can set
seed = 1
numCellType <- 5
gamma_sd <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
nSample <- 20
nTimepoints <- 2

A = 3 
################################################################################
############################## Useful Functions ################################
################################################################################

generate_theta <- function(nSample, nTimepoints = 2, nCellType, mean, sd, num_Cov = 1){
  Simplex <- nCellType - 1
  sample_ids <- rep(1:nSample, each = nTimepoints) 
  Timepoint <- rep(c(0, 1), nSample)
  Response <- rep(c(0, 1), each = (nSample * nTimepoints) / 2)
  # interaction Term
  Timepoint_Response <- Timepoint * Response
  
  # Create design matrix for each sample/time combination
  X <- data.frame(Timepoint = Timepoint,
                  Response = Response,
                  Timepoint_Response = Timepoint_Response)
  rownames(X) <- paste(paste0("Sample", sample_ids), ifelse(Response == 0, "nonResponse", "Response"),
                       ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  
  
  # Generate Effect Size Matrix
  gamma <- vector(mode = "list")
  for (i in 1:Simplex) {
    gamma[[i]] <- diag(mvrnorm(n = num_Cov, mu = mean[,i], Sigma = sd[[i]]))
  }
  gamma <- do.call(cbind, gamma)
  # gamma <- rbind(1, gamma)
  colnames(gamma) <- c(paste0("K", 1:ncol(gamma)))
  rownames(gamma) <- c("Timepoint", "Response", "Timepoint:Response")
  
  mu <- t(t(as.matrix(gamma))  %*% t(as.matrix(X)))
  
  # simulate patient random effect
  psi <- mvrnorm(n = nSample, mu = rep(0, Simplex),  Sigma= diag(rep(1, Simplex)))
  
  # adding random effect eta = mu + psi
  for (i in 1:nSample) {
    idx_t1 <- 2 * (i - 1) + 1
    idx_t2 <- 2 * i
    
    mu[idx_t1, ] <- mu[idx_t1, ] + psi[i, ]
    mu[idx_t2, ] <- mu[idx_t2, ] + psi[i, ]
  }
  
  # transfer the proportion into simplex.
  # adding 0 to every sample/time for identifiability in the model
  eta <- cbind(mu, 0)
  colnames(eta) <- paste0("K", 1:ncol(eta))
  
  map_to_simplx <- function(x) {
    exp(x - log(sum(exp(x))))
  }
  
  theta <- t(apply(eta, 1, map_to_simplx))
  theta <- list(t1 = theta[grep("t1", rownames(theta)), ], t2 = theta[grep("t2", rownames(theta)), ])
  # return(list(theta = theta, gamma = gamma, psi = psi))
  return(list(theta = theta, gamma = gamma))
}

param_dat <- expand.grid(numCellType, gamma_sd)
colnames(param_dat) <- c("numCellType", "gamma_sd")

# #############################################################
true_param <- vector(mode = "list")
for (i in 1:nrow(param_dat)){
  nCellType <- param_dat[i,1]
  interaction_effect <- param_dat[i,2]
  simplex <- nCellType- 1 
  mean = matrix(rep(0, simplex*A), nrow = A)
  sd <- replicate(simplex, diag(c(rep(0, A -1), interaction_effect), A), simplify = FALSE)
  true_param[[i]] <- generate_theta(nSample = nSample, nTimepoints = 2, nCellType = nCellType, mean = mean, sd = sd, num_Cov = A)
}

# For each parameter combination, we want to calculate the relative change in cell type proportion
# Then group the relative change by response group
proportion_change <- lapply(true_param, function(x){
  RC <- (x$theta$t2-x$theta$t1)/x$theta$t1 %>% 
    as.data.frame()
  RC$Response =str_extract(rownames(RC), "(?<=_)[^_]+(?=_)")
  RC <- RC %>%
    group_by(Response) %>%
    summarise(across(starts_with("K"), ~mean(.x, na.rm = TRUE)))
})
names(proportion_change) <- gamma_sd
###############################################################################
############################### Plot ##########################################
###############################################################################
# row bind the list and then pivot to long format
proportion_change <- bind_rows(proportion_change, .id = "gamma") %>%
  mutate(gamma = as.numeric(gamma)) %>%
  pivot_longer(cols = starts_with("K"), names_to = "K", values_to = "value")

ggplot(proportion_change, aes(x = K, y = value, color = Response, group = Response)) +
  geom_line() +
  facet_wrap(~ gamma, scales = "free_y") +  # Facet by each K variable
  labs(x = "Gamma", y = "Value", color = "Response") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12), legend.position = "top")
 