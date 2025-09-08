#' Run scSTMseq Simulation
#'
#' This function generates synthetic single-cell RNA sequencing data with structured topic modeling features,
#' including support for cell-level covariates and sample-level random variation.
#'
#' @param output_dir Directory to save simulated datasets.
#' @param seed Random seed for reproducibility.
#' @param nSample Number of samples (subjects).
#' @param nTimepoints Number of time points per subject.
#' @param nGenes Number of genes.
#' @param de.prob Probability of differential expression.
#' @param de.facLoc Log-fold change for differential expression.
#' @param numCellType_set Vector of cell type counts to simulate.
#' @param gamma_sd_set Vector of standard deviations for the gamma effect.
#' @param A Number of covariates.
#' @param cancerCellGroup Index of the cancer cell group (or NULL).
#' @param batch.rmEffect Logical; whether to remove batch effects.
#'
#' @return Saves RDS files containing simulated `SingleCellExperiment` objects to disk.
#' @export
#'
#' @importFrom scuttle logNormCounts
#' @importFrom splatter newSplatParams
#' @importFrom MASS mvrnorm
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix Matrix
#' @importFrom stats rnorm
#'
SimPairedData <- function(
    seed = 1234,
    nSample = 10,
    nTimepoints = 2,
    nGenes = 3000,
    nCell = c(250),
    de.prob = 0.3,
    de.facLoc = 0.5,
    numCellType_set = c(10),
    gamma_sd_set = c(0, 0.3),
    A = 1,
    cancerCellGroup = 2,
    batch.rmEffect = FALSE
) {
  batchCells <- rep(c(nCell, nCell), each = nSample)
  batch.facLoc <- runif(nSample, min = 0, max = 0.5)
  batch.facLoc <- rep(batch.facLoc, times = 2)

  save_batch <- ifelse(batch.rmEffect, "noBatch", "Batch")
  save_cancer <- ifelse(is.null(cancerCellGroup), "StromalCell", "CancerCell")

  param_dat <- expand.grid(numCellType_set, gamma_sd_set)
  colnames(param_dat) <- c("numCellType", "gamma_sd")

  for (i in 1:nrow(param_dat)) {
    nCellType <- param_dat[i, "numCellType"]
    gamma_sd_tmp <- param_dat[i, "gamma_sd"]
    simplex <- nCellType - 1
    mean <- matrix(rep(0, simplex), nrow = 1)
    sd <- replicate(simplex, diag(gamma_sd_tmp, A), simplify = FALSE)
    type <- if (gamma_sd_tmp == 0) "NullModel" else paste0("HighVar", gamma_sd_tmp)

    true_param <- generate_theta(
      nSample = nSample, nTimepoints = nTimepoints,
      nCellType = nCellType, mean = mean, sd = sd
    )

    sims <- delta_sim(
      true_param = true_param, seed = seed, nSample = nSample,
      nGenes = nGenes, nCellType = nCellType,
      de.prob = de.prob, de.facLoc = de.facLoc,
      batchCells = batchCells, batch.facLoc = batch.facLoc,
      cancerCellGroup = cancerCellGroup, batch.rmEffect = batch.rmEffect
    )

    dir_path <- file.path(output_dir, paste0("nSample", nSample,
                                             "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "sims"))
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

    saveRDS(sims, file = file.path(dir_path, paste0("sims_", seed, "_", type, ".rds")))
    rm(sims)

    message("Generated simulation with ", nCellType, " cell types and gamma SD = ", gamma_sd_tmp)
  }
}

#' @keywords internal
generate_theta <- function(nSample, nTimepoints = 2, nCellType, mean, sd){
  Simplex <- nCellType - 1
  Timepoint <- rep(c(0, 1), nSample)
  X <- data.frame(Timepoint = Timepoint)

  if(nSample > 1){
    sample_ids <- rep(1:nSample, each = nTimepoints)
    rownames(X) <- paste(paste0("Sample", sample_ids), ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  } else {
    rownames(X) <- paste(ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  }

  map_to_simplx <- function(x) {
    exp(x - log(sum(exp(x))))
  }

  gamma <- vector(mode = "list")
  for (i in 1:Simplex) {
    gamma[[i]] <- rnorm(n = 1, mean = mean[i], sd = sd[[i]])
  }
  gamma <- do.call(cbind, gamma)
  colnames(gamma) <- paste0("K", 1:ncol(gamma))
  rownames(gamma) <- c("Timepoint")

  mu <- t(t(as.matrix(gamma)) %*% t(as.matrix(X)))

  if(nSample > 1){
    psi <- MASS::mvrnorm(n = nSample, mu = rep(0, Simplex), Sigma = diag(rep(1, Simplex)))
    for (i in 1:nSample) {
      mu[paste0("Sample", i, "_t1"), ] <- mu[paste0("Sample", i, "_t1"), ] + psi[i, ]
      mu[paste0("Sample", i, "_t2"), ] <- mu[paste0("Sample", i, "_t2"), ] + psi[i, ]
    }
  } else {
    psi <- NULL
  }

  eta <- cbind(mu, 0)
  colnames(eta) <- paste0("K", 1:ncol(eta))

  theta <- t(apply(eta, 1, map_to_simplx))
  theta <- list(t1 = theta[grep("t1", rownames(theta)), ], t2 = theta[grep("t2", rownames(theta)), ])

  return(list(theta = theta, psi = psi))
}

#' @keywords internal
delta_sim <- function(true_param, seed, nSample, nGenes, nCellType,
                      de.prob, de.facLoc, batchCells, batch.facLoc,
                      cancerCellGroup = NULL, batch.rmEffect = TRUE) {
  params <- newSTMParams()
  params <- setParams(params, nGenes = nGenes,
                      group.prob = true_param$theta,
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      batchCells = batchCells, batch.facLoc = batch.facLoc,
                      seed = seed)
  sims <- scSTMseqSimulate(params, method = "groups",
                        verbose = FALSE, batch.rmEffect = batch.rmEffect,
                        cancerCellGroup = cancerCellGroup,
                        true_param = true_param)
  return(sims)
}

