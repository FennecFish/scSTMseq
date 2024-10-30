# This script is to apply manova with repeated measurement
library(MASS)
library(MANOVA.RM)
# library(CompDTUReg)
library(stats)
library(dplyr)
library(SingleCellExperiment)
library(tibble)
# library(MCMCglmm)
# This function is to calculate pValue directly copied from compDTU
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, lm_model_fit = NULL, q, nsamp, df_residual = NA){
  if(is.null(lm_model_fit) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in lm_model_fit to extract df.residual from")
  }
  
  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)
  # Etilde <- SigmaTildeAlt
  # Htilde <- SigmaTildeNull - SigmaTildeAlt
  
  vals <- diag((Htilde %*%solve(Etilde + Htilde)))
  pill_stat <- sum(vals)
  
  
  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas
  #This proved to be the easiest way to calculate the statistic and are confirmed to match R's anova.mlm function
  #v <- nsamp - ncol(stats::model.matrix(lm_model_fit))
  
  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- lm_model_fit$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit
  
  
  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)
  
  s <- min(p,q)
  
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)
  
  #Formulas come from the SAS help file Multivariate ANOVA Testing and are confirmed to match R's anova.mlm function
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))
  
  if(fstat_pillai < 0){
    return(list(FStat = NA, NumDF = NA, DenomDF = NA, pval_pillai = NA))
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-stats::pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  return(list(FStat = fstat_pillai, NumDF = numdf_pillai, DenomDF = denomdf_pillai, pval_pillai = pval_pillai))
}

ThetaManova <- function(model, nsims = 100){
  
  #########################################
  ## Step 1: simulate eta 
  #########################################
  SimEta <- etaPosterior(model, nsims = nsims)
  metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, Sample, Time, Response)
  #########################################
  ## Step 2: get topic proportion 
  ## then collapse by time and sample
  ## renormalize theta to proportion
  #########################################
  collaposed.theta <- lapply(SimEta, function(x){
    x <- cbind(x,0)
    theta.old <- exp(x - log(rowSums(exp(x))))
    rownames(theta.old) <- model$DocName
    theta.old <- theta.old %>%
      as.data.frame() %>%
      rownames_to_column("Cell") %>%
      left_join(metadata, by = "Cell")
    
    theta.collapsed <- theta.old %>%
      group_by(Sample, Time) %>%
      summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
      ungroup()
    
    theta.new <- theta.collapsed %>%
      rowwise() %>%
      mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
      ungroup()
    
    return(theta.new)
  })
  
  
  # #########################################
  # ## Step 3a: Directly using CompDTU
  # #########################################
  # # calculate irl for each replicate
  # irl.InfRep <- lapply(1:length(collaposed.theta), function(i){
  #   y <- collaposed.theta[[i]] %>% dplyr::mutate(InfRep = paste0("infRep", i))
  #   Y <- y %>% 
  #     dplyr::select(starts_with("V"))
  #   Y <- compositions::ilr(Y)
  #   res <- y %>% 
  #     dplyr::select(-starts_with("V"))
  #   res <- cbind(res, Y)
  #   return(res)
  # })
  # irl.InfRep <- do.call(rbind, irl.InfRep)
  # 
  # # calculate within sample covariance for each sample
  # WithinSampleCov <- vector(mode = "list")
  # unique.sample <- unique(irl.InfRep$Sample)
  # Z <- stats::model.matrix(~unique.sample)
  # for (sample in unique.sample){
  #   irl.sample <- irl.InfRep %>% 
  #     dplyr::filter(Sample == sample) %>%
  #     dplyr::select(starts_with("V"))
  #   ilrCov <- stats::cov(irl.sample)
  #   WithinSampleCov[[sample]] <- ilrCov
  # }
  # # calculate mean across samples
  # WithinSampleCov <- Reduce(`+`, WithinSampleCov) / length(WithinSampleCov)
  # TimeEffect <- irl.InfRep$Time
  # ncond <- length(unique(TimeEffect))
  # nsamp <- length(unique.sample)
  # # the following code are adapted from CompDTU
  # 
  # # alternative model
  # XAlt <- stats::model.matrix(~TimeEffect)
  # YInfRep <- irl.InfRep %>% dplyr::select(starts_with("V")) %>% as.matrix()
  # ns <- nrow(YInfRep)
  # XAltT <- t(XAlt)
  # bhatalt <- solve(crossprod(XAlt)) %*% (XAltT %*% YInfRep)
  # pie1 <- YInfRep - (XAlt%*%bhatalt)
  # SigmaTildeAltNewModeling <- (crossprod(pie1))/ns
  # 
  # # null model
  # XNull <- stats::model.matrix(~1, data = irl.InfRep)
  # XNullT <- t(XNull)
  # bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% YInfRep)
  # pie2 <- YInfRep - (XNull%*%bhatnull)
  # SigmaTildeNullNewModeling <- (crossprod(pie2))/ns
  # 
  # # Between_Sample_covariance = covariance - Within Sample Covariance
  # UpdatedCovAlt <- SigmaTildeAltNewModeling - WithinSampleCov
  # UpdatedCovNull <- SigmaTildeNullNewModeling - WithinSampleCov
  # 
  # #If there is a negative variance term, the pvalue for CompDTUme will be undefined
  # statement1 <- sum(diag(UpdatedCovAlt)<=0) !=0
  # statement2 <- sum(diag(UpdatedCovNull)<=0) !=0
  # ret.mixed <- data.frame(pval_CompDTUme = NA, FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)
  # ret <- data.frame(pval_CompDTUme = NA, FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)
  # 
  # if((statement1==TRUE | statement2==TRUE)){
  #   ret.mixed$pval_CompDTUme <- NA
  #   ret$pval_CompDTUme <- NA
  # }else{
  #   qq <- ncond - 1
  #   CompDTUmeRes.mixed <- tryCatch(calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
  #                                           lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
  #                            error = function(x){})
  #   
  #   if(is.null(CompDTUmeRes.mixed)==FALSE){
  #     ret.mixed$pval_CompDTUme <- CompDTUmeRes.mixed$pval_pillai
  #     ret.mixed$FStat <- CompDTUmeRes.mixed$FStat
  #     ret.mixed$NumDF <- CompDTUmeRes.mixed$NumDF
  #     ret.mixed$DenomDF <- CompDTUmeRes.mixed$DenomDF
  #   }
  #   CompDTUmeRes <- tryCatch(calcPillaiPval(SigmaTildeNull = SigmaTildeNullNewModeling, SigmaTildeAlt = SigmaTildeAltNewModeling,
  #                                                 lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
  #                                  error = function(x){})
  #   if(is.null(CompDTUmeRes)==FALSE){
  #     ret$pval_CompDTUme <- CompDTUmeRes$pval_pillai
  #     ret$FStat <- CompDTUmeRes$FStat
  #     ret$NumDF <- CompDTUmeRes$NumDF
  #     ret$DenomDF <- CompDTUmeRes$DenomDF
  #   }
  # }
  
  #########################################
  ## Step 3b: MANOVA + Repeated Measure
  #########################################
  # fit a manova model
  manova.fit <- lapply(collaposed.theta, function(x){
    Y <- x %>% 
      dplyr::select(starts_with("V"))
    Y <- compositions::ilr(Y)
    x <- x %>%
      dplyr::select(-starts_with("V")) %>%
      dplyr::mutate(Time = as.factor(Time),
             Sample = as.factor(Sample))
    x <- cbind(x, Y)
    response_vars <- grep("^V", names(x), value = TRUE)
    fit.temp <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
                       data = x, subject = "Sample", within = "Time", iter = 1000)
  })
 

  # #########################################
  # ## Step 3c: multivariate mixed linear effect
  # #########################################
  # # fit a manova model
  # mcmcglmm.fit <- lapply(collaposed.theta, function(x){
  #   Y <- x %>% 
  #     dplyr::select(starts_with("V"))
  #   Y <- compositions::ilr(Y)
  #   x <- x %>%
  #     dplyr::select(-starts_with("V")) %>%
  #     dplyr::mutate(Time = as.factor(Time),
  #                   Sample = as.factor(Sample))
  #   x <- cbind(x, Y)
  #   response_vars <- grep("^V", names(x), value = TRUE)
  #   fixed_formula <- as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time "))
  #   ilr_dim <- ncol(Y)
  #   priors = list(R = list(V = diag(ilr_dim), nu = 0.002),
  #                 G = list(G1 = list(V = diag(ilr_dim), nu = 2,
  #                                    alpha.mu = rep(0, ilr_dim),
  #                                    alpha.V = diag(1, ilr_dim, ilr_dim))))
  #   mcmc.fit <- tryCatch(
  #     MCMCglmm(fixed_formula, random = ~us(trait):Sample, rcov = ~us(trait):units, data = x,
  #       family = rep("gaussian", ilr_dim), verbose = FALSE, prior = priors),
  #     error = function(e) NA)
  # })
  # 
  toreturn <- list(SimEta = SimEta, collaposed.theta = collaposed.theta, 
                   manova.fit = manova.fit)
                   # mcmcglmm.fit = mcmcglmm.fit)
                   # compDTU.mixed = ret.mixed,
                   # compDTU = ret)
  return(toreturn)
}

