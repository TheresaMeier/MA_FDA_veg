# The only change to the original: output all available eigenvalues, not just the first M ones

MFPCA_2 <- function(mFData, M, uniExpansions, weights = rep(1, length(mFData)), fit = FALSE, approx.eigen = FALSE,
                  bootstrap = FALSE, nBootstrap = NULL, bootstrapAlpha = 0.05, bootstrapStrat = NULL, 
                  verbose = options()$verbose)
{
  if(! inherits(mFData, "multiFunData"))
    stop("Parameter 'mFData' must be passed as a multiFunData object.")
  
  # number of components
  p <- length(mFData)
  # number of observations
  N <- nObs(mFData)
  
  if(!all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a number > 0.")
  
  if(!(is.list(uniExpansions) & length(uniExpansions) == p))
    stop("Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")
  
  if(!(is.numeric(weights) & length(weights) == p))
    stop("Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")
  
  if(!is.logical(fit))
    stop("Parameter 'fit' must be passed as a logical.")
  
  if(!is.logical(approx.eigen))
    stop("Parameter 'approx.eigen' must be passed as a logical.")
  
  if(!is.logical(bootstrap))
    stop("Parameter 'bootstrap' must be passed as a logical.")
  
  if(bootstrap)
  {
    if(is.null(nBootstrap))
      stop("Specify number of bootstrap iterations.")
    
    if(any(!(0 < bootstrapAlpha & bootstrapAlpha < 1)))
      stop("Significance level for bootstrap confidence bands must be in (0,1).")
    
    if(!is.null(bootstrapStrat))
    {
      if(!is.factor(bootstrapStrat))
        stop("bootstrapStrat must be either NULL or a factor.")
      
      if(length(bootstrapStrat) != nObs(mFData))
        stop("bootstrapStrat must have the same length as the number of observations in the mFData object.")
    }
  }
  
  if(!is.logical(verbose))
    stop("Parameter 'verbose' must be passed as a logical.")
  
  # dimension for each component
  dimSupp <- dimSupp(mFData)
  
  # get type of univariate expansions
  type <- vapply(uniExpansions, function(l){l$type}, FUN.VALUE = "")
  
  # de-mean functions -> coefficients are also de-meaned!
  # do not de-mean in uFPCA, as PACE gives a smooth estimate of the mean (see below)
  m <- meanFunction(mFData, na.rm = TRUE) # ignore NAs in data
  for(j in seq_len(p))
  { 
    if(type[j] != "uFPCA")
      mFData[[j]] <- mFData[[j]] - m[[j]]
  }
  
  if(verbose)
    cat("Calculating univariate basis expansions (", format(Sys.time(), "%T"), ")\n", sep = "")
  
  # calculate univariate basis expansion for all components
  uniBasis <- mapply(function(expansion, data){do.call(univDecomp, c(list(funDataObject = data), expansion))},
                     expansion = uniExpansions, data = mFData, SIMPLIFY = FALSE)
  
  # for uFPCA: replace estimated mean in m
  for(j in seq_len(p))
  {
    if(type[j] == "uFPCA")
      m[[j]] <- uniBasis[[j]]$meanFunction
  }
  
  # Multivariate FPCA
  npc <- vapply(uniBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions
  
  if(M > sum(npc))
  {
    M <- sum(npc)
    warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
  } 
  
  # check if non-orthonormal basis functions used
  if(all(foreach::foreach(j = seq_len(p), .combine = "c")%do%{uniBasis[[j]]$ortho}))
    Bchol = NULL
  else
  {
    # Cholesky decomposition of B = block diagonal of Cholesky decompositions
    Bchol <- Matrix::bdiag(lapply(uniBasis, function(l){
      if(l$ortho)
        res <- Matrix::Diagonal(n = ncol(l$scores))
      else
        res <- Matrix::chol(l$B)
      
      return(res)}))
  }
  
  if(verbose)
    cat("Calculating MFPCA (", format(Sys.time(), "%T"), ")\n", sep = "")
  
  mArgvals <- if (utils::packageVersion("funData") <= "1.2") {
    getArgvals(mFData)
  } else {
    funData::argvals(mFData)
  }
  
  res <- calcMFPCA_2(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                   npc = npc, argvals = mArgvals, uniBasis = uniBasis, fit = fit, approx.eigen = approx.eigen)
  
  res$meanFunction <- m # return mean function, too
  
  names(res$functions) <- names(mFData)
  
  if(fit)
  {
    res$fit <- m + res$fit # add mean function to fits
    names(res$fit) <- names(mFData)
  } 
  
  # give correct names
  namesList <- lapply(mFData, names)
  if(!all(vapply(namesList, FUN = is.null, FUN.VALUE = TRUE))) 
  {
    if(length(unique(namesList)) != 1)
      warning("Elements have different curve names. Use names of the first element for the results.")
    
    row.names(res$scores) <- namesList[[1]]
    
    if(fit)
      for(i in seq_len(p))
        names(res$fit[[i]]) <- namesList[[1]]
  }
  
  # bootstrap for eigenfunctions
  if(bootstrap)
  {
    if(verbose)
      cat("Bootstrapping results:\n")
    
    booteFuns <- vector("list", p)
    
    for(j in seq_len(p))
      booteFuns[[j]] <- array(NA, dim  = c(nBootstrap, M, vapply(mFData[[j]]@argvals, FUN = length, FUN.VALUE = 0)))
    
    booteVals <- matrix(NA, nrow = nBootstrap, ncol = M)
    
    for(n in seq_len(nBootstrap))
    {
      if(verbose)
      {
        if(n %% 10 == 0)
          cat("\t n = ", n, " (", format(Sys.time(), "%T"), ")\n", sep = "")
      }
      
      if(is.null(bootstrapStrat))
        bootObs <- sample(N, replace = TRUE)
      else
        bootObs <- stratSample(bootstrapStrat)
      
      bootBasis <- vector("list", p)
      
      for(j in seq_len(p))
      {
        if(!is.null(uniBasis[[j]]$functions)) # re-estimate scores AND functions
        {
          bootBasis[[j]] <- do.call(univDecomp, c(list(funDataObject = (mFData[bootObs])[[j]]), uniExpansions[[j]]))
          
          # recalculate Bchol if necessary
          if(!bootBasis[[j]]$ortho)
            Bchol[[j]] <- Matrix::chol(bootBasis[[j]]$B)
        } 
        else # resample scores (functions are given and scores can simply be resampled)
          bootBasis[[j]] <- list(scores = uniBasis[[j]]$scores[bootObs, ], B = uniBasis[[j]]$B, ortho = uniBasis[[j]]$ortho, functions = uniBasis[[j]]$functions,
                                 settings = uniBasis[[j]]$settings)
      }
      
      npcBoot <- vapply(bootBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions
      
      if(M > sum(npcBoot))
        stop("Function MFPCA (bootstrap): total number of univariate basis functions must be greater or equal M!")
      
      # calculate MFPCA for bootstrap sample (Bchol has been updated for UMPCA)
      bootMFPCA <- calcMFPCA_2(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                             npc = npcBoot, argvals = mArgvals, uniBasis = bootBasis, fit = FALSE, approx.eigen = approx.eigen)
      
      # save eigenvalues
      booteVals[n,] <- bootMFPCA$values
      
      # flip bootstrap estimates if necessary
      tmpFuns <- flipFuns(res$functions, bootMFPCA$functions)
      
      # save in booteFuns
      for(j in seq_len(p))
      {
        if(dimSupp[j] == 1)
          booteFuns[[j]][n,,] <- tmpFuns[[j]]@X
        if(dimSupp[j] == 2)
          booteFuns[[j]][n,,,] <- tmpFuns[[j]]@X
        if(dimSupp[j] == 3)
          booteFuns[[j]][n,,,,] <- tmpFuns[[j]]@X
      }
    }
    
    CIvalues <- vector("list", length(bootstrapAlpha))
    CI <- vector("list", length(bootstrapAlpha))
    
    for(alpha in seq_len(length(bootstrapAlpha)))
    {
      if(verbose)
        cat("Calculating bootstrap quantiles for alpha = ", bootstrapAlpha[alpha], " (", format(Sys.time(), "%T"), ")\n", sep = "")
      
      CIvalues[[alpha]]$lower <- apply(booteVals, 2, quantile, bootstrapAlpha[alpha]/2)
      CIvalues[[alpha]]$upper <-  apply(booteVals, 2, quantile, 1 - bootstrapAlpha[alpha]/2)
      names(CIvalues)[alpha] <- paste("alpha", bootstrapAlpha[alpha], sep = "_")
      
      bootCI_lower <- bootCI_upper <-  vector("list", p)
      for(j in seq_len(p))
      {
        bootCI_lower[[j]] <- funData(mFData[[j]]@argvals, apply(booteFuns[[j]], 2:length(dim(booteFuns[[j]])),
                                                                quantile, bootstrapAlpha[alpha]/2))
        bootCI_upper[[j]] <- funData(mFData[[j]]@argvals, apply(booteFuns[[j]],  2:length(dim(booteFuns[[j]])),
                                                                quantile, 1 - bootstrapAlpha[alpha]/2))
      }
      
      CI[[alpha]]$lower <- multiFunData(bootCI_lower)
      CI[[alpha]]$upper <- multiFunData(bootCI_upper)
      
      names(CI)[alpha] <- paste("alpha", bootstrapAlpha[alpha], sep = "_")
    }
    
    res$CIvalues <- CIvalues
    
    res$CI <- CI
  }
  
  class(res) <- "MFPCAfit"
  
  return(res)
}


calcMFPCA_2 <- function(N, p, Bchol, M, type, weights, npc, argvals, uniBasis, fit = FALSE, approx.eigen = FALSE)
{
  # combine all scores
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}
  
  # block vector of weights
  allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}
  
  Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)
  
  # check if approximation is appropriate (cf. irlba)
  if(approx.eigen & (M > min(N, sum(npc))/2))
  {
    warning("Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
    approx.eigen = FALSE
  }
  
  # check if non-orthonormal basis functions used and calculate PCA on scores
  if(is.null(Bchol))
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Z), nv = M)
      
      vectors <- tmpSVD$v
      values <- tmpSVD$d^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
      
      e <- eigen(stats::cov(allScores) * outer(allWeights, allWeights, "*"))
      
      values <- e$values
      vectors <- e$vectors[,seq_len(M)]
    }
  }
  else
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Matrix::tcrossprod(Z, Bchol)), nv = M)
      
      vectors <- Matrix::crossprod(Bchol, tmpSVD$v)
      values <- tmpSVD$d^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
      
      e <- eigen(Matrix::crossprod(Bchol) %*% (stats::cov(allScores) * outer(allWeights, allWeights, "*")))
      
      values <- Re(e$values)
      vectors <- Re(e$vectors[,seq_len(M)])
    }
  }
  
  # normalization factors
  normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z %*% vectors))))
  
  ### Calculate scores
  scores <- Z %*% vectors * sqrt(N-1) # see defintion of Z above!
  scores <- as.matrix(scores %*% diag(sqrt(values[seq_len(M)]) * normFactors, nrow = M, ncol = M)) # normalization
  
  ### Calculate eigenfunctions (incl. normalization)
  npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)
  
  tmpWeights <- as.matrix(Matrix::crossprod(Z, Z %*%vectors))
  eFunctions <- foreach::foreach(j = seq_len(p)) %do% {
    univExpansion(type = type[j],
                  scores = 1/sqrt(weights[j] * values[seq_len(M)]) * normFactors * t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                  argvals = argvals[[j]],
                  functions = uniBasis[[j]]$functions,
                  params = uniBasis[[j]]$settings)
  }
  
  res <- list(values = values,
              functions = multiFunData(eFunctions),
              scores = scores,
              vectors = vectors,
              normFactors = normFactors)
  
  # calculate truncated Karhunen-Loeve representation (no mean here)
  if(fit)
    res$fit <- multivExpansion(multiFuns = res$functions, scores = scores)
  
  return(res)
}


multivExpansion <- function(multiFuns, scores)
{
  if(nObs(multiFuns) != NCOL(scores))
    stop("Number of scores does not match number of eigenfunctions.")
  
  # calculate linear combination of multivariate basis functions
  univExp <- foreach::foreach(j = seq_len(length(multiFuns))) %do% { # %do% might require extra loading
    univExpansion(type = "default", 
                  scores = scores,
                  functions = multiFuns[[j]])
  }
  
  # return as multiFunData object
  return(multiFunData(univExp))
}


plot.MFPCAfit_2 <- function(x, plotPCs = seq_len(nObs(x$functions)), stretchFactor = NULL, combined = FALSE, cols, main, cex.main, cex.lab, ...)
{
  if(!inherits(x, "MFPCAfit"))
    stop("Argument is not of class 'MFPCAfit'.")
  
  if(!(is.numeric(plotPCs) & 0 < length(plotPCs) & length(plotPCs) <= length(x$values) & all(0 < plotPCs, plotPCs <= length(x$values))))
    stop("Parameter 'plotPCs' must be a vector with values between 1 and ", length(x$values), ".")
  
  if(!(is.null(stretchFactor) | all(is.numeric(stretchFactor), length(stretchFactor) == 1, stretchFactor > 0)))
    stop("Parameter 'stretchFactor' must be either NULL or a positive number.")
  
  if(!is.logical(combined))
    stop("Parameter 'combined' must be passed as a logical.")
  
  # check dimensions
  dims <- funData::dimSupp(x$functions)
  
  if(any(dims > 2))
    stop("Cannot plot principal components having a 3- or higher dimensional domain.")
  
  # Wait for user input for each new eigenfunction
  oldPar <- graphics::par(no.readonly = TRUE)
  graphics::par(ask = TRUE)
  
  # set number of rows:
  # 1: all dimensions from left to right, "+" and "-" in one plot
  # 2: all dimensions from left to right, upper row for "+", lower for "-"
  nRows <- if(combined == TRUE)
  {
    if(all(dims == 1)) {1} else {
      warning("Cannot combine plots for two-dimensional elements. Will use separate plots (combined = FALSE)")
      2
    } 
  } else{2}
  
  graphics::par(mfrow = c(nRows, length(x$functions)))
  
  
  for(ord in plotPCs) # for each order
  {
    # calculate stretch factor if not given
    if(is.null(stretchFactor))
      stretchFactor <- stats::median(abs(x$scores[,ord]))
    
    PCplus <- x$meanFunction + stretchFactor * x$functions
    PCminus <- x$meanFunction - stretchFactor * x$functions
    
    for(rows in seq_len(nRows))
    {
      for(i in seq_len(length(x$functions))) # for each element
      {
        yRange <- range(PCplus[[i]]@X, PCminus[[i]]@X)
        #main <- paste("PC", ord, "(explains", round(x$values[ord]/sum(x$values)*100, 2), "% of total variability)")
        
        if(dims[i] == 1)
        {
          funData::plot(x$meanFunction[[i]], lwd = 2, col = "black", 
                        main = main[i], ylim = yRange, cex.main = cex.main, cex.lab = cex.lab, ...)
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, cex = 2,
                          add = TRUE, type = "p", pch = "+", col = "darkred", ...)
          if(rows == 2 | combined == TRUE)
            funData::plot( PCminus[[i]], obs = ord, cex = 2,
                           add = TRUE, type = "p", pch = "-", col = "cornflowerblue", ...)
        }  
        else # dims[i] == 2 (higher dimensional domains are not supported)
        {
          if(rows == 1)
            funData::plot(PCplus[[i]], obs = ord, main = main, ylim = yRange, ...)
          else
            funData::plot(PCminus[[i]], obs = ord, main = main, ylim = yRange, ...)
          
        }
      }
    }
  }
  
  graphics::par(oldPar)
  
  invisible(NULL)
}
