#' Linear Regression Gibbs Sampler
#'
#' Estimates the posterior distribution of a Hierarchical Linear Model.
#' Places Half-Cauchy priors on the regression coefficient prior variance and the regression variance.
#' Introduces auxiliary scale parameters in order to sample a Half Cauchy variable
#' via sampling an inverse gamma - inverse gamma mixture using Gibbs sampling.
#'
#' @param X Design Matrix
#' @param testX Design Matrix
#' @param Y Outcome
#' @param testY Outcome
#' @param regVarPrior Scale hyperparameter for the Half-Cauchy prior for the regression variance.
#' @param lambdaSqPrior Scale hyperparameter for the Half-Cauchy prior for the regression coefficient prior variance.
#' @param numEpochs Number of epochs to run the Gibbs Sampler
#' @param numDiscard Number of epochs to throw away as burn in
#' @return List of posterior samples
linRegGibbs <- function(X,testX,Y,testY,numEpochs,regVarPrior,lambdaSqPrior) {

  # Initialize Markov Chain And Compute Sufficient Statistics
  n <- length(Y)
  numWeights <- ncol(X)
  beta <- matrix(0.0, nrow = numWeights, ncol = numEpochs + 1)
  svdRes <- svd(X, nu = nrow(X), nv = ncol(X))
  U <- svdRes$u
  D <- svdRes$d
  V <- svdRes$v
  tV <- t(svdRes$v)
  Zdiag <- rep(0.0, ncol(X))
  if (ncol(X) > nrow(X)) {
    D <- c(D, rep(0.0,ncol(X) - nrow(X)))
  }
  lambdaSqScale <- rep(0.0,numEpochs + 1)
  lambdaSqScale[1] <- 1.0
  lambdaSq <- rep(0.0,numEpochs + 1)
  lambdaSq[1] <- 1.0/numWeights

  regVarScale <- rep(0.0,numEpochs + 1)
  regVarScale[1] <- 1.0
  regVar <- rep(0.0,numEpochs + 1)
  regVar[1] <- 1.0

  # Run Markov Chain
  for (epoch in 1:(numEpochs)) {
    # Sample regression coefficients
    Zdiag <- (1.0/regVar[epoch]*D^2.0 + (1.0/regVar[epoch])*(1.0/lambdaSq[epoch]))^(-1.0)
    evBeta <- (1.0/regVar[epoch])*(V %*% (diag(Zdiag) %*% (tV %*% (t(X) %*% Y))))
    beta[,epoch+1] <- evBeta + V %*% (sqrt(Zdiag) * rnorm(numWeights))

    # Sample global shrinkage parameter
    lambdaSqScaleAlpha <- 1.0
    lambdaSqScaleBeta <- 1.0/lambdaSq[epoch] + lambdaSqPrior^(-2.0)
    lambdaSqScale[epoch + 1] <- 1.0/rgamma(1,lambdaSqScaleAlpha,lambdaSqScaleBeta)
    lambdaSqAlpha <- 0.5*(numWeights + 1.0)
    lambdaSqBeta <- 0.5*(1.0/regVar[epoch])*sum(beta[,epoch+1] * beta[,epoch+1]) + (1.0/lambdaSqScale[epoch+1])
    lambdaSq[epoch + 1] <- 1.0/rgamma(1,lambdaSqAlpha,lambdaSqBeta)

    # Sample regression variance
    trainResiduals <- Y - X %*% beta[,epoch + 1]
    regVarScaleAlpha <- 1.0
    regVarScaleBeta <- 1.0/regVar[epoch] + regVarPrior^(-2.0)
    regVarScale[epoch + 1] <- 1.0/rgamma(1,regVarScaleAlpha,regVarScaleBeta)
    regVarAlpha <- 0.5*(n + numWeights + 1.0)
    regVarBeta <- 1.0/regVarScale[epoch+1] + 0.5*sum(trainResiduals*trainResiduals) + 0.5*(1.0/lambdaSq[epoch+1])*sum(beta[,epoch+1] * beta[,epoch+1])
    regVar[epoch + 1] <-1.0/rgamma(1,regVarAlpha,regVarBeta)
  }

  return(list(coefBeta = beta,
              lambdaSqScale = lambdaSqScale,
              lambdaSq = lambdaSq,
              regVarScale = regVarScale,
              regVar = regVar))
}

#' Calculate Rhat Convergence Diagnostics From MCMC Output
#'
#' @param samples The number of samples in each chain.
#' @param numChains The number of chains.
#' @return Vector Of Rhat Convergence Diagnostics For Model Parameters.
#' @examples
Rhat <- function(samples,numChains) {
  n <- length(samples)/numChains
  si2 <- rep(NA,numChains)
  mui <- rep(NA,numChains)
  for (i in 1:numChains) {

    # Calculate chain mean
    chainMean <- samples[(i-1)*n + 1]
    for (j in 2:n) {
      chainMean <- chainMean + samples[(i-1)*n+j]
    }

    chainMean <- chainMean/n

    # Calculate chain sample variance
    chainVar <- (samples[(i-1)*n + 1] - chainMean)^2
    for (j in 2:n) {
      chainVar <- chainVar + (samples[(i-1)*n+j] - chainMean)^2
    }

    chainVar <- chainVar/(n - 1)
    mui[i] <- chainMean
    si2[i] <- chainVar
  }

  # Calculate overall mean
  mu <- mui[1]
  for (i in 2:numChains) {
    mu <- mu + mui[i]
  }
  mu <- mu/numChains

  # Calculate average chain sample variance
  s2 <- si2[1]
  for (i in 2:numChains) {
    s2 <- s2 + si2[i]
  }
  s2 <- s2/numChains

  # Calculate scaled between chain sample variance
  B <- (mui[1] - mu)^2
  for (i in 2:numChains) {
    B <- B + (mui[i] - mu)^2
  }
  B <- (n/(numChains - 1))*B

  #Calculate Rhat
  sigmaHat2 <- ((n-1)/n)*s2 + B/n
  Rhat <- sqrt(sigmaHat2/s2)
  return(Rhat)
}

#' Gibbs Sampling With Multiple MCMC Chains Running In Parallel
#'
#' @param X Design Matrix
#' @param testX Design Matrix
#' @param Y Outcome
#' @param testY Outcome
#' @param regVarPrior Scale hyperparameter for the Half-Cauchy prior for the regression variance.
#' @param lambdaSqPrior Scale hyperparameter for the Half-Cauchy prior for the regression coefficient prior variance.
#' @param numEpochs Number of epochs to run the Gibbs Sampler
#' @param numDiscard Number of epochs to throw away as burn in
#' @param numChains Number of MCMC chains to run
#' @return List of posterior samples, posterior means, and Rhat criteria
linRegGibbsSampler <- function(X,testX,Y,testY,
                               regVarPrior,lambdaSqPrior,
                               numEpochs,numDiscard,numChains) {

  library(parallel)
  library(foreach)

  ### Parallel Setup
  n.cores <- parallel::detectCores() - 1
  #create the cluster
  my.cluster <- parallel::makeCluster(n.cores,type = "PSOCK")

  #check cluster definition (optional)
  print(my.cluster)

  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  #check if it is registered (optional)
  foreach::getDoParRegistered()

  #how many workers are available? (optional)
  foreach::getDoParWorkers()

  ### Run Markov Chains In Parallel
  resMat <- foreach(k = 1:numChains, .combine = "cbind") %dopar% {
    res <- linRegGibbs(X = X,
                       testX = testX,
                       Y = Y,
                       testY = testY,
                       numEpochs = numEpochs,
                       regVarPrior = regVarPrior,
                       lambdaSqPrior = lambdaSqPrior)

    postBeta <- res$coefBeta[,(numDiscard+2):(numEpochs+1)]
    postLambdaSq <- res$lambdaSq[(numDiscard+2):(numEpochs+1)]
    postLambdaSqScale <- res$lambdaSqScale[(numDiscard+2):(numEpochs+1)]
    postRegVar <- res$regVar[(numDiscard+2):(numEpochs+1)]
    postRegVarScale <- res$regVarScale[(numDiscard+2):(numEpochs+1)]

    samples <- rbind(postBeta,postLambdaSq,postLambdaSqScale,
                     postRegVar,postRegVarScale)
    samples
  }

  RhatBeta <- rep(NA,ncol(X))
  for (i in 1:ncol(X)) {
    RhatBeta[i] <- Rhat(resMat[i,],numChains)
  }
  RhatLambdaSq <- Rhat(resMat[ncol(X)+1,],numChains)
  RhatLambdaSqScale <- Rhat(resMat[ncol(X)+2,],numChains)
  RhatRegVar <- Rhat(resMat[ncol(X)+3,],numChains)
  RhatRegVarScale <- Rhat(resMat[ncol(X)+4,],numChains)

  RhatList <- list(RhatBeta = RhatBeta,
                   RhatLambdaSq = RhatLambdaSq,
                   RhatLambdaSqScale = RhatLambdaSqScale,
                   RhatRegVar = RhatRegVar,
                   RhatRegVarScale = RhatRegVarScale)

  postMean <- rowMeans(resMat[,(numDiscard+2):ncol(resMat)])
  postMeanList <- list(beta = as.numeric(postMean[1:ncol(X)]),
                       lambdaSq = as.numeric(postMean[ncol(X)+1]),
                       lambdaSqScale = as.numeric(postMean[ncol(X)+2]),
                       regVar = as.numeric(postMean[ncol(X)+3]),
                       regVarScale = as.numeric(postMean[ncol(X)+4]))

  return(list(resMat = resMat,
              postMeanList = postMeanList,
              RhatList = RhatList))
}















