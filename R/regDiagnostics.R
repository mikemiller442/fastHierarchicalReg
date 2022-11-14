library(Rcpp)
library(RcppArmadillo)
library(foreach)

# Calculate Rhat

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

# Define Gibbs Sampler

linRegGibbsSampler <- function(X,testX,Y,testY,
                               regVarPrior,lambdaSqPrior,
                               numEpochs,numDiscard,numChains) {

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
  resMat <- foreach(k = 1:numChains, .combine = "cbind", .packages = c("Rcpp","RcppArmadillo")) %dopar% {
    sourceCpp("linRegGS.cpp")
    res <- linRegGibbs(X = X,
                       testX = testX,
                       Y = Y,
                       testY = testY,
                       numEpochs = numEpochs,
                       regVarPrior = regVarPrior,
                       lambdaSqPrior = lambdaSqPrior)

    postIntercept <- res$intercept[,(numDiscard+2):(numEpochs+1)]
    postBeta <- res$coefBeta[,(numDiscard+2):(numEpochs+1)]
    postLambdaSq <- res$lambdaSq[(numDiscard+2):(numEpochs+1)]
    postLambdaSqScale <- res$lambdaSqScale[(numDiscard+2):(numEpochs+1)]
    postRegVar <- res$regVar[(numDiscard+2):(numEpochs+1)]
    postRegVarScale <- res$regVarScale[(numDiscard+2):(numEpochs+1)]

    samples <- rbind(postIntercept,postBeta,postLambdaSq,postLambdaSqScale,
                     postRegVar,postRegVarScale)
    samples
  }

  RhatIntercept <- Rhat(resMat[1,],numChains)
  RhatBeta <- rep(NA,ncol(X))
  for (i in 1:ncol(X)) {
    RhatBeta[i] <- Rhat(resMat[1+i,],numChains)
  }
  RhatLambdaSq <- Rhat(resMat[ncol(X)+2,],numChains)
  RhatLambdaSqScale <- Rhat(resMat[ncol(X)+3,],numChains)
  RhatRegVar <- Rhat(resMat[ncol(X)+4,],numChains)
  RhatRegVarScale <- Rhat(resMat[ncol(X)+5,],numChains)

  RhatList <- list(RhatIntercept = RhatIntercept,
                   RhatBeta = RhatBeta,
                   RhatLambdaSq = RhatLambdaSq,
                   RhatLambdaSqScale = RhatLambdaSqScale,
                   RhatRegVar = RhatRegVar,
                   RhatRegVarScale = RhatRegVarScale)

  return(list(resMat = resMat,
              postMean = rowMeans(resMat[,(numDiscard+2):(numEpochs+1)]),
              RhatList = RhatList))
}














