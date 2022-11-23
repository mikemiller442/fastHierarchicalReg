fastHierarchicalReg
================
Michael Miller
2022-11-22

## Overview

fastHierarchicalReg provides functionality for estimating Bayesian
hierarchical linear models (HLMs). Bayesian HLMs are very useful in
applied statistics because the model learns the prior variance of the
model parameters, which functions as a model regularizer, from the data
rather than from cross validation (CV). This is critical from an
inferential point of view as we can then use the resulting credible
intervals to perform inference without having to perform CV for the
regularization parameter in addition to calculating bootstrapped
confidence intervals.

The primary benchmark used in this package is the HLM presented in the
fastBayesReg repository hosted by KangJian2016. The speed of these
implementations follows in part from pre-computing the SVD of the design
matrix to obtain sufficient statistics for the conditional posterior
distributions. This allows us to reduce the time complexity of each
iteration to $O(p^2)$ in the number of model covariates rather than
$O(p^3)$ required for the matrix inversion in the naive blocked Gibbs
sampler. Implementing the blocked Gibbs sampler is critical because it
drastically reduces the random walk behavior that Gibbs sampling
exhibits for high dimensional posterior distributions. The functionality
of this package includes:

- Parallel computing of the Markov chains using the foreach package
- Efficient Blocked Gibbs sampler implementation in R
- Calculation of Rhat convergence diagnostics for model parameters

## Installation

``` r
library(devtools)
devtools::install_github("mikemiller442/fastHierarchicalReg")
library(fastHierarchicalReg)
```

## Usage

``` r
# Simulate data
n <- 10000   # subjects
numBeta <- 5   # covariates
betaSD <- 0.75   # standard deviation of the randomly sampled covariates
XSD <- 0.5   # standard deviation of the generated features
errorSD <- 2.0   # regression standard deviation

e <- rnorm(n, mean = 0, sd = errorSD)
beta <- rnorm(numBeta, mean = 0, sd = betaSD*errorSD)
Z <- matrix(NA, nrow = n, ncol = numBeta)
for (i in 1:ncol(Z)) {
  Z[,i] <- rnorm(n, mean = 0, sd = XSD)
}

y <- Z %*% beta + e
output <- y

X <- Z
testX <- Z
resp <- output
testResp <- output

numEpochs <- 10000
numDiscard <- 2000
numChains <- 4
numCores <- 8
lambdaSqPrior <- 1.0
regVarPrior <- 1.0

# Call Gibbs Sampler Directly For A Single Chain
res <- fastHierarchicalReg::linRegGibbs(X = X,
                                        testX = testX,
                                        Y = resp,
                                        testY = testResp,
                                        numEpochs = numEpochs,
                                        regVarPrior = regVarPrior,
                                        lambdaSqPrior = lambdaSqPrior)

# Call R Function To Compute Multiple Chains In Parallel
res <- fastHierarchicalReg::linRegGibbsProcessed(X = X,
                                                 testX = testX,
                                                 Y = resp,
                                                 testY = testResp,
                                                 lambdaSqPrior = lambdaSqPrior,
                                                 regVarPrior = regVarPrior,
                                                 numEpochs = numEpochs,
                                                 numDiscard = numDiscard,
                                                 numChains = numChains,
                                                 numCores = numCores)
```

    ## socket cluster with 8 nodes on host 'localhost'

``` r
# Compare true beta values to the posterior means
resDF <- data.frame(betaTrue = beta,
                    betaHat = res$postMeanList$beta)
knitr::kable(resDF)
```

|   betaTrue |    betaHat |
|-----------:|-----------:|
| -1.2124422 | -1.2235755 |
|  0.4166135 |  0.4690173 |
|  2.3911720 |  2.3808407 |
| -0.7461658 | -0.7664415 |
| -1.4909645 | -1.4654161 |
