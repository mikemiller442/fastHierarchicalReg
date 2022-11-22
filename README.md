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
```

## Usage

``` r
# Simulate data
n <- 10000
numBeta <- 5
betaSD <- 0.75
XSD <- 0.5
errorSD <- 2.0
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
res <- fastHierarchicalReg::linRegGibbsSampler(X = X,
                                               testX = testX,
                                               Y = resp,
                                               testY = testResp,
                                               regVarPrior = regVarPrior,
                                               lambdaSqPrior = lambdaSqPrior,
                                               numEpochs = numEpochs,
                                               numDiscard = numDiscard,
                                               numChains = numChains)
```

    ## socket cluster with 15 nodes on host 'localhost'

``` r
# Compare true beta values to the posterior means
resDF <- data.frame(betaTrue = beta,
                    betaHat = res$postMeanList$beta)
knitr::kable(resDF)
```

|   betaTrue |    betaHat |
|-----------:|-----------:|
| -3.4491918 | -3.4421431 |
| -1.3283198 | -1.3133284 |
| -0.2653254 | -0.2608867 |
|  2.5677206 |  2.6069016 |
|  1.3268278 |  1.3802128 |
