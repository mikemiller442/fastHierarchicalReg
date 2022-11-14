#include <RcppArmadillo.h>;
#include <math.h>
#include <chrono>
#include <thread>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

class FeatureSet {
public:
  arma::Mat<double> features;
  arma::Mat<double> testFeatures;
  arma::Mat<double> beta;
  arma::Mat<double> priorMat;
  arma::Mat<double> covMat;
  int numWeights;
  int numEpochs;
  double interceptPrior;
  
  void construct(int n, int testN, int ne) {
    features = arma::ones(n,1);
    testFeatures = arma::ones(testN,1);
    numWeights = features.n_cols;
    numEpochs = ne;
    beta = arma::zeros(numWeights,numEpochs + 1);
    beta.col(0) = arma::randu(numWeights,arma::distr_param(-1.0,1.0));
    interceptPrior = 100.0;
  }
};

class FeatureSetHP {
public:
  arma::Mat<double> features;
  arma::Mat<double> testFeatures;
  arma::Mat<double> beta;
  arma::Col<double> lambdaSqScale;
  arma::Col<double> lambdaSq;
  int numWeights;
  int numEpochs;
  arma::Mat<double> U;
  arma::Col<double> D;
  arma::Mat<double> V;
  arma::Mat<double> tV;
  arma::Mat<double> Z;
  arma::Col<double> Zdiag;
  double lambdaSqScaleAlpha;
  double lambdaSqScaleBeta;
  double lambdaSqAlpha;
  double lambdaSqBeta;
  double lambdaSqPrior;
  
  void construct(arma::Mat<double> X, arma::Mat<double> testX,
                 double lsp, int ne) {
    features = X;
    testFeatures = testX;
    numWeights = features.n_cols;
    numEpochs = ne;
    beta = arma::zeros(numWeights,numEpochs + 1);
    arma::svd(U,D,V,features);
    tV = arma::trans(V);
    Z = arma::eye(X.n_cols, X.n_cols);
    Zdiag = arma::vec(X.n_cols);
    if (X.n_cols > X.n_rows) {
      D = join_cols(D, arma::zeros(X.n_cols - X.n_rows));
    }
    lambdaSqPrior = lsp;
    lambdaSqScale = arma::zeros(numEpochs + 1);
    lambdaSqScale(0) = 1.0;
    lambdaSq = arma::zeros(numEpochs + 1);
    lambdaSq(0) = 1.0/(float)numWeights;
  }
};

class LinRegModel {
public:
  arma::Col<double> Y;
  arma::Col<double> testY;
  arma::Col<double> regVarScale;
  arma::Col<double> regVar;
  int n;
  int testN;
  int numEpochs;
  double regVarScaleAlpha;
  double regVarScaleBeta;
  double regVarAlpha;
  double regVarBeta;
  double regVarPrior;
  FeatureSet feat;
  FeatureSetHP featHP;
  
  void construct(arma::Mat<double> X, arma::Mat<double> testX,
                 arma::Col<double> resp, arma::Col<double> testResp, int ne,
                 double rvp, double lambdaSqPrior) {
    Y = resp;
    testY = testResp;
    n = Y.n_elem;
    numEpochs = ne;
    regVarPrior = rvp;
    regVarScale = arma::zeros(numEpochs + 1);
    regVarScale(0) = 1.0;
    regVar = arma::zeros(numEpochs + 1);
    regVar(0) = 1.0;
    feat.construct(X.n_rows,testX.n_rows,ne);
    featHP.construct(X,testX,lambdaSqPrior,ne);
  }
};

// [[Rcpp::export]]
// Linear Regression Gibbs Sampler
Rcpp::List linRegGibbs(NumericMatrix X, NumericMatrix testX, NumericVector Y, NumericVector testY,
                       int numEpochs, double regVarPrior, double lambdaSqPrior) {
  
  LinRegModel mod;
  mod.construct(as<arma::mat>(X), as<arma::mat>(testX),
                as<arma::vec>(Y), as<arma::vec>(testY), numEpochs,
                regVarPrior, lambdaSqPrior);
  
  arma::Col<double> trainResiduals;
  arma::Col<double> evBeta;
  arma::Col<double> stdNormal;
  
  for (int epoch = 0; epoch < numEpochs; epoch++) {
    
    // Sample intercept
    trainResiduals = mod.Y - mod.feat.features*mod.feat.beta.col(epoch)
                   - mod.featHP.features*mod.featHP.beta.col(epoch);
    mod.feat.priorMat = arma::eye(mod.feat.numWeights,mod.feat.numWeights)
                      *mod.feat.interceptPrior;
    mod.feat.covMat = arma::inv(arma::trans(mod.feat.features)*mod.feat.features + mod.feat.priorMat);
    evBeta = mod.feat.covMat*(arma::trans(mod.feat.features)*
             (trainResiduals + mod.feat.features*mod.feat.beta.col(epoch)));
    mod.feat.beta.col(epoch + 1) = (arma::mvnrnd(evBeta,mod.feat.covMat));
    
    // Sample regression coefficients
    trainResiduals = mod.Y - mod.feat.features*mod.feat.beta.col(epoch + 1)
                   - mod.featHP.features*mod.featHP.beta.col(epoch);
    mod.featHP.Zdiag = arma::pow((1.0/mod.regVar(epoch))*arma::pow(mod.featHP.D, 2.0)
                     + (1.0/mod.regVar(epoch))*(1.0/mod.featHP.lambdaSq(epoch)),-1.0);
    mod.featHP.Z.diag() = mod.featHP.Zdiag;
    evBeta = (1.0/mod.regVar(epoch))*mod.featHP.V*(mod.featHP.Z
           *(mod.featHP.tV*(arma::trans(mod.featHP.features)
           *(trainResiduals + mod.featHP.features*mod.featHP.beta.col(epoch)))));
    stdNormal = arma::randn(mod.featHP.numWeights,arma::distr_param(0.0,1.0));
    mod.featHP.beta.col(epoch + 1) = evBeta + mod.featHP.V*
                                     (arma::pow(mod.featHP.Zdiag, 0.5) % stdNormal);
    
    // Sample global shrinkage parameters
    mod.featHP.lambdaSqScaleAlpha = 1.0;
    mod.featHP.lambdaSqScaleBeta = 1.0/mod.featHP.lambdaSq(epoch)
                                 + pow(mod.featHP.lambdaSqPrior, -2.0);
    mod.featHP.lambdaSqScale(epoch + 1) = 1.0/arma::randg(arma::distr_param(mod.featHP.lambdaSqScaleAlpha,
                                          1.0/mod.featHP.lambdaSqScaleBeta));
    mod.featHP.lambdaSqAlpha = 0.5*(mod.featHP.numWeights + 1.0);
    mod.featHP.lambdaSqBeta = 0.5*(1.0/mod.regVar(epoch))*arma::dot(mod.featHP.beta.col(epoch+1),
                              mod.featHP.beta.col(epoch+1)) + (1.0/mod.featHP.lambdaSqScale(epoch+1));
    mod.featHP.lambdaSq(epoch + 1) = 1.0/randg(arma::distr_param(mod.featHP.lambdaSqAlpha,
                                     1.0/mod.featHP.lambdaSqBeta));
    
    // Sample regression variance
    trainResiduals = mod.Y - mod.feat.features*mod.feat.beta.col(epoch + 1)
                   - mod.featHP.features*mod.featHP.beta.col(epoch + 1);
    mod.regVarScaleAlpha = 1.0;
    mod.regVarScaleBeta = 1.0/mod.regVar(epoch) + pow(mod.regVarPrior, -2.0);
    mod.regVarScale(epoch + 1) = 1.0/arma::randg(arma::distr_param(mod.regVarScaleAlpha,
                                 1.0/mod.regVarScaleBeta));
    mod.regVarAlpha = 0.5*(mod.n + mod.featHP.numWeights + 1.0);
    mod.regVarBeta = 1.0/mod.regVarScale(epoch+1) + 0.5*arma::dot(trainResiduals,trainResiduals)
                   + 0.5*(1.0/mod.featHP.lambdaSq(epoch+1))
                   *arma::dot(mod.featHP.beta.col(epoch+1),mod.featHP.beta.col(epoch+1));
    mod.regVar(epoch + 1) = 1.0/arma::randg(arma::distr_param(mod.regVarAlpha,1.0/mod.regVarBeta));
  }
  
  return Rcpp::List::create(Named("intercept") = mod.feat.beta,
                            Named("coefBeta") = mod.featHP.beta,
                            Named("regVarScale") = mod.regVarScale,
                            Named("regVar") = mod.regVar,
                            Named("lambdaSqScale") = mod.featHP.lambdaSqScale,
                            Named("lambdaSq") = mod.featHP.lambdaSq);
}
















