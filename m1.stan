data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of covariates
  matrix[N, K] X;  // design matrix
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector[K] b;  // population-level effects
  real<lower=0> lambda;  // regression coefficient prior variance
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(Intercept | 0, 100.0);
  lprior += normal_lpdf(b | 0, sigma*lambda);
  lprior += student_t_lpdf(lambda | 1, 0, 1) - 1 * student_t_lccdf(0 | 1, 0, 1);
  lprior += student_t_lpdf(sigma | 1, 0, 1) - 1 * student_t_lccdf(0 | 1, 0, 1);
}
model {
  // initialize linear predictor term
  vector[N] mu = rep_vector(0.0, N);
  mu += Intercept;
  target += normal_id_glm_lpdf(Y | X, mu, b, sigma);
  target += lprior;
}










