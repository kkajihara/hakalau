data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> n_obs;
  vector[J] x[n_obs];
  int<lower=0> y1[n_obs];
  int<lower=0> y2[n_obs];
}
parameters {
  matrix[K, J] beta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
  vector[K] y[n_obs];
}
model {
  vector[K] mu[n_obs];
  matrix[K, K] L_Sigma;

  for (n in 1:n_obs)
    mu[n] = beta * x[n];

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  to_vector(beta) ~ normal(0, 5);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  y ~ multi_normal_cholesky(mu, L_Sigma);
  
  for (n in 1:n_obs) {
    y1[n] ~ poisson(exp(y[n, 1]));
    y2[n] ~ poisson(exp(y[n, 2]));
  }

}
generated quantities {
  corr_matrix[K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

