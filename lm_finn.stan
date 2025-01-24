
data {
  int<lower=1> N; // Number of data
  int<lower=1> M; // Number of covariates
  matrix[N, M] X;
  real y[N];
  vector[N] y_sd;
  real<lower=0> scale;
}

// slab_scale = 5, slab_df = 25 -> 8 divergences

transformed data {
  real m0 = 4;           // Expected number of large slopes
  real slab_scale = scale;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df =  25;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;
}

parameters {
  vector[M] beta_tilde;
  real<lower=0,upper=1> gamma[M];
  vector<lower=0>[M] lambda;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
  real alpha;
  real<lower=0> sigma;
  real y_raw[N];
}

transformed parameters {
  vector[M] beta;
  vector[N] y_lat;
  vector[N] y_lat_st;
  {
    real tau0 = (m0 / (M - m0)) * (sigma / sqrt(1.0 * N));
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[M] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    // beta ~ normal(0, tau * lambda_tilde)
    beta = tau * lambda_tilde .* beta_tilde;
  }
  
  y_lat = alpha + X*beta;
  for (i in 1:N) {
    y_lat_st[i] = y_lat[i] + y_raw[i]*sigma;
  }
}

model {
  beta_tilde ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y_raw ~ normal(0, 1);

  y ~ normal(y_lat_st, y_sd);
}
