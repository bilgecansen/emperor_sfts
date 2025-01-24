
data {
  int<lower=0> N;
  int<lower=0> L;
  int<lower=0> S;
  vector[N] X1;
  vector[L] X_pred1;
  vector[N] X2;
  vector[L] X_pred2;
  vector[N] y;
  vector[N] y_sd;
  vector[S] S_pred;
}

parameters {
  real alpha;
  real beta1;
  real beta2;
  real<lower=0> sigma;
  real y_raw[N];
}

transformed parameters {
  real y_lat[N];
  
  for (i in 1:N) {
    y_lat[i] = X1[i]*beta1 + X2[i]*beta2 + alpha + y_raw[i]*sigma;
  }
  
}

model {
  sigma ~ normal(0, 1);
  y_raw ~ normal(0, 1);
  y ~ normal(y_lat, y_sd);
}

generated quantities {
  real Rsq;
  vector[N] res;
  vector[N] error;
  vector[L] mu_pred1;
  vector[L] mu_pred2;
  vector[N] mu_full;
  vector[N] mu1;
  vector[N] mu2;
  vector[S] mu_s_pred;
  
  mu_full = X1*beta1 + X2*beta2 + alpha;
  mu1 = X1*beta1 + alpha;
  mu2 = X2*beta2 + alpha;
  mu_pred1 = X_pred1*beta1 + alpha;
  mu_pred2 = X_pred2*beta2 + alpha;
  mu_s_pred = S_pred*beta2 + alpha;
  
  for (i in 1:N) {
    res[i] = (mu_full[i]- mean(y))^2;
    error[i] = (y[i] - mu_full[i])^2;
  }
  
  Rsq = (sum(res)/(N-1))/((sum(res)/(N-1)) + (sum(error)/(N-1)));
}
