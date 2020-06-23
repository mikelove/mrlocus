data {
  int n; 
  vector[n] beta_hat_a;
  vector[n] beta_hat_b;
  vector[n] sd_a;
  vector[n] sd_b;
  real alpha_sd;
  real sigma_sd;
}
parameters {
  real alpha;
  real<lower=0> sigma;
  vector[n] beta_a;
  vector[n] beta_b;
}
model {
  beta_hat_a ~ normal(beta_a, sd_a);
  beta_hat_b ~ normal(beta_b, sd_b);
  beta_b ~ normal(alpha * beta_a, sigma);
  alpha ~ normal(0, alpha_sd);
  sigma ~ normal(0, sigma_sd);
}
