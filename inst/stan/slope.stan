data {
  int n; 
  vector[n] beta_hat_a;
  vector[n] beta_hat_b;
  vector[n] sd_a;
  vector[n] sd_b;
  real sd_beta;
  real mu_alpha;
  real sd_alpha;
  real sd_sigma;
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
  beta_a ~ normal(0, sd_beta);
  beta_b ~ normal(alpha * beta_a, sigma);
  alpha ~ normal(mu_alpha, sd_alpha);
  sigma ~ normal(0, sd_sigma);
}
