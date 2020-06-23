data {
  int n;
  vector[n] beta_hat_a;
  vector[n] beta_hat_b;
  vector[n] se_a;
  vector[n] se_b;
  matrix[n,n] Sigma_a;
  matrix[n,n] Sigma_b;
}
parameters {
  vector[n] beta_a;
  vector[n] beta_b;
  vector<lower=0>[n] lambda;
  real<lower=0> tau;
}
model {
  tau ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  beta_hat_a ~ normal(Sigma_a * beta_a, se_a);
  beta_hat_b ~ normal(Sigma_b * beta_b, se_b);
  for (i in 1:n) {
    beta_a[i] ~ normal(0, lambda[i] * tau);
    beta_b[i] ~ normal(0, lambda[i] * tau);
  }
}
