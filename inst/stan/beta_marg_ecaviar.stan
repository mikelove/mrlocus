data {
  int nsnp;
  vector[nsnp] beta_hat_a;
  vector[nsnp] beta_hat_b;
  vector[nsnp] se_a;
  vector[nsnp] se_b;
  matrix[nsnp,nsnp] Sigma_a;
  matrix[nsnp,nsnp] Sigma_b;
}
parameters {
  vector[nsnp] beta_a;
  vector[nsnp] beta_b;
  vector<lower=0>[nsnp] lambda;
  real<lower=0> tau;
}
model {
  tau ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  beta_hat_a ~ normal(Sigma_a * beta_a, se_a);
  beta_hat_b ~ normal(Sigma_b * beta_b, se_b);
  for (i in 1:nsnp) {
    beta_a[i] ~ normal(0, lambda[i] * tau);
    beta_b[i] ~ normal(0, lambda[i] * tau);
  }
}
