data {
  int nsnp;
  int ncond;
  matrix[nsnp,ncond] beta_a;
  matrix[nsnp,ncond] beta_b;
}
parameters {
  real<lower=0, upper=1> theta;
  real mu;
  real gamma;
  real<lower=0> sigma_1b;
}
model {
  for (j in 1:ncond) {
    for (i in 1:nsnp) {
      target += log_mix(theta,
                        normal_lpdf(beta_a[i,j] | 0, 0.5),
			normal_lpdf(beta_a[i,j] | mu, 2));
      target += log_mix(theta,
                        normal_lpdf(beta_b[i,j] | 0, 0.5),
			normal_lpdf(beta_b[i,j] | gamma * beta_a[i,j], sigma_1b));
    }
  }
  theta ~ beta(1, 1);
  gamma ~ normal(0, 1);
  mu ~ normal(8, 2);
  sigma_1b ~ normal(0, 5);
}
