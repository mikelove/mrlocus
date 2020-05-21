data {
  int nsnp;
  int ncond;
  int n[ncond]; // first n[j] SNPs used
  matrix[nsnp,ncond] beta_hat_a;
  matrix[nsnp,ncond] beta_hat_b;
  matrix[nsnp,ncond] sd_a;
  matrix[nsnp,ncond] sd_b;
}
parameters {
  real<lower=0, upper=1> theta;
  real mu;
  real gamma;
  real<lower=0> sigma_1b;
  matrix[nsnp,ncond] beta_a;
  matrix[nsnp,ncond] beta_b;
}
model {
  for (j in 1:ncond) {
    segment(beta_hat_a[,j],1,n[j]) ~
      normal(head(beta_a[,j],n[j]), head(sd_a[,j],n[j]));
    segment(beta_hat_b[,j],1,n[j]) ~
      normal(head(beta_b[,j],n[j]), head(sd_b[,j],n[j]));
    for (i in 1:n[j]) {
      target += log_mix(theta,
                        normal_lpdf(beta_a[i,j] | 0, 0.5),
			normal_lpdf(beta_a[i,j] | mu, 2));
      target += log_mix(theta,
                        normal_lpdf(beta_b[i,j] | 0, 0.5),
			normal_lpdf(beta_b[i,j] | gamma * beta_a[i,j], sigma_1b));
    }
    if (n[j] < nsnp) {
      for (i in (n[j]+1):nsnp) {
        beta_a[i,j] ~ normal(0, 1);
	beta_b[i,j] ~ normal(0, 1);
      }
    }
  }
  theta ~ beta(1, 1);
  gamma ~ normal(0, 1);
  mu ~ normal(8, 2);
  sigma_1b ~ normal(0, 5);
}
