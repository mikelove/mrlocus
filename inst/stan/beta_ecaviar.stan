data {
  int nsnp;
  int ncond;
  matrix[nsnp,ncond] beta_hat_a;
  matrix[nsnp,ncond] beta_hat_b;
  real Sigma_a[nsnp,nsnp,ncond];
  real Sigma_b[nsnp,nsnp,ncond];
  real se2;
}
parameters {
  matrix[nsnp,ncond] beta_a;
  matrix[nsnp,ncond] beta_b;
}
model {
  for (j in 1:ncond) {
    beta_hat_a[,j] ~ multi_normal(to_matrix(Sigma_a[,,j]) * beta_a[,j], se2 * to_matrix(Sigma_a[,,j]));
    beta_hat_b[,j] ~ multi_normal(to_matrix(Sigma_b[,,j]) * beta_b[,j], se2 * to_matrix(Sigma_b[,,j]));
    for (i in 1:nsnp) {
      beta_a[i,j] ~ normal(0, 10);
      beta_b[i,j] ~ normal(0, 10);
    }
  }
}
