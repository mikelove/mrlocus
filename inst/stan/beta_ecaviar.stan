data {
  int nsnp;
  int ncond;
  int n[ncond]; // first n[j] SNPs used
  matrix[nsnp,ncond] beta_hat_a;
  matrix[nsnp,ncond] beta_hat_b;
  matrix[nsnp,ncond] se_a;
  matrix[nsnp,ncond] se_b;
  real Sigma_a[nsnp,nsnp,ncond];
  real Sigma_b[nsnp,nsnp,ncond];
}
parameters {
  matrix[nsnp,ncond] beta_a;
  matrix[nsnp,ncond] beta_b;
}
model {
  for (j in 1:ncond) {
    segment(beta_hat_a[,j],1,n[j]) ~
      multi_normal(block(to_matrix(Sigma_a[,,j]),1,1,n[j],n[j]) *
        head(beta_a[,j],n[j]),
    	quad_form_diag(block(to_matrix(Sigma_a[,,j]),1,1,n[j],n[j]),
	head(se_a[,j],n[j])));
    segment(beta_hat_b[,j],1,n[j]) ~
      multi_normal(block(to_matrix(Sigma_b[,,j]),1,1,n[j],n[j]) *
        head(beta_b[,j],n[j]),
        quad_form_diag(block(to_matrix(Sigma_b[,,j]),1,1,n[j],n[j]),
        head(se_b[,j],n[j])));
    for (i in 1:n[j]) {
      beta_a[i,j] ~ normal(0, 10);
      beta_b[i,j] ~ normal(0, 10);
    }
    if (n[j] < nsnp) {
      for (i in (n[j]+1):nsnp) {
        beta_a[i,j] ~ normal(0, 1);
	beta_b[i,j] ~ normal(0, 1);
      }
    }
  }
}
