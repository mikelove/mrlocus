data {
  int tot; // sum of nsnp across clusters
  vector[tot] beta_hat_a;
  vector[tot] beta_hat_b;
  vector[tot] sd_a;
  vector[tot] sd_b;
  real sigma_0a;
  real sigma_0b;
  real sigma_1a;
  real mu_loc;
  real mu_sd;
  real sigma1b_sd;
}
parameters {
  real<lower=0, upper=1> theta;
  real mu;
  real gamma;
  real<lower=0> sigma_1b;
  vector[tot] beta_a;
  vector[tot] beta_b;
}
model {
  beta_hat_a ~ normal(beta_a, sd_a);
  beta_hat_b ~ normal(beta_b, sd_b);
  for (i in 1:tot) {
    target += log_mix(theta,
    	              normal_lpdf(beta_a[i] | 0, sigma_0a),
		      normal_lpdf(beta_a[i] | mu, sigma_1a));
    target += log_mix(theta,
                      normal_lpdf(beta_b[i] | 0, sigma_0b),
	              normal_lpdf(beta_b[i] | gamma * beta_a[i], sigma_1b));
  }
  theta ~ beta(1, 1);
  gamma ~ normal(0, 1);
  mu ~ normal(mu_loc, mu_sd);
  sigma_1b ~ normal(0, sigma1b_sd);
}
