context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  library(MASS)
  
  set.seed(1)
  nsnp <- 29
  ncond <- 4
  Sigma_a <- Sigma_b <- array(NA, dim=c(nsnp,nsnp,ncond))
  for (j in 1:ncond) {
    Sigma_a[,,j] <- diag(nsnp) # A will be eQTL
    Sigma_b[,,j] <- diag(nsnp) # B will be GWAS
    idx <- ceiling(nsnp/2) + -5:5
    Sigma_a[idx,idx,j] <- ifelse(Sigma_a[idx,idx,j] == 0, .5, 1)
    Sigma_b[idx,idx,j] <- ifelse(Sigma_b[idx,idx,j] == 0, .5, 1)
  }
  x <- (nsnp - 1)/2
  beta <- sapply(1:ncond, function(j) rep(c(0,6 + j,0),c(x,1,x)))
  beta_hat_a <- beta
  beta_hat_b <- beta
  se_a <- se_b <- matrix(0.2, nsnp, ncond)
  gamma <- 0.5 
  sigma <- 0.1
  theta <- 1 - 1/nsnp
  for (j in 1:ncond) {
    beta_hat_a[,j] <- mvrnorm(1,
                              mu=Sigma_a[,,j] %*% beta[,j],
                              diag(se_a[,j]) * Sigma_a[,,j] * diag(se_a[,j]))
    beta_hat_b[,j] <- mvrnorm(1,
                              mu=Sigma_b[,,j] %*% (gamma * beta[,j] + rnorm(nsnp,0,sigma)),
                              diag(se_b[,j]) * Sigma_b[,,j] * diag(se_b[,j]))
  }
  data <- list(nsnp=nsnp,
               ncond=ncond,
               beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               se_a=se_a,
               se_b=se_b,
               Sigma_a=Sigma_a,
               Sigma_b=Sigma_b)

  options(mc.cores=4)
  fit1 <- fitBetaEcaviar(data)

  print(fit1, pars=c("beta_a[1,1]","beta_b[1,1]"), digits=3)
  rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp,",1]"))
  rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp,",1]"))

  # here dropping the inferential uncertainty
  coefs1 <- rstan::extract(fit1)
  beta_clean_a <- sapply(1:ncond, function(j) colMeans(coefs1$beta_a[,,j]))
  beta_clean_b <- sapply(1:ncond, function(j) colMeans(coefs1$beta_b[,,j]))

  data <- list(nsnp=nsnp,
               ncond=ncond,
               beta_a=beta_clean_a,
               beta_b=beta_clean_b)

  options(mc.cores=4)
  fit2 <- fitBetaMixture(data)

  print(fit2, pars=c("theta","mu","gamma","sigma_1b"), digits=3)
  rstan::stan_plot(fit2, pars="theta")
  rstan::stan_plot(fit2, pars="mu")
  rstan::stan_plot(fit2, pars="gamma")
  rstan::stan_plot(fit2, pars="sigma_1b")

  coefs2 <- rstan::extract(fit2)
  
  plot(beta_clean_a, beta_clean_b)
  abline(0, mean(coefs2$gamma), lwd=2)
  abline(mean(coefs2$sigma_1b), mean(coefs2$gamma), col="blue")
  abline(-mean(coefs2$sigma_1b), mean(coefs2$gamma), col="blue")
  
})
