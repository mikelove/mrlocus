context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  library(MASS)
  
  set.seed(1)
  nsnp <- 15
  ncond <- 3
  n <- c(5,10,15)
  Sigma_a <- Sigma_b <- array(NA, dim=c(nsnp,nsnp,ncond))
  for (j in 1:ncond) {
    Sigma_a[,,j] <- diag(nsnp) # A will be eQTL
    Sigma_b[,,j] <- diag(nsnp) # B will be GWAS
    idx <- ceiling(nsnp/2) + -2:2
    Sigma_a[idx,idx,j] <- ifelse(Sigma_a[idx,idx,j] == 0, .5, 1)
    Sigma_b[idx,idx,j] <- ifelse(Sigma_b[idx,idx,j] == 0, .5, 1)
  }
  Sigma_a[6:15,6:15,1] <- 0
  Sigma_a[11:15,11:15,2] <- 0
  Sigma_b[6:15,6:15,1] <- 0
  Sigma_b[11:15,11:15,2] <- 0
  x <- (nsnp - 1)/2
  beta <- sapply(1:ncond, function(j) rep(c(0,6 + j,0),c(x,1,x)))
  beta_hat_a <- beta
  beta_hat_b <- beta
  se_a <- se_b <- matrix(0.2, nsnp, ncond)
  gamma <- 0.5 
  sigma <- 0.1
  theta <- 1 - 1/nsnp
  mu <- mean(beta[x+1,])
  for (j in 1:ncond) {
    beta_a_j <- beta[,j]
    beta_hat_a[,j] <- mvrnorm(1,
                              mu=Sigma_a[,,j] %*% beta_a_j,
                              diag(se_a[,j]) %*% Sigma_a[,,j] %*% diag(se_a[,j]))
    beta_b_j <- gamma * beta_a_j + ifelse(beta_a_j==0,0,rnorm(nsnp,0,sigma))
    beta_hat_b[,j] <- mvrnorm(1,
                              mu=Sigma_b[,,j] %*% beta_b_j,
                              diag(se_b[,j]) %*% Sigma_b[,,j] %*% diag(se_b[,j]))
  }
  beta_hat_a[6:15,1] <- 0
  beta_hat_a[11:15,2] <- 0
  beta_hat_b[6:15,1] <- 0
  beta_hat_b[11:15,2] <- 0
  se_a[6:15,1] <- 0
  se_a[11:15,2] <- 0
  se_b[6:15,1] <- 0
  se_b[11:15,2] <- 0
  data <- list(nsnp=nsnp,
               ncond=ncond,
               n=n,
               beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               se_a=se_a,
               se_b=se_b,
               Sigma_a=Sigma_a,
               Sigma_b=Sigma_b)

  options(mc.cores=2)
  fit1 <- fitBetaEcaviar(data)

  print(fit1, pars=c("beta_a[1,1]","beta_b[1,1]"), digits=3)
  
  rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp,",1]"))
  rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp,",1]"))
  rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp,",2]"))
  rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp,",2]"))
  rstan::stan_plot(fit1, pars=paste0("beta_a[",1:nsnp,",3]"))
  rstan::stan_plot(fit1, pars=paste0("beta_b[",1:nsnp,",3]"))

  library(matrixStats)
  coefs1 <- rstan::extract(fit1)
  beta_clean_a <- sapply(1:ncond, function(j) colMeans(coefs1$beta_a[,,j]))
  beta_clean_b <- sapply(1:ncond, function(j) colMeans(coefs1$beta_b[,,j]))
  beta_sd_a <- sapply(1:ncond, function(j) colSds(coefs1$beta_a[,,j]))
  beta_sd_b <- sapply(1:ncond, function(j) colSds(coefs1$beta_b[,,j]))
  
  data <- list(nsnp=nsnp,
               ncond=ncond,
               beta_hat_a=beta_clean_a,
               beta_hat_b=beta_clean_b,
               sd_a=beta_sd_a,
               sd_b=beta_sd_b)

  options(mc.cores=2)
  fit2 <- fitBetaMixture(data)

  # naive
  a <- beta_clean_a[x+1,]
  b <- beta_clean_b[x+1,]
  fit <- lm(b ~ 0 + a)

  print(fit2, pars=c("theta","mu","gamma","sigma_1b"), digits=3)
  data.frame(par=c("theta","mu","gamma","sigma_1b"),
             value=c(theta, mu, gamma, sigma))
  data.frame(par=c("mu","gamma","sigma_1b"),
             value=c(mean(a), coef(fit), summary(fit)$sigma))
  
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