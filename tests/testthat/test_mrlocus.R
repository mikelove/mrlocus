context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  library(MASS)
  
  set.seed(5)
  ncond <- 3
  nsnp <- c(10,15,20)
  Sigma_a <- Sigma_b <- list()
  for (j in 1:ncond) {
    Sigma_a[[j]] <- diag(nsnp[j]) # A will be eQTL
    Sigma_b[[j]] <- diag(nsnp[j]) # B will be GWAS
    idx <- 8 + -2:2
    Sigma_a[[j]][idx,idx] <- ifelse(Sigma_a[[j]][idx,idx] == 0, .5, 1)
    Sigma_b[[j]][idx,idx] <- ifelse(Sigma_b[[j]][idx,idx] == 0, .5, 1)
  }
  x <- 7
  y <- 12
  beta <- lapply(1:ncond, function(j) (rep(c(0,6 + j,0),c(x,1,y)))[1:nsnp[j]])
  beta_hat_a <- beta_hat_b <- beta
  se_a <- se_b <- lapply(1:ncond, function(j) rep(0.2, nsnp[j]))
  gamma <- 0.5 
  sigma <- 0.1
  theta <- 1 - ncond/sum(nsnp)
  mu <- mean(sapply(beta, `[`, x+1))
  for (j in 1:ncond) {
    beta_a_j <- beta[[j]]
    beta_hat_a[[j]] <- mvrnorm(1,
                         mu=Sigma_a[[j]] %*% beta_a_j,
                         diag(se_a[[j]]) %*% Sigma_a[[j]] %*% diag(se_a[[j]]))
    beta_b_j <- gamma * beta_a_j + ifelse(beta_a_j==0,0,rnorm(nsnp[j],0,sigma))
    beta_hat_b[[j]] <- mvrnorm(1,
                         mu=Sigma_b[[j]] %*% beta_b_j,
                         diag(se_b[[j]]) %*% Sigma_b[[j]] %*% diag(se_b[[j]]))
  }

  fit1 <- list()
  options(mc.cores=2)
  for (j in 1:ncond) {
    print(paste("\n\n-----",j,"-----\n\n"))
    fit1[[j]] <- fitBetaEcaviar(nsnp=nsnp[j],
                                beta_hat_a=beta_hat_a[[j]],
                                beta_hat_b=beta_hat_b[[j]],
                                se_a=se_a[[j]],
                                se_b=se_b[[j]],
                                Sigma_a=Sigma_a[[j]],
                                Sigma_b=Sigma_b[[j]])
  }

  j <- 1
  print(fit1[[j]], pars=c(paste0("beta_a[",1:nsnp[j],"]")), digits=3)
  rstan::stan_plot(fit1[[j]], pars=paste0("beta_a[",1:nsnp[j],"]"))
  rstan::stan_plot(fit1[[j]], pars=paste0("beta_b[",1:nsnp[j],"]"))

  library(matrixStats)
  coefs1 <- lapply(fit1, function(x) rstan::extract(x))
  beta_clean_a <- do.call(c, lapply(coefs1, function(x) colMeans(x$beta_a)))
  beta_clean_b <- do.call(c, lapply(coefs1, function(x) colMeans(x$beta_b)))
  beta_sd_a <- do.call(c, lapply(coefs1, function(x) colSds(x$beta_a)))
  beta_sd_b <- do.call(c, lapply(coefs1, function(x) colSds(x$beta_b)))

  options(mc.cores=2)
  fit2 <- fitBetaMixture(nsnp=nsnp,
                         beta_hat_a=beta_clean_a,
                         beta_hat_b=beta_clean_b,
                         sd_a=beta_sd_a,
                         sd_b=beta_sd_b, 
                         sigma_0a=.1,
                         sigma_0b=.1)
  
  # naive
  a <- beta_clean_a[beta_clean_a > 4]
  b <- beta_clean_b[beta_clean_b > 2]
  plot(a,b)
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

  rstan::stan_plot(fit2, pars=paste0("beta_a[",1:sum(nsnp),"]"))
  rstan::stan_plot(fit2, pars=paste0("beta_b[",1:sum(nsnp),"]"))
  
})
