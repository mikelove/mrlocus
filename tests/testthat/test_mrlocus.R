context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  ### simulate data using multivariate normal ###
  
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
  alpha <- 0.5 
  sigma <- 0.1
  mu <- mean(sapply(beta, `[`, x+1))
  for (j in 1:ncond) {
    beta_a_j <- beta[[j]]
    beta_hat_a[[j]] <- mvrnorm(1,
                         mu=Sigma_a[[j]] %*% beta_a_j,
                         diag(se_a[[j]]) %*% Sigma_a[[j]] %*% diag(se_a[[j]]))
    beta_b_j <- alpha * beta_a_j + ifelse(beta_a_j==0,0,rnorm(nsnp[j],0,sigma))
    beta_hat_b[[j]] <- mvrnorm(1,
                         mu=Sigma_b[[j]] %*% beta_b_j,
                         diag(se_b[[j]]) %*% Sigma_b[[j]] %*% diag(se_b[[j]]))
  }

  ### end simulation ###

  # colocalization:
  
  fit1 <- list()
  options(mc.cores=2)
  for (j in 1:ncond) {
    print(paste("\n\n-----",j,"-----\n\n"))
    fit1[[j]] <- fitBetaColoc(nsnp=nsnp[j],
                              beta_hat_a=beta_hat_a[[j]],
                              beta_hat_b=beta_hat_b[[j]],
                              se_a=se_a[[j]],
                              se_b=se_b[[j]],
                              Sigma_a=Sigma_a[[j]],
                              Sigma_b=Sigma_b[[j]])
  }

  j <- 1
  print(fit1[[j]]$stanfit, pars=c(paste0("beta_a[",1:nsnp[j],"]")), digits=3)
  rstan::stan_plot(fit1[[j]]$stanfit, pars=paste0("beta_a[",1:nsnp[j],"]"))
  rstan::stan_plot(fit1[[j]]$stanfit, pars=paste0("beta_b[",1:nsnp[j],"]"))

  # extract results from colocalization and slope fitting:
  
  res <- list(beta_hat_a=lapply(fit1, `[[`, "beta_hat_a"),
              beta_hat_b=lapply(fit1, `[[`, "beta_hat_b"),
              sd_a=lapply(fit1, `[[`, "sd_a"),
              sd_b=lapply(fit1, `[[`, "sd_b"))

  res <- extractForSlope(res, plot=TRUE)
  res <- fitSlope(res, iter=10000)

  # print the estimated slope
  
  print(res$stanfit, pars=c("alpha","sigma"), digits=3)
  
})
