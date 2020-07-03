context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  set.seed(1)
  out <- makeSimDataForMrlocus()
  plotInitEstimates(out)

  # colocalization:
  fit <- list()
  nsnp <- lengths(out$beta_hat_a)
  nclust <- length(nsnp)
  for (j in 1:nclust) {
    print(paste("----",j,"----"))
    fit[[j]] <- with(out, 
                fitBetaColoc(nsnp=nsnp[j],
                beta_hat_a=beta_hat_a[[j]], beta_hat_b=beta_hat_b[[j]],
                se_a=se_a[[j]], se_b=se_b[[j]],
                Sigma_a=Sigma_a[[j]], Sigma_b=Sigma_b[[j]],
                verbose=FALSE, open_progress=FALSE,
                show_messages=FALSE, refresh=-1))
  }
  j <- 1
  rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_a[",1:nsnp[j],"]"))
  rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_b[",1:nsnp[j],"]"))
  # extract results from colocalization for slope fitting:
  res <- list(beta_hat_a=lapply(fit, `[[`, "beta_hat_a"),
              beta_hat_b=lapply(fit, `[[`, "beta_hat_b"),
              sd_a=out$se_a,
              sd_b=out$se_b)
  res <- extractForSlope(res, plot=TRUE)  
  res <- fitSlope(res, iter=10000)

  # print the estimated slope
  library(rstan)
  print(res$stanfit, pars=c("alpha","sigma"), digits=3)

  # plot
  plotMrlocus(res, main="mrlocus")
  
})
