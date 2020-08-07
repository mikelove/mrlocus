context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  set.seed(1)
  out <- makeSimDataForMrlocus(n_mult=10)
  #set.seed(2)
  #out <- makeSimDataForMrlocus(n_mult=10,sigma=.5)
  #set.seed(2)
  #out <- makeSimDataForMrlocus(alpha=0,n_mult=10,sigma=.5)
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
  #j <- 1
  #rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_a[",1:nsnp[j],"]"))
  #rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_b[",1:nsnp[j],"]"))
  # extract results from colocalization for slope fitting:
  res <- list(beta_hat_a=lapply(fit, `[[`, "beta_hat_a"),
              beta_hat_b=lapply(fit, `[[`, "beta_hat_b"),
              sd_a=out$se_a,
              sd_b=out$se_b)
  res <- extractForSlope(res, plot=TRUE)  
  res <- fitSlope(res, iter=10000)

  # print the estimated slope
  library(rstan)
  print(res$stanfit, pars=c("alpha","sigma"), probs=c(.1,.9), digits=3)

  # plot
  plotMrlocus(res)

  ## png(file="~/Desktop/sim.png", width=1500, height=500, res=150)
  ## par(mfrow=c(1,3))
  ## load("~/Desktop/sim_lowdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(mediation with low dispersion)",
  ##             label="Effect size of", legend=FALSE, pointers=TRUE, ylim=c(-2.5,2.5))
  ## load("~/Desktop/sim_highdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(mediation with high dispersion)",
  ##             label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
  ## load("~/Desktop/sim_alpha0_highdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(no mediation)",
  ##             label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
  ## dev.off()

})
