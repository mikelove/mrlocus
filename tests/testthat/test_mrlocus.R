context("mrlocus")
library(mrlocus)

test_that("mrlocus works on simple sim data", {

  set.seed(1)
  data <- makeSimDataForMrlocus(n_mult=10)
  #set.seed(2)
  #data <- makeSimDataForMrlocus(n_mult=10,sigma=.5)
  #set.seed(2)
  #data <- makeSimDataForMrlocus(alpha=0,n_mult=10,sigma=.5)
  plotInitEstimates(data)
  dev.off()

  # step 1) colocalization:
  coloc_fit <- list()
  nsnp <- lengths(data$beta_hat_a)
  nclust <- length(nsnp)
  for (j in 1:nclust) {
    print(paste("----",j,"----"))
    coloc_fit[[j]] <- suppressWarnings({
      with(data, 
           fitBetaColoc(
             beta_hat_a=beta_hat_a[[j]], beta_hat_b=beta_hat_b[[j]],
             se_a=se_a[[j]], se_b=se_b[[j]],
             Sigma_a=Sigma_a[[j]], Sigma_b=Sigma_b[[j]],
             verbose=FALSE, open_progress=FALSE,
             show_messages=FALSE, refresh=-1
           ))
    })
  }
  #j <- 1
  #rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_a[",1:nsnp[j],"]"))
  #rstan::stan_plot(fit[[j]]$stanfit, pars=paste0("beta_b[",1:nsnp[j],"]"))
  # extract results from colocalization for slope fitting:
  res <- list(beta_hat_a=lapply(coloc_fit, `[[`, "beta_hat_a"),
              beta_hat_b=lapply(coloc_fit, `[[`, "beta_hat_b"),
              sd_a=data$se_a,
              sd_b=data$se_b)

  res0 <- res # for another experiment below

  # step 2) slope fitting:
  res <- extractForSlope(res)
  dev.off()

  suppressWarnings({
    res <- fitSlope(res, iter=10000)
  })

  # print the estimated slope
  library(rstan)
  print(res$stanfit, pars=c("alpha","sigma"), probs=c(.1,.9), digits=3)

  # basic prior check
  priorCheck(res)
  dev.off()
  
  # plot
  plotMrlocus(res)
  dev.off()

  plotMrlocus(res, sigma_mult=1)
  dev.off()

  ### FIGURE 1 ###

  ## png(file="fig1_sim.png", width=1500, height=500, res=150)
  ## par(mfrow=c(1,3))
  ## load("sim_lowdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(mediation with low dispersion)",
  ##             label="Effect size of", legend=FALSE, pointers=TRUE, ylim=c(-2.5,2.5))
  ## load("sim_highdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(mediation with high dispersion)",
  ##             label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
  ## load("sim_alpha0_highdisp.rda")
  ## plotMrlocus(res, main="SNPs → gene → trait\n(no mediation)",
  ##             label="Effect size of", legend=FALSE, ylim=c(-2.5,2.5))
  ## dev.off()

  
  # test slope estimation with just a single cluster
  res <- res0
  res <- extractForSlope(res, plot=FALSE)
  res <- lapply(res, `[[`, 4) # just the last
  res <- fitSlope(res, iter=10000)
  names(res)
  res$est

  # extract potentially more than one SNP per cluster with Mclust
  res <- res0
  res <- extractForSlope(res, niter=3, plot=TRUE)
  dev.off()
  
  suppressWarnings({
    res <- fitSlope(res, iter=10000)
  })
  
  print(res$stanfit, pars=c("alpha","sigma"), probs=c(.1,.9), digits=3)
  
})
