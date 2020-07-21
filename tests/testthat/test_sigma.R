context("sigma")
library(mrlocus)

test_that("benchmark sigma estimation", {

  if (FALSE) {

    #ptm <- proc.time()
    library(pbapply)
    sigmas <- rep(1:8 * .1, each=3)

    ests <- pblapply(seq_along(sigmas), function(i) {
      set.seed(i)
      out <- makeSimDataForMrlocus(nsnp=rep(10,5),betas=1:5,sigma=sigmas[i],se=.05)
      fit <- list()
      nsnp <- lengths(out$beta_hat_a)
      nclust <- length(nsnp)
      for (j in 1:nclust) {
        cap.out <- capture.output({
          fit[[j]] <- with(out, 
                           fitBetaColoc(nsnp=nsnp[j],
                                        beta_hat_a=beta_hat_a[[j]],
                                        beta_hat_b=beta_hat_b[[j]],
                                        se_a=se_a[[j]], se_b=se_b[[j]],
                                        Sigma_a=Sigma_a[[j]], Sigma_b=Sigma_b[[j]],
                                        verbose=FALSE, open_progress=FALSE,
                                        show_messages=FALSE, refresh=-1))
        })
      }
      res <- list(beta_hat_a=lapply(fit, `[[`, "beta_hat_a"),
                  beta_hat_b=lapply(fit, `[[`, "beta_hat_b"),
                  sd_a=out$se_a, sd_b=out$se_b)
      res <- extractForSlope(res, plot=FALSE)
      cap.out <- capture.output({
        res <- fitSlope(res, iter=10000)
      })
      rstan::summary(res$stanfit, pars="sigma")$summary["sigma","mean"]
    }, cl=3)

    ests <- unlist(ests)
    dat <- data.frame(sigmas, ests)
    format(dat, digits=2)
  }
})
