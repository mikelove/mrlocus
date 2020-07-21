context("sigma")
library(mrlocus)

test_that("benchmark sigma estimation", {

  if (FALSE) {

    library(pbapply)
    library(ggplot2)
    n <- 4
    sigmas <- rep(1:4 * .25, each=n)

    ests <- pblapply(seq_along(sigmas), function(i) {
      set.seed(i)
      out <- makeSimDataForMrlocus(nsnp=rep(7,6),
                                   idx=4,
                                   alpha=1,
                                   sigma=sigmas[i],
                                   betas=1:6,
                                   se=0.1)
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
      rstan::summary(res$stanfit, pars="sigma",
                     probs=c(.1,.9))$summary["sigma",c("mean","10%","90%")]
    }, cl=4)

    dat <- as.data.frame(do.call(rbind, ests))
    names(dat) <- c("estimate","ymin","ymax")
    dat$rep <- factor(rep(1:n,times=nrow(dat)/n))
    dat$sigma <- sigmas
    dat$cover <- factor(dat$sigma > dat$ymin & dat$sigma < dat$ymax, c("FALSE","TRUE"))

    png(file="~/Desktop/sigma_est.png", width=800, height=800, res=125)
    ggplot(dat, aes(sigma, estimate, ymin=ymin, ymax=ymax, group=rep, col=cover)) +
      geom_pointrange(position=position_dodge(width=.1)) +
      scale_color_manual(values=c("FALSE"="red","TRUE"="black")) +
      geom_abline(slope=1,intercept=0,alpha=0.25)
    dev.off()
    
  }
})
