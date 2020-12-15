context("sigma")
library(mrlocus)

test_that("benchmark sigma estimation", {

  if (FALSE) {

    library(pbapply)
    library(ggplot2)

    n <- 20
    sigmas <- rep(c(1,2,4,6) * .25, each=3 * n)
    nclusts <- rep(rep(c(4,6,8), each=n), times=4)

    ests <- pblapply(seq_along(sigmas), function(i) {
      set.seed(i)
      nclust <- nclusts[i]
      nsnp <- rep(5,nclust)
      out <- makeSimDataForMrlocus(nsnp=nsnp,
                                   idx=3,
                                   alpha=1,
                                   sigma=sigmas[i],
                                   betas=1:nclust,
                                   se=0.1)
      fit <- list()
      for (j in 1:nclust) {
        cap.out <- capture.output({
          suppressWarnings({
            out2 <- with(out, 
                         fitBetaColoc(beta_hat_a=beta_hat_a[[j]],
                                      beta_hat_b=beta_hat_b[[j]],
                                      se_a=se_a[[j]], se_b=se_b[[j]],
                                      Sigma_a=Sigma_a[[j]], Sigma_b=Sigma_b[[j]],
                                      iter=1000))
          })
        })
        fit[[j]] <- out2[c("beta_hat_a","beta_hat_b")]
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

    #save(sigmas, nclusts, ests, file="~/proj/mrlocus/mrlocusPaper/supp/figs/sigma_bench.rda")
    
    dat <- as.data.frame(do.call(rbind, ests))
    names(dat) <- c("estimate","ymin","ymax")
    dat$rep <- factor(rep(1:n,times=nrow(dat)/n))
    dat$sigma <- sigmas
    dat$nclust <- nclusts
    dat$cover <- factor(dat$sigma > dat$ymin & dat$sigma < dat$ymax, c("FALSE","TRUE"))

    #png(file="~/proj/mrlocus/mrlocusPaper/supp/figs/sigma_est.png", width=1200, height=1200, res=150)
    ggplot(dat, aes(sigma, estimate, ymin=ymin, ymax=ymax, group=rep, col=cover)) +
      geom_pointrange(position=position_dodge(width=.2)) +
      scale_color_manual(values=c("FALSE"="red","TRUE"="black")) +
      geom_abline(slope=1,intercept=0,alpha=0.25) +
      facet_wrap(~nclust, labeller=label_both, ncol=1) +
      scale_x_continuous(breaks=unique(dat$sigma))
    #dev.off()
    
  }
})
