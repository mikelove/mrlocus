#' Fit the beta colocalization per signal cluster
#'
#' @param beta_hat_a vector of length nsnp, estimated coefficients for A
#' @param beta_hat_b " " for B
#' @param se_a vector of length nsnp, standard errors for beta_hat_a
#' @param se_b " " for beta_hat_b
#' @param Sigma_a correlation matrix of SNPs for A, dimension should be nsnp x nsnp
#' @param Sigma_b " " for B (this could be different for different LD structures)
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitBetaColoc <- function(nsnp, beta_hat_a, beta_hat_b,
                         se_a, se_b, Sigma_a, Sigma_b, ...) {
  n <- length(beta_hat_a)
  stopifnot(length(beta_hat_b) == n)
  stopifnot(length(se_a) == n)
  stopifnot(length(se_b) == n)
  stopifnot(dim(Sigma_a) == c(n,n))
  stopifnot(dim(Sigma_b) == c(n,n))
  # pick out largest z score in 'a'
  idx <- which.max(beta_hat_a/se_a)
  # scale 'a' to be around 1
  scale_a <- 1/abs(beta_hat_a)[idx]
  scaled_beta_hat_a <- scale_a * beta_hat_a
  scaled_se_a <- scale_a * se_a
  # scale 'b' to be around 1
  scale_b <- 1/abs(beta_hat_b)[idx]
  scaled_beta_hat_b <- scale_b * beta_hat_b
  scaled_se_b <- scale_b * se_b
  data <- list(n=n, beta_hat_a=scaled_beta_hat_a,
               beta_hat_b=scaled_beta_hat_b,
               se_a=scaled_se_a, se_b=scaled_se_b,
               Sigma_a=Sigma_a, Sigma_b=Sigma_b)
  stanfit <- rstan::sampling(stanmodels$beta_coloc, data, ...)
  coefs <- rstan::extract(stanfit)
  # outgoing estimates
  beta_hat_a <- colMeans(coefs$beta_a) / scale_a
  beta_hat_b <- colMeans(coefs$beta_b) / scale_b
  sd_a <- matrixStats::colSds(coefs$beta_a) / scale_a
  sd_b <- matrixStats::colSds(coefs$beta_b) / scale_b
  list(stanfit=stanfit, 
       beta_hat_a=beta_hat_a, beta_hat_b=beta_hat_b,
       sd_a=sd_a, sd_b=sd_b,
       scale_a=scale_a, scale_b=scale_b)
}

#' Fit the beta slope model: effect of A on B
#'
#' @param res list with the following named elements:
#' \itemize{
#' \item beta_hat_a - vector of length sum(nsnp), first step point estimates of beta for A
#' \item beta_hat_b - " " for B
#' \item sd_a - vector of length sum(nsnp), first step posterior SD for beta for A
#' \item sd_b - " " for B
#' }
#' @param alpha_mu prior mean for alpha
#' @param alpha_sd prior SD for alpha
#' @param sigma_sd prior SD for sigma
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitSlope <- function(res,
                     alpha_mu=NULL,
                     alpha_sd=NULL,
                     sigma_sd=1,
                     ...) {
  n <- length(res$beta_hat_a)
  stopifnot(length(res$beta_hat_b) == n)
  stopifnot(length(res$sd_a) == n)
  stopifnot(length(res$sd_b) == n)
  stopifnot(alpha_sd > 0)
  stopifnot(sigma_sd > 0)

  # specify prior for alpha
  naive <- unname(coef(lm(res$beta_hat_b ~ res$beta_hat_a + 0)))
  if (is.null(alpha_mu)) {
    alpha_mu=naive[1]
  }
  if (is.null(alpha_sd)) {
    alpha_sd=abs(naive[1])/2
  }

  if (n > 1) {
    data <- list(n=n,
                 beta_hat_a=res$beta_hat_a,
                 beta_hat_b=res$beta_hat_b,
                 sd_a=res$sd_a,
                 sd_b=res$sd_b,
                 alpha_mu=alpha_mu,
                 alpha_sd=alpha_sd,
                 sigma_sd=sigma_sd)
    stanfit <- rstan::sampling(stanmodels$slope, data, ...)
    out <- list(stanfit=stanfit)
  } else {
    # parametric simulation if just one SNP
    m <- 1e5
    slope <- rnorm(m, beta_hat_b, sd_b)/rnorm(m, beta_hat_a, sd_a)
    est <- c(median(slope), mad(slope))
    out <- list(est=est)
  }
  out <- c(list(beta_hat_a=res$beta_hat_a,
                beta_hat_b=res$beta_hat_b,
                sd_a=res$sd_a,
                sd_b=res$sd_b),
           out)

  return(out)
}
