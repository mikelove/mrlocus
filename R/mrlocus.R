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
#' @param sd_beta prior SD for beta A
#' @param mu_alpha prior mean for alpha
#' @param sd_alpha prior SD for alpha
#' @param sd_sigma prior SD for sigma
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitSlope <- function(res,
                     sd_beta=NULL,
                     mu_alpha=NULL,
                     sd_alpha=NULL,
                     sd_sigma=1,
                     ...) {
  
  n <- length(res$beta_hat_a)
  stopifnot(length(res$beta_hat_b) == n)
  stopifnot(length(res$sd_a) == n)
  stopifnot(length(res$sd_b) == n)
  stopifnot(sd_beta > 0)
  stopifnot(sd_alpha > 0)
  stopifnot(sd_sigma > 0)

  # specify SD for prior for beta A
  if (is.null(sd_beta)) {
    sd_beta <- 2 * max(abs(res$beta_hat_a))
  }
  
  # specify prior for alpha
  lmfit <- lm(res$beta_hat_b ~ res$beta_hat_a + 0)
  lmsum <- summary(lmfit)$coefficients
  if (is.null(mu_alpha)) {
    mu_alpha=lmsum[1,1]
  }
  if (is.null(sd_alpha)) {
    sd_alpha=2*abs(lmsum[1,1])
  }

  if (n > 1) {
    data <- list(n=n,
                 beta_hat_a=res$beta_hat_a,
                 beta_hat_b=res$beta_hat_b,
                 sd_a=res$sd_a,
                 sd_b=res$sd_b,
                 sd_beta=sd_beta,
                 mu_alpha=mu_alpha,
                 sd_alpha=sd_alpha,
                 sd_sigma=sd_sigma)
    stanfit <- rstan::sampling(stanmodels$slope, data, ...)
    out <- list(stanfit=stanfit)
  } else {
    # parametric simulation if just one SNP
    m <- 1e5
    slope <- rnorm(m, res$beta_hat_b, res$sd_b) /
      rnorm(m, res$beta_hat_a, res$sd_a)
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
