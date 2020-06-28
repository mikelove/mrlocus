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
  out <- rstan::sampling(stanmodels$beta_coloc, data, ...)
  list(stan=out, scale_a=scale_a, scale_b=scale_b)
}

#' Fit the beta slope model: effect of A on B
#'
#' @param beta_hat_a vector of length sum(nsnp), first step point estimates of beta for A
#' @param beta_hat_b " " for B
#' @param sd_a vector of length sum(nsnp), first step posterior SD for beta for A
#' @param sd_b " " for B
#' @param alpha_sd prior SD for alpha
#' @param sigma_sd prior SD for sigma
#' @param gamma_sd prior SD for gamma
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitSlope <- function(beta_hat_a, beta_hat_b, sd_a, sd_b,
                     alpha_sd=1, sigma_sd=1,
                     gamma_sd=1e-3, ...) {
  n <- length(beta_hat_a)
  stopifnot(length(beta_hat_b) == n)
  stopifnot(length(sd_a) == n)
  stopifnot(length(sd_b) == n)
  if (n > 1) {
    data <- list(n=n,
                 beta_hat_a=beta_hat_a,
                 beta_hat_b=beta_hat_b,
                 sd_a=sd_a, sd_b=sd_b,
                 alpha_sd=alpha_sd,
                 sigma_sd=sigma_sd,
                 gamma_sd=gamma_sd)
    rstan::sampling(stanmodels$slope, data, ...)
  } else {
    m <- 1e5
    slope <- rnorm(m, beta_hat_b, sd_b)/rnorm(m, beta_hat_a, sd_a)
    c(median(slope), mad(slope))
  }
}
