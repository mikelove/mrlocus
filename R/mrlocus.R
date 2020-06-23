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
  # scale 'b' to match the scale of 'a'
  scale_b <- max(abs(beta_hat_a))/max(abs(beta_hat_b))
  scaled_beta_hat_b <- scale_b * beta_hat_b
  scaled_se_b <- scale_b * se_b
  data <- list(n=n, beta_hat_a=beta_hat_a,
               beta_hat_b=scaled_beta_hat_b,
               se_a=se_a, se_b=scaled_se_b,
               Sigma_a=Sigma_a, Sigma_b=Sigma_b)
  out <- rstan::sampling(stanmodels$beta_coloc, data, ...)
  list(stan=out, scale_b=scale_b)
}

#' Fit the mixture model for null and large effects
#'
#' @param beta_hat_a vector of length sum(nsnp), first step point estimates of beta for A
#' @param beta_hat_b " " for B
#' @param sd_a vector of length sum(nsnp), first step posterior SD for beta for A
#' @param sd_b " " for B
#' @param sigma_0a prior SD of the null component for A
#' @param sigma_0b " " for B
#' @param sigma_1a prior SD of the non-null component for beta for A (SD for B is fitted parameter)
#' @param alpha_sd prior SD for alpha
#' @param mu_loc center of prior for mu
#' @param mu_sd prior SD for mu
#' @param sigma1b_sd prior SD for sigma_1b
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitMixture <- function(beta_hat_a, beta_hat_b,
                       sd_a, sd_b,
                       sigma_0a, sigma_0b, sigma_1a,
                       alpha_sd, mu_loc, mu_sd,
                       sigma1b_sd, ...) {
  n <- length(beta_hat_a)
  stopifnot(length(beta_hat_b) == n)
  stopifnot(length(sd_a) == n)
  stopifnot(length(sd_b) == n) 
  data <- list(n=n,
               beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               sd_a=sd_a, sd_b=sd_b,
               sigma_0a=sigma_0a,
               sigma_0b=sigma_0b,
               sigma_1a=sigma_1a,
               alpha_sd=alpha_sd,
               mu_loc=mu_loc, mu_sd=mu_sd,
               sigma1b_sd=sigma1b_sd)
  rstan::sampling(stanmodels$mixture, data, ...)
}

#' Fit the beta slope model: effect of A on B
#'
#' @param beta_hat_a vector of length sum(nsnp), first step point estimates of beta for A
#' @param beta_hat_b " " for B
#' @param sd_a vector of length sum(nsnp), first step posterior SD for beta for A
#' @param sd_b " " for B
#' @param alpha_sd prior SD for alpha
#' @param sigma_sd prior SD for sigma
#' @param ... further arguments passed to rstan::sampling
#' 
#' @export
fitSlope <- function(beta_hat_a, beta_hat_b, sd_a, sd_b,
                     alpha_sd=1, sigma_sd=1, ...) {
  n <- length(beta_hat_a)
  stopifnot(length(beta_hat_b) == n)
  stopifnot(length(sd_a) == n)
  stopifnot(length(sd_b) == n) 
  data <- list(n=n,
               beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               sd_a=sd_a, sd_b=sd_b,
               alpha_sd=alpha_sd,
               sigma_sd=sigma_sd)
  rstan::sampling(stanmodels$slope, data, ...)
}
