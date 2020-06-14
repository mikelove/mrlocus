#' Fit the beta eCAVIAR model for a signal cluster
#'
#' @param nsnp number of SNPs
#' @param beta_hat_a vector of length nsnp, estimated coefficients for A
#' @param beta_hat_b " " for B
#' @param se_a vector of length nsnp, standard errors for beta_hat_a
#' @param se_b " " for beta_hat_b
#' @param Sigma_a correlation matrix of SNPs for A, dimension should be nsnp x nsnp
#' @param Sigma_b " " for B (this could be different for different LD structures)
#' 
#' @export
fitBetaEcaviar <- function(nsnp, beta_hat_a, beta_hat_b,
                           se_a, se_b, Sigma_a, Sigma_b) {
  stopifnot(length(beta_hat_a) == nsnp)
  stopifnot(length(beta_hat_b) == nsnp)
  stopifnot(length(se_a) == nsnp)
  stopifnot(length(se_b) == nsnp)
  stopifnot(dim(Sigma_a) == c(nsnp,nsnp))
  stopifnot(dim(Sigma_b) == c(nsnp,nsnp))
  data <- list(nsnp=nsnp, beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               se_a=se_a, se_b=se_b,
               Sigma_a=Sigma_a, Sigma_b=Sigma_b)
  rstan::sampling(stanmodels$beta_marg_ecaviar, data)
}

#' Fit the beta mixture model
#'
#' @param nsnp vector, number of SNPs per signal cluster
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
#' 
#' @export
fitBetaMixture <- function(nsnp,
                           beta_hat_a, beta_hat_b, sd_a, sd_b,
                           sigma_0a=0.5, sigma_0b=0.5, sigma_1a=2,
                           alpha_sd=1, mu_loc=8, mu_sd=2, sigma1b_sd=3) {
  tot <- sum(nsnp)
  stopifnot(length(beta_hat_a) == tot)
  stopifnot(length(beta_hat_b) == tot)
  stopifnot(length(sd_a) == tot)
  stopifnot(length(sd_b) == tot) 
  data <- list(tot=tot,
               beta_hat_a=beta_hat_a,
               beta_hat_b=beta_hat_b,
               sd_a=sd_a, sd_b=sd_b,
               sigma_0a=sigma_0a,
               sigma_0b=sigma_0b,
               sigma_1a=sigma_1a,
               mu_loc=mu_loc, mu_sd=mu_sd,
               sigma1b_sd=sigma1b_sd)
  rstan::sampling(stanmodels$beta_mixture, data)
}
