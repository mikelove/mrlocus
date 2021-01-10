#' Fit the colocalization model per signal cluster
#'
#' Implements the MRLocus colocalization step, in which summary
#' statistics from A and B studies are used in a generative
#' model, with a horseshoe prior on the latent true effect sizes.
#' Posterior effect sizes and posterior SD are returned,
#' although in pratice we use the original SE in the subsequent
#' \code{\link{fitSlope}} model.
#' See Supplementary Methods of the MRLocus manuscript.
#' 
#' @param beta_hat_a vector of estimated coefficients for A
#' @param beta_hat_b " " for B
#' @param se_a vector of standard errors for beta_hat_a
#' @param se_b " " for beta_hat_b
#' @param Sigma_a correlation matrix of SNPs for A
#' @param Sigma_b " " for B (this could be different for different LD matrix)
#' @param ... further arguments passed to rstan::sampling
#'
#' @importFrom matrixStats colSds
#' @importFrom stats as.dist cutree hclust kmeans lm mad median rnorm runif sd
#' @importFrom graphics abline arrows par points polygon segments text
#' @importFrom grDevices rgb
#' @importFrom utils packageVersion
#' @importFrom rstan summary
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{stanfit} object
#' \item posterior means for estimated coefficients for A and B
#' \item posterior standard deviations for A and B
#' \item scaling factors for A and B
#' }
#' Two important notes: (1) in MRLocus manuscript, original
#' SE are used instead of posterior SD in the slope fitting step,
#' (2) the posterior means and SD for estimated coefficients 
#' are appropriately scaled, while the results from the
#' \code{stanfit} object are not scaled. In order to scale the
#' results from the stanfit object, \code{scale_a} and \code{scale_b}
#' should be divided out from both coefficients and SDs
#' (see Supplementary Methods).
#'
#' @references
#'
#' Anqi Zhu*, Nana Matoba*, Emmaleigh Wilson, Amanda L. Tapia, Yun Li,
#' Joseph G. Ibrahim, Jason L. Stein, Michael I. Love.
#' MRLocus: identifying causal genes mediating a trait through
#' Bayesian estimation of allelic heterogeneity. (2020) bioRxiv
#' \url{https://doi.org/10.1101/2020.08.14.250720}
#' 
#' @export
fitBetaColoc <- function(beta_hat_a, beta_hat_b,
                         se_a, se_b, Sigma_a, Sigma_b, ...) {
  n <- length(beta_hat_a)
  if (n == 1) {
    stop("colocalization not needed with n=1 SNP in cluster
  pass the effect sizes and SE's directly to slope estimation")
  }
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

#' Fit the mediation slope model: effect of A on B
#'
#' Implements the MRLocus slope fitting step, in which the estimated
#' coefficients and their original SE are used to determine
#' the mediation slope (\code{alpha}), and the dispersion of individual
#' signal clusters around the slope (\code{sigma}). This function
#' follows the colocalization step \code{\link{fitBetaColoc}}
#' and \code{\link{extractForSlope}}.
#' The output \code{fitSlope} can be visualized with
#' \code{\link{plotMrlocus}}. For details on the model,
#' see Supplementary Methods of the MRLocus manuscript.
#' See vignette for example of model interpretation.
#'
#' Note that if summary statistics for only one SNP are provided
#' a warning will be printed (this is not a recommended use
#' of MRLocus) and a parametric simulation is used to estimate the
#' slope, instead of the Bayesian model.
#' 
#' @param res list with the following named elements:
#' \itemize{
#' \item \code{beta_hat_a} - point estimates of coefficients for A from colocalization
#' \item \code{beta_hat_b} - " " for B
#' \item \code{sd_a} - sampling SD for \code{beta_hat_a} (in practice original
#' SE are provided here)
#' \item \code{sd_b} - " " for \code{beta_hat_b} " "
#' \item {alleles} (optional) data.frame with allele information
#' }
#' @param sd_beta prior SD for beta A (default value will be derived from data)
#' @param mu_alpha prior mean for \code{alpha} (default value will be derived from data)
#' @param sd_alpha prior SD for \code{alpha} (default value will be derived from data)
#' @param sd_sigma prior SD for \code{sigma} (default value of 1)
#' @param ... further arguments passed to \code{rstan::sampling}
#'
#' @return list with the following elements: \code{stanfit} object,
#' original estimated coefficients and standard deviations,
#' as well as the \code{alleles} data.frame (if it was provided)
#' 
#' @references
#'
#' Anqi Zhu*, Nana Matoba*, Emmaleigh Wilson, Amanda L. Tapia, Yun Li,
#' Joseph G. Ibrahim, Jason L. Stein, Michael I. Love.
#' MRLocus: identifying causal genes mediating a trait through
#' Bayesian estimation of allelic heterogeneity. (2020) bioRxiv
#' \url{https://doi.org/10.1101/2020.08.14.250720}
#' 
#' @export
fitSlope <- function(res,
                     sd_beta=NULL,
                     mu_alpha=NULL,
                     sd_alpha=NULL,
                     sd_sigma=NULL,
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
    mu_alpha <- lmsum[1,1]
  }
  if (is.null(sd_alpha)) {
    sd_alpha <- 2*abs(lmsum[1,1])
  }

  # specify SD for prior for sigma
  if (is.null(sd_sigma)) {
    # a wide prior based on SD of B coefs,
    # or max B coef, whichever is larger
    sd_sigma <- max( 2*sd(res$beta_hat_b), max(abs(res$beta_hat_b)) )
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
    priors <- list(sd_beta=sd_beta,
                   mu_alpha=mu_alpha,
                   sd_alpha=sd_alpha,
                   sd_sigma=sd_sigma)
    out <- list(stanfit=stanfit, priors=priors)
  } else {
    # parametric simulation if just one SNP
    warning("It is recommended to input more than one signal cluster to MRLocus,
  but only one signal cluster was provided to `fitSlope`.
  Parametric simulation will be used to estimate slope")
    m <- 1e5
    slope <- rnorm(m, res$beta_hat_b, res$sd_b) /
      rnorm(m, res$beta_hat_a, res$sd_a)
    est <- c(median(slope), mad(slope))
    out <- list(est=est)
  }
  out <- c(out,
           list(beta_hat_a=res$beta_hat_a,
                beta_hat_b=res$beta_hat_b,
                sd_a=res$sd_a,
                sd_b=res$sd_b))
  if ("alleles" %in% names(res)) {
    out$alleles <- res$alleles
  }
  out$mrlocus_version <- packageVersion("mrlocus")
  return(out)
}
