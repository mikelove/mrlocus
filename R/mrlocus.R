#' Fit the beta eCAVIAR model
#'
#' @param data a list with the following named elements:
#' \itemize{
#' \item{nsnp - number of SNPs}
#' \item{beta_hat_a - vector of length nsnp, estimated coefficients for A}
#' \item{beta_hat_b - " " for B}
#' \item{se_a - vector of length nsnp, standard errors for beta_hat_a}
#' \item{se_b - " " for beta_hat_b}
#' \item{Sigma_a - correlation matrix of SNPs for A, dimension should be nsnp x nsnp}
#' \item{Sigma_b - " " for B (this could be different for different LD structures)}
#' }
#' 
#' @export
fitBetaEcaviar <- function(data) {
  checkData(data, 1)
  rstan::sampling(stanmodels$beta_marg_ecaviar, data)
}

#' Fit the beta mixture model
#'
#' @param data a list with the following named elements:
#' \itemize{
#' \item{nsnp - vector, number of SNPs per cluster}
#' \item{beta_hat_a - vector of length sum(nsnp), first step point estimates of beta for A}
#' \item{beta_hat_b - " " for B}
#' \item{sd_a - vector of length sum(nsnp), first step posterior SD for beta for A}
#' \item{sd_b - " " for B}
#' }
#' 
#' @export
fitBetaMixture <- function(data) {
  checkData(data, 2)
  data$tot <- sum(data$nsnp)
  data$nsnp <- NULL
  rstan::sampling(stanmodels$beta_mixture, data)
}

checkData <- function(data, fit) {
  if (fit == 1) {
    exp.args <- c("nsnp",
                  "beta_hat_a","beta_hat_b",
                  "se_a","se_b",
                  "Sigma_a","Sigma_b")
  } else if (fit == 2) {
    exp.args <- c("nsnp",
                  "beta_hat_a", "beta_hat_b",
                  "sd_a", "sd_b")
  }
  if (!all(exp.args %in% names(data))) {
    stop(
      paste("missing argument(s) in 'data':",
            paste(exp.args[!exp.args %in% names(data)], collapse=" "))
    )
  }
  if (fit == 1) {
    stopifnot(length(data$beta_hat_a) == data$nsnp)
    stopifnot(length(data$beta_hat_b) == data$nsnp)
    stopifnot(length(data$se_a) == data$nsnp)
    stopifnot(length(data$se_b) == data$nsnp)
    stopifnot(dim(data$Sigma_a) == c(data$nsnp,data$nsnp))
    stopifnot(dim(data$Sigma_b) == c(data$nsnp,data$nsnp))
  } else if (fit == 2) {
    tot <- sum(data$nsnp)
    stopifnot(length(data$beta_hat_a) == tot)
    stopifnot(length(data$beta_hat_b) == tot)
    stopifnot(length(data$sd_a) == tot)
    stopifnot(length(data$sd_b) == tot) 
  }
}
