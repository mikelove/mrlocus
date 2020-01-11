#' Fit the beta eCAVIAR model
#'
#' @param data a list with the following named elements:
#' \itemize{
#' \item{nsnp - number of SNPs}
#' \item{ncond - number of conditional signal clusters}
#' \item{beta_hat_a - matrix nsnp x ncond, estimated coefficients for A}
#' \item{beta_hat_b - " " for B}
#' \item{se_a - matrix nsnp x ncond, standard errors for beta_hat_a}
#' \item{se_b - " " for beta_hat_b}
#' \item{Sigma_a - correlation matrix of SNPs for A}
#' \item{Sigma_b - " " for B}
#' }
#' 
#' @export
fitBetaEcaviar <- function(data) {
  checkData(data, 1)
  rstan::sampling(stanmodels$beta_ecaviar, data)
}

#' Fit the beta mixture model
#'
#' @param data a list with the following named elements:
#' \itemize{
#' \item{nsnp - number of SNPs}
#' \item{ncond - number of conditional signal clusters}
#' \item{beta_hat_a - matrix nsnp x ncond, first step point estimates of beta for A}
#' \item{beta_hat_b - " " for B}
#' }
#' 
#' @export
fitBetaMixture <- function(data) {
  checkData(data, 2)
  rstan::sampling(stanmodels$beta_mixture, data)
}

checkData <- function(data, fit) {
  if (fit == 1) {
    exp.args <- c("nsnp","ncond",
                  "beta_hat_a","beta_hat_b",
                  "se_a","se_b",
                  "Sigma_a","Sigma_b")
  } else if (fit == 2) {
    exp.args <- c("nsnp","ncond", "beta_a", "beta_b")
  }
  if (!all(exp.args %in% names(data))) {
    stop(
      paste("missing argument(s) in 'data':",
            paste(exp.args[!exp.args %in% names(data)], collapse=" "))
    )
  }
  if (fit == 1) {
    stopifnot(dim(data$beta_hat_a) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$beta_hat_b) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$se_a) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$se_b) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$Sigma_a) == c(data$nsnp,data$nsnp,data$ncond))
    stopifnot(dim(data$Sigma_b) == c(data$nsnp,data$nsnp,data$ncond))
  } else if (fit == 2) {
    stopifnot(dim(data$beta_a) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$beta_b) == c(data$nsnp,data$ncond))
  }
}
