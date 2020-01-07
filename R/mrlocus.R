#' Fit the beta eCAVIAR model
#'
#' @export
fitBetaEcaviar <- function(data) {
  checkData(data, 1)
  rstan::sampling(stanmodels$beta_ecaviar, data)
}

#' Fit the beta mixture model
#'
#' @export
fitBetaMixture <- function(data) {
  checkData(data, 2)
  rstan::sampling(stanmodels$beta_mixture, data)
}

checkData <- function(data, fit) {
  if (fit == 1) {
    exp.args <- c("nsnp","ncond","beta_hat_a","beta_hat_b",
                  "Sigma_a","Sigma_b","se2")
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
    stopifnot(dim(data$Sigma_a) == c(data$nsnp,data$nsnp,data$ncond))
    stopifnot(dim(data$Sigma_b) == c(data$nsnp,data$nsnp,data$ncond))
  } else if (fit == 2) {
    stopifnot(dim(data$beta_a) == c(data$nsnp,data$ncond))
    stopifnot(dim(data$beta_b) == c(data$nsnp,data$ncond))
  }
}
