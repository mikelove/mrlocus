% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrlocus.R
\name{fitSlope}
\alias{fitSlope}
\title{Fit the mediation slope model: effect of A on B}
\usage{
fitSlope(
  res,
  sd_beta = NULL,
  mu_alpha = NULL,
  sd_alpha = NULL,
  sd_sigma = NULL,
  ...
)
}
\arguments{
\item{res}{list with the following named elements:
\itemize{
\item \code{beta_hat_a} - point estimates of coefficients for A from colocalization
\item \code{beta_hat_b} - " " for B
\item \code{sd_a} - sampling SD for \code{beta_hat_a} (in practice original
SE are provided here)
\item \code{sd_b} - " " for \code{beta_hat_b} " "
\item {alleles} (optional) data.frame with allele information
}}

\item{sd_beta}{prior SD for beta A (default value will be derived from data)}

\item{mu_alpha}{prior mean for \code{alpha} (default value will be derived from data)}

\item{sd_alpha}{prior SD for \code{alpha} (default value will be derived from data)}

\item{sd_sigma}{prior SD for \code{sigma} (default value of 1)}

\item{...}{further arguments passed to \code{rstan::sampling}}
}
\value{
list with the following elements: \code{stanfit} object,
original estimated coefficients and standard deviations,
as well as the \code{alleles} data.frame (if it was provided)
}
\description{
Implements the MRLocus slope fitting step, in which the estimated
coefficients and their original SE are used to determine
the mediation slope (\code{alpha}), and the dispersion of individual
signal clusters around the slope (\code{sigma}). This function
follows the colocalization step \code{\link{fitBetaColoc}}
and \code{\link{extractForSlope}}.
The output \code{fitSlope} can be visualized with
\code{\link{plotMrlocus}}. For details on the model,
see Supplementary Methods of the MRLocus manuscript.
See vignette for example of model interpretation.
}
\details{
Note that if summary statistics for only one SNP are provided
a warning will be printed (this is not a recommended use
of MRLocus) and a parametric simulation is used to estimate the
slope, instead of the Bayesian model.
}
\references{
Anqi Zhu*, Nana Matoba*, Emma P. Wilson, Amanda L. Tapia, Yun Li,
Joseph G. Ibrahim, Jason L. Stein, Michael I. Love.
MRLocus: identifying causal genes mediating a trait through
Bayesian estimation of allelic heterogeneity. (2021) PLOS Genetics
\url{https://doi.org/10.1371/journal.pgen.1009455}
}
