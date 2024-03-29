% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{plotMrlocus}
\alias{plotMrlocus}
\title{Plot estimates from MRLocus slope fitting step}
\usage{
plotMrlocus(
  res,
  q = c(0.1, 0.9),
  sigma_mult = 1.28,
  label = "Effect size of",
  a = "eQTL",
  b = "GWAS",
  xlim = NULL,
  ylim = NULL,
  legend = TRUE,
  digits = 3,
  col_slope = "blue",
  col_band = rgb(0, 0, 1, 0.1),
  col_dashed = rgb(0, 0, 1, 0.5),
  ...
)
}
\arguments{
\item{res}{the output from \code{\link{fitSlope}}}

\item{q}{the quantiles of the posterior
to use for drawing the uncertainty on the slope.
The default is an 80 percent interval}

\item{sigma_mult}{multiplier on estimate of sigma
for drawing the dispersion band
(e.g. \code{qnorm(1 - .2/2) ~= 1.28} should include
80 percent of coefficient pairs)}

\item{label}{what preceeds \code{a} and \code{b} in
the x- and y-axis labels}

\item{a}{name of A experiment}

\item{b}{name of B experiment}

\item{xlim}{xlim (if NULL will be set automatically)}

\item{ylim}{ylim (if NULL will be set automatically)}

\item{legend}{logical, whether to show a legend}

\item{digits}{number of digits to show in legend}

\item{col_slope}{the color of the slope (alpha)}

\item{col_band}{the color of the band}

\item{col_dashed}{the color of the dashed lines}

\item{...}{arguments passed to \code{plot}}
}
\description{
Plot estimates from MRLocus slope fitting step
}
